#!/usr/bin/env nextflow

params.input = '/lustre/scratch/users/falko.hofmann/isoseq/test/*.bam'
params.outdir = '/lustre/scratch/users/falko.hofmann/isoseq/test/results'
params.primers = '/lustre/scratch/users/falko.hofmann/pipelines/isoseq3/primers.fasta'
params.genome = 'tair10'
params.annotation = params.genome ? params.genomes[ params.genome ].annotation ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.index = params.genome ? params.genomes[ params.genome ].index ?: false : false
params.intron_max = params.genome ? params.genomes[ params.genome ].intron_max ?: false : false
params.transcript_max = params.genome ? params.genomes[ params.genome ].transcript_max ?: false : false


params.intron_max  = 6000
params.transcript_max = 60000

Channel
    .fromPath(params.input)
    .ifEmpty { error "Cannot find any matching bam files: $params.input" }
    .into {input_ccs; input_polish}

Channel
    .fromPath(params.input + '.pbi')
    .ifEmpty { error "Cannot find matching bam.pbi files: $params.input" }
    .set {input_pbi}

Channel
    .fromPath(params.primers)
    .ifEmpty { error "Cannot find primer file: $params.primers" }
    .set {primers_file}


process make_paths{
    input:
    var path from params.input

    output:

    """
    """
    
}


process run_ccs{

        tag "ccs: $name"

        publishDir "$params.outdir/ccs", mode: 'copy'

        input:
        file ccs_input from input_ccs

        output:
        file '$name.ccs.*'
        file 'ccs_report.txt'
        file "$name.ccs.bam" into ccs_out

        """
        time ccs $ccs_input $name.ccs.bam --noPolish --minPasses 1
        """
}


process run_lima{

    tag "lima: $name"

    publishDir "$params.outdir/lima", mode: 'copy'

    input:
    file ccs_bam from ccs_out
    file primers from primers_file

    output:
    file '$name.demux.ccs.*'
    file '$name.demux.ccs.primer_5p--primer_3p.bam' into lima_out
    
    """
    time lima $ccs_bam $primers $name.demux.ccs.bam --isoseq --no-pbi --dump-clips --dump-removed
    """

}


process cluster_reads{

    tag "clustering : $name"

    publishDir "$params.outdir/cluster", mode: 'copy'

    input:
    file lima_demux from lima_out
    
    output:
    file '$name.unpolished.*'
    file '$name.unpolished.bam' into cluster_out

    """
    time isoseq3 cluster $lima_demux $name.unpolished.bam --require-polya --verbose 
    """
}


process polish_reads{
    
    tag "polishing : $name"

    publishDir "$params.outdir/polish", mode: 'copy'

    input:
    file pbi from input_pbi
    file cluster_bam from cluster_out
    file all_reads_bam from input_polish

    output:
    file '$name.polished.*'
    file '$name.polished.hq.fastq.gz' into polish_out

    """
    time isoseq3 polish $cluster_bam $all_reads_bam polished.bam
    """

}


process build_index{

    tag "STARlong index: $params.index"

    storeDir "$params.index"

    input:
    file annotation from params.annotation
    file fasta from params.fasta

    output:
    file "${params.genome}.*" into star_index
    file 'Log.out' into log

    """
    mkdir -p $params.index
    STARlong --runThreadN ${task.cpus} \
        --runMode genomeGenerate \
        --genomeDir ${params.genome} \ 
        --genomeFastaFiles $fasta\
        --sjdbGTFfile $annotation
    """

}


// TODO: implement process for alignment with STARlong
// TODO: check which files are beeing generated and are needed
process align_reads{

    tag "aligning: $name"

    publishDir "$params.outdir/alignment", mode: 'copy'

    input:
    var intron_max from params.intron_max
    var transcript_max from params.transcript_max
    file index from star_index
    file hq_fastq from polish_out

    output:
    file "*.*" into star_out
    file "*.bam*" into bam_files

    """"
    time STARlong --readFilesIn ${hq_fastq} --genomeDir $index \
        --runMode alignReads \
        --runThreadN ${task.cpus}
        --outSAMattributes NH HI NM MD \
        --readNameSeparator space \
        --outFilterMultimapScoreRange 1 \
        --outFilterMismatchNmax 1000 \
        --alignIntronMax $intron_max \
        --alignMatesGapMax $transcript_max \
        --limitBAMsortRAM ${task.memory} \
        --genomeLoad NoSharedMemory \
        --scoreGapNoncan -20 \
        --scoreGapGCAG -4 \
        --scoreGapATAC -8 \
        --scoreDelOpen -1 \
        --scoreDelBase -1 \
        --scoreInsOpen -1 \
        --scoreInsBase -1 \
        --alignEndsType Local \
        --seedSearchStartLmax 50 \
        --seedPerReadNmax 100000 \
        --seedPerWindowNmax 1000 \
        --alignTranscriptsPerReadNmax 100000 \
        --alignTranscriptsPerWindowNmax 10000 \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix $name
    
    samtools index $name
    """
}


process bam_to_bed{
    
    tag "bamToBed: $name"

    publishDir "$params.outdir/bed", mode: 'copy'

    input:
    file '$name.bam' from bam_files

    output:
    file '$name.bed' into bed

    """   
    bedtools bamtobed -bed12 -i name.bam > $name.bed
    """
}