#!/usr/bin/env nextflow
params.index_dir =  params.genome ? params.genomes[ params.genome ].index_dir ?: false : false
params.star_index = params.genome ? params.genomes[ params.genome ].star_index ?: false : false
params.annotation = params.genome ? params.genomes[ params.genome ].annotation ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.intron_max = params.genome ? params.genomes[ params.genome ].intron_max ?: false : false
params.transcript_max = params.genome ? params.genomes[ params.genome ].transcript_max ?: false : false

log.info "IsoSeq3 NF  ~  version 0.1"
log.info "====================================="
log.info "input paths: ${params.input}"
log.info "genome: ${params.genome}"
log.info "genome annotation: ${params.annotation}"
log.info "genome sequence: ${params.fasta}"
log.info "index location: ${params.index_dir}/${params.star_index}"
log.info "intron max length: ${params.intron_max}"
log.info "transcript max length: ${params.transcript_max}"
log.info "\n"

// input channels
Channel
    .fromFilePairs(params.input + '*.{bam,bam.pbi}') { file -> file.name.replaceAll(/.bam|.pbi$/,'') }
    .ifEmpty { error "Cannot find matching bam and pbi files: $params.input." }
    .into {input_ccs; input_polish}
// see https://github.com/nextflow-io/patterns/blob/926d8bdf1080c05de406499fb3b5a0b1ce716fcb/process-per-file-pairs/main2.nf

Channel
    .fromPath(params.input + '/results')
    .into {outdir_ccs; outdir_lima; outdir_cluster; outdir_polish; outdir_alignment; outdir_bed}

Channel
    .fromPath(params.primers)
    .ifEmpty { error "Cannot find primer file: $params.primers" }
    .set {primers_file}

// Channels for index building
Channel
    .fromPath(params.fasta)
    .ifEmpty { error "Cannot find fasta files: $params.fasta" }
    .collectFile(name: 'merged.fa', newLine: true, sort: true)
    .set {fasta_files}

Channel
    .fromPath(params.annotation)
    .ifEmpty { error "Cannot find annotation file: $params.annotation" }
    .set {annotation_file}


process run_ccs{

        tag "ccs: $name"

        publishDir "$params.outdir/ccs", mode: 'copy'

        input:
        set name, file(bam) from input_ccs

        output:
        file "${name}.ccs.*"
        file 'ccs_report.txt'
        file "${name}.ccs.bam" into ccs_out
        val name into sample_id_ccs
    
        //TODO increase number of passes. Using passes 1 for debugging reasons only.
        """
        time ccs ${name}.bam ${name}.ccs.bam --noPolish --minPasses 1
        """
}


process run_lima{

    tag "lima: $name"

    publishDir "$params.outdir/lima", mode: 'copy'

    input:
    val name from sample_id_ccs
    val outdir from outdir_lima
    file ccs_bam from ccs_out
    file primers from primers_file

    output:
    file "${name}.demux.ccs.*"
    file "${name}.demux.ccs.primer_5p--primer_3p.bam" into lima_out
    val name into sample_id_lima

    
    """
    time lima $ccs_bam $primers ${name}.demux.ccs.bam --isoseq --no-pbi --dump-clips --dump-removed
    """

}


process cluster_reads{

    tag "clustering : $name"
    publishDir "$params.outdir/cluster", mode: 'copy'

    input:
    val name from sample_id_lima
    file lima_demux from lima_out
    
    output:
    file "${name}.unpolished.*"
    file "${name}.unpolished.bam" into cluster_out
    val name into sample_id_cluster


    """
    time isoseq3 cluster $lima_demux ${name}.unpolished.bam --require-polya --verbose 
    """
}


process polish_reads{
    
    tag "polishing : $name"

    publishDir "$params.outdir/polish", mode: 'copy'

    input:
    set name, file(bam) from input_polish
    file cluster_bam from cluster_out
    // file pbi from input_pbi
    // file all_reads_bam from input_polish
 
    output:
    file "${name}.polished.*"
    file "${name}.polished.hq.fastq.gz" into polish_out
    val name into sample_id_polish
    
    """
    time isoseq3 polish $cluster_bam ${name}.bam ${name}.polished.bam
    """

}

// TODO make run conditional
process build_index{

    tag "STARlong index: $params.star_index"

    storeDir "$params.index_dir"

    input:
    file annotation from annotation_file
    file fasta from fasta_files
    val genome_dir from params.star_index

    output:
    file "${genome_dir}" into sl_index
    file 'Log.out' into log

    """
    mkdir -p ${genome_dir}
    STARlong --runThreadN ${task.cpus} \
        --runMode genomeGenerate \
        --genomeDir ${genome_dir} \
        --genomeFastaFiles ${fasta} \
        --sjdbGTFfile ${annotation}
    """

}


// TODO: implement process for alignment with STARlong
// TODO: check which files are beeing generated and are needed
process align_reads{

    tag "aligning: $name"

    publishDir "$params.outdir/alignment", mode: 'copy'

    input:
    val name from sample_id_polish
    val intron_max from params.intron_max
    val transcript_max from params.transcript_max
    file index from sl_index
    file hq_fastq from polish_out


    output:
    file "${name}.*" into star_out
    file "${name}.Aligned.sortedByCoord.out.{bam,bam.bai}" into bam_files
    val "${name}.Aligned.sortedByCoord.out"  into sample_id_align


    """
    time STARlong --readFilesIn ${hq_fastq} --genomeDir $index \
	--readFilesCommand zcat \
        --runMode alignReads \
        --runThreadN ${task.cpus} \
        --outSAMattributes NH HI NM MD \
        --readNameSeparator space \
        --outFilterMultimapScoreRange 1 \
        --outFilterMismatchNmax 1000 \
        --alignIntronMax $intron_max \
        --alignMatesGapMax $transcript_max \
        --limitBAMsortRAM 0 \
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
        --outFileNamePrefix ${name}.

    samtools index ${name}.Aligned.sortedByCoord.out.bam   
 
    """
}


process bam_to_bed{
    
    tag "bamToBed: $name"

    publishDir "$params.outdir/bed", mode: 'copy'

    input:
    val name from sample_id_align
    set file(bam), file(bam_index) from bam_files
    val outdir from outdir_bed
    // file '$name.bam' from bam_files

    output:
    file "${name}.bed" into bed

    """   
    bedtools bamtobed -bed12 -i ${name}.bam > ${name}.bed
    """
}

// TODO: add processes for staging in and out.
