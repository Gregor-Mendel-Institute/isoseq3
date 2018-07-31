#!/usr/bin/env nextflow

params.input = '/lustre/scratch/users/falko.hofmann/isoseq/test/*.bam'
params.outdir = '/lustre/scratch/users/falko.hofmann/isoseq/test/results'
params.primers = '/lustre/scratch/users/falko.hofmann/pipelines/isoseq3/primers.fasta'
params.annotation = 'tair10'
params.intron_max  = 6000
params.transcript_max = 60000

Channel.fromPath(params.input).into {input_ccs; input_polish}
primers_file = file(params.primers)


process run_ccs{

        publishDir "$params.outdir/ccs", mode: 'copy'

        input:
        file ccs_input from input_ccs

        output:
        file 'ccs.*'
	file 'ccs_report.txt'
        file "ccs.bam" into ccs_out

        """
        time ccs $ccs_input ccs.bam --noPolish --minPasses 1
        """
}



process run_lima{

    publishDir "$params.outdir/lima", mode: 'copy'

    input:
    file ccs_bam from ccs_out
    file primers from primers_file

    output:
    // file 'demux.css.json'
    // file 'demux.css.lima.clips'
    // file 'demux.css.lima.counts'
    // file 'demux.css.lima.report'
    // file 'demux.css.lima.summary'
    // file 'demux.css.removed.bam'
    // file 'demux.css.removed.bam.pbi'
    // file 'demux.css.removed.subreadset.xml'
    // file 'demux.css.primer_5p--primer_3p.subreadset.xml'
    file 'demux.ccs.*'
    file 'demux.ccs.primer_5p--primer_3p.bam' into lima_out
    
    """
    time lima $ccs_bam $primers demux.ccs.bam --isoseq --no-pbi --dump-clips --dump-removed
    """

}


process cluster_reads{

    publishDir "$params.outdir/cluster", mode: 'copy'

    input:
    file lima_demux from lima_out
    
    output:
    file 'unpolished.*'
    file 'unpolished.bam' into cluster_out

    """
    time isoseq3 cluster $lima_demux unpolished.bam --require-polya --verbose 
    """
}


process polish_reads{

    publishDir "$params.outdir/polish", mode: 'copy'

    input:
    file cluster_bam from cluster_out
    file all_reads_bam from input_polish

    output:
    file 'polished.*'
    file 'polished.hq.fastq.gz' into polish_out

    """
    time isoseq3 polish $cluster_bam $all_reads_bam polished.bam
    """

}

// process align_reads{
//     tag "$name"

//     publishDir "$params.outdir/star", mode: 'copy'

//     input:
//     var intron_max from
//     var transcript_max 

//     file polished_bam from polish_out
//     file genome from genome_file
//     file index from starlong_index


//     output:
//     file 

//     """"
//     time STARlong --readFilesIn $polished_bam --genomeDir $genome \
//         --runMode alignReads \
//         --outSAMattributes NH HI NM MD \
//         --readNameSeparator space \
//         --outFilterMultimapScoreRange 1 \
//         --outFilterMismatchNmax 1000 \
//         --alignIntronMax 6000 \
//         --alignMatesGapMax 6000 \
//         --limitBAMsortRAM 10000000000\
//         --genomeLoad NoSharedMemory \
//         --scoreGapNoncan -20 \
//         --scoreGapGCAG -4 \
//         --scoreGapATAC -8 \
//         --scoreDelOpen -1 \
//         --scoreDelBase -1 \
//         --scoreInsOpen -1 \
//         --scoreInsBase -1 \
//         --alignEndsType Local \
//         --seedSearchStartLmax 50 \
//         --seedPerReadNmax 100000 \
//         --seedPerWindowNmax 1000 \
//         --alignTranscriptsPerReadNmax 100000 \
//         --alignTranscriptsPerWindowNmax 10000 \
//         --outSAMtype BAM SortedByCoordinate \
//         --outFileNamePrefix BAM/$sample.
//     """
