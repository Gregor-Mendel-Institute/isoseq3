#!/usr/bin/env nextflow


params.merge = true
params.align = true

params.ref_fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.intron_max = params.genome ? params.genomes[ params.genome ].intron_max ?: false : false
params.transcript_max = params.genome ? params.genomes[ params.genome ].transcript_max ?: false : false

log.info "IsoSeq3 NF  ~  version 0.1"
log.info "====================================="
log.info "input paths: ${params.input}"
log.info "output paths: ${params.output}"
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
    .fromPath(params.primers)
    .ifEmpty { error "Cannot find primer file: $params.primers" }
    .set {primers_remove; primers_refine}

Channel
    .fromPath(params.ref_fasta)
    .ifEmpty { error "Cannot find primer file: $params.primers" }
    .set {ref_fasta}




process ccs_calling{

        tag "circular consensus sequence calling: $name"
        
        publishDir "$params.output/$name/ccs", mode: 'copy'

        input:
        set name, file(bam) from input_ccs
        val outdir from outdir_ch

        output:
        file "*"
        file "${name}.ccs.bam" into ccs_out
        val name into name_ccs
    
        //TODO make minPasses param as parameter
        """
        ccs ${name}.bam ${name}.ccs.bam --noPolish --minPasses 1
        """
}


process primers_rm{

    tag "primer removal: $name"

    publishDir "$params.output/$name/lima", mode: 'copy'

    input:
    val name from name_ccs
    file ccs_bam from ccs_out
    file primers from primers_remove

    output:
    file "*"
    file "${name}.fl.primer_5p--primer_3p.bam" into primers_removed
    val name into name_primers_rm

    """
    lima $ccs_bam $primers ${name}.fl.bam --isoseq --no-pbi
    """
}


process run_refine{

    tag "refining : $name"
    publishDir "$params.output/$name/refine", mode: 'copy'

    input:
    val name from name_primers_rm
    file p_rm_bam from primers_removed
    file primers from primers_refine
    
    output:
    file "*"
    file "${name}.flnc.bam" into refine_out
    file "${name}.flnc.bam" into refine_merge_out
    val name into name_refine

    //TODO update input & output channels
    """
    isoseq3 refine $p_rm_bam $primers ${name}.flnc.bam --require-polya
    """

}


process merge_samples{

    tag "merging SMRT cells"

    publishDir "$params.input/merged", mode: 'copy'

    input:
    file '*.bam' from refine_merge_out.collect()

    output:
    file "merged.flnc.xml" into merge_out
    val "merged" into name_merge_out

    when:
    params.merge

    """
    dataset create --type TranscriptSet merged.flnc.xml *.bam
    """
}


process cluster_reads{

    tag "clustering : $name"
    publishDir "$params.output/$name/cluster", mode: 'copy'

    input:
    file refined from refine_out.mix(merge_out)
    val name from name_refine.mix(name_merge_out)

    output:
    file "*"
    file "${name}.unpolished.bam" into cluster_out
    val name into name_cluster

    """
    isoseq3 cluster $refined ${name}.unpolished.bam --verbose 
    """
}


process polish_reads{
    
    tag "polishing : $name"

    publishDir "$params.output/$name/polish", mode: 'copy'

    input:
    set name, file(bam) from input_polish
    file cluster_bam from cluster_out

    output:
    file "*"
    file "${name}.polished.hq.fastq.gz" into polish_out
    val name into name_polish
    
    """
    isoseq3 polish $cluster_bam ${name}.bam ${name}.polished.bam
    """

}

process align_reads{

    tag "mapping : $name"

    publishDir "$params.output/$name/minimap2", mode: 'copy'

    input:
    file fasta from ref_fasta
    file sample from polish_out
    val name from name_polish

    output:
    file "*"

    when:
    params.align

    """
    minimap2 $fasta $sample \
        -G $params.intron_max \
        -H \
        -ax splice \
        -C 5 \
        -u f \
        -p 0.9 \
        -t ${task.cpus} > $name.paf
    """



}


// process build_index{

//     tag "STARlong index: $params.star_index"

//     storeDir "$params.index_dir"

//     input:
//     file annotation from annotation_file
//     file fasta from fasta_files
//     val genome_dir from params.star_index

//     output:
//     file "${genome_dir}" into sl_index
//     file 'Log.out' into log

//     """
//     mkdir -p ${genome_dir}
//     STARlong --runThreadN ${task.cpus} \
//         --runMode genomeGenerate \
//         --genomeDir ${genome_dir} \
//         --genomeFastaFiles ${fasta} \
//         --sjdbGTFfile ${annotation}
//     """

// }


// process align_reads{

//     tag "aligning: $name"

//     publishDir "$params.input/alignment", mode: 'copy'

//     input:
//     val name from name_polish
//     val intron_max from params.intron_max
//     val transcript_max from params.transcript_max
//     file index from sl_index
//     file hq_fastq from polish_out


//     output:
//     file "*"
//     file "${name}.Aligned.sortedByCoord.out.{bam,bam.bai}" into bam_files
//     val "${name}.Aligned.sortedByCoord.out" into name_align


//     """
//     time STARlong --readFilesIn ${hq_fastq} --genomeDir $index \
// 	--readFilesCommand zcat \
//         --runMode alignReads \
//         --runThreadN ${task.cpus} \
//         --outSAMattributes NH HI NM MD \
//         --readNameSeparator space \
//         --outFilterMultimapScoreRange 1 \
//         --outFilterMismatchNmax 1000 \
//         --alignIntronMax $intron_max \
//         --alignMatesGapMax $transcript_max \
//         --limitBAMsortRAM 0 \
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
//         --outFileNamePrefix ${name}.

//     samtools index ${name}.Aligned.sortedByCoord.out.bam   
 
//     """
// }


process bam_to_bed{
    
    tag "bamToBed: $name"

    publishDir "$params.input/bed", mode: 'copy'

    input:
    val name from name_align
    set file(bam), file(bam_index) from bam_files
    
    output:
    file "${name}.bed" into bed

    """   
    bedtools bamtobed -bed12 -i ${name}.bam > ${name}.bed
    """
}

// TODO: add processes for staging in and out.
