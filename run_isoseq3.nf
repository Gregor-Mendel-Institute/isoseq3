#!/usr/bin/env nextflow

params.merge = true
params.align = true
params.bookend = false


//params.ref_fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
//params.intron_max = params.genome ? params.genomes[ params.genome ].intron_max ?: false : false
//params.transcript_max = params.genome ? params.genomes[ params.genome ].transcript_max ?: false : false

log.info "IsoSeq3 NF  ~  version 3.1"
log.info "====================================="
log.info "input paths: ${params.input}"
log.info "output paths: ${params.output}"
log.info "merge smrt cells: ${params.merge}"
log.info "align reads: ${params.align}"
log.info "genome: ${params.genome}"
log.info "genome sequence: ${params.ref_fasta}"
log.info "intron max length: ${params.intron_max}"
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
    .into {primers_remove; primers_refine}

Channel
    .fromPath(params.ref_fasta)
    .ifEmpty { error "Cannot find primer file: $params.primers" }
    .set {ref_fasta}

process ccs_calling{

        tag "circular consensus sequence calling: $name"
        
        publishDir "$params.output/$name/ccs", mode: 'copy'

        input:
        set name, file(bam) from input_ccs

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

    publishDir "$params.output/merged", mode: 'copy'

    input:
    file bam from refine_merge_out.collect()

    output:
    file "merged.flnc.xml" into merge_out
    val "merged" into name_merge_out

    when:
    params.merge

    """
    dataset create --type TranscriptSet merged.flnc.xml $bam
    """
}


process cluster_reads{

    tag "clustering : $name"
    publishDir "$params.output/$name/cluster", mode: 'copy'

    input:
    file refined from refine_out.concat(merge_out)
    val name from name_refine.concat(name_merge_out)

    output:
    file "*"
    file "${name}.unpolished.bam" into cluster_out
    val name into name_cluster

    """
    isoseq3 cluster ${refined} ${name}.unpolished.bam --verbose 
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
    isoseq3 polish ${cluster_bam} ${name}.bam ${name}.polished.bam
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
    file "*.{bam,bed,log}"

    when:
    params.align

    """
    minimap2 $fasta ${sample} \
        -G $params.intron_max \
        -H \
        -ax splice \
        -C 5 \
        -u f \
        -p 0.9 \
        -t ${task.cpus} > ${name}.aln.sam \
        2> ${name}.log

    samtools view -Sb ${name}.aln.sam > ${name}.aln.bam

    bedtools bamtobed -bed12 -i ${name}.aln.bam > ${name}.aln.bed
    """
}

// process bam_to_bed{
    
//     tag "bamToBed: $name"

//     publishDir "$params.input/bed", mode: 'copy'

//     input:
//     val name from name_align
//     set file(bam), file(bam_index) from bam_files
    
//     output:
//     file "${name}.bed" into bed

//     """   
//     bedtools bamtobed -bed12 -i ${name}.bam > ${name}.bed
//     """
// }

// TODO: add processes for staging in and out.
