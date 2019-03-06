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
    .ifEmpty { error "Cannot find matching bam and pbi files: $params.input. Make 
     sure your bam files are pb indexed." }
    .set { ccs_in }
// see https://github.com/nextflow-io/patterns/blob/926d8bdf1080c05de406499fb3b5a0b1ce716fcb/process-per-file-pairs/main2.nf

Channel
    .fromPath(params.input + '*.bam')
    .ifEmpty { error "Cannot find matching bam files: $params.input." }
    .tap { merge_subread_in }
    .map{ file -> tuple(file.name.replaceAll(/.bam$/,''), file) }
    .tap { polish_in }

Channel
    .fromPath(params.primers)
    .ifEmpty { error "Cannot find primer file: $params.primers" }
    .into { primers_remove; primers_refine }

Channel
    .fromPath(params.ref_fasta)
    .ifEmpty { error "Cannot find primer file: $params.primers" }
    .set {ref_fasta}


process ccs_calling{

        tag "circular consensus sequence calling: $name"
        
        publishDir "$params.output/$name/ccs", mode: 'copy'

        input:
        set name, file(bam) from ccs_in.//dump(tag: 'input')

        output:
        file "*"
        set val(name), file("${name}.ccs.bam") into ccs_out
     
        //TODO make minPasses param as parameter
        """
        ccs ${name}.bam ${name}.ccs.bam --noPolish --minPasses 1
        """
}


process remove_primers{

    tag "primer removal: $name"

    publishDir "$params.output/$name/lima", mode: 'copy'

    input:
    set name, file(bam) from ccs_out//.dump(tag: 'ccs_name')
    file primers from primers_remove.collect()
    
    output:
    file "*"
    set val(name), file("${name}.fl.primer_5p--primer_3p.bam") into primers_removed_out
 
    """
    lima $bam $primers ${name}.fl.bam --isoseq --no-pbi
    """
}


process run_refine{

    tag "refining : $name"
    publishDir "$params.output/$name/refine", mode: 'copy'

    input:
    set name, file(bam) from primers_removed_out//.dump(tag: 'primers_removed')
    file primers from primers_refine.collect()
    
    output:
    file "*"
    file("${name}.flnc.bam") into refine_merge_out
    set val(name), file("${name}.flnc.bam") into refine_out
 
    //TODO update input & output channels
    """
    isoseq3 refine $bam $primers ${name}.flnc.bam --require-polya
    """

}


process merge_transcripts{

    tag "merging transcript sets ${bam}"

    publishDir "$params.output/merged", mode: 'copy'

    input:
    file(bam) from refine_merge_out.collect().dump(tag: 'merge transcripts')

    output:
    set val("merged"), file("merged.flnc.xml") into cluster_in

    when:
    params.merge

    """
    dataset create --type TranscriptSet merged.flnc.xml ${bam}

    """
}


process merge_subreads{

    tag "merging subreads ${bam}"

    publishDir "$params.output/merged", mode: 'copy'

    input:
    file(bam), from merge_subread_in.collect().dump(tag: 'merge subreads')

    output:
    set val("merged"), file("merged.subreadset.xml") into merged_subreads

    when:
    params.merge

    """
    dataset create --type SubreadSet merged.subreadset.xml ${bam}
    """
}


process cluster_reads{

    tag "clustering : $name"
    publishDir "$params.output/$name/cluster", mode: 'copy'

    input:
    set name, file(refined) from refine_out.concat(cluster_in).dump(tag: 'cluster')

    output:
    file "*"
    set val(name), file( "${name}.unpolished.bam") into cluster_out

    """
    isoseq3 cluster ${refined} ${name}.unpolished.bam --verbose 
    """
}


polish_in.concat(merge_subreads)
polish_in.join(cluster_out)

process polish_reads{
    
    tag "polishing : $name"

    publishDir "$params.output/$name/polish", mode: 'copy'

    input:
    set name, file(unpolished_bam), file(subreads_bam) from polish_in.dump(tag: 'polish')

    output:
    file "*"
    set name, file("${name}.polished.hq.fastq.gz") into polish_out
    
    """
    isoseq3 polish ${unpolished_bam} ${subreads_bam} ${name}.polished.bam
    """

}

process align_reads{

    tag "mapping : $name"

    publishDir "$params.output/$name/minimap2", mode: 'copy'

    input:
    set name, file(sample) from polish_out.dump(tag: 'align')
    file fasta from ref_fasta.collect()

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