#!/usr/bin/env nextflow

Channel
    .fromPath(params.input)
    .splitCsv(header:true, sep:'\t', quote:'"')
    .map{ row-> tuple(row.sampleId, file(row.read1), file(row.read2)) }
    .into { samples_fastqc_ch; samples_trim_ch; }

reference_ch = Channel.fromPath(params.reference)

def summary = [:]
summary['Pipeline Name']  = 'depth-summary'
summary['Input']          = params.input
summary['Reference']      = params.reference
summary['Output dir']     = params.outdir
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile


summary.each{ k, v -> println "${k}: ${v}" }


/*
 * FastQC
 */
process fastqc_pre_trim {
    tag "$sampleId"
    cpus 4
    conda 'fastqc=0.11.8'
    publishDir "${params.outdir}/${sampleId}/fastqc_pre_trim", mode: 'copy', pattern: "*fastqc*"

    input:
	set sampleId, file(read1), file(read2) from samples_fastqc_ch

    output:
	file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $read1 $read2
    """
}

/*
 * trim with trim_galore, run fastqc on trimmed files.
 */
process trim {
    tag "$sampleId"
    cpus 4
    memory '1 GB'
    conda 'trim-galore=0.5.0'
    publishDir "${params.outdir}/${sampleId}/fastqc_post_trim", mode: 'copy', pattern: "*fastqc*"
    input:
        set val(sampleId), file(read1), file(read2) from samples_trim_ch
    output:
        set val("$sampleId"), file("*_val_1.fq"), file("*_val_2.fq") into samples_trimmed_ch
        file "*fastqc*"
    script:
    """
    trim_galore --fastqc --paired $read1 $read2
    """
}

/*
 * bwa mem & samtools sort
 */
process bwa_align {
    cpus 8
    memory '2 GB'
    tag "$sampleId"
    conda 'bwa=0.7.17 samtools=1.9'
    publishDir "${params.outdir}/${sampleId}/alignment", mode: 'copy'

    input:
	set val(sampleId), file(read1), file(read2) from samples_trimmed_ch
	file reference from reference_ch
    output:
        set val("$sampleId"), file("*.sort.bam"), file("*.sort.bam.bai") into bwa_align_results
    script:
    """
    bwa index $reference
    bwa mem -t 8 $reference $read1 $read2 | samtools sort -@8 -o '$sampleId'.sort.bam -
    samtools index '$sampleId'.sort.bam
    """
}



