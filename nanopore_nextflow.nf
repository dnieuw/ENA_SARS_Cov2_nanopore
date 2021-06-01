#!/usr/bin/env nextflow

/*
 * Define some parameters in order to specify the refence genomes
 * barcodes to be checked in the demultiplexing process and the coverage
 * thresholds set for the consensus calling
 */

params.COVERAGE = 30
params.PIPELINE_STARTED = "pipeline_started"
params.PIPELINE_FINISHED = "pipeline_finished"
params.PIPELINE_FIELD = "pipeline_analysis"
params.EXPORT_STARTED = "export_started"
params.EXPORT_FINISHED = "export_finished"
params.EXPORT_FIELD = "export_to_ena"
params.OUTDIR = "results"

/*
 * Report to mongodb that pipeline started
 */
 process report_pipeline_started {
     cpus 1
     memory '1 GB'
     container 'alexeyebi/ena-sars-cov2-nanopore'

     input:
     val run_id from params.RUN
     val status from params.PIPELINE_STARTED
     val field from params.PIPELINE_FIELD

     output:
     path 'pipeline_started.log' into pipeline_started_ch

     script:
     """
     touch pipeline_started.log
     """
 }

/*
 * Trim 30 nucleotides of each end of the reads using cutadapt to ensure that primer derived sequences are not used to generate a consensus sequence
 */  
process cut_adapters {
    cpus 10 /* more is better, parallelizes quiet well*/
    memory '20 GB'
    container 'kfdrc/cutadapt'
    
    input:
    path input_file from params.INPUT
    path pipeline_started from pipeline_started_ch
    
    output:
    path 'trimmed.fastq' into trimmed_ch
    
    script:
    """
    cutadapt -u 30 -u -30 -o trimmed.fastq ${input_file} -m 75 -j ${task.cpus} --quiet
    """
}

/*
 * Map reads to SARS-CoV-2 reference genome
 */ 
process map_to_reference {
    publishDir params.OUTDIR, mode:'copy'

    cpus 10 /* more is better, parallelizes very well*/
    memory '20 GB'
    container 'alexeyebi/ena-sars-cov2-nanopore'
    
    input:
    path trimmed from trimmed_ch
    path ref from params.REFERENCE
    val run_id from params.RUN
    
    output:
    path "${run_id}.bam" into mapped_ch1, mapped_ch2, mapped_ch3
    path("${run_id}.bam")
    
    script:
    """
    minimap2 -Y -t ${task.cpus} -x map-ont -a ${ref} ${trimmed} | samtools view -bF 4 - | samtools sort -@ ${task.cpus} - > ${run_id}.bam
    """
}


process bam_to_vcf {
    tag '$run_id'
    publishDir params.OUTDIR, mode:'copy'
    cpus 10
    container 'rmwthorne/ena-sars-cov2-nanopore'

    input:
    path bam from mapped_ch1
    path ref from params.REFERENCE
    val run_id from params.RUN
    
    output:
    path "${run_id}.vcf" into vcf_ch

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    bam_to_vcf.py -b ${bam} -r ${ref} --mindepth 30 --minAF 0.1 -c ${task.cpus} -o ${run_id}.vcf
    """
}

