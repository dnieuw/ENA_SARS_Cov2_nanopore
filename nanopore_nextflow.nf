#!/usr/bin/env nextflow

/*
 * Define some parameters in order to specify the refence genomes
 * barcodes to be checked in the demultiplexing process and the coverage
 * thresholds set for the consensus calling
 */

params.COVERAGE = 30
params.STARTED = "pipeline started"
params.FINISHED = "pipeline finished"

/*
 * Report to mongodb that pipeline started
 */
 process report_pipeline_started {
     cpus 1
     memory '1 GB'
     container 'alexeyebi/ena-sars-cov2-nanopore'

     input:
     val name from params.NAME
     val status from params.STARTED

     output:
     path 'pipeline_started.log' into pipeline_started_ch

     script:
     """
     update_samples_status.py ${name} ${status}
     touch pipeline_started.log
     """
 }

/*
 * Trim 30 nucleotides of each end of the reads using cutadapt to ensure that primer derived sequences are not used to generate a consensus sequence
 */  
process cut_adapters {
    
    cpus 19 /* more is better, parallelizes quiet well*/
    memory '1 GB'
    container 'kfdrc/cutadapt'
    
    input:
    path input_file from params.INPUT
    path pipeline_started from pipeline_started_ch
    
    output:
    path 'trimmed.fastq' into trimmed
    
    script:
    """
    cutadapt -u 30 -u -30 -o trimmed.fastq ${input_file} -m 75 -j ${task.cpus} --quiet
    """
}

/*
 * Map reads to SARS-CoV-2 reference genome
 */ 
process map_to_reference {

    cpus 19 /* more is better, parallelizes very well*/
    memory '1 GB'
    container 'alexeyebi/ena-sars-cov2-nanopore'
    
    input:
    path trimmed
    path ref from params.REFERENCE
    
    output:
    path 'mapped.fastq' into mapped
    
    script:
    """
    minimap2 -Y -t ${task.cpus} -x map-ont -a ${ref} ${trimmed} | samtools view -bF 4 - | samtools sort -@ ${task.cpus} - > mapped.fastq
    """
}

/*
 * Generate consensus genome from mapped reads
 */ 
process create_consensus {

    cpus 1
    memory '1 GB'
    container 'alexeyebi/ena-sars-cov2-nanopore'
    
    input:
    path mapped
    val cov from params.COVERAGE
    
    output:
    path 'consensus.fasta' into consensus
    
    script:
    """
    samtools index -@ ${task.cpus} ${mapped}
    bam2consensus.py -i ${mapped} -o consensus.fasta -d ${cov} -g 1
    """
}

/*
 * Align consensus with reference to keep coordinates the same
 */ 
process align_consensus {

    publishDir params.outdir, mode:'copy'
    cpus 1 /* doesn't benefit from more cores*/
    memory '10 GB'
    container 'alexeyebi/ena-sars-cov2-nanopore'
    
    input:
    path consensus
    path ref from params.REFERENCE
    val name from params.NAME
    val status from params.FINISHED

    output:
    path('results.fasta')
    
    script:
    """
    align_to_ref.py -i ${consensus} -o results.fasta -r ${ref} -n ${name}
    update_samples_status.py ${name} ${status}
    """
}
