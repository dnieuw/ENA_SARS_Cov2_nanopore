#!/usr/bin/env nextflow

/*
 * Define some parameters in order to specify the refence genomes
 * barcodes to be checked in the demultiplexing process and the coverage
 * thresholds set for the consensus calling
 */
 
params.REFERENCE = '/home/dnieuwenhuijse/Projects/SARS-CoV-2/nanopore_nextflow/reference.fasta'
params.COVERAGE = 30
params.NAME = "consensus_test"

/*
 * Trim 30 nucleotides of each end of the reads using cutadapt to ensure that primer derived sequences are not used to generate a consensus sequence
 */  
process cut_adapters {
    
    cpus 1 /* more is better, parallelizes quiet well*/
    memory '1 GB' 
    /*container '###'*/
    
    input:
    path input_file from '/home/dnieuwenhuijse/Projects/SARS-CoV-2/nanopore_nextflow/input.fastq'
    
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

    cpus 1 /* more is better, parallelizes very well*/
    memory '1 GB'
    /*container '###'*/
    
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
    /*container '###'*/
    
    input:
    path mapped
    val cov from params.COVERAGE
    
    output:
    path 'consensus.fasta' into consensus
    
    script:
    """
    samtools index -@ ${task.cpus} ${mapped}
    bam2consensus.py -i ${mapped} -o consensus.fasta -d ${cov}
    """
}

/*
 * Align consensus with reference to keep coordinates the same
 */ 
process align_consensus {

    cpus 1 /* doesn't benefit from more cores*/
    memory '1 GB'
    /*container '###'*/
    
    input:
    path consensus
    path ref from params.REFERENCE
    val name from params.NAME
    
    script:
    """
    align_to_ref.py -i ${consensus} -o consensus.fasta -r ${ref} -n 'consensus'
    """
}