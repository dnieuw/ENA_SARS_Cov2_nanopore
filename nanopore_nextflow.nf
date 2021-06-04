#!/usr/bin/env nextflow

/*
 * Define some parameters in order to specify the refence genomes
 * barcodes to be checked in the demultiplexing process and the coverage
 * thresholds set for the consensus calling
 */

params.COVERAGE = 30
params.OUTDIR = "results"


/*
 * Trim 30 nucleotides of each end of the reads using cutadapt to ensure that primer derived sequences are not used to generate a consensus sequence
 */  
process cut_adapters {
    cpus 19
    memory '30 GB'
    container 'kfdrc/cutadapt'
    
    input:
    path input_file from params.INPUT
    
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

    cpus 19 /* more is better, parallelizes very well*/
    memory '30 GB'
    container 'alexeyebi/ena-sars-cov2-nanopore'
    
    input:
    path trimmed from trimmed_ch
    path ref from params.REFERENCE
    val run_id from params.RUN_ID
    
    output:
    path "${run_id}.bam" into mapped_ch
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
    path bam from mapped_ch
    path ref from params.REFERENCE
    val run_id from params.RUN_ID
    
    output:
    path "${run_id}.vcf" into vcf_ch

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    bam_to_vcf.py -b ${bam} -r ${ref} --mindepth 30 --minAF 0.1 -c ${task.cpus} -o ${run_id}.vcf
    """
}

process annotate_snps {
    publishDir params.OUTDIR, mode:'copy'
    cpus 1
    memory '30 GB'
    container 'alexeyebi/snpeff'

    input:
    path vcf from vcf_ch
    val run_id from params.RUN_ID

    output:
    path("${run_id}.annot.vcf")

    script:
    """
    zcat ${vcf} | sed "s/^NC_045512.2/NC_045512/" > \
    ${run_id}.newchr.vcf
    java -Xmx4g -jar /data/tools/snpEff/snpEff.jar -v -s ${run_id}.snpEff_summary.html sars.cov.2 \
    ${run_id}.newchr.vcf > ${run_id}.annot.vcf
    """
}
