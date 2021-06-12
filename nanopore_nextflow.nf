#!/usr/bin/env nextflow

/*
 * Define some parameters in order to specify the refence genomes
 * barcodes to be checked in the demultiplexing process and the coverage
 * thresholds set for the consensus calling
 */

params.OUTDIR = "results"
params.SARS2_FA = "/data/ref/sars2/fa/NC_045512.2.fa"

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
    path ref from params.SARS2_FA
    val run_id from params.RUN_ID
    
    output:
    path "${run_id}.bam" into sars2_aligned_reads_ch, sars2_aligned_reads_ch2
    path("${run_id}.bam")
    
    script:
    """
    minimap2 -Y -t ${task.cpus} -x map-ont -a ${ref} ${trimmed} | samtools view -bF 4 - | samtools sort -@ ${task.cpus} - > ${run_id}.bam
    """
}

process check_coverage {
    publishDir params.OUTDIR, mode:'copy'
    cpus 1
    memory '30 GB'
    container 'alexeyebi/bowtie2_samtools'

    input:
    path bam from sars2_aligned_reads_ch
    val run_id from params.RUN_ID
    path sars2_fasta from params.SARS2_FA

    output:
    path "${run_id}.pileup" into check_coverage_ch
    path("${run_id}.pileup")

    script:
    """
    samtools mpileup -a -A -Q 0 -d 8000 -f ${sars2_fasta} ${bam} > \
    ${run_id}.pileup
    """
}

process make_small_file_with_coverage {
    publishDir params.OUTDIR, mode:'copy'
    cpus 1
    memory '30 GB'
    container 'alexeyebi/bowtie2_samtools'

    input:
    path pileup from check_coverage_ch
    val run_id from params.RUN_ID

    output:
    path("${run_id}.coverage")

    script:
    """
    cat ${pileup} | awk '{print \$2,","\$3,","\$4}' > ${run_id}.coverage
    """
}


process bam_to_vcf {
    tag '$run_id'
    publishDir params.OUTDIR, mode:'copy'
    cpus 10
    memory '30 GB'
    container 'alexeyebi/ena-sars-cov2-nanopore'

    input:
    path bam from sars2_aligned_reads_ch2
    path ref from params.SARS2_FA
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
    cat ${vcf} | sed "s/^NC_045512.2/NC_045512/" > \
    ${run_id}.newchr.vcf
    java -Xmx4g -jar /data/tools/snpEff/snpEff.jar -q -no-downstream -no-upstream -noStats sars.cov.2 ${run_id}.newchr.vcf > ${run_id}.annot.vcf
    """
}
