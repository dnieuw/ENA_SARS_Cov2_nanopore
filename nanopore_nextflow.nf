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
     update_samples_status.py ${run_id} ${status} ${field}
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
    memory '10 GB'
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

/*
 * Generate consensus genome from mapped reads
 */ 
process create_consensus {
    cpus 1
    memory '1 GB'
    container 'alexeyebi/ena-sars-cov2-nanopore'
    
    input:
    path bam from mapped_ch1
    val cov from params.COVERAGE
    
    output:
    path 'consensus.fasta' into consensus_ch
    
    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    bam2consensus.py -i ${bam} -o consensus.fasta -d ${cov} -g 1
    """
}

/*
 * Align consensus with reference to keep coordinates the same
 */ 
process align_consensus {
    publishDir params.OUTDIR, mode:'copy'

    cpus 1 /* doesn't benefit from more cores*/
    memory '10 GB'
    container 'alexeyebi/ena-sars-cov2-nanopore'
    
    input:
    path consensus from consensus_ch
    path ref from params.REFERENCE
    val name from params.NAME
    val status from params.PIPELINE_FINISHED
    val run_id from params.RUN
    val field from params.PIPELINE_FIELD

    output:
    path "${run_id}.fasta.gz" into align_consensus_ch
    path("${run_id}.fasta")
    
    script:
    """
    align_to_ref.py -i ${consensus} -o ${run_id}.fasta -r ${ref} -n ${name}
    gzip ${run_id}.fasta
    update_samples_status.py ${run_id} ${status} ${field}
    """
}


/*
 * Upload BAM file to ENA ftp
 */
process upload_files_to_ena {
    cpus 1
    memory '1 GB'
    container 'alexeyebi/ena-sars-cov2-nanopore'

    input:
    path consensus from align_consensus_ch
    path bam from mapped_ch2
    val run_id from params.RUN
    val status from params.EXPORT_STARTED
    val user from params.USER
    val password from params.PASSWORD
    val field from params.EXPORT_FIELD

    output:
    path 'files_uploaded_to_ena.log' into files_uploaded_to_ena_ch
    path "${consensus}" into gzip_consensus_ch

    script:
    """
    update_samples_status.py ${run_id} ${status} ${field}
    curl -T ${bam} ftp://webin.ebi.ac.uk --user ${user}:${password}
    curl -T ${consensus} ftp://webin.ebi.ac.uk --user ${user}:${password}
    touch files_uploaded_to_ena.log
    """
}


/*
 * Create analysis xml file, required for ENA submission
 */
process create_analysis_xml {
    cpus 1
    memory '1 GB'
    container 'alexeyebi/ena-sars-cov2-nanopore'

    input:
    path files_uploaded_to_ena from files_uploaded_to_ena_ch
    path bam from mapped_ch3
    path consensus from gzip_consensus_ch
    val run_id from params.RUN
    val user from params.USER
    val password from params.PASSWORD

    output:
    path "${run_id}_alignment_analysis.xml" into create_alignment_analysis_xml_ch1, create_alignment_analysis_xml_ch2
    path "${run_id}_consensus_sequence_analysis.xml" into create_consensus_sequence_analysis_xml_ch1, create_consensus_sequence_analysis_xml_ch2

    script:
    """
    create_analysis_xml.py ${bam} ${user} ${password}
    create_analysis_xml.py ${consensus} ${user} ${password}
    """
}

/*
 * Create submission xml file, required for ENA submission
 */
 process create_submission_xml {
     cpus 1
     memory '1 GB'
     container 'alexeyebi/ena-sars-cov2-nanopore'

     input:
     val run_id from params.RUN
     path alignment_analysis_xml from create_alignment_analysis_xml_ch1
     path consensus_sequence_analysis_xml from create_consensus_sequence_analysis_xml_ch1

     output:
     path "${run_id}_alignment_submission.xml" into create_alignment_submission_xml_ch
     path "${run_id}_consensus_sequence_submission.xml" into create_consensus_sequence_submission_xml_ch

     script:
     """
     create_submission_xml.py ${alignment_analysis_xml}
     create_submission_xml.py ${consensus_sequence_analysis_xml}
     """
 }

 /*
  * Submit xml files to ena
  */
  process submit_xml_files_to_ena {
      cpus 1
      memory '1 GB'
      container 'alexeyebi/ena-sars-cov2-nanopore'

      input:
      val run_id from params.RUN
      val status from params.EXPORT_FINISHED
      path alignment_analysis_xml from create_alignment_analysis_xml_ch2
      path alignment_submission_xml from create_alignment_submission_xml_ch
      path consensus_sequence_analysis_xml from create_consensus_sequence_analysis_xml_ch2
      path consensus_sequence_submission_xml from create_consensus_sequence_submission_xml_ch
      val user from params.USER
      val password from params.PASSWORD
      val field from params.EXPORT_FIELD


      script:
      """
      submit_xml_files_to_ena.py ${alignment_submission_xml} ${alignment_analysis_xml} ${user} ${password}
      submit_xml_files_to_ena.py ${consensus_sequence_submission_xml} ${consensus_sequence_analysis_xml} ${user} ${password}
      update_samples_status.py ${run_id} ${status} ${field}
      """
  }


