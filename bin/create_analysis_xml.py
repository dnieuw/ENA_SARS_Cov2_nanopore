#!/usr/bin/env python

import sys
import os
import subprocess

from pymongo import MongoClient

from analysis_xml import AnalysisXML


def main():
    sample = DB.samples.find_one({'id': FILENAME.split(".")[0]})
    md5 = check_md5_values()
    create_analysis_xml(sample, md5)


def check_md5_values():
    md5_uploaded_command = f"curl -s ftp://webin.ebi.ac.uk/" \
                           f"{os.path.basename(FILENAME)} " \
                           f"--user {USER}:{PASSWORD} | md5sum | cut -f1 -d ' '"
    md5_original_command = f"md5sum {FILENAME} | cut -f1 -d ' '"
    completed_process_original = subprocess.run(md5_original_command,
                                                shell=True,
                                                capture_output=True,
                                                timeout=300)
    completed_process_uploaded = subprocess.run(md5_uploaded_command,
                                                shell=True,
                                                capture_output=True,
                                                timeout=300)
    md5_original = completed_process_original.stdout.decode('utf-8').rstrip()
    md5_uploaded = completed_process_uploaded.stdout.decode('utf-8').rstrip()
    if md5_original != md5_uploaded:
        sys.exit(1)
    else:
        return md5_original


def create_analysis_xml(sample, md5):
    if 'bam' in FILENAME:
        run_id = FILENAME.split('.bam')[0]
        alias = f"{run_id}_alignment"
        title = f"{run_id} alignment"
        description = f"{run_id} alignment"
        analysis_type = 'REFERENCE_ALIGNMENT'
        analysis_xml_file = f"{run_id}_alignment_analysis.xml"
    elif 'fasta.gz' in FILENAME:
        run_id = FILENAME.split('.fasta.gz')[0]
        alias = f"{run_id}_consensus_sequence"
        title = f"{run_id} consensus sequence"
        description = f"{run_id} consensus sequence"
        analysis_type = 'PATHOGEN_ANALYSIS'
        analysis_xml_file = f"{run_id}_consensus_sequence_analysis.xml"
    else:
        print(f"File format is not supported: {FILENAME}")
        sys.exit(1)
    sample_accession = sample['sample_id']
    study_accession = sample['study_id']
    analysis_xml_obj = AnalysisXML(
        alias=alias, title=title, description=description,
        sample_accession=sample_accession, study_accession=study_accession,
        run_accession=run_id, checksum=md5, analysis_type=analysis_type,
        analysis_xml_file=analysis_xml_file)
    analysis_xml_obj.build_xml()


if __name__ == "__main__":
    FILENAME = sys.argv[1]
    USER = sys.argv[2]
    PASSWORD = sys.argv[3]

    # Getting access to MongoDB
    CLIENT = MongoClient('mongodb://samples-logs-db-svc')
    DB = CLIENT.samples

    main()
