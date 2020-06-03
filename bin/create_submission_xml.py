#!/usr/bin/env python

import sys

from submission_xml import SubmissionXML


def main():
    if 'alignment' in FILENAME:
        run_id = FILENAME.split('_alignment_analysis.xml')[0]
        alias = f"sub_{run_id}_alignment"
        submission_xml_file = f"{run_id}_alignment_submission.xml"
    elif 'consensus_sequence' in FILENAME:
        run_id = FILENAME.split('_consensus_sequence_analysis.xml')[0]
        alias = f"sub_{run_id}_consensus_sequence"
        submission_xml_file = f"{run_id}_consensus_sequence_submission.xml"
    else:
        print(f"File format is not supported: {FILENAME}")
        sys.exit(1)
    action = 'ADD'
    schema = 'analysis'
    submission_xml_obj = SubmissionXML(
        alias=alias, action=action,
        submission_xml_file=submission_xml_file, source_xml=FILENAME,
        schema=schema)
    submission_xml_obj.build_xml()


if __name__ == "__main__":
    FILENAME = sys.argv[1]
    main()
