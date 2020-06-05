#!/usr/bin/env python

import sys
import subprocess

from lxml import etree


def main():
    command = f'curl -k  -F "SUBMISSION=@{SUBMISSION_XML}" ' \
              f'-F "ANALYSIS=@{ANALYSIS_XML}" ' \
              f'"{ANALYSIS_SUBMISSION_URL_DEV}%20{USER}%20{PASSWORD}"'
    completed_process_command = subprocess.run(command, shell=True,
                                               capture_output=True)
    returned_xml = completed_process_command.stdout
    submission_error_messages = list()
    if returned_xml:
        root = etree.XML(returned_xml)
        root = etree.fromstring(returned_xml)
        for messages in root.findall('MESSAGES'):
            for mess in messages.findall('ERROR'):
                submission_error_messages.append('ERROR:' + mess.text)
    if len(submission_error_messages) > 0:
        for message in submission_error_messages:
            print(message)
        sys.exit(1)


if __name__ == "__main__":
    SUBMISSION_XML = sys.argv[1]
    ANALYSIS_XML = sys.argv[2]
    ANALYSIS_SUBMISSION_URL_DEV = "https://www-test.ebi.ac.uk/ena/submit/" \
                                  "drop-box/submit/?auth=ENA"
    USER = sys.argv[3]
    PASSWORD = sys.argv[4]
    main()
