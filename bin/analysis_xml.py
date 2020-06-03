#!/usr/bin/env python

"""
This class will generate Analysis XML that could be send to ENA
"""

from lxml import etree


class AnalysisXML:
    def __init__(self, alias, title, description, sample_accession,
                 study_accession, run_accession, checksum, analysis_type,
                 analysis_xml_file):
        self.alias = alias
        self.title = title
        self.description = description
        self.sample_accession = sample_accession
        self.study_accession = study_accession
        self.run_accession = run_accession
        self.checksum = checksum
        self.analysis_xml_file = analysis_xml_file
        self.analysis_type = analysis_type

    def build_xml(self):
        analysis_set = etree.Element('ANALYSIS_SET')
        analysis_xml = etree.ElementTree(analysis_set)
        analysis_elt = etree.SubElement(analysis_set, 'ANALYSIS',
                                        alias=self.alias,)
        title = etree.SubElement(analysis_elt, 'TITLE')
        title.text = self.title
        description = etree.SubElement(analysis_elt, 'DESCRIPTION')
        description.text = self.description
        etree.SubElement(analysis_elt, 'STUDY_REF',
                         accession=self.study_accession)
        etree.SubElement(analysis_elt, 'SAMPLE_REF',
                         accession=self.sample_accession)
        etree.SubElement(analysis_elt, 'RUN_REF',
                         accession=self.run_accession)
        analysis_type_elt = etree.SubElement(analysis_elt, 'ANALYSIS_TYPE')
        if self.analysis_type == 'REFERENCE_ALIGNMENT':
            reference_alignment = etree.SubElement(analysis_type_elt,
                                                   'REFERENCE_ALIGNMENT')
            assembly = etree.SubElement(reference_alignment, 'ASSEMBLY')
            etree.SubElement(assembly, 'STANDARD', accession='GCA_009858895.3')
            etree.SubElement(reference_alignment, 'SEQUENCE',
                             accession='NC_045512.2')
            files = etree.SubElement(analysis_elt, 'FILES')
            etree.SubElement(files, 'FILE',
                             filename=f"{self.run_accession}.bam",
                             filetype="bam", checksum_method="MD5",
                             checksum=self.checksum)
        elif self.analysis_type == 'PATHOGEN_ANALYSIS':
            etree.SubElement(analysis_type_elt, 'PATHOGEN_ANALYSIS')
            files = etree.SubElement(analysis_elt, 'FILES')
            etree.SubElement(files, 'FILE',
                             filename=f"{self.run_accession}.fasta.gz",
                             filetype="other", checksum_method="MD5",
                             checksum=self.checksum)

        # Writing xml file
        analysis_xml.write(self.analysis_xml_file, pretty_print=True,
                           xml_declaration=True, encoding='UTF-8')
