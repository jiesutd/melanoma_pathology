# -*- coding: utf-8 -*-
# @Author: Jie Yang
# @Date:   2019-08-20 14:52:25
# @Last Modified by:   Jie Yang,     Contact: jieynlp@gmail.com
# @Last Modified time: 2020-05-07 02:17:18

import datetime 
import os
from os.path import isfile, join
import matplotlib.pyplot as plt
import datetime
import csv
import re
from datetime import date, timedelta
from dateutil.relativedelta import relativedelta
import pandas as pd
import numpy as np

import scipy.stats


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h

the_count = 0
general_count = 0

note_pair_dict = {}
sample_dict = []

regression_dict = {
	"PRESENT (SEE COMMENT)":"PRESENT",
	"PRESENT (EXTENSIVE)":"PRESENT",
	"PRESENT, EXTENSIVE":"PRESENT",
	"PRESENT, PROMINENT":"PRESENT",
	"PRESENT":"PRESENT",
	"PRESENT (MARKED)":"PRESENT",
	"INFLAMMATORY REGRESSION PRESENT":"PRESENT",
	"PRESENT (EXTENSIVE); EXTENDING TO A SIDE MARGIN":"PRESENT",
	"PRESENT, INFLAMMATORY":"PRESENT",
	"INFLAMMATORY REGRESSION":"PRESENT",
	"PRESENT (INFLAMMATORY)":"PRESENT",
	"PRESENT (EXTENSIVE AND EXPANSILE)":"PRESENT",
	"PRESENT, EXTENSIVE":"PRESENT",
	"PRESENT (PER REPORT)":"PRESENT",
	"POSSIBLE (SEE COMMENT)":"PRESENT",
	"SCARRED DERMIS":"PRESENT",
	"PRESENT, EARLY":"PRESENT",
	"POSSIBLE PARTIAL":"PRESENT",
	"FOCAL REGRESSION":"PRESENT",
	"PROBABLE/PARTIAL":"PRESENT",
	"MELANOPHAGIC":"PRESENT",
	"PRESENT FOCALLY":"PRESENT",
	"FOCAL, SUPERFICIAL":"PRESENT",
	"POSSIBLE, PARTIAL":"PRESENT",
	"PROBABLE PARTIAL":"PRESENT",
	"PRESENT; FOCAL":"PRESENT",
	"FOCAL, PARTIAL REGRESSION":"PRESENT",
	"FOCAL DERMAL REGRESSION":"PRESENT",
	"FOCAL PARTIAL":"PRESENT",
	"PRESENT (FOCAL)":"PRESENT",
	"PRESENT, FOCALLY":"PRESENT",
	"PROBABLE (INFLAMMATORY REGRESSION)":"PRESENT",
	"PRESENT, FOCAL":"PRESENT",
	"POSSIBLE":"PRESENT",
	"FOCALLY PRESENT":"PRESENT",
	"PRESENT (PARTIAL)":"PRESENT",
	"PARTIAL":"PRESENT",
	"PROBABLE":"PRESENT",
	"PARTIAL, PRESENT":"PRESENT",
	"PRESENT/PARTIAL":"PRESENT",
	"PRESENT, MULTIFOCAL":"PRESENT",
	"FOCAL":"PRESENT",
	"PRESENT, PARTIAL":"PRESENT",
	"FOCAL INFLAMMATORY":"PRESENT",
	"PARTIAL/INFLAMMATORY":"PRESENT",
	"EARLY INFLAMMATORY":"PRESENT",
	"FOCAL REGRESSION-LIKE STROMAL RESPONSE":"PRESENT",
	"PRESENT, PARTIAL INFLAMMATORY":"PRESENT",
	"PROBABLE PARTIAL INFLAMMATORY":"PRESENT",
	"PROBABLE PARTIAL, INFLAMMATORY":"PRESENT",
	"PROBABLE, INFLAMMATORY":"PRESENT",
	"INFLAMMATORY, PARTIAL":"PRESENT",
	"POSSIBLE PARTIAL INFLAMMATORY":"PRESENT",
	"PROBABLE (INFLAMMATORY REGRESSION)":"PRESENT",
	"PARTIAL INFLAMMATORY":"PRESENT",
	"FOCAL INFLAMMATORY REGRESSION PRESENT":"PRESENT",
	"POSSIBLE PARTIAL (INFLAMMATORY)":"PRESENT",
	"POSSIBLE PARTIAL":"PRESENT",
	"PROBABLE PARTIAL":"PRESENT",
	"FOCAL DERMAL":"PRESENT",
	"PRESENT, FOCAL":"PRESENT",
	"FOCAL, SUPERFICIAL":"PRESENT",
	"POSSIBLE FOCAL":"PRESENT",
	"PROBABLE PARTIAL, INFLAMMATORY":"PRESENT",
	"PARTIAL REGRESSION PRESENT":"PRESENT",
	"PRESENT, PARTIAL INFLAMMATORY":"PRESENT",
	"PRESENT (FOCAL)":"PRESENT",
	"PARTIAL":"PRESENT",
	"PROBABLE":"PRESENT",
	"PRESENT; SEE COMMENT":"PRESENT",
	"PROBABLE (SEE COMMENT)":"PRESENT",
	"FOCAL":"PRESENT",
	"PRESENT, PARTIAL":"PRESENT",
	"ABSENT (SEE COMMENT)":"ABSENT",
	"NO DEFINITE":"ABSENT",
	"NEGATIVE":"ABSENT",
	"NOT DEFINITE":"ABSENT",
	"ABSENT":"ABSENT",
	"NONE":"ABSENT",
	"NOT IDENTIFIED":"ABSENT",
	"ABSENT (PARTIAL SAMPLE)":"ABSENT",
	"NONE SEEN":"ABSENT",
	"NOT SEEN":"ABSENT",
	"NOT PRESENT":"ABSENT",
	"ABSENT":"ABSENT",
	"NONE":"ABSENT",
	"NOT IDENTIFIED":"ABSENT",
	"NOT IDENTIFIED (PER PRIOR BIOPSY; SEE PART 2)":"ABSENT",
	"FOCALLY PRESENT (SEE COMMENT)":"PRESENT",
	"NOT SEEN":"ABSENT",
	"SEE COMMENT":"UNKNOWN",
	"NOT APPLICABLE":"UNKNOWN",
	"NOT EVALUABLE":"UNKNOWN",
	"UNABLE TO EVALUATE":"UNKNOWN",
	"CANNOT BE DETERMINED":"UNKNOWN",
	"CANNOT BE EVALUATED":"UNKNOWN",
	"UNABLE TO FULLY EVALUATE":"UNKNOWN",
	"CHANGES CONSISTENT WITH INFLAMMATORY REGRESSION RESPONSE":"UNKNOWN",
	"PRESENT, PARTIAL INFLAMMATORY IS IDENTIFIED":"PRESENT",
	"NOT FULLY EVALUABLE":"UNKNOWN",
	"DIFFICULT TO ASSESS DUE TO PRIOR BIOPSY":"UNKNOWN",
	"POSSIBLE PARTIAL, INFLAMMATORY":"PRESENT",
	"SEE NOTE":"SEE NOTE",
	"CANNOT EVALUATE":"UNKNOWN",
	"CANNOT BE ASSESSED DUE TO PRIOR":"UNKNOWN",
	"CANNOT EVALUATE":"UNKNOWN",
	"PRESENT (VERTICAL GROWTH PHASE MELANOPHAGIC TYPE)":"PRESENT",
	"PRESENT, NOT IDENTIFIED":"PRESENT",
	"PRESENT (SEE PART )":"PRESENT",
	"PRESENT (INFLAMMATORY REGRESSION)":"PRESENT",
	"STROMAL CHANGES SUGGESTIVE OF REGRESSION":"PRESENT",
	"NO EVALUABLE":"UNKNOWN",
	"N/A":"UNKNOWN",
	"NULL":"UNKNOWN",
	'PRESENT/EXTENSIVE': "PRESENT",
	'PRESENT, EXTENSIVE (SEE COMMENT)':"PRESENT",
}


TIL_dict = {
	"PRESENT (SEE COMMENT)":"PRESENT_NON-BRISK",
	"PRESENT":"PRESENT_NON-BRISK",
	"PRESENT, FOCAL BRISK":"PRESENT_BRISK",
	"PRESENT, BRISK":"PRESENT_BRISK",
	"PRESENT; FOCAL BRISK":"PRESENT_BRISK",
	"PRESENT; BRISK":"PRESENT_BRISK",
	"PRESENT/BRISK":"PRESENT_BRISK",
	"AT LEAST FOCALLY BRISK":"PRESENT_BRISK",
	"FOCAL BRISK":"PRESENT_BRISK",
	"BRISK":"PRESENT_BRISK",
	"BRISK HALO-LIKE TUMOR-INFILTRATING":"PRESENT_BRISK",
	"FOCALLY BRISK":"PRESENT_BRISK",
	"PRESENT, FOCALLY BRISK":"PRESENT_BRISK",
	"PRESENT; FOCALLY BRISK":"PRESENT_BRISK",
	"PRESENT, AT LEAST BRISK":"PRESENT_BRISK",
	"FOCAL, BRISK":"PRESENT_BRISK",
	"PRESENT (BRISK)":"PRESENT_BRISK",
	"PRESENT, BRISK":"PRESENT_BRISK",
	"PRESENT; BRISK":"PRESENT_BRISK",
	"FOCAL BRISK":"PRESENT_BRISK",
	"BRISK":"PRESENT_BRISK",
	"PRESENT, AT LEAST FOCAL BRISK":"PRESENT_BRISK",
	"FOCALLY BRISK":"PRESENT_BRISK",
	"PRESENT, FOCALLY BRISK":"PRESENT_BRISK",
	"FOCAL, BRISK":"PRESENT_BRISK", 
	"PRESENT (NON-BRISK)":"PRESENT_NON-BRISK",
	"NON-BRSK":"PRESENT_NON-BRISK",
	"FOCALLY NON-BRISK":"PRESENT_NON-BRISK",
	"AT LEAST NON-BRISK":"PRESENT_NON-BRISK",
	"FOCAL, NON-BRISK":"PRESENT_NON-BRISK",
	"PRESENT; NON-BRISK LYMPHOID AGGREGATES":"PRESENT_NON-BRISK",
	"NON-BRISK/PATCHY":"PRESENT_NON-BRISK",
	"PRESENT NOT-BRISK":"PRESENT_NON-BRISK",
	"PRESENT; NON-BRISK":"PRESENT_NON-BRISK",
	"NON-BRISK/INFILTRATIVE":"PRESENT_NON-BRISK",
	"PRESENT, AT LEAST NON-BRISK":"PRESENT_NON-BRISK",
	"NON-BRISK":"PRESENT_NON-BRISK",
	"NON-BRISK/NOT FULLY EVALUABLE":"PRESENT_NON-BRISK",
	"PRESENT, NON-BRISK":"PRESENT_NON-BRISK",
	"FOCAL PRESENT, NON-BRISK":"PRESENT_NON-BRISK",
	"PRESENT/NON-BRISK":"PRESENT_NON-BRISK",
	"FOCAL, AT LEAST NON-BRISK":"PRESENT_NON-BRISK",
	"PRESENT, NOT BRISK":"PRESENT_NON-BRISK",
	"PRESENT, NONE BRISK":"PRESENT_NON-BRISK",
	"PRESENT; NON-BRISK":"PRESENT_NON-BRISK",
	"PRESENT, AT LEAST NON-BRISK":"PRESENT_NON-BRISK",
	"NON-BRISK":"PRESENT_NON-BRISK",
	"NON BRISK":"PRESENT_NON-BRISK",
	"PRESENT, NON-BRISK":"PRESENT_NON-BRISK",
	"FOCAL AGGREGATES":"PRESENT_NON-BRISK",
	"PRESENT; MINIMAL":"PRESENT_NON-BRISK",
	"PRESENT; FOCAL":"PRESENT_NON-BRISK",
	"PATCHY LYMPHOID AGGREGATES":"PRESENT_NON-BRISK",
	"VERY FOCALLY PRESENT, NON-BRISK":"PRESENT_NON-BRISK",
	"FOCALLY PRESENT, NON-BRISK":"PRESENT_NON-BRISK",
	"FOCAL; NON-BRISK":"PRESENT_NON-BRISK",
	"FOCAL":"PRESENT_NON-BRISK",
	"PRESENT, MINIMAL":"PRESENT_NON-BRISK",
	"FOCAL, NON-BRISK":"PRESENT_NON-BRISK",
	"NON-BRISK, FOCAL":"PRESENT_NON-BRISK",
	"FOCAL; NON-BRISK":"PRESENT_NON-BRISK",
	"MINIMAL":"ABSENT",
	"ABSENT":"ABSENT",
	"NOT SEEN":"ABSENT",
	"ABSENT IN THIS SAMPLE":"ABSENT",
	"NOT PRESENT":"ABSENT",
	"MINIMAL":"ABSENT",
	"ABSENT":"ABSENT",
	"NOT IDENTIFIED":"ABSENT",
	"CANNOT BE DETERMINED":"NON-EVALUABLE",
	"CANNOT BE FULLY ASSESSED":"NON-EVALUABLE",
	"UNABLE TO FULLY EVALUATE":"NON-EVALUABLE",
	"UNABLE TO EVALUATE":"NON-EVALUABLE",
	"NOT EVALUABLE (VGP TOO SMALL)":"NON-EVALUABLE",
	"NOT FULLY EVALUABLE":"NON-EVALUABLE",
	"CANNOT ASSESS":"NON-EVALUABLE",
	"NOT EVALUABLE":"NON-EVALUABLE",
	"NO EVALUABLE":"NON-EVALUABLE",
	"CANNOT FULLY EVALUATE":"NON-EVALUABLE",
	"NOT FULLY EVALUABLE (AT LEAST NON-BRISK)":"NON-EVALUABLE",
	"NOT EVALUATED":"NON-EVALUABLE",
	"NOT ASSESSABLE":"NON-EVALUABLE",
	"NOT EVALUABLE DUE TO EARLY STAGE OF VERTICAL GROWTH PHASE":"NON-EVALUABLE",
	"NOT FULLY EVALUABLE":"NON-EVALUABLE",
	"NOT APPLICABLE":"NON-EVALUABLE",
	"NOT EVALUABLE":"NON-EVALUABLE",
	"NOT OBSERVED":"NON-EVALUABLE",
	"NOT SEEN":"NON-EVALUABLE",
	"NOT ASSESSABLE":"NON-EVALUABLE",
	"PRESENT, NON-BRISK VASCULAR/LYMPHATIC INVASION NOT IDENTIFIED":"PRESENT_NON-BRISK",
	"PRESENT, NON-BRISK PERINEURAL INVASION":"PRESENT_NON-BRISK",
	"PRESENT, NON-BRISK VASCULAR/LYMPHATIC INVASION HIGHLY SUSPICIOUS":"PRESENT_NON-BRISK",
	"PRESENT, NON-BRISK VASCULAR/LYMPHATIC INVASION":"PRESENT_NON-BRISK",
	"ABSENT VASCULAR/LYMPHATIC INVASION":"ABSENT",
	"REFER TO INITIAL BIOPSY":"UNKNOWN",
	"SEE COMMENT":"SEE COMMENT",
	"SPARSE PRESENT":"PRESENT_NON-BRISK",
	"NON-BRISK":"PRESENT_NON-BRISK",
	"NON-BRISK, SPARSE":"PRESENT_NON-BRISK",
	"NON-EVALUABLE":"UNKNOWN",
	"NULL":"UNKNOWN",
	"PRESENT, MILD": "PRESENT_NON-BRISK",
	"NON-BRISK, INFILTRATIVE": "PRESENT_NON-BRISK",
	"PRESENT AT LEAST NON-BRISK":"PRESENT_NON-BRISK"
	## add new
	# "NON-BRISK" :"PRESENT_NON-BRISK",         
	# "PRESENT, NON-BRISK" :"PRESENT_NON-BRISK",            
	# "PRESENT; NON-BRISK":"PRESENT_NON-BRISK",     
	# "BRISK":"PRESENT_BRISK",
	# "FOCAL, BRISK": "PRESENT_FOCAL",
	# "PRESENT; BRISK" :"PRESENT_BRISK",
	# "FOCAL BRISK":"PRESENT_FOCAL"
}


ulceration_dict = {
	"PRESENT, EXTENSIVE":"PRESENT",
	"PRESENT, EXTENSIVE":"PRESENT",
	"PRESENT, FOCAL":"PRESENT",
	"PRESENT (SEE COMMENT)":"PRESENT",
	"NOT IDENTIFIED":"ABSENT",
	"NONE":"ABSENT",
	"PRESENT, 3.5 MM (PRIOR BIOPSY AT THIS SITE MAY HAVE RESULTED IN ULCERATION)":"PRESENT",
	"PRESENT IN ORIGINAL BIOPSY, PER REPORT":"PRESENT",
	"PRESENT (FOCAL)":"PRESENT",
	"PRESENT; SEE COMMENT":"PRESENT",
	"ABSENT (SMALL EROSION CONSISTENT WITH TRAUMA PRESENT)":"ABSENT",
	"ABSENT (HEALED BIOPSY SITE; SEE PREVIOUS BIOPSY BS-09-W21559)":"ABSENT",
	"ABSENT (SEE COMMENT)":"ABSENT",
	"PRESENT (EXTENSIVE)":"PRESENT",
	"PRESENT EXTENSIVELY":"PRESENT",
	"NOT SEEN":"ABSENT",
	"PRESENT (ON BIOPSY SPECIMEN)":"PRESENT",
	'FOCALLY PRESENT':"PRESENT",
	"SEE COMMENT":"SEE COMMENT",
	"PROBABLE (SMALL SAMPLE)":"PRESENT", 
	"SURFACE EROSION":"UNKNOWN",
	"POSSIBLE":"PRESENT",
	"NOT DEFINITE":"ABSENT",
	"SEE NOTE":"UNKNOWN",
	"NO DEFINITE":"ABSENT",
	"THERE IS A FOCAL AREA CONSISTENT":"PRESENT",
	"MINIMAL, FOCAL":"PRESENT",
	"PROBABLE (SEE COMMENT)":"PRESENT", 
	"MINUTE FOCUS OF UNCERTAIN PROGNOSTIC SIGNIFICANCE":"UNKNOWN",
	"NOT EVALUABLE (S/P BIOPSY)', 'NOT EVALUABLE (SEE NOTE)":"UNKNOWN", 
	"EQUIVOCAL":"PRESENT", 
	"SEE COMMENT":"UNKNOWN",
	"NULL":"UNKNOWN",
	"(SEE COMMENT)":"UNKNOWN", 
	"ABSENT (PARTIAL SAMPLE)":"ABSENT", 
	"ABSENT, SEE NOTE":"ABSENT",
	"PRESENT (ULCER MAY RELATE TO PRIOR BIOPSY)":"PRESENT",
	"":"?",
	"EQUIVOCAL (SEE COMMENT)":"?",
	"PROBABLE FOCAL":"?",
	"FOCAL, POSSIBLY TRAUMA-INDUCED":"?",
	"PRESENT (INCIPIENT)":"PRESENT",
	"REFER TO INITIAL BIOPSY (SEE COMMENT)":"?",
	"PRESENT (PREVIOUS BIOPSY)":"PRESENT",
	"BIOPSY SITE CHANGES":"?",
	"FAVOR PRESENT":"PRESENT",
	"EARLY INCIPIENT ULCERATION POSSIBLE (MICROSCOPIC FOCUS)":"?",
	"NO UNEQUIVOCAL":"?",
	"PRESENT FOCALLY":"PRESENT",
	"DIFFICULT TO ASSESS DUE TO FRAGMENTATION":"?",
	"PRESENT (HOWEVER, PREVIOUSLY BIOPSIED)":"PRESENT",
	"NO DEFINITE ULCERATION":"?",
	"FOCAL EQUIVOCAL":"?",
	"NOT IDENTIFIED (SEE COMMENT)":"?",
	"POSSIBLE MINUTE FOCUS":"?",
	"NOT SEEN (EPIDERMAL EROSION)":"?",
	"EARLY, FOCAL":"?",
	"PRESENT (EARLY ULCERATION)":"PRESENT",
	"NO IDENTIFIED":"?",
	"CANNOT FULLY EVALUATE":"?",
	"NOT EVALUABLE (SEE NOTE)":"?",
	"NOT SEEN (SEE NOTE)":"?",
	"FOCAL BUT POSSIBLY TRAUMA-INDUCED":"?",
	"NOT DEFINITE (SURFACE DISRUPTED)":"?",
	"PRESENT, BUT POSSIBLY DUE TO PRIOR BIOPSY":"PRESENT",
	"MINUTE FOCUS, POSSIBLY TRAUMA RELATED":"?",
	"PRESENT (SEE PART \"B\")":"PRESENT",
	"PRESENT (SUPERFICIAL TRAUMA NOT EXCLUDED)":"PRESENT",
	"PRESENT (SEE BS-11-E31033)":"PRESENT",
	"FOCAL (POSSIBLY RELATED TO PREVIOUS BIOPSY)":"?",
	"EXTENSIVE":"?",
	"NOT EVALUABLE (PRIOR PROCEDURE)":"?",
	"PRESENT (PREVIOUSLY BIOPSIED)":"PRESENT",
	"PRESENT (NOTE PREVIOUSLY BIOPSIED)":"PRESENT",
	"NOT IDENTIFIED; SEE COMMENT":"?",
	"NOT EVALUABLE":"?",
	"PRESENT (IN PRIOR BIOPSY)":"PRESENT",
	"PROBABLE, INCIPIENT":"?",
	"EQUIVOCAL (NO DEFINITE)":"?",
	"PROBABLE":"?",
	"FOCAL, PROBABLE":"?",
	"NOT EVALUABLE DUE TO PRIOR SURGERY":"?",
	"NOT DETECTED":"?",
	"NEGATIVE":"?",
	"PRESENT, BUT POSSIBLY TRAUMA-INDUCED":"PRESENT",
	"NOT EVALUABLE (S/P BIOPSY)":"?",
	"FOCAL":"?",
	"NONE SEEN":"?",
	"CANNOT FULLY ASSESS DUE TO INCOMPLETE SURFACE REPRESENTATION":"?",
	"PRESENT (S/P PRIOR SHAVE BIOPSY)":"PRESENT",
	"PRIOR SURGICAL SITE":"?",
	"NOT PRESENT":"ABSENT",
	"PRESENT (SEE NOTE)":"PRESENT",
	"NOT IDENTIFIED (PER PRIOR BIOPSY; SEE PART 2)":"ABSENT",
	"DEFINITIVE ULCERATION NOT IDENTIFIED":"ABSENT",
	"DEFINITIVE ULCERATION IS NOT SEEN IN THE MATERIAL PROVIDED":"ABSENT",
	"NO DEFINITIVE ULCERATION IDENTIFIED":"ABSENT"
}



# feature_list = ["subtype", "intraepidermal component", "vertical growth phase", "ulceration", "regression","mitotic rate", "tumor-infiltrating lymphocytes", 
	# "perineural invasion", "vascular/lymphatic invasion", "microscopic satellites", "precursor"]


def load_standlize_label_excel(excel_file):
	all_dict = {}
	print(excel_file)
	the_data = pd.read_excel(excel_file,sheet_name=1)
	column_num = len(the_data.columns)
	for idx in range(0, column_num, 2):
		original_term = the_data.columns[idx]
		standard_term = the_data.columns[idx+1]
		new_df = the_data[[original_term,standard_term]][1:]
		the_dict = {}
		for index, row in new_df.iterrows():
			ori_result = row[original_term]
			if pd.isna(ori_result):
				continue 
			stand_result = row[standard_term]
			if type(stand_result) == str:
				stand_result = stand_result.upper()
			the_dict[ori_result] = stand_result
		if original_term.lower() != 'cell type':
			all_dict[original_term.lower()] = the_dict
	return all_dict


all_dict = load_standlize_label_excel("Standard_result-final.xlsx")






class PathologyReport:
	def __init__(self, report_text):
		self.patient_name = "NULL"
		self.date_of_received = "NULL"
		self.assertion_number = "NULL"
		self.date_of_birth = "NULL"
		self.date_of_dead = "NULL"
		self.sex = "NULL"
		self.mrn = "NULL"
		self.phone = "NULL"
		self.address = "NULL"
		self.referring_physician = "NULL"
		self.specimen = "NULL"
		self.section_num = 0
		self.section_texts = []
		self.section_results = []
		self.invasive_depth = None
		self.malignant_depth = False
		self.sentinel_lymph_node = False
		self.original_text = report_text
		self.load_pr_text(report_text)
		self.patient_encode = self.mrn
		self.patient_second_encode = self.patient_name+"/"+self.date_of_birth+"/"+self.sex

	def load_pr_text(self, report_text):
		sections, front = split_each_report_to_section(report_text)
		demographic = sections[0]
		# print("demographic:", demographic)
		self.parse_demographic(demographic)
		# print(self.patient_name, self.sex, self.date_of_birth, self.specimen)
		section_num = len(sections) - 1
		for idx in range(section_num):
			section_text = sections[idx+1]
			section_result = SectionResult(section_text, self.date_of_received)
			print(self.assertion_number, section_result.note_pair)
			self.section_texts.append(section_text)
			if section_result.invasive_depth:
				if self.invasive_depth:
					if self.malignant_depth and section_result.malignant_depth:
						self.invasive_depth = "MULTIPLE"
					elif section_result.malignant_depth:
						self.invasive_depth = section_result.invasive_depth
						self.malignant_depth = True
				else:
					self.invasive_depth = section_result.invasive_depth
					self.malignant_depth = section_result.malignant_depth
			self.section_results.append(section_result)
			# print("note pair:",section_result.note_pair)

		self.section_num = len(self.section_results)


	def parse_demographic(self, demographic_text):
		if '\r\n' in demographic_text:
			lines = demographic_text.split("\r\n")
		else:
			lines = demographic_text.split("\n")
		demo_dict = {}
		for line in lines:
			split_list = line.strip("\"\"\n").split("	")
			len_list = len(split_list)
			if len_list > 0:
				if len_list%2 == 0:
					for idx in range(int(len_list/2)):
						demo_dict[split_list[2*idx].strip(":\"")] = split_list[2*idx+1].strip("\"")
				elif len_list == 1:
					if ":\"" in split_list[0]:
						new_list = split_list[0].split(":\"")
						if len(new_list) == 3:
							middle_split = new_list[1].split("\" \"")
							demo_dict[new_list[0].strip(" \"")] = middle_split[0].strip(" \"")
							demo_dict[middle_split[1].strip(" \"")] = new_list[2].strip(" \"")

					
		# print("demo dict:",demo_dict)
		# exit(0)
		the_string = "Patient"
		if the_string in demo_dict:
			self.patient_name = demo_dict[the_string]
		the_string = "Date Received"
		if the_string in demo_dict:
			self.date_of_received = demo_dict[the_string]
		the_string = "Date of Birth"
		if the_string in demo_dict:
			self.date_of_birth = demo_dict[the_string]
		the_string = "Sex"
		if the_string in demo_dict:
			self.sex = demo_dict[the_string]
		the_string = "Accession No"
		if the_string in demo_dict:
			self.assertion_number = demo_dict[the_string]
		the_string = "Referring Physician"
		if the_string in demo_dict:
			self.referring_physician = demo_dict[the_string]
		the_string = "Phone"
		if the_string in demo_dict:
			self.phone = demo_dict[the_string]
		the_string = "MRN"
		if the_string in demo_dict:
			self.mrn = demo_dict[the_string]
		the_string = "Address"
		if the_string in demo_dict:
			self.address = demo_dict[the_string]
		if '''"Specimen(s):"''' in lines:
			the_id = lines.index('''"Specimen(s):"''')
			self.specimen = lines[the_id+1]
		

class SectionResult:
	def __init__(self, section_text, report_date):
		self.date = report_date
		self.consult = False
		self.body_part = "NULL"
		self.content = "NULL" ## e.g. invasive to a depth of 2.5 mm, anatomic level III/early IV; mrargins appear negative
		self.note_pair = {} ## pairwise note, e.g. {Subtype:Unclassified, Regression:Not Seen}
		self.comment = "NULL"
		self.invasive_depth = None
		self.malignant_depth = False ## if the invasive depth is right after "malignant manaloma"
		self.subsection_results = [] ##
		self.load_section_text(section_text)
		self.section_text = section_text


	def load_section_text(self, section_text):
		## TODO:  Extract the section text
		title, num_subpart, note_list, the_date, invasive_depth, malignant_depth = parse_each_diagnosis(section_text)
		if invasive_depth:
			self.invasive_depth = invasive_depth
			self.malignant_depth = malignant_depth
		if the_date != "NULL":
			self.date = the_date
		if note_list[-1] and len(note_list[-1][-1]) > 0:
			for a in note_list[-1][-1]:
				pair = a.split("->")
				key = pair[0].lower()
				if key == "tumor infiltrating lymphocytes":
					key = "tumor-infiltrating lymphocytes"
				self.note_pair[key] = pair[1]

		# return section_result_list
		# for each_section in section_result_list:
		
		

		

class Patient:
	def __init__(self, each_report=None):
		self.last_name = "NULL"
		self.first_name = "NULL"
		self.name = "NULL"
		self.date_of_birth = "NULL"
		self.date_of_dead = "NULL"
		self.date_of_lastvisit = "NULL"
		self.sex = "NULL"
		self.mrn = "NULL"
		self.race = "NULL"
		self.ethnicity = "NULL"
		self.marital_status = "NULL"
		self.education_level = "NULL"
		self.veteran_status = "NULL"
		self.patient_encode = "NULL"
		self.patient_second_encode = "NULL"
		self.invasive_depth = []
		self.pathology_report_list = []
		self.pathology_results = []
		if each_report:
			self.load_report(each_report)


	def load_report(self, each_report):
		self.name = each_report.patient_name
		self.date_of_birth = each_report.date_of_birth
		self.date_of_dead = each_report.date_of_dead
		self.sex = each_report.sex 
		self.mrn = each_report.mrn 
		self.patient_encode = each_report.patient_encode
		self.patient_second_encode = each_report.patient_second_encode
		self.invasive_depth.append(each_report.invasive_depth)
		self.pathology_report_list.append(each_report)
		self.pathology_results.append(each_report.section_results)

	def add_report(self, each_report):
		self.pathology_report_list.append(each_report)

	def config_from_dict(self, config_dict):
		if "last_name" in config_dict:
			self.last_name = config_dict['last_name']
		if "first_name" in config_dict:
			self.first_name = config_dict['first_name']
		if "name" in config_dict:
			self.name = config_dict['name']
		if "date_of_birth" in config_dict:
			self.date_of_birth = config_dict['date_of_birth']
		if "date_of_dead" in config_dict:
			self.date_of_dead = config_dict['date_of_dead']
		if "date_of_lastvisit" in config_dict:
			self.date_of_lastvisit = config_dict['date_of_lastvisit']
		if "sex" in config_dict:
			self.sex = config_dict['sex']
		if "mrn" in config_dict:
			self.mrn = config_dict['mrn']
		if "race" in config_dict:
			self.race = config_dict['race']
		if "ethnicity" in config_dict:
			self.ethnicity = config_dict['ethnicity']
		if "marital_status" in config_dict:
			self.marital_status = config_dict['marital_status']
		if "education_level" in config_dict:
			self.education_level = config_dict['education_level']
		if "veteran_status" in config_dict:
			self.veteran_status = config_dict['veteran_status']


	def get_dod_etc(self, mrn_dead, other_dead):
		the_target = None
		if self.patient_encode in mrn_dead:
			the_target = mrn_dead[self.patient_encode]
			# self.date_of_dead = mrn_dead[self.patient_encode].date_of_dead
			# self.date_of_lastvisit = mrn_dead[self.patient_encode].date_of_lastvisit
		elif self.patient_second_encode in other_dead:
			the_target = other_dead[self.patient_second_encode]
		if the_target:
			self.date_of_dead = the_target.date_of_dead
			self.date_of_lastvisit = the_target.date_of_lastvisit
			self.race = the_target.race
			self.ethnicity = the_target.ethnicity
			self.marital_status = the_target.marital_status
			self.education_level = the_target.education_level
			self.veteran_status = the_target.veteran_status


	

def split_all_reports_to_list(all_reports):
	## split the all reports into report list, based on the occurence of "Accession No:".
	print("Split reports...")
	report_list = []
	with open(all_reports,'r') as infile:
		lines = infile.readlines()
		report_text = ""
		for line in lines:
			if line.startswith('''"Accession No:"'''):
				if report_text != "":
					report_list.append(report_text)
				report_text = line
			else:
				if report_text != "":
					report_text += line
	if report_text != "":
		report_list.append(report_text)
	print("%s reports founded."%len(report_list))
	return report_list



def split_each_report_to_section(report_text):
	## split each report into section, including the [demographic information, Section A, Section B,..]
	## each section (A, B, C) is corresponding to the words ("A. SKIN...", "B. CONSULT ...") of the report.
	# extract the demographic information
	if ('''"Results:"	"PATHOLOGIC DIAGNOSIS:"''' not in report_text) \
	and ('''"Results:"	"REPORT NOT FINAL"''' not in report_text) \
	and ('''"Results:"	"FINAL ANATOMIC DIAGNOSIS:"''' not in report_text)\
	and ('''"Results:"	"FINAL DIAGNOSIS:"''' not in report_text)\
	and ('''"Results:" "PATHOLOGIC DIAGNOSIS:"''' not in report_text)\
	and ('''"Results:" "PATHOLOGIC DIAGNOSIS:"''' not in report_text)\
	and ('''""Results:""	""PATHOLOGIC DIAGNOSIS:""''' not in report_text)\
	and ('''"Results:"	"RESULT:"''' not in report_text):
		print("Error, not result found.")
		print(report_text)
		exit(0)
	if '''"Results:"	"PATHOLOGIC DIAGNOSIS:"''' in report_text:
		demographic, result = report_text.split('''"Results:"	"PATHOLOGIC DIAGNOSIS:"''')
	elif '''"Results:"	"REPORT NOT FINAL"''' in report_text:
		demographic, result = report_text.split('''"Results:"	"REPORT NOT FINAL"''')
	elif '''"Results:"	"RESULT:"''' in report_text:
		demographic, result = report_text.split('''"Results:"	"RESULT:"''')
	elif '''"Results:"	"FINAL DIAGNOSIS:"''' in report_text:
		demographic, result = report_text.split('''"Results:"	"FINAL DIAGNOSIS:"''')
	elif '''"Results:" "PATHOLOGIC DIAGNOSIS:"''' in report_text:
		demographic, result = report_text.split('''"Results:" "PATHOLOGIC DIAGNOSIS:"''')
	elif '''""Results:""	""PATHOLOGIC DIAGNOSIS:""''' in report_text:
		demographic, result = report_text.split('''""Results:""	""PATHOLOGIC DIAGNOSIS:""''')
	else:
		demographic, result = report_text.split('''"Results:"	"FINAL ANATOMIC DIAGNOSIS:"''')
	front, section_list = section_split(result)
	## return a list of [demographic, Section A, section B,...]
	return [demographic]+ section_list, front


def section_split(input_text):
	## TODO: split the result into sections (A, B, C), return a list of [demographic, Section A, section B,...]
	## Input: String, the whole text without the demographic information
	## Output: List, a list of string, each string is the text of each section
	## Each section should have a date.
	diagnosis_line = input_text.split("\"HISTORY:\"")[0].strip("\r\"\"\n").split("\"CLINICAL DATA:\"\r\n")[0]
	front = "NULL"
	section_list = re.sub( r"((?<!Dr.)(?<!Drs.)(?<!Pathologist:) [A-Z](\.|\))( |\t))", r"$$$", diagnosis_line).split("$$$")
	if  section_list[0].startswith("A."): ## if A. in the diagnosis
		front = "NULL"
		section_list[0] = section_list[0][2:].strip()
	elif section_list[0].startswith("CONSULT SLIDES") or section_list[0].startswith("CONSULTATION RECEIVED"):
		try:
			front = section_list[0][0:section_list[0].find(':')]
			section_list[0]=section_list[0][section_list[0].index(":") + len(":"):]
		#eil[0:neil.find('r.')]n
		except:
			front = section_list[0]
			del section_list[0]
	return front, section_list


def parse_accession_line(accession_line):
	## parse the line with accession no or patient
	item_list = accession_line.strip("\"\"\n").split("\"	\"")
	assert(len(item_list) == 3)
	new_item_list = item_list[:-1] + item_list[-1].split("\"\t")
	return new_item_list


def reports2patients(all_reports):
	report_list = split_all_reports_to_list(all_reports)
	patient_dict = {}
	count = 1
	for each_report in report_list:
		count += 1
		the_p_report = PathologyReport(each_report)
		# print(the_p_report.patient_encode)
		if len(the_p_report.patient_encode) == 0:
			if the_p_report.patient_second_encode not in patient_dict:
				patient_dict[the_p_report.patient_second_encode] = new_patient
			else:
				patient_dict[the_p_report.patient_second_encode].add_report(the_p_report)
		if the_p_report.patient_encode not in patient_dict:
			new_patient = Patient(the_p_report)
			patient_dict[the_p_report.patient_encode] = new_patient
		else:
			patient_dict[the_p_report.patient_encode].add_report(the_p_report)
		# if count > 10000:
		# 	exit(0)
	print("Report list number:", len(report_list))
	print("Patient num:", len(patient_dict))
	print("Avg reports per patient:", (len(report_list)+0.)/len(patient_dict))
	# d_view = [ (v,k) for k,v in patient_dict.iteritems() ]
	# d_view.sort(reverse=True) # natively sort tuples by first element
	# for v,k in d_view:
	# 	if len(v.pathology_report_list) > 10:
	# 		print("%s: %d" % (k,len(v.pathology_report_list)))
	# the_patient  = patient_dict["CHASE, DANIEL/08/04/1948/M"]
	# for each_report in the_patient.pathology_report_list:
	# 	print(each_report.date_of_received, each_report.referring_physician, each_report.specimen, string2date(each_report.date_of_received).year)
	return patient_dict


def extracted_all_reports(all_reports):
	report_list = split_all_reports_to_list(all_reports)
	structured_report_list = []
	for each_report in report_list:
		
		the_p_report = PathologyReport(each_report)
		structured_report_list.append(the_p_report)
	print("Report list number:", len(structured_report_list))
	return structured_report_list


def string2date(date_string):
	time_list = date_string.split('/')
	if len(time_list)!= 3:
		return datetime.date(2000,1,1)
	return datetime.date(int(time_list[2]), int(time_list[0]), int(time_list[1]))

def load_motality_csv(input_file):
	## load motality csv and extract the name/dateofbirth/dateofdead
	encode2patient_dict = {}
	mrn2patient_dict = {}
	with open(input_file, 'r') as csvfile:
		readlines = csv.reader(csvfile)
		first_line = True
		
		for row in readlines:
			if first_line:
				first_line = False
				continue
			inf_dict = {}
			name = row[0]
			inf_dict["name"] = name
			first_name, last_name = name.split(",")
			sex = row[2]
			mrn = row[3]
			date_of_birth = row[5]
			race = row[6]
			ethnicity = row[7]
			marital_status = row[8]
			education = row[9]
			veteran = row[10]
			date_of_dead = row[12]
			inf_dict["name"] = name
			inf_dict["first_name"] = first_name
			inf_dict["last_name"] = last_name
			inf_dict["sex"] = sex
			inf_dict["mrn"] = mrn
			def fullstring2date(full_string):
				date_string = full_string.split(" ")[0]
				year, month, day = date_string.split('-')
				return date(int(year), int(month), int(day))
			inf_dict["date_of_birth"] = fullstring2date(date_of_birth)
			inf_dict["race"] = race
			inf_dict["ethnicity"] = ethnicity
			inf_dict["marital_status"] = marital_status
			inf_dict["education"] = education
			inf_dict["veteran"] = veteran
			inf_dict["date_of_dead"] = fullstring2date(date_of_dead)
			patient = Patient()		
			patient.config_from_dict(inf_dict)	
			dob_string = inf_dict["date_of_birth"].strftime('%m/%d/%Y')
			key = name.upper()+"/"+dob_string+"/"+sex
			encode2patient_dict[key] = patient
			mrn2patient_dict[mrn] = patient
			dead_age = patient.date_of_dead-patient.date_of_birth
			# print(key, patient.date_of_birth, patient.date_of_dead, dead_age.days/365)
	return encode2patient_dict, mrn2patient_dict


def load_motality_lastvisit_csv(input_file):
	## load motality csv and extract the name/dateofbirth/dateofdead
	encode2patient_dict = {}
	mrn2patient_dict = {}
	with open(input_file, 'r') as csvfile:
		readlines = csv.reader(csvfile)
		first_line = True
		
		for row in readlines:
			if first_line:
				first_line = False
				continue
			inf_dict = {}
			name = row[4]
			inf_dict["name"] = name
			the_names = name.strip(',').split(",")
			first_name, last_name = the_names[0], the_names[1]
			sex = row[6]
			mrn = row[0]
			date_of_birth = row[5]
			if date_of_birth == "1885-12-12 00:00:00.0000000": ## correct one error in data
				date_of_birth = "1942-11-16  12:00:00 AM"
			race = row[8]
			ethnicity = row[9]
			marital_status = row[10]
			education = row[11]
			veteran = row[12]
			date_of_dead = row[14]
			date_of_lastvisit = row[15]
			inf_dict["name"] = name
			inf_dict["first_name"] = first_name
			inf_dict["last_name"] = last_name
			inf_dict["sex"] = sex
			inf_dict["mrn"] = mrn
			def fullstring2date(full_string):
				if full_string == "NULL":
					return "NULL"
				date_string = full_string.split(" ")[0].replace('/','-')
				a, b, c = date_string.split('-')
				if int(a) > 1000:
					year, month, day  = a, b, c
				else:
					month, day, year = a, b, c
				return date(int(year), int(month), int(day))

			if date_of_birth == "NULL" and row[3] != "NULL":
				month, day, year = row[3].split('/')
				year = int(year)
				if year > 2020:
					year = year -100
				inf_dict["date_of_birth"] = date(year, int(month), int(day))
			else:
				inf_dict["date_of_birth"] = fullstring2date(date_of_birth)
			inf_dict["race"] = race
			inf_dict["ethnicity"] = ethnicity
			inf_dict["marital_status"] = marital_status
			inf_dict["education"] = education
			inf_dict["veteran_status"] = veteran
			# if date_of_dead == "NULL" and row[4] != "NULL":
			# 	month, day, year = row[4].split('/')
			# 	inf_dict["date_of_dead"] = date(int(year), int(month), int(day))
			# else:
			inf_dict["date_of_dead"] = fullstring2date(date_of_dead)
			inf_dict['date_of_lastvisit'] = fullstring2date(date_of_lastvisit)
			patient = Patient()		
			patient.config_from_dict(inf_dict)
			if inf_dict["date_of_birth"] != "NULL":
				dob_string = inf_dict["date_of_birth"].strftime('%m/%d/%Y')
			else:
				dob_string = "NULL"
			key = name.upper()+"/"+dob_string+"/"+sex
			encode2patient_dict[key] = patient
			mrn2patient_dict[mrn] = patient
			if patient.date_of_dead != "NULL":
				dead_age = patient.date_of_dead-patient.date_of_birth
			# print(key, patient.date_of_birth, patient.date_of_dead, dead_age.days/365)
	return encode2patient_dict, mrn2patient_dict


def parse_each_diagnosis(each_diagnosis_line):
	## parse each diagnosis within A/B/C. num_subpart counts the number of sub part within A/B.. i.e. A. 1. 2.
	## return [title, num_subpart, [[sub_title,content, comment, note_pair_list],[...]]]
	## add date if possible
	global the_count, sample_dict
	sentinel_lymph_node = None
	sentinel_line = each_diagnosis_line.lower().replace("non-sentinel", "sent-no-tinel").strip("\r\n").replace("\t", " ").replace(", by report","")
	sentinel_line = re.sub(r'\s+', " ", sentinel_line)
	if " note:" in sentinel_line:
		sentinel_line = sentinel_line.split(" note:")[0]
	if " comment:" in sentinel_line:
		sentinel_line = sentinel_line.split(" comment:")[0]

	search_result = re.search(r'\(\d[:|/]\d\)', sentinel_line)
	sentinel_line = sentinel_line.replace(" is ", " ")
	first_part = sentinel_line.split(":")[0]
	if "sentinel lymph node" in first_part or "sentinel node" in first_part:
		if "no tumor seen" in sentinel_line or "negative for tumor" in sentinel_line \
		or "with no evidence of metastasis" in sentinel_line or\
		 "no tumor present" in sentinel_line or "no metastatic tumor" in sentinel_line \
		or "negative for metastasis" in sentinel_line or "negative for melanoma" in sentinel_line \
		or "negative for metastatic" in sentinel_line:
			sentinel_lymph_node = 0 
		elif search_result:
			result = search_result.group().strip("()")[0]
			sentinel_lymph_node = result
		else:
			# print(the_count, each_diagnosis_line)
			sample_dict.append(sentinel_line)
		# the_count Pres
	title = ""
	note_list = []
	the_date = "NULL"
	num_subpart = 0
	each_diagnosis_line = re.sub(r' +\t+ +', ' ', each_diagnosis_line)
	each_diagnosis_line = re.sub(r'\t+ +', ' ', each_diagnosis_line)
	each_diagnosis_line = re.sub(r'\t+', ' ', each_diagnosis_line)
	fit_invasive_text = each_diagnosis_line.lower().replace("microinvasive", "invasive")
	invasive_depth = re.search(r"invasive to((?!mm).)+\d mm", fit_invasive_text)
	malignant_depth = False
	if "malignant melanoma, invasive to" in fit_invasive_text:
		malignant_depth = True
	global general_count
	if invasive_depth:
		general_count = general_count + 1
		invasive_depth = re.findall("([0-9]+[,.]+[0-9]+)", invasive_depth.group())


	# if "CONSULT SLIDES FROM " in each_diagnosis_line:
	# print(each_diagnosis_line)
	the_date = re.findall(r"\d{1,2}/\d{1,2}/\d{2,4}\)",each_diagnosis_line)
	# print(each_diagnosis_line)
	# print(the_date)
	
	if len(the_date) >= 1:
		the_date = the_date[0]
		pair = the_date.strip(')').split('/')
		month = int(pair[0])
		day = int(pair[1])
		year = int(pair[2])
		if year > 1000:
			the_date = date(year, month, day)
		elif year < 19:
			the_date = date(2000+year, month, day)
		else:
			the_date = date(1900+year, month, day)
	elif len(the_date) < 1:
		the_date = "NULL"
		# exit(0)
	if " 1. " in each_diagnosis_line:
		diagnosis_num_split_list = re.split(r'( \d\. )',each_diagnosis_line)
		title = diagnosis_num_split_list[0]
		del diagnosis_num_split_list[0]
		split_num = len(diagnosis_num_split_list)
		assert(split_num%2==0)
		num_subpart = split_num/2
		
		## for each 1. or 2. and so on
		for idx in range(0, split_num-1, 2):
			sub_title = ""
			# print(diagnosis_num_split_list[idx+1])
			splited_items = re.split('\s{2,}', diagnosis_num_split_list[idx+1], 1)
			if len(splited_items) == 2:
				sub_title, other_context = splited_items
			else:
				sub_title = "INVAILD"
				other_context = splited_items[0]
			[content, comment, note_pair_list] = parse_content_note_comment(other_context)
			note_list.append([sub_title, content, comment, note_pair_list])
	else:
		split_title_context= re.split('\s{2,}',each_diagnosis_line, 1)
		if len(split_title_context) == 2:
			title, other_context = split_title_context
		else:
			return ["INVAILD", 0, [["INVAILD","INVAILD", "INVAILD",[]]], the_date, invasive_depth, malignant_depth]
		[content, comment, note_pair_list] = parse_content_note_comment(other_context)
		note_list = [[title, content, comment, note_pair_list]]
	return [title, num_subpart, note_list, the_date, invasive_depth, malignant_depth]


def cut_phrase(input_text, cut_text):
	if cut_text in input_text:
		input_text = " ".join(input_text.split(cut_text))
	return input_text


def parse_content_note_comment(input_line):
	## parse the content with notes and comment, return content/comment/node_pair_list
	# print("input:", input_line)
	content = ""
	comment = ""
	note_pair_list = []
	note_text = ""

	input_line = cut_phrase(input_line, "Anatomic (Clark) Level:  IV")
	input_line = cut_phrase(input_line, "Anatomic (Clark) Level:  III")
	input_line = cut_phrase(input_line, "Anatomic (Clark) Level:  II")
	input_line = cut_phrase(input_line, "Anatomic (Clark) Level:  I")
	input_line = cut_phrase(input_line, "Anatomic (Clark) Level: At least level III")
	input_line = cut_phrase(input_line, "Anatomic (Clark) Level: At least level IV")
	input_line = cut_phrase(input_line, "Anatomic (Clark) Level: At least level II")
	input_line = cut_phrase(input_line, "Anatomic (Clark) Level: At least level I")
	if "NOTE:  Other attributes include:" in input_line or "(see NOTE) anatomic level IV; margins involved." in input_line or "Maximum Tumor Thickness:" in input_line:
		if "NOTE:  Other attributes include:" in input_line:
			first_split = input_line.split("NOTE:  Other attributes include:")
		elif "(see NOTE) anatomic level IV; margins involved." in input_line:
			first_split = input_line.split("(see NOTE) anatomic level IV; margins involved.")
		elif "Maximum Tumor Thickness:" in input_line:
			first_split = input_line.split("Maximum Tumor Thickness:")
			
		content = first_split[0].strip()

		if " Comment: " in first_split[-1]:
			second_split = first_split[-1].split(" Comment: ")
			if "Immunohistochemistry performed at" in second_split[0]:
				content = content + " " + second_split[0].split("Immunohistochemistry performed at")[-1].strip()
				second_split = [second_split[0].split("Immunohistochemistry performed at")[0].strip(), second_split[-1]]
		elif " COMMENT: " in first_split[-1]:
			second_split = first_split[-1].split(" COMMENT: ")
		elif " Comment/Recommendation: " in first_split[-1]:
			second_split = first_split[-1].split(" Comment/Recommendation: ")
		elif "Immunohistochemistry performed at" in first_split[-1]:
			second_split = first_split[-1].split("Immunohistochemistry performed at")
			second_split[-1] = "Immunohistochemistry performed at" + second_split[-1].strip()
		elif "The immunoperoxidase, immunofluorescence and" in first_split[-1]:
			second_split = first_split[-1].split("The immunoperoxidase, immunofluorescence and")
			second_split[-1] = "The immunoperoxidase, immunofluorescence and" + second_split[-1].strip()
		elif "Provided Melan-A/D2-40" in first_split[-1]:
			second_split = first_split[-1].split("Provided Melan-A/D2-40")
			second_split[-1] = "Provided Melan-A/D2-40" + second_split[-1].strip()
		elif "Case reviewed at Dermatopathology Division Staff Conference." in first_split[-1]:
			second_split = first_split[-1].split("Case reviewed at Dermatopathology Division Staff Conference.")
			second_split[-1] = "Case reviewed at Dermatopathology Division Staff Conference." + second_split[-1].strip()
		elif "This case has been reviewed at the" in first_split[-1]:
			second_split = first_split[-1].split("This case has been reviewed at the ")
			second_split[-1] = "This case has been reviewed at the " + second_split[-1].strip()
		elif "Incidental predominantly dermal nevus; margins appear negative." in first_split[-1]:
			second_split = first_split[-1].split("Incidental predominantly dermal nevus; margins appear negative.")
			second_split[-1] = "Incidental predominantly dermal nevus; margins appear negative." + second_split[-1].strip()
		else:
			second_split = [first_split[-1], ""]
		note_text = second_split[0].strip()
		comment = second_split[-1].strip()
		note_pair_list = parse_note_table_text(note_text)
	else:
		content = input_line
	return [content, comment, note_pair_list]





def parse_note_table_text(input_note_table_line):
	##corner case
	if "Minute focus, possibly trauma    related" in input_note_table_line:
		input_note_table_line = input_note_table_line.replace("Minute focus, possibly trauma    related", "Minute focus, possibly trauma related")
	if "Perineural invasion is not    identified" in input_note_table_line:
		input_note_table_line = input_note_table_line.replace("Perineural invasion is not    identified", "Perineural invasion is not identified")
	if "Perineural invasion is not    evaluable" in input_note_table_line:
		input_note_table_line = input_note_table_line.replace("Perineural invasion is not    evaluable", "Perineural invasion is not evaluable")
	if "No definitive dermal mitosis identified." in input_note_table_line:
		input_note_table_line = input_note_table_line.replace("No definitive dermal mitosis identified.", "No definitive dermal mitosis  identified.")
	input_note_table_line = input_note_table_line.replace('	', '')
	# print(input_note_table_line)
	note_pair = re.split('\s{2,}', input_note_table_line)
	# print(note_pair)
	# exit(0)
	# print note_pair
	if len(note_pair)%2 != 0:
		if note_pair[-1].lower() == "see comment":
			del note_pair[-1]
		else:
			note_pair.append("INVAILD:"+input_note_table_line)
		# print input_note_table_line, len(note_pair)
	assert(len(note_pair)%2==0)
	note_list = []
	for idx in range(0, len(note_pair)-1, 2):
		note_list.append(note_pair[idx] + "->"+note_pair[idx+1])
		if note_pair[idx] not in note_pair_dict:
			note_pair_dict[note_pair[idx]] = 1
		else:
			note_pair_dict[note_pair[idx]] += 1
	return note_list


def patient_analysis_overall_dead(patient_dict, index_string, focus_result_dict = {}):
	reformat_dict = {}
	if index_string == "Tumor-infiltrating lymphocytes":
		reformat_dict = TIL_dict
	elif index_string == "Regression":
		reformat_dict = regression_dict
	dead_dict = {}
	unknown_dict = {}
	for k, patient in patient_dict.iteritems():
		if patient.date_of_dead != "NULL":
			# print(patient.date_of_dead.strftime('%m/%d/%Y'))
			place_holder = dead_dict
		else: 
			place_holder = unknown_dict
		matched = False
		for each_report_results in patient.pathology_results:
			for each_section_result in each_report_results:
				if matched:
					continue
				if index_string in each_section_result.note_pair:
					result = each_section_result.note_pair[index_string].strip(',:;-_. *').upper()
					if result in reformat_dict:
						result = reformat_dict[result]
					if len(focus_result_dict) > 0:
						if  result in focus_result_dict:
							if result in place_holder:
								place_holder[result] += 1 
							else:
								place_holder[result] = 1 
					else:
						if result in place_holder:
							place_holder[result] += 1 
						else:
							place_holder[result] = 1 
					matched = True
	
	dead_count = 0 
	for k,v in dead_dict.iteritems():
		dead_count += v
	unknown_count = 0 
	for k,v in unknown_dict.iteritems():
		unknown_count += v
	print("INDEX: "+index_string)
	print("Dead Patient, total number:",dead_count)
	print("%s\t%s\t%s\t%s\t%s"%("Result","Unknown(#"+str(unknown_count)+")","Dead(#"+str(dead_count)+")", "Test","OneTailP-value")).expandtabs(20)
	for k,v in  sorted(dead_dict.items(), key=lambda x: x[1], reverse=True):
		if k in unknown_dict:
			unknown_rate = (unknown_dict[k] +0.)/unknown_count
			dead_rate = (v+0.)/dead_count

			if unknown_rate > dead_rate:
				significance = calculate_two_sample_p_value(dead_rate, unknown_rate, dead_count, unknown_count)
				
				test = "P(unknown)>P(dead)"
			else:
				significance = calculate_two_sample_p_value(unknown_rate, dead_rate, unknown_count, dead_count)
				test = "P(dead)>P(unknown)"
			print("%s\t%s\t%s\t%s\t%s"%(k, round(unknown_rate,4), round(dead_rate,4), test, significance)).expandtabs(20)



	# for k,v in  sorted(dead_dict.items(), key=lambda x: x[1], reverse=True):
	# 	print("   Result: " + str(k)+";    Num:" +str(v)+";    Rate:" + str((v+0.)/dead_count))
	# print("Unknown Patient,  total number:",unknown_count)
	# for k,v in  sorted(unknown_dict.items(), key=lambda x: x[1], reverse=True):
	# 	print ("   Result: " + str(k)+";    Num:" +str(v)+";    Rate:" + str((v+0.)/unknown_count))



def new_string2date(input_string):
	## month/dat/year -> date
	if not isinstance(input_string, str):
		return input_string
	if input_string == "NULL":
		return input_string
	pair = input_string.split(' ')[0].split("/")
	return date(year=int(pair[2]), month = int(pair[0]), day = int(pair[1]))


def calculate_two_sample_p_value(p1, p2, n1, n2):
	import numpy as np
	import scipy.special as scsp
	import scipy.stats as st
	dif_p = (p1*n1+p2*n2)/(n1+n2)
	down = dif_p*(1-dif_p)*(1.0/n1+1.0/n2)
	z_score = (p1-p2)/(down)**(0.5)
	return 0.5 * (1 + scsp.erf(z_score / np.sqrt(2)))


def patient_analysis_period_dead(patient_dict, index_string, investigate_months=12, focus_result_dict = {}):
	
	reformat_dict = {}
	if index_string == "tumor-infiltrating lymphocytes":
		reformat_dict = TIL_dict
	elif index_string == "regression":
		reformat_dict = regression_dict
	elif index_string == "ulceration":
		reformat_dict = ulceration_dict
	dead_dict = {}
	live_dict = {}
	for k, patient in patient_dict.iteritems():
		matched = False
		pat_dod = patient.date_of_dead 
		pat_lastvisit = patient.date_of_lastvisit
		for each_report_results in patient.pathology_results:
			for each_section_result in each_report_results:
				if matched:
					continue
				if index_string in each_section_result.note_pair:
					date_of_examined = new_string2date(each_section_result.date)
					live_months = "NULL"
					at_least_live_months = "NULL"
					if pat_dod != "NULL":
						live_months = relativedelta(pat_dod, date_of_examined).years*12 + relativedelta(pat_dod, date_of_examined).months
					elif pat_lastvisit != "NULL":
						at_least_live_months = relativedelta(pat_lastvisit, date_of_examined).years*12 + relativedelta(pat_lastvisit, date_of_examined).months
					if live_months != "NULL":
						if live_months> investigate_months:
							place_holder = live_dict 
						else:
							place_holder = dead_dict
					elif at_least_live_months != "NULL" and at_least_live_months > investigate_months:
						place_holder = live_dict 
					else:
						continue

					result = each_section_result.note_pair[index_string].strip(',:;-_. *').upper()
					if result in reformat_dict:
						result = reformat_dict[result]
					if len(focus_result_dict) > 0:
						if  result in focus_result_dict:
							if result in place_holder:
								place_holder[result] += 1 
							else:
								place_holder[result] = 1 
					else:
						if result in place_holder:
							place_holder[result] += 1 
						else:
							place_holder[result] = 1 
					matched = True
	dead_count = 0 
	for k,v in dead_dict.iteritems():
		dead_count += v
	live_count = 0 
	for k,v in live_dict.iteritems():
		live_count += v
	print("INDEX: "+index_string)
	print("Focus result:", ";  ".join(list(focus_result_dict)))
	print("Period live analysis, period: %s months"%investigate_months)
	print("Total types of features:", set(dead_dict.keys()).union(set(live_dict.keys())))
	print("%s\t%s\t%s\t%s\t%s"%("Result","Live(#"+str(live_count)+")","Dead(#"+str(dead_count)+")", "Test","OneTailP-value")).expandtabs(20)
	for k,v in  sorted(dead_dict.items(), key=lambda x: x[1], reverse=True):
		if k in live_dict:
			living_rate = (live_dict[k] +0.)/live_count
			dead_rate = (v+0.)/dead_count

			if living_rate > dead_rate:
				significance = calculate_two_sample_p_value(dead_rate, living_rate, dead_count, live_count)
				
				test = "P(live) > P(dead)"
			else:
				significance = calculate_two_sample_p_value(living_rate, dead_rate, live_count, dead_count)
				test = "P(dead) > P(live)"
			print("%s\t%s\t%s\t%s\t%s"%(k, round(living_rate,4), round(dead_rate,4), test, significance)).expandtabs(20)
	# print("Live Patient,  total number:",live_count)
	# for k,v in  sorted(live_dict.items(), key=lambda x: x[1], reverse=True):
	# 	print "   Result: " + str(k)+";    Num:" +str(v)+";    Rate:" + str((v+0.)/live_count)


class GeneReport:
	def __init__(self, report_text):
		self.patient_name = "NULL"
		self.date_of_received = "NULL"
		self.assertion_number = "NULL"
		self.date_of_birth = "NULL"
		self.date_of_dead = "NULL"
		self.sex = "NULL"
		self.mrn = "NULL"
		self.phone = "NULL"
		self.address = "NULL"
		self.referring_physician = "NULL"
		self.specimen = "NULL"
		self.load_gene_report(report_text)
		self.patient_encode = self.mrn
		self.patient_second_encode = self.patient_name+"/"+self.date_of_birth+"/"+self.sex

	def load_gene_report(self, report_text):
		demographic_text, sections = self.split_gene_report_to_sectons(report_text)
		self.parse_demographic(demographic_text)


	def split_gene_report_to_sectons(self, report_text):
		if "\"Results:\"	\"RESULT:\"" in report_text:
			demographic_text, section_text = report_text.split("\"Results:\"	\"RESULT:\"")
		elif "\"DNA VARIANTS:\"" in report_text:
			demographic_text, section_text = report_text.split("\"DNA VARIANTS:\"")
		else:
			print("ERROR format")
			print(report_text)
			exit(0)
		return demographic_text, section_text




	def parse_demographic(self, demographic_text):
		lines = demographic_text.split("\r\n")
		demo_dict = {}
		for line in lines:
			split_list = line.strip("\"\"\n").split("	")
			len_list = len(split_list)
			if len_list > 0 and len_list%2 == 0:
				for idx in range(int(len_list/2)):
					demo_dict[split_list[2*idx].strip(":\"")] = split_list[2*idx+1].strip("\"")
		the_string = "Patient"
		if the_string in demo_dict:
			self.patient_name = demo_dict[the_string]
		the_string = "Date Received"
		if the_string in demo_dict:
			self.date_of_received = demo_dict[the_string]
		the_string = "Date of Birth"
		if the_string in demo_dict:
			self.date_of_birth = demo_dict[the_string]
		the_string = "Sex"
		if the_string in demo_dict:
			self.sex = demo_dict[the_string]
		the_string = "Accession No"
		if the_string in demo_dict:
			self.assertion_number = demo_dict[the_string]
		the_string = "Referring Physician"
		if the_string in demo_dict:
			self.referring_physician = demo_dict[the_string]
		the_string = "Phone"
		if the_string in demo_dict:
			self.phone = demo_dict[the_string]
		the_string = "MRN"
		if the_string in demo_dict:
			self.mrn = demo_dict[the_string]
		the_string = "Address"
		if the_string in demo_dict:
			self.address = demo_dict[the_string]
		if '''"Specimen(s):"''' in lines:
			the_id = lines.index('''"Specimen(s):"''')
			self.specimen = lines[the_id+1]



def read_gene_report(gene_report_file):
	reports = split_all_reports_to_list(gene_report_file)
	gene_report_list = []
	for each_report in reports:
		the_report = GeneReport(each_report)
		gene_report_list.append(the_report)
		# print(the_report.patient_second_encode)
	return gene_report_list


def gene_with_patient_analysis(gene_report_list, patient_dict):
	mrn_dict ={}
	second_encode_patient_dict = {}
	for code, patient in patient_dict.iteritems():
		if patient.patient_second_encode not in second_encode_patient_dict:
			second_encode_patient_dict[patient.patient_second_encode] = patient
	for each_report in gene_report_list:
		mrn_dict[each_report.mrn] = each_report
	match_pathology = 0 
	match_dead_pathology = 0
	for mrn in list(set(mrn_dict.keys())):
		if mrn in patient_dict:
			match_pathology += 1 
			if patient_dict[mrn].date_of_dead != "NULL":
				match_dead_pathology += 1 
				print(mrn, patient_dict[mrn].date_of_birth, mrn_dict[mrn].date_of_received, patient_dict[mrn].date_of_dead)
	print("gene patient num:", len(set(mrn_dict.keys())))
	print(match_pathology, match_dead_pathology)


def report_analysis_feature_with_invasion(patient_dict, index_string, target_results):
	import numpy as np
	index_string = index_string.lower()
	reformat_dict = {}
	print(index_string)
	if index_string == "tumor-infiltrating lymphocytes":
		reformat_dict = TIL_dict
	elif index_string == "regression":
		reformat_dict = regression_dict
	elif index_string == "ulceration":
		reformat_dict = ulceration_dict
	result_pairs = {}
	for key, patient in patient_dict.iteritems():
		for report in patient.pathology_report_list:
			for section in report.section_results:
				if section.invasive_depth:
					the_depth = float(section.invasive_depth[0])
					# if the_depth >20 :
					# 	print(section.section_text)
					if index_string in section.note_pair:
						result = section.note_pair[index_string].strip(',:;-_. *').upper()
						if result in reformat_dict:
							result = reformat_dict[result]
						if target_results:
							if result in target_results:
								if result not in result_pairs:
									result_pairs[result] = [the_depth]
								else:
									result_pairs[result].append(the_depth)
						else:
							if result not in result_pairs:
								result_pairs[result] = [the_depth]
							else:
								result_pairs[result].append(the_depth)
	for key, value in result_pairs.iteritems():
		# if key == "PRESENT_BRISK":
		# 	print(sorted(value))
		value = np.asarray(value)
		print(key, len(value),  np.mean(value), np.std(value))



def patient_cohort_analysis(patient_dict):
	patients = list(patient_dict.values())
	# feature_list = ["Subtype", 
	# "Intraepidermal component", 
	# "Vertical growth phase",
	# "Ulceration",
	# "Regression",
	# "mitotic rate",
	# "Tumor-infiltrating lymphocytes",
	# "perineural invasion",
	# "vascular/lymphatic invasion",
	# "microscopic satellites",
	# "cell type",
	# "precursor"
	# ]
	feature_list = [
	"Ulceration"
	]
	# feature_result_convert_list = [subtype_dict,
	# intraepidermal_component_dict,
	# vertical_growth_phase_dict,
	# ulceration_dict,
	# regression_dict,
	# mitotic_rate_dict,
	# TIL_dict,
	# perineural_invasion_dict,
	# vascular_lymphatic_invasion_dict,
	# microscopic_satellites_dict,
	# cell_type,
	# precursor
	# ]
	feature_result_convert_list = [{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{}
	]
	feature_list = [a.lower() for a in feature_list]
	feature_num = len(feature_list)
	feature_distribution_list = [{} for a in range(feature_num)]
	feature_patient_list = [0]*feature_num
	comb_patient = 0
	error_count = 0
	report_num  = 0
	for patient in patients:
		feature_count_list = [0]*feature_num
		comb_feature_count = 0
		for report in patient.pathology_report_list:
			for section in report.section_results:
				# if len(section.note_pair) > 0:# and len(section.note_pair) < 5 :#and (feature1 not in section.note_pair or feature2 not in section.note_pair or "ulceration" not in section.note_pair) :
				# 	print(report.assertion_number, section.note_pair)
				# 	print('\n')
				# 	exit(0)
				# 	error_count += 1
				for idx, each_feature in enumerate(feature_list):
					if each_feature in section.note_pair:
						feature_count_list[idx] = 1
						result = section.note_pair[each_feature].strip(',:;-_. *').upper()
						## normorize
						if result in feature_result_convert_list[idx]:
							result = feature_result_convert_list[idx][result]
						if result in feature_distribution_list[idx]:
							feature_distribution_list[idx][result] += 1 
						else:
							feature_distribution_list[idx][result] = 1 
		if 	sum(feature_count_list) == feature_num:
			comb_feature_count = 1
		for f_id in range(feature_num):
			feature_patient_list[f_id] += feature_count_list[f_id]
		comb_patient += comb_feature_count
		if comb_feature_count:
			report_num += len(patient.pathology_report_list)
	for f_id in range(feature_num):
		print(feature_list[f_id], feature_patient_list[f_id])
	print("Combine", comb_patient)
	print("report num:", report_num)

	for a in range(feature_num):
		print(feature_list[a], feature_distribution_list[a])
		print("----"*20)
	with open("../name_ulceration.csv", 'wb') as fout:
		csv_writer = csv.writer(fout)
		new_feature_list = []
		for a in feature_list:
			new_feature_list.append(a)
			new_feature_list.append("Standard "+a)
		csv_writer.writerow(new_feature_list)
		result_lists = []
		max_num = 0
		for a in range(feature_num):
			result_lists.append(feature_distribution_list[a].keys())
			if len(feature_distribution_list[a].keys()) > max_num:
				max_num = len(feature_distribution_list[a].keys())
		for a in range(feature_num):
			if len(result_lists[a]) < max_num:
				result_lists[a] += [""]*(max_num-len(result_lists[a]))
		for a in range(max_num):
			new_list = [result_lists[idx][a] for idx in range(feature_num)]
			out_new_list = []
			for b in range(feature_num):
				out_new_list.append(new_list[b])
				out_new_list.append("")
			csv_writer.writerow(out_new_list)



def all_reports_export_to_excel(all_reports, csv_outfile, with_death_info=False):
	structured_report_list = extracted_all_reports(all_reports)
	result_dict = {"Assertion":[], "Patient":[], 'MRN':[], 'Sex':[],  'Date of Birth:':[], 'Date Received:':[]}
	for report in structured_report_list:
		result_dict['Assertion'].append(report.assertion_number)
		result_dict['Patient'].append(report.patient_name)
		result_dict['MRN'].append(report.mrn)
		result_dict['Sex'].append(report.sex)
		result_dict['Date of Birth:'].append(report.date_of_birth)
		result_dict['Date Received:'].append(report.date_of_received)
	for k,v in result_dict.iteritems():
		print(k, len(v))
	pd_data = pd.DataFrame.from_dict(result_dict)
	pd_data.to_csv(csv_outfile, header=True, index=False)





def patient_cohort_export_to_excel(patient_dict, csv_outfile):
	vertical_data = pd.read_excel("vertical_check.xlsx")
	assertion_vertical_dict = dict(zip(vertical_data.Assertion_Num, vertical_data.validate))

	patients = list(patient_dict.values())
	result_dict = {"name":[], 'mrn':[], 'sex':[], 'race':[], 'ethnicity':[], 'marital':[], 'education':[], 'veteran':[], 'age':[], 'age.type':[], 'live_time':[], 'death':[], 'invasion':[],'T.stage':[], 'dob':[], 'Assertion_Num':[], 'Need_manual_check':[]}
	feature_list = list(all_dict.keys())
	all_report_count = 0 
	TIL_report_count = 0
	for each_feature in feature_list:
		result_dict[each_feature] = []
	for patient in patients:
		valid_section = None
		valid_report = None
		valid_date = datetime.date(1800,1,1)
		for report in patient.pathology_report_list:
			for section in report.section_results:
				if len(section.note_pair) > 0:
					if new_string2date(section.date) > valid_date:
						valid_section = section
						valid_report = report
						valid_date = new_string2date(section.date)
		if  valid_section and (not (patient.date_of_dead == "NULL" and patient.date_of_lastvisit == "NULL")):
			result_dict['name'].append(patient.name)
			result_dict['mrn'].append(patient.mrn)
			result_dict['sex'].append(patient.sex)
			result_dict['race'].append(patient.race)
			result_dict['ethnicity'].append(patient.ethnicity)
			if patient.marital_status == "PARTNERED" or patient.marital_status == "MARRIED":
				patient.marital_status = "MARRIED/PARTNERED"
			result_dict['marital'].append(patient.marital_status)
			result_dict['education'].append(patient.education_level)
			result_dict['veteran'].append(patient.veteran_status)
			result_dict['dob'].append(patient.date_of_birth)
			result_dict['Assertion_Num'].append(valid_report.assertion_number)
			for feature in feature_list:
				if feature in valid_section.note_pair:
					the_result = valid_section.note_pair[feature].strip(',:;-_. *').upper()
					the_dict = all_dict[feature]
					if feature == "tumor-infiltrating lymphocytes" and the_result in TIL_dict:
						the_result = TIL_dict[the_result]
					elif feature == "regression" and the_result in regression_dict:
						the_result = regression_dict[the_result]
					if the_result in the_dict:
						the_result = the_dict[the_result]
					if the_result == "N/A":
						the_result = "UNKNOWN"
				else:
					the_result = "UNKNOWN"
				if feature == 'mitotic rate' and the_result == "UNKNOWN":
					the_result = -1
				if feature == 'vertical growth phase' and valid_report.assertion_number in assertion_vertical_dict:
					old_result = the_result
					the_result = assertion_vertical_dict[valid_report.assertion_number].upper()
					print(old_result, the_result, valid_report.assertion_number)
				result_dict[feature].append(the_result)
			if len(valid_section.note_pair) > 0 and "tumor-infiltrating lymphocytes" not in valid_section.note_pair:
				need_manual = True 
			else:
				need_manual = False 
			result_dict['Need_manual_check'].append(need_manual)
			

			
			# print(patient.mrn, section.date, patient.date_of_birth , patient.date_of_lastvisit, patient.date_of_dead)
			dob = new_string2date(patient.date_of_birth)
			dod = new_string2date(patient.date_of_dead)
			do_lastvisit = new_string2date(patient.date_of_lastvisit)
			recev_date = new_string2date(section.date)
			age = relativedelta(recev_date, dob).years + (relativedelta(recev_date, dob).months+0.)/12
			
			if patient.date_of_dead != "NULL":
				live_time = relativedelta(dod,recev_date).years + (relativedelta(dod,recev_date).months+0.)/12 + (relativedelta(dod,recev_date).days+0.)/365
				death = 1
			else:
				live_time = relativedelta(do_lastvisit,recev_date).years + (relativedelta(do_lastvisit,recev_date).months+0.)/12 + (relativedelta(do_lastvisit,recev_date).days+0.)/365
				death = 0
			# print(age,live_time, death)
			result_dict['age'].append(age)
			if age > 70:
				result_dict['age.type'].append(">70")
			elif age > 50:
				result_dict['age.type'].append("50-70")
			else:
				result_dict['age.type'].append("<=50")
			result_dict['live_time'].append(live_time)
			result_dict['death'].append(death)
			if valid_section.invasive_depth:
				max_invasive_depth = max([float(a) for a in valid_section.invasive_depth])
			else:
				max_invasive_depth =-1
			result_dict['invasion'].append(max_invasive_depth)
			## add T.stage based on invasion depth
			if max_invasive_depth < 0:
				t_stage = "UNKNOWN"
			elif max_invasive_depth <= 1:
				t_stage = "T1"
			elif max_invasive_depth <= 2:
				t_stage = "T2"
			elif max_invasive_depth < 4:
				t_stage = "T3"
			else:
				t_stage = "T4"
			result_dict['T.stage'].append(t_stage)
	print(len(result_dict['age']))
	print(result_dict.keys())
	for k,v in result_dict.items():
		print(k, len(v))
	pd_data = pd.DataFrame.from_dict(result_dict)
	## filter data:
	pd_data = pd_data.drop(pd_data[pd_data[u"tumor-infiltrating lymphocytes"] == 'UNKNOWN'].index)
	pd_data = pd_data.drop(pd_data[pd_data[u"live_time"] < 0].index)
	pd_data = pd_data.drop(pd_data[pd_data[u"vertical growth phase"] != 'PRESENT'].index)
	pd_data.to_csv(csv_outfile, header=True, index=False)
	print("File %s generated."%csv_outfile, "output shape:",pd_data.shape)
	# for a in result_dict[u'microscopic satellites']:
	# 	if a != "PRESENT" and a != "ABSENT" and a!="UNKNOWN":
	# 		print(a)
	

def analysis_all_patient(patient_dict, other_dead, mrn_dead):
	report_num_list = []
	sex_dict = {}
	ethnicity_dict = {}
	race_dict = {}
	for k, v in patient_dict.iteritems():
		report_num_list.append(len(v.pathology_report_list))
		print(v.patient_encode, v.patient_second_encode)
		v.get_dod_etc(mrn_dead, other_dead)
		if v.sex not in sex_dict:
			sex_dict[v.sex] = 1 
		else:
			sex_dict[v.sex] += 1
		if v.ethnicity not in ethnicity_dict:
			ethnicity_dict[v.ethnicity] = 1 
		else:
			ethnicity_dict[v.ethnicity] += 1
		if v.race not in race_dict:
			race_dict[v.race] = 1 
		else:
			race_dict[v.race] += 1
	mean, lower, upper = mean_confidence_interval(report_num_list)
	print("Mean: %s, 95 CI (%s, %s)"%(mean, lower, upper))
	print("Sex:", sex_dict)
	print("Ethnicity: ", ethnicity_dict)
	print("Race:", race_dict)

def random_select_reports(patient_dict, rand_num, out_file):
	import random
	random.seed(2020)
	feature_list = list(all_dict.keys())
	all_list = []
	valid_date = datetime.date(1800,1,1)
	for k, v in patient_dict.items():
		for report in v.pathology_report_list:
			the_list = [report.assertion_number]
			valid_section = None
			for section in report.section_results:
				if len(section.note_pair) > 0:
					if new_string2date(section.date) > valid_date:
						valid_section = section	
			if valid_section:
				valid_til = True
				for feature in feature_list:
					if feature in valid_section.note_pair:
						the_result = valid_section.note_pair[feature].strip(',:;-_. *').upper()
						the_dict = all_dict[feature]
						if feature == "tumor-infiltrating lymphocytes":
							if the_result in TIL_dict:
								the_result = TIL_dict[the_result]
						elif feature == "regression":
							if the_result in regression_dict:
								the_result = regression_dict[the_result]
						if the_result in the_dict:
							the_result = the_dict[the_result]
						if the_result == "N/A":
							the_result = "UNKNOWN"
					else:
						the_result = "UNKNOWN"
					if feature == 'mitotic rate' and the_result == "UNKNOWN":
						the_result = -1
					the_list.append(the_result)
			else:
				the_list += ['UNKNOWN']*len(feature_list)
			if valid_section and len(valid_section.note_pair) > 0 and "tumor-infiltrating lymphocytes" not in valid_section.note_pair:
				need_manual = True 
			else:
				need_manual = False
			the_list.append(need_manual)
			the_list.append(report.original_text)
			all_list.append(the_list)

	random.shuffle(all_list)
	column_list = ['Assertion'] + feature_list + ['Need_manual', 'Report']
	df = pd.DataFrame(all_list[:rand_num],columns=column_list)
	df.to_excel(out_file, encoding='utf-8', index=False)





def random_select(original_file):
	
	patient_dict = reports2patients(original_file)
	print("patient num:",len(patient_dict))
	
	random_select_reports(patient_dict, 200, "random_validate_200.xlsx")


if __name__ == '__main__':	
	

	# disable##other_dead, mrn_dead = load_motality_csv("../Motality_Pathology.csv")
	other_dead, mrn_dead = load_motality_lastvisit_csv("../motality_lastvisit.csv")
	# print(len(mrn_dead))
	all_reports = "../BS-2004.06-2019.12.ttx"
	patient_dict = reports2patients(all_reports)
	print("patient num:",len(patient_dict))
		#all_reports_export_to_excel(all_reports, "initial_results_20200406.csv")
	# analysis_all_patient(patient_dict, other_dead, mrn_dead)
	# exit(0)
	
	# for k, v in patient_dict.iteritems():
	# 	if v.mrn == '12071189':
	# 		for report in v.pathology_report_list:
	# 			print(report.date_of_received)
	# 			for sec in report.section_results:
	# 				print(sec.date, sec.section_text)
	# exit(0)
	#patient_cohort_analysis(patient_dict)
	
	# exit(0)
	# reports_num = 0
	# multi_invasive_num = 0
	# invasive_num = 0 
	# hard= 0
	# for key, value in patient_dict.iteritems():
	# 	for each_report in value.pathology_report_list:
	# 		reports_num += 1 
	# 		if each_report.invasive_depth:
	# 			invasive_num += 1 
	# 		if each_report.invasive_depth == "MULTIPLE":
	# 			multi_invasive_num += 1
	# 			print(each_report.assertion_number)
	# 		if each_report.malignant_depth:
	# 			hard += 1
	# print(reports_num, invasive_num, multi_invasive_num, hard)

	# exit(0)
	# the_list = []
	# for k,v in patient_dict.iteritems():
	# 	the_list.append(v.sex)
	# print(len(the_list))
	# for a in set(the_list):
	# 	print(a, the_list.count(a))

	# exit(0)
	for k,v in patient_dict.items():
		v.get_dod_etc(mrn_dead, other_dead)
	count = 0
	for k,v in patient_dict.items():
		if v.date_of_dead is not "NULL":
			count += 1 
	print("Total #dead patient for sure:", count)
	patient_cohort_export_to_excel(patient_dict, "../melanoma.final.all.csv")
	exit(0)
	# gene_with_patient_analysis(gene_list, patient_dict)
	# exit(0)

	index_string = "Ulceration"#"Regression"#"Tumor-infiltrating lymphocytes"#"Ulceration"#"Regression"#
	index_string = index_string.lower()
	if index_string == "tumor-infiltrating lymphocytes":
		target_results = {"PRESENT_NON-BRISK", "ABSENT", "PRESENT_BRISK"}
	elif index_string == "regression":
		target_results = {"PARTIAL_PRESENT", "PRESENT", "ABSENT"}
	elif index_string == "ulceration":
		target_results = {"PRESENT", "ABSENT"}
	else:
		target_results = {}
	patient_analysis_period_dead(patient_dict, index_string, 120)#, target_results)
	# report_analysis_feature_with_invasion(patient_dict, index_string,  target_results)
	# patient_analysis_overall_dead(patient_dict,index_string, target_results)



	