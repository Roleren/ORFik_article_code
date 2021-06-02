# Update annotation path below to make it work for you

#!/usr/bin/env python
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

from builtins import str, zip, map, range
# -*- coding:UTF-8 -*-

import RiboCode
from collections import OrderedDict, defaultdict

from RiboCode.detectORF import start_check, percentage_format
from RiboCode.orf_finder import orf_finder
from RiboCode.prepare_transcripts import *
from RiboCode.loadconfig import LoadConfig
from RiboCode.parsing_opts import parsing_ribo

from datetime import datetime
import os
import sys

def find_orfs(transcript_dict, annot_dir, START_CODON, ALTERNATIVE_START_CODON_LIST, STOP_CODON_LIST, MIN_AA_LENGTH):
   # Minimal version of ORF sequence finder in Ribo-code, from RiboCode:Ribocode.main
	transcript_seq = GenomeSeq(os.path.join(annot_dir,"transcripts_sequence.fa"))

	orf_results = []
	tid_num = len(transcript_dict)

	for i,tid in enumerate(transcript_dict.keys()):
		if i % 1000 == 0:
			sys.stderr.write(percentage_format(i / tid_num, decimal=0) + " has finished! \r")

		tobj = transcript_dict[tid]

		tseq = transcript_seq.get_seq(tobj.transcript_id)
		orf_dict = orf_finder(tseq,START_CODON,ALTERNATIVE_START_CODON_LIST,STOP_CODON_LIST,MIN_AA_LENGTH)
		orf_results.append(orf_dict)
	return orf_results

def main():
	"""
	Master function to call different functionalities of RiboCode
	"""

	#args = parsing_ribo()

	# read the config file
	#configIn = LoadConfig(args.config_file)
	now = datetime.now()
	current_time = now.strftime("%H:%M:%S")
	print("Current Time =", current_time)
	# Change to your path
	annot_dir = "/media/roler/S/data/Bio_data/references/zebrafish/Ribocode_anot"
	if not os.path.exists(annot_dir):
		sys.stderr.write("Error, the annotation directory not exists, pls run prepare_transcript.py first!\n")
		sys.exit()
	else:
		gene_dict, transcript_dict = load_transcripts_pickle(os.path.join(annot_dir,"transcripts.pickle"))

	print("finished loading")
	now = datetime.now()

	current_time = now.strftime("%H:%M:%S")
	print("Current Time =", current_time)

	#if args.start_codon:
	START_CODON = "ATG"
	ALTERNATIVE_START_CODON_LIST = None

	stops = "TAA,TAG,TGA"
	#if args.stop_codon:
	STOP_CODON_LIST = stops.strip().split(",")

	res = find_orfs(transcript_dict=transcript_dict, annot_dir = annot_dir,
	          START_CODON=START_CODON, ALTERNATIVE_START_CODON_LIST=ALTERNATIVE_START_CODON_LIST,
			  STOP_CODON_LIST=STOP_CODON_LIST, MIN_AA_LENGTH=30)
	print("finished ORFs")
	now = datetime.now()
	current_time = now.strftime("%H:%M:%S")
	print("Current Time =", current_time)
	return res

final = main()
#print(final)
