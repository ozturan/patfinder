#!/usr/bin/python
#
# Patfinder
# Mutation Pattern Finder 
# github.com/ozturan/patfinder
#
# Dogancan Ozturan
#

import os
import sys
import subprocess
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import time

# Methods

def get_seqs_ids(file_name):
	fasta_sequences = SeqIO.parse(file_name,'fasta')
	sequences = []
	ids = []
	for fasta in fasta_sequences:
		sequences.append(str(fasta.seq))
		ids.append(str(fasta.id))
	fasta_sequences.close()
	return sequences, ids


def get_positions():
	print "Please give the positions of the pattern. To quit type 'q'."
	print "First pattern positions: "
	position = raw_input('')
	if position == 'q':
		print "Closed."
		exit(1)
	print "Second pattern positions: " 
	position2 = raw_input('')
	if position2 == 'q':
		print "Closed."
		exit(1)
	positions = [position.split(), position2.split()]
	if positions[0] == [] or positions[1] == []:
		get_positions()
	if position > position2:
		print "It is impossible!"
		print "Hint: Second positions must be higher than the first positions"
		get_positions()
	return positions

def get_pattern(p, seq, id_):
	sub_pattern = seq[int(p[0][0])-1:int(p[1][-1])]
	match_list = []
	if len(p[0]) == 2:
		first = int(p[0][1]) - int(p[0][0])
		second = int(p[1][0]) - int(p[0][1])
		for b in range(len(seq)+1-len(sub_pattern)):
			e = b + len(sub_pattern)
			current = seq[b:e]
			for b2 in range(b+1, len(seq)+1-len(sub_pattern)):
				e2 = b2 + len(sub_pattern)
				current2 = seq[b2:e2]
				if current[0] == current2[0] and current[-1] == current2[-1] and current[first] == current2[first] and current[first+second] == current2[first+second]:
					match_list.append([id_, current, current2])
	return match_list

def initial():
	global sequences
	global ids
	seqs_ids = list(get_seqs_ids(file_name))
	sequences = list(seqs_ids[0])
	ids = seqs_ids[1]
	print "|-------------------------------------------|"
	print "|    Patfinder - Mutation Pattern Finder    |"
	print "|       github.com/ozturan/patfinder        |"
	print "|-------------------------------------------|"
	if len(sequences) > 1:
		print "File opened: %s (total %s sequences)" %(f, len(sequences))
	else:
		print "File opened: %s (total %s sequence)" %(f, len(sequences))

def msa():
	flags = ["./clustalo", "-i", "{0}".format(f), "-o", "{0}".format(filename)+".aln", "--force", "-v", "--outfmt=fasta"]
	print "Current command: %s" %(' '.join(map(str, flags)))
	flag_choice = raw_input("Please specify if you use any flags in clustalo: ")
	if flag_choice != '':
		flags.append(flag_choice)
	else:
		pass
	os.environ['f'] = f
	os.environ['filename'] = filename
	subprocess.call(flags)

def blast(sequence):
	print "BLASTP - Non-Redundant database (nr)."
	print "http://blast.ncbi.nlm.nih.gov"
	print "Query sequence: %s" %sequence[0]
	hits = raw_input("Please give the number of hits you wanted (Default is 50): ")
	if hits == None or hits == 0:
		hits = 50
	print "Running..."
	now = time.time()
	blast_handle = NCBIWWW.qblast("blastp", "nr", sequence[0], hitlist_size=hits).read()
	end = time.time()
	elapsed = int(end-now)/60
	elapsed_s = int(int(end-now)/60.0 - int(end-now)/60)*60
	print "It took %s minutes %s seconds." %(elapsed, elapsed_s)
	save_file = open("patfinder_blast.xml", "w")
	save_file.write(blast_handle)
	save_file.close()
	blast_handle.close()
	for record in NCBIXML.parse(open("patfinder_blast.xml")):
		print record.alignments.title

# Routine

try:
	f = sys.argv[1] 
	file_name = open(f,'r')  # open the file
	filename, file_extension = os.path.splitext(f)  # to get the file extension
	filename
except:
	print "Please provide a FASTA or an aligment (.aln) file as an argument!"
	exit(1)


if file_extension == '.fasta':
	initial()
	if len(sequences) > 1:
		print "BLAST search is skipped.."
		choice = raw_input("Please specify if multiple sequence alingment is needed. y/n: ")
		if choice == "n" or choice == "N":
			print "Multiple sequence alingment is skipped.. Pattern search is initiated.."
			positions = get_positions()
			result_ = []
			for seq, id_ in zip(sequences, ids):
				result_.append(get_pattern(positions, seq, id_))
			if result_ == None or len(result_) == 0:
				print "No pattern has been found between the sequences."
			else:
				print "|---------|"
				print "| Results |"
				print "|-------------------------|"
				cleaned_result = filter(None, result_)
				cleaned_result.insert(0, [["| Name", "Seq I", "Seq II |"], 
										  ["|-------------------------|"]])
				file = open("patfinder_results.txt", "w")
				for i in cleaned_result:
					for row in i:
						a = '    '.join(row)
						print '    '.join(row)
						file.write(a+"\n")
				file.close()
				print "Results are written in patfinder_results.txt file. Enjoy!"
		elif choice == "y" or choice == "Y":
			msa()
			print "Multiple sequence alingment is done.. Pattern search is initiated.."
			file_name_msa = open("{0}".format(filename)+".aln",'r') 
			seqs_ids_msa = list(get_seqs_ids(file_name_msa))
			sequences_msa = list(seqs_ids_msa[0])
			ids_msa = seqs_ids_msa[1]
			positions = get_positions()
			sub_pattern = int(int(positions[-1][-1]) - int(positions[0][0]) +1)
			result_ = []
			for seq, id_ in zip(sequences_msa, ids_msa):
				result_.append(get_pattern(positions, seq, id_))
			if result_ == None or len(result_) == 0:
				print "No pattern has been found between the sequences."
			else:
				print "|---------|"
				print "| Results |"
				print "|-------------------------|"
				thing = str('-'*sub_pattern)
				cleaned_result = filter(None, result_)
				for i in cleaned_result:
					if thing in i:
						cleaned_result.remove(thing)
					for j in i:
						if thing in j:
							i.remove(thing)
				print cleaned_result
				file = open("patfinder_results.txt", "w")
				for i in cleaned_result:
					file.write(str(i))
				file.close()
				print "Results are written to patfinder_results.txt ~ Enjoy!"
		else:
			print "Please provide an answer!"
			print "Closed."
			exit(1)
	else:
		print "BLAST is initiated."
		blast_results = blast(sequences)
		sequences_blast = blast_results[0]
		ids_blast = blast_results[1]
		positions = get_positions()
		sub_pattern = int(int(positions[-1][-1]) - int(positions[0][0]) +1)
		result_ = []
		for seq, id_ in zip(sequences_blast, ids_blast):
			result_.append(get_pattern(positions, seq, id_))
		if result_ == None or len(result_) == 0:
			print "No pattern has been found between the sequences."
		else:
			print "|---------|"
			print "| Results |"
			print "|-------------------------|"
			thing = str('-'*sub_pattern)
			cleaned_result = filter(None, result_)
			for i in cleaned_result:
				if thing in i:
					cleaned_result.remove(thing)
				for j in i:
					if thing in j:
						i.remove(thing)
			print cleaned_result
			file = open("patfinder_results.txt", "w")
			for i in cleaned_result:
				file.write(str(i))
			file.close()
			print "Results are written to patfinder_results.txt ~ Enjoy!"

elif file_extension == '.aln':
	initial()
	file_name_aln = open(f,'r')
	seqs_ids_aln = list(get_seqs_ids(file_name_aln))
	sequences_aln = list(seqs_ids_aln[0])
	ids_aln = list(seqs_ids_aln[1])
	positions = get_positions()
	result_ = []
	for seq, id_ in zip(sequences_aln, ids_aln):
		result_.append(get_pattern(positions, seq, id_))
		if result_ == None or len(result_) == 0:
			print "No pattern has been found between the sequences."
		else:
			print "|---------|"
			print "| Results |"
			print "|-------------------------|"
			cleaned_result = filter(None, result_)
			thing = '-'*sub_pattern
			file = open("patfinder_results.txt", "w")
			for i in cleaned_result:
				file.write(str(i))
			file.close()
			print "Results are written in patfinder_results.txt file. Enjoy!"
else:
	print "Please provide a FASTA or an aligment (.aln) file file as an argument!"
	exit(1)
