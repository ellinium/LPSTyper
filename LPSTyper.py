import re
import os
import sys
import csv
from Bio import SeqIO

LPS_types = ['R1', 'R2', 'R3', 'R4', 'K-12']

##Primers
######R1#############
# 842 bp 1520 waaW->waaV
R1_pr = ('TTCTTGCTGGAATGATGTGG', 'TGCGATATGTATGCGGTTC')
R1_1_pr = ('GAACCGCATACATATCGCA', 'CCACATCATTCCAGCAAGAA')

R1_2_pr = ('CTCTTGCTGGAATGATGTGG', 'TGCGATATGTATGCGGTTC')
R1_3_pr = ('GAACCGCATACATATCGCA', 'CCACATCATTCCAGCAAGAG')

R1_4_pr = ('TTCTTGCTGGAATGATGTGG', 'TGCGATATGTATGCGATTC')
R1_5_pr = ('GAATCGCATACATATCGCA', 'CCACATCATTCCAGCAAGAA')

R1_6_pr = ('TTCTTGCTGGAATGATGTGG', 'TGTGATATGTATGCGGTTC')
R1_7_1_pr = ('GAACCGCATACATATCGCA', 'CCACATCATTCCAGCAAGAA')

R1_8_pr = ('TTCTTGCTGGAATGATGTTG', 'TGCGATATGTATGCGGTTC')
R1_9_pr = ('GAACCGCATACATATCGCA', 'CAACATCATTCCAGCAAGAA')

######R2#############
# 609 bp waaZ - waaK
R2_pr = ('TACGTTCATTTTACCGCAGAG', 'CATGGATTTATCAGGTCGCAA')
R2_1_pr = ('TTGCGACCTGATAAATCCATG', 'CTCTGCGGTAAAATGAACGTA')

R2_2_pr = ('TATGTTCATTTTACCGCAGAG', 'CATGGATTTATCAGGTCGCAA')
R2_3_pr = ('TTGCGACCTGATAAATCCATG', 'CTCTGCGGTAAAATGAACATA')

R2_4_pr = ('TACGTTCATTTTACCGCAGAG', 'CATGGATTTATCATGTCGCAA')
R2_5_pr = ('TTGCGACATGATAAATCCATG', 'CTCTGCGGTAAAATGAACGTA')

R2_6_pr = ('TACGTTCATTTTACCGCGGAG', 'CATGGATTTATCAGGTCGCAA')
R2_7_pr = ('TTGCGACCTGATAAATCCATG', 'CTCCGCGGTAAAATGAACGTA')

######R3#############
# 1021 bp waaJ - waaD
R3_pr = ('AATGTCTCATGGGGTATTGATGA', 'GCTGTTGAAACGTGGATGTA')
R3_1_pr = ('TACATCCACGTTTCAACAGC', 'TCATCAATACCCCATGAGACATT')

R3_2_pr = ('AATGTCTCGTGGGGTATTGATGA', 'GCTGTTGAAACGTGGATGTA')
R3_3_pr = ('TACATCCACGTTTCAACAGC', 'TCATCAATACCCCACGAGACATT')

R3_4_pr = ('AATGTCTCATGGGGGATTGATGA', 'GCTGTTGAAACGTGGATGTA')
R3_5_pr = ('TACATCCACGTTTCAACAGC', 'TCATCAATCCCCCATGAGACATT')

######R4#############
# 380 bp waaW-waaX
R4_pr = ('TTCTTGCTGGAATGATGTGG', 'CGCAGAAATGACTGATGGA')
R4_1_pr = ('TCCATCAGTCATTTCTGCG', 'CCACATCATTCCAGCAAGAA')

R4_2_pr = ('TTCTTGCTGGAATGATGTGG', 'CGCAGAAATGACTGAAGGA')
R4_3_pr = ('TCCTTCAGTCATTTCTGCG', 'CCACATCATTCCAGCAAGAA')

######K-12#############
#wwaU-waaY
K12_pr = ('AACTGTAAATCCCGCAGATGC', 'TCAATACGATCTTTCGCTTTACG',)
K12_1_pr = ('CGTAAAGCGAAAGATCGTATTGA', 'GCATCTGCGGGATTTACAGTT')

K12_2_pr = ('CAGCGTAAAGCGAAAGATCG', 'CTTACTCTATGTAAGCATCTGCG')
K12_3_pr = ('CGCAGATGCTTACATAGAGTAAG', 'CGATCTTTCGCTTTACGCTG')

# List of primers to process
R1_primers = [R1_pr, R1_1_pr, R1_2_pr, R1_3_pr, R1_4_pr, R1_5_pr, R1_6_pr, R1_7_1_pr, R1_8_pr, R1_9_pr]
R2_primers = [R2_pr, R2_1_pr, R2_2_pr, R2_3_pr, R2_4_pr, R2_5_pr, R2_6_pr, R2_7_pr]
R3_primers = [R3_pr, R3_1_pr, R3_2_pr, R3_3_pr, R3_4_pr, R3_5_pr]
R4_primers = [R4_pr, R4_1_pr, R4_2_pr, R4_3_pr]
K12_primers = [K12_pr, K12_1_pr, K12_2_pr, K12_3_pr]

# Dictionary with all the primers
prim_dic = {"R1": R1_primers, "R2": R2_primers, "R3": R3_primers,
			"R4": R4_primers, "K-12": K12_primers}

columns_lst = ["File", "Strain", "LPS Type", "Product size"]


# Get a directory parameter
def get_dir():
	dir_name = sys.argv[1]
	if dir_name is None:
		raise Exception('Please provide a directory with fasta file/files')
	return dir_name


# Identify LPS outer core type for a fasta file
def identify_LPS(filename):
	print(filename)
	records = list(SeqIO.parse(filename, "fasta"))

	match = ''
	for rec in records:
		match, match_len = primer_seq(rec.seq.upper())
		if len(match) != 0:
			break

	if len(match) > 0:
		for idx, m in enumerate(match):
			# print(str(m.group()))
			print(str(m))
			print('product size = ' + str(match_len[idx]))
			if m in LPS_types:
				res_row_lst = [filename, records[0].description, str(m), match_len[idx]]
	else:
		res_row_lst = [filename, records[0].description, "Not found"]

	return res_row_lst


def main():
	fasta_dir = get_dir()

	files2process = list()
	valid_extensions = {'.FASTA', '.fasta', '.fna', '.FNA'}
	exclude_keyword = 'cds_from_genomic'

	files2process = [
		os.path.join(dirpath, file)
		for dirpath, _, filenames in os.walk(fasta_dir)
		for file in filenames
		if any(file.endswith(ext) for ext in valid_extensions)
		   and exclude_keyword not in file
	]

	if len(files2process) == 0:
		raise Exception("No files to process are found in " + fasta_dir)
	print("Processing " + str(len(files2process)) + " files")

	results_lst = [identify_LPS(file) for file in files2process]
	# csv header
	results_lst.insert(0, columns_lst)

	not_found_num_lst = [1 if "Not found" in res else 0 for res in results_lst]
	not_found_num = sum(not_found_num_lst)

	print("Number of matches found: " + str(len(files2process) - not_found_num))
	print("Number of files " + str(len(files2process)))

	with open(fasta_dir+'/LPSTyper_results.csv', 'w', newline='') as file:
		writer = csv.writer(file, delimiter=',')
		writer.writerows(results_lst)

	print("The results are saved to LPSTyper_results.csv")


# Search forward and reverse primers in a file or contig
def search_pattern(fwd, rev, sequence):
	pattern = re.compile(r'' + fwd + r'.*' + rev)
	matches = re.search(pattern, sequence)
	return matches


# Run primers in a file or contig
def primer_seq(seq):
	fnl_matches = list()
	fnl_len_matches = list()
	for lps_type in prim_dic:
		primers = prim_dic[lps_type]

		matches = [search_pattern(fwd, rev, str(seq)) for [fwd, rev] in primers]
		final_match = next((item for item in matches if item is not None), None)
		if final_match is not None:
			fnl_matches.append(lps_type)
			# appending PCR product size
			fnl_len_matches.append(final_match.span()[1] - final_match.span()[0])

	return fnl_matches, fnl_len_matches


if __name__ == '__main__':
	main()
