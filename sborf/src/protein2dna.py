#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  protein2dna.py
#  
#  Copyright 2019 Gabriele Orlando <orlando.gabriele89@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
import os,pickle,multiprocessing,string,random
def leggifasta(database): #legge un file fasta e lo converte in un dizionario
		f=open(database)
		uniprot=f.readlines()
		f.close()
		dizio={}
		for i in uniprot:
			#c=re.match(">(\w+)",i)  4
		
			if i[0]=='>':
					if '|' in i:
						uniprotid=i.strip('>\n').split('|')[1]
					else:
						uniprotid=i.strip('>\n').split(' ')[0]
					dizio[uniprotid]=''
			else:
				dizio[uniprotid]=dizio[uniprotid]+i.strip('\n')
		return dizio
		
def leggifasta_dna(database): #legge un file fasta e lo converte in un dizionario
		f=open(database)
		uniprot=f.readlines()
		f.close()
		dizio={}
		for i in uniprot:
			#c=re.match(">(\w+)",i)  4
		
			if i[0]=='>':

					uniprotid=i.strip('>\n').split(' ')[0]
					dizio[uniprotid]=''
			else:
				dizio[uniprotid]=dizio[uniprotid]+i.strip('\n')
		return dizio


def run_blast_protein2dna(inp):
	tmp_folder = "blast_tmp/"
	if not os.path.isdir(tmp_folder):
		os.system("mkdir "+tmp_folder)
	seq,blastDB = inp
	name = tmp_folder+"/"+''.join(random.choice(string.ascii_uppercase) for i in range(8))

	f = open(name+".fasta", "w")
	f.write(">" + "tmp" + "\n" + seq)
	f.close()

	os.system(
		"tblastn -db " + blastDB + " -query "+name+".fasta -evalue 1e-80 -outfmt \'6 sseqid sseq\' | awk \'BEGIN{FS=\"\\t\"; OFS=\"\\n\"}{gsub(/-/, \"\", $2); print \">\"$1}\' > "+name+".blast")

	# os.system("tblastn -db blast_DNA_DB/mouse/GCF_000001635.27_top_level -query tmp.fasta -evalue 1e-50 -outfmt \'6 sseqid sseq\' | awk \'BEGIN{FS=\"\\t\"; OFS=\"\\n\"}{gsub(/-/, \"\", $2); print \">\"$1}\' > tmp.blast")
	m = open(name+".blast", "r").readlines()
	if len(m) > 0:
		res = open(name+".blast", "r").readlines()[0].strip().strip(">")
	else:
		print("salto")
		res = False
	os.system("rm " + name + ".fasta")
	os.system("rm "+name+".blast")
	#os.system("rmdir "+name)

	print(name,res)
	return res


def protein2dna(seqs,blastDB="dataset/blast_databases/fungi/cds_from_genomic.fna",write_dna_file=None,cpus=multiprocessing.cpu_count()):  # dna file from https://www.ncbi.nlm.nih.gov/datasets/taxonomy/4932/
	# to build the blast DB, use the .pl file from the blast+ package. for coli, do makeblastdb -in dnaColiSequences.fna -parse_seqids -blastdb_version 5 -title "coli_db" -dbtype nucl

	dna = leggifasta_dna(blastDB)
	mapping = {}
	todo = list(seqs.keys())[:]
	#cpus = 1
	if cpus > 1 :

		inps = [(seqs[i],blastDB) for i in todo]
		with multiprocessing.Pool() as pool:
			res = pool.map(run_blast_protein2dna, inps)
		for k in range(len(res)):
			if res[k]:
				mapping[todo[k]] = res[k]
	else:
		for i in list(seqs.keys())[:]:
				f = open("tmp.fasta", "w")
				f.write(">" + i + "\n" + seqs[i])
				f.close()

				os.system("tblastn -db " + blastDB + " -query tmp.fasta -evalue 1e-80 -outfmt \'6 sseqid sseq\' | awk \'BEGIN{FS=\"\\t\"; OFS=\"\\n\"}{gsub(/-/, \"\", $2); print \">\"$1}\' > tmp.blast")

				# os.system("tblastn -db blast_DNA_DB/mouse/GCF_000001635.27_top_level -query tmp.fasta -evalue 1e-50 -outfmt \'6 sseqid sseq\' | awk \'BEGIN{FS=\"\\t\"; OFS=\"\\n\"}{gsub(/-/, \"\", $2); print \">\"$1}\' > tmp.blast")
				m = open("tmp.blast", "r").readlines()
				if len(m) > 0:
					mapping[i] = open("tmp.blast", "r").readlines()[0].strip().strip(">")
				else:
					print("salto", i)
					continue
				print(i, mapping[i], "asd")
	# asd
	final={}
	for i in mapping.keys():
		if mapping[i] in dna:
			final[i] = dna[mapping[i]]
		elif 'lcl|'+mapping[i] in dna: # some data has a lcl| at the beginning. God only knows why
			final[i] = dna['lcl|'+mapping[i]]
		else:
			print("missing seq in the uniprot mapping")

	if write_dna_file is not None:
		f= open(write_dna_file,"w")
		for k in final.keys():
			f.write(">"+k+"\n"+final[k]+"\n")
		f.close()

	return final

def get_dna_sequences(paxDB_folder="dataset/sacaromices/",blastDB="blast_databases/fungi/cds_from_genomic.fna",sequence_folder="dna_sequences/"):
	seqs = get_sequences(paxDB_folder)
	print("got", len(seqs),"sequences")
	dna_seqs = protein2dna(seqs,blastDB=blastDB)
	#os.system("mkdir -p "+sequence_folder)
	f=open(sequence_folder,"w")
	for i in dna_seqs.keys():
		f.write(">"+i+"\n"+dna_seqs[i]+"\n")
	f.close()
	return dna_seqs

def get_sequences(fol="dataset/sacaromices/"):
	paxdb_seqs = readPaxDBsequences() # protein sequences
	final_list = {}

	for fil in os.listdir(fol):
		for l in open(fol + fil, "r").readlines():
			if l[0] == "#":
				continue
			if l.strip() == '':
				continue
			a = l.split()
			paxID = a[0]
			if not paxID in paxdb_seqs:
				continue
			final_list[paxID] = paxdb_seqs[paxID]
	return 	final_list

def readPaxDBsequences(folder = "dataset/paxdb-protein-sequences-v5.0/"):
	fin={}
	for k in os.listdir(folder):
		a = leggifasta(folder+k)
		fin.update(a)
	return fin
