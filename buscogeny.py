#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
import numpy as np
import subprocess
import shutil
from version import __version__
from Bio import SeqIO
from alive_progress import alive_bar
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from argparse import RawTextHelpFormatter

class Buscogeny():
	def __init__(self,input, odb,threads):
		self.input=input
		self.odb=odb
		self.threads=threads

	def BUSCO(self,input,odb, threads):
		if not os.path.exists("Buscogeny_out/BUSCO"):
			os.mkdir("Buscogeny_out/BUSCO")
		for file in os.listdir(input):
			subprocess.run(['busco', '-i', input+file, '-o', os.path.basename(file)+'_busco_out','-c',threads,'-m', 'genome', '-l', odb])
			shutil.move(os.path.basename(file)+'_busco_out',"Buscogeny_out/BUSCO/"+os.path.basename(file)+'_busco_out')

	def database_parsing(self,odb):
		if os.path.exists(os.path.join(odb,"links_to_ODB10.txt")):
			df=pd.read_csv(os.path.join(odb,"links_to_ODB10.txt"), sep="\t", header=None)
			targets = df[0].tolist()
			return(targets)

	def busco_targs(self,targets,odb):
		if not os.path.exists("Buscogeny_out/targets"):
			os.mkdir("Buscogeny_out/targets")
		n_isolates=len(os.listdir("Buscogeny_out/BUSCO/"))
		n_targs=len(targets)
		with alive_bar(n_targs, title="Generating target multifastas") as bar:
			for target in targets:
				with open(f"Buscogeny_out/targets/{target}_all.faa","w") as out_h:
					for dir in os.listdir("Buscogeny_out/BUSCO/"):
						db=os.path.basename(odb)
						if os.path.exists(os.path.join("Buscogeny_out/BUSCO/",dir,"run_"+db,"busco_sequences","single_copy_busco_sequences",target+".faa")):
							for seqrecord in SeqIO.parse(os.path.join("Buscogeny_out/BUSCO/",dir,"run_"+db,"busco_sequences","single_copy_busco_sequences",target+".faa"), "fasta"):
								seqrecord.id=f"{dir}_{target}"
								SeqIO.write(seqrecord, out_h, "fasta")
					bar()

	def target_alignments(self,input,threads):
		if not os.path.exists("Buscogeny_out/alignments"):
			os.mkdir("Buscogeny_out/alignments")
		n_targs=len(os.listdir("Buscogeny_out/targets"))
		n_isolates=len(os.listdir(input))
		with alive_bar(n_targs, title="Alignment completion status") as bar:
			for file in os.listdir("Buscogeny_out/targets"):
				(base,ext)=os.path.splitext(file)
				n=0
				for seqrecord in SeqIO.parse(os.path.join("Buscogeny_out/targets",file), "fasta"):
					n+=1
					if n==n_isolates:
						with open(f"Buscogeny_out/alignments/{base}.aln","w") as out_h:
							subprocess.run(['mafft','--adjustdirection','--thread',threads,os.path.join("Buscogeny_out/targets",file)], stdout=out_h,stderr=subprocess.DEVNULL)
					else:
						pass
				bar()

	def alignment_concatenate(self):
		if not os.path.exists("Buscogeny_out/supermatrix"):
			os.mkdir("Buscogeny_out/supermatrix")
		if not os.path.exists("Buscogeny_out/trimmed_alignments"):
			os.mkdir("Buscogeny_out/trimmed_alignments")
		i=0
		matrix=[]
		n_alns=len(os.listdir("Buscogeny_out/alignments"))
		for file in os.listdir("Buscogeny_out/alignments"):
			alignment = AlignIO.read(os.path.join("Buscogeny_out/alignments",file), "fasta")
			alignment.sort()
			for seqrecord in alignment:
				seqrecord.id=seqrecord.id.split(".")[0]
			if i==0:
				cat_algn = alignment
			else:
				cat_algn += alignment
			i += 1
		print("Concatenated:")
		print("Alignment of length %i" % cat_algn.get_alignment_length())
		outfh = open("Buscogeny_out/supermatrix/superaln.clstl", "w")
		AlignIO.write(cat_algn, outfh, "clustal")
		outfh.close()
		subprocess.run(['clipkit',"Buscogeny_out/supermatrix/superaln.clstl",'-m','gappy','-o',"Buscogeny_out/supermatrix/superaln.degapped.clstl"])

	def ML(self):
		if not os.path.exists("Buscogeny_out/iqtree/"):
			os.mkdir("Buscogeny_out/iqtree/")
		subprocess.run(['iqtree', '-s', "Buscogeny_out/supermatrix/superaln.degapped.clstl", "--keep-ident",'-T', '60', '-B', '1000', '--prefix', 'Buscogeny_out/iqtree/iqtree'])


	def run(self):
		if not os.path.exists("Buscogeny_out"):
			os.mkdir("Buscogeny_out")
		self.BUSCO(self.input,self.odb,self.threads)
		targets=self.database_parsing(self.odb)
		self.busco_targs(targets, self.odb)
		self.target_alignments(self.input,self.threads)
		self.alignment_concatenate()
		self.ML()


def main():
	description=f"""
	Create Phylogenies from BUSCO output. \n
	Version: {__version__}"""

	parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter, add_help=False)

	required = parser.add_argument_group('Required Arguments')
	required.add_argument('-i', '--input', type=str, required=True, help="Input folder of Genomes")
	required.add_argument('-d', '--db', type=str, required=True, help="Location of odb database")


	optional = parser.add_argument_group('Optional Arguments')
	optional.add_argument("-h", "--help", action="help", help="show this help message and exit")
	optional.add_argument("-t", "--threads", type=str, required=False, default=8, help="Number of threads to use. Default 8")

	args=parser.parse_args()
	input=args.input
	db=args.db
	threads=args.threads

	job = Buscogeny(input, db,threads)
	job.run()


if __name__ == '__main__':
	main()
