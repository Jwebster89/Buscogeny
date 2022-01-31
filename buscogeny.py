#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
import numpy as np
import subprocess
from Bio import SeqIO
from alive_progress import alive_bar
from Bio import AlignIO


class Buscogeny():
	def __init__(self,dirs, odb):
		self.dirs=dirs
		self.odb=odb

	def database_parsing(self,odb):
		if os.path.exists(os.path.join(odb,"links_to_ODB10.txt")):
			df=pd.read_csv(os.path.join(odb,"links_to_ODB10.txt"), sep="\t", header=None)
			targets = df[0].tolist()
			return(targets)

	def busco_targs(self,targets,dirs,odb):
		if not os.path.exists("Buscogeny_out/targets"):
			os.mkdir("Buscogeny_out/targets")
		n_isolates=len(os.listdir(dirs))
		n_targs=len(targets)
		with alive_bar(n_targs, title="Generating target multifastas") as bar:
			for target in targets:
				with open(f"Buscogeny_out/targets/{target}_all.faa","w") as out_h:
					for dir in os.listdir(dirs):
						db=os.path.basename(odb)
						if os.path.exists(os.path.join(dirs,dir,"run_"+db,"busco_sequences","single_copy_busco_sequences",target+".faa")):
							for seqrecord in SeqIO.parse(os.path.join(dirs,dir,"run_"+db,"busco_sequences","single_copy_busco_sequences",target+".faa"), "fasta"):
								seqrecord.id=f"{dir}_{target}"
								SeqIO.write(seqrecord, out_h, "fasta")
					bar()

	def target_alignments(self,dirs):
		if not os.path.exists("Buscogeny_out/alignments"):
			os.mkdir("Buscogeny_out/alignments")
		n_targs=len(os.listdir("Buscogeny_out/targets"))
		n_isolates=len(os.listdir(dirs))
		with alive_bar(n_targs, title="Alignment completion status") as bar:
			for file in os.listdir("Buscogeny_out/targets"):
				(base,ext)=os.path.splitext(file)
				n=0
				for seqrecord in SeqIO.parse(os.path.join("Buscogeny_out/targets",file), "fasta"):
					n+=1
					if n==n_isolates:
						with open(f"Buscogeny_out/alignments/{base}.aln","w") as out_h:
							subprocess.run(['mafft','--adjustdirection','--thread','60',os.path.join("Buscogeny_out/targets",file)], stdout=out_h,stderr=subprocess.DEVNULL)
					else:
						pass
				bar()

	def alignment_concatenate(self):
		if not os.path.exists("Buscogeny_out/supermatrix"):
			os.mkdir("Buscogeny_out/supermatrix")
		i=0
		for file in os.listdir("Buscogeny_out/alignments"):
			alignment = AlignIO.read(os.path.join("Buscogeny_out/alignments",file), "fasta")
			alignment.sort()
			for seqrecord in alignment:
				seqrecord.id=seqrecord.id.split(".")[0]
			# print("Alignment of length %i" % alignment.get_alignment_length())
			if i==0:
				cat_algn = alignment
			else:
				cat_algn += alignment
			i += 1
		print("Concatenated:")
		print("Alignment of length %i" % cat_algn.get_alignment_length())
		outfh = open("Buscogeny_out/supermatrix/superalgn.clstl", "w")
		AlignIO.write(cat_algn, outfh, "clustal")
		outfh.close()

	def ML(self):
		if not os.path.exists("Buscogeny_out/iqtree/"):
			os.mkdir("Buscogeny_out/iqtree/")
		# subprocess.run(['raxmlHPC', '-m', 'PROTGAMMAWAG', '-p', '1337', '-s', 'Buscogeny_out/supermatrix/superalgn.clstl', '-#', '1', '-n', 'output'])
		subprocess.run(['iqtree', '-s', 'Buscogeny_out/supermatrix/superalgn.clstl', '-T', '60', '-B', '1000', '--prefix', 'Buscogeny_out/iqtree/iqtree'],stderr=subprocess.DEVNULL)


	def run(self):
		if not os.path.exists("Buscogeny_out"):
			os.mkdir("Buscogeny_out")
		targets=self.database_parsing(self.odb)
		self.busco_targs(targets,self.dirs, self.odb)
		self.target_alignments(self.dirs)
		self.alignment_concatenate()
		self.ML()


def main():
	parser = argparse.ArgumentParser(description="Create Phylogenies from BUSCO output", add_help=False)

	required = parser.add_argument_group('Required Arguments')
	required.add_argument('-i', '--input', type=str, required=True, help="input folder of Busco output directories")
	required.add_argument('-d', '--db', type=str, required=True, help="Location of odb database")


	optional = parser.add_argument_group('Optional Arguments')
	optional.add_argument("-h", "--help", action="help", help="show this help message and exit")

	args=parser.parse_args()
	input=args.input
	db=args.db

	job = Buscogeny(input, db)
	job.run()


if __name__ == '__main__':
	main()
