#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
import numpy as np
import subprocess
import shutil
from version import __version__
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from alive_progress import alive_bar
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.io as pio

class Buscogeny():
	def __init__(self,input, odb,threads, gap_threshold,exclude_threshold, seq_type):
		self.input=input
		self.odb=odb
		self.threads=threads
		self.gap_threshold=gap_threshold
		self.exclude_threshold=exclude_threshold
		if seq_type == "prot":
			self.seq_ext="faa"
		else:
			self.seq_ext="fna"

	def BUSCO(self,input,odb, threads):
		if not os.path.exists("Buscogeny_out/BUSCO"):
			os.mkdir("Buscogeny_out/BUSCO")
		for file in os.listdir(input):
			if not os.path.exists("Buscogeny_out/BUSCO/"+os.path.basename(file)+'_busco_out'):
				subprocess.run(['busco', '-i', input+file, '-o', os.path.basename(file)+'_busco_out','-c',threads,'-m', 'genome', '-l', odb])
				shutil.move(os.path.basename(file)+'_busco_out',"Buscogeny_out/BUSCO/"+os.path.basename(file)+'_busco_out')
			else:
				print(f"{file} has already been processed by BUSCO, skipping..")

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
		if odb.endswith('/'):
			path = odb[:-1]
			db=os.path.basename(path)
		else:
			db=os.path.basename(odb)
		with alive_bar(n_targs, title="Generating target multifastas") as bar:
			for target in targets:
				with open(f"Buscogeny_out/targets/{target}_all.{self.seq_ext}","w") as out_h:
					for dir in os.listdir("Buscogeny_out/BUSCO/"):
						# db=os.path.basename(odb)
						if os.path.exists(os.path.join("Buscogeny_out/BUSCO/",dir,"run_"+db,"busco_sequences","single_copy_busco_sequences",target+f".{self.seq_ext}")):
							for seqrecord in SeqIO.parse(os.path.join("Buscogeny_out/BUSCO/",dir,"run_"+db,"busco_sequences","single_copy_busco_sequences",target+f".{self.seq_ext}"), "fasta"):
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
				for seqrecord in SeqIO.parse(os.path.join("Buscogeny_out/targets",file), "fasta"):
					alignment_path = os.path.join("Buscogeny_out/alignments", f"{base}.aln")
					if not os.path.exists(alignment_path):
						with open(f"Buscogeny_out/alignments/{base}.aln","w") as out_h:
							subprocess.run(['mafft','--adjustdirection','--thread',threads,os.path.join("Buscogeny_out/targets",file)], stdout=out_h,stderr=subprocess.DEVNULL)
						alignment = AlignIO.read(alignment_path, "fasta")
						for record in alignment:
							if record.id.startswith("_R_"):
								record.id = record.id.replace("_R_", "")
							record.description = ""
						AlignIO.write(alignment, alignment_path, "fasta")
				bar()

	def identify_incomplete_isolates(self, all_isolates, alignments_directory):
		print("Identifying missing isolates")
		missing_counts = {iso: 0 for iso in all_isolates}

		for file in os.listdir(alignments_directory):
			alignment = AlignIO.read(os.path.join(alignments_directory, file), "fasta")
			present_isolates = {record.id.split("_busco_out")[0] for record in alignment} 

			for iso in all_isolates:
				if iso not in present_isolates:
					missing_counts[iso] += 1

		return missing_counts

	def remove_threshold(self, missing_counts):
		threshold = len(os.listdir("Buscogeny_out/alignments")) * float(self.exclude_threshold)
		to_remove = set()
		for iso, count in missing_counts.items():
			if count >= threshold:
				to_remove.add(iso)
		if len(to_remove) > 0:
			print(f"WARNING: Removing {len(to_remove)} isolate(s) from alignments")
			print(', '.join(to_remove))  
		else:
			print(f"Removing 0 isolates from alignments")
	
		return to_remove

	def alignment_concatenate(self):
		if not os.path.exists("Buscogeny_out/supermatrix"):
			os.mkdir("Buscogeny_out/supermatrix")

		all_isolates = set()
		for file in os.listdir("Buscogeny_out/alignments"):
			alignment = AlignIO.read(os.path.join("Buscogeny_out/alignments",file), "fasta")
			for seq in alignment:
				seq.id=seq.id.split("_busco_out")[0]
				all_isolates.add(seq.id)
		
		missing_counts=self.identify_incomplete_isolates(all_isolates,"Buscogeny_out/alignments")
		# print(missing_counts)
		to_remove=self.remove_threshold(missing_counts)
		

		print("Adding missing orthologs as gaps in alignments")
		n_to_cat=len(os.listdir("Buscogeny_out/alignments"))
		with alive_bar(n_to_cat, title="Concatenation completion status") as bar:
			for file in os.listdir("Buscogeny_out/alignments"):
				alignment = AlignIO.read(os.path.join("Buscogeny_out/alignments",file), "fasta")
				for iso in all_isolates:
					iso_present = False
					for seqrecord in alignment:
						seqrecord.id=seqrecord.id.split("_busco_out")[0]
						if iso in seqrecord.id:
							iso_present = True
							break
					if not iso_present:
						gap_seq = SeqRecord(Seq("-" * alignment.get_alignment_length()), id=iso, description="")
						alignment.append(gap_seq)
						
				# Remove isolates that are missing in more than X% of the alignments
				if len(to_remove) > 0:
					alignment = MultipleSeqAlignment([seq for seq in alignment if seq.id not in to_remove])
					alignment.sort()

				alignment.sort(key=lambda x: x.id)

				if file == os.listdir("Buscogeny_out/alignments")[0]:
					cat_algn = alignment
				else:
					cat_algn += alignment
				bar()

		print("Concatenated!")
		print("Concatenated alignment of length %i" % cat_algn.get_alignment_length())

		outfh = open("Buscogeny_out/supermatrix/superaln.fasta", "w")
		AlignIO.write(cat_algn, outfh, "fasta")
		outfh.close()

	def run_clipkit(self):
		print("Running clipkit")
		subprocess.run(['clipkit', "Buscogeny_out/supermatrix/superaln.fasta", '-m', 'gappy', '-g', self.gap_threshold, '-o', "Buscogeny_out/supermatrix/superaln.degapped.fasta"])
		print("clipkit finished")

	def ML(self,threads):
		if not os.path.exists("Buscogeny_out/iqtree/"):
			os.mkdir("Buscogeny_out/iqtree/")
		subprocess.run(['iqtree', '-s', "Buscogeny_out/supermatrix/superaln.degapped.fasta", "--keep-ident","-T", threads, '-B', '1000', '--prefix', 'Buscogeny_out/iqtree/iqtree'])

	def count_ortho_in_alignments(self):
		ortho_counts = {}
		alignments_directory="Buscogeny_out/alignments"
		for alignment_file in os.listdir(alignments_directory):
			if alignment_file.endswith(".aln"): 
				alignment = AlignIO.read(os.path.join(alignments_directory, alignment_file), "fasta")
				for record in alignment:
					genome_id = record.id.split("_busco_out")[0]
					ortho_counts[genome_id] = ortho_counts.get(genome_id, 0) + 1

		return ortho_counts

	def plot_ortho_counts(self,ortho_counts):
		names = list(ortho_counts.keys())
		values = list(ortho_counts.values())
		# print(names)

		output_file="Buscogeny_out/Genome_ortholog_counts.png"
		fig_size = max(10, len(names)* 0.35)  # Adjust the multiplier as needed
		plt.figure(figsize=(fig_size, 6))

		plt.bar(range(len(ortho_counts)), values, tick_label=names)
		plt.xlabel('Genomes')
		plt.ylabel('Ortholog Counts')
		plt.title('Ortholog Counts in Each Genome')
		plt.xticks(rotation=90)
		plt.savefig(output_file, format='png',bbox_inches="tight")
		plt.close()

	def plotly_ortho_counts(self, ortho_counts):
		names = list(ortho_counts.keys())
		values = list(ortho_counts.values())

		output_file = "Buscogeny_out/Genome_ortholog_counts.html"

		# Create a Plotly bar chart
		fig = go.Figure(data=[go.Bar(x=names, y=values)])
		fig.update_layout(
			title='Ortholog Counts in Each Genome',
			xaxis=dict(title='Genomes'),
			yaxis=dict(title='Ortholog Counts'),
		)

		# Save the plot as an interactive HTML file
		pio.write_html(fig, file=output_file, auto_open=False)

	def run(self):
		if not os.path.exists("Buscogeny_out"):
			os.mkdir("Buscogeny_out")
		self.BUSCO(self.input,self.odb,self.threads)
		targets=self.database_parsing(self.odb)
		self.busco_targs(targets, self.odb)
		self.target_alignments(self.input,self.threads)
		ortho_counts=self.count_ortho_in_alignments()
		self.plot_ortho_counts(ortho_counts)
		self.plotly_ortho_counts(ortho_counts)
		self.alignment_concatenate()
		self.run_clipkit()
		self.ML(self.threads)


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
	optional.add_argument("-t", "--threads", type=str, required=False, default='8', help="Number of threads to use. Default 8")
	optional.add_argument("-g", "--gappy_threshold", type=str, required=False, default='0.1', help="Specifies gaps threshold used by clipkit. Default 0.1")
	optional.add_argument("-e", "--exclude_threshold", type=str, required=False, default='0.2', help="Specifies proportion of alignments an isolate is allowed to be missing from.")
	optional.add_argument("-s", "--seq_type", type=str, choices=['prot', 'nucl'], default='prot', help="Alignment using: 'prot' for protein, 'nucl' for nucleotide. Default is 'prot'.")

	args=parser.parse_args()
	input=args.input
	db=args.db
	threads=args.threads
	gap_threshold=args.gappy_threshold
	exclude_threshold=args.exclude_threshold
	seq_type=args.seq_type

	job = Buscogeny(input, db, threads, gap_threshold, exclude_threshold, seq_type)
	job.run()


if __name__ == '__main__':
	main()