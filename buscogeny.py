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


from argparse import RawTextHelpFormatter

class Buscogeny():
	def __init__(self,input, prefix, odb,threads, gap_threshold, exclude_threshold, seq_type, rc_filt):
		self.input=input
		self.output=f"{prefix}_Buscogeny_out"
		self.odb=odb
		self.threads=threads
		self.gap_threshold=gap_threshold
		self.exclude_threshold=exclude_threshold
		self.rc_filt=rc_filt
		if seq_type == "prot":
			self.seq_ext="faa"
		else:
			self.seq_ext="fna"

	def BUSCO(self,input,odb, threads):
		if not os.path.exists(f"{self.output}/BUSCO"):
			os.mkdir(f"{self.output}/BUSCO")
		for file in os.listdir(input):
			if not os.path.exists(f"{self.output}/BUSCO/"+os.path.basename(file)+'_busco_out'):
				subprocess.run(['busco', '-i', f"{input}/{file}", '-o', os.path.basename(file)+'_busco_out','-c',threads,'-m', 'genome', '-l', odb, "--offline"])
				shutil.move(os.path.basename(file)+'_busco_out',f"{self.output}/BUSCO/"+os.path.basename(file)+'_busco_out')
			else:
				print(f"{file} has already been processed by BUSCO, skipping..")

	def database_parsing(self,odb):
		if os.path.exists(os.path.join(odb,"scores_cutoff")):
			df=pd.read_csv(os.path.join(odb,"scores_cutoff"), sep="\t", header=None)
			targets = df[0].tolist()
			return(targets)

	def busco_targs(self,targets,odb):
		if not os.path.exists(f"{self.output}/targets"):
			os.mkdir(f"{self.output}/targets")
		n_isolates=len(os.listdir(f"{self.output}/BUSCO/"))
		n_targs=len(targets)
		if odb.endswith('/'):
			path = odb[:-1]
			db=os.path.basename(path)
		else:
			db=os.path.basename(odb)
		with alive_bar(n_targs, title="Generating target multifastas") as bar:
			for target in targets:
				with open(f"{self.output}/targets/{target}_all.{self.seq_ext}","w") as out_h:
					for dir in os.listdir(f"{self.output}/BUSCO/"):
						if os.path.exists(os.path.join(f"{self.output}/BUSCO/",dir,"run_"+db,"busco_sequences","single_copy_busco_sequences",target+f".{self.seq_ext}")):
							for seqrecord in SeqIO.parse(os.path.join(f"{self.output}/BUSCO/",dir,"run_"+db,"busco_sequences","single_copy_busco_sequences",target+f".{self.seq_ext}"), "fasta"):
								seqrecord.id=f"{dir}_{target}"
								SeqIO.write(seqrecord, out_h, "fasta")
					bar()

	def target_alignments(self,input,threads):
		if not os.path.exists(f"{self.output}/alignments"):
			os.mkdir(f"{self.output}/alignments")
		n_targs=len(os.listdir(f"{self.output}/targets"))
		n_isolates=len(os.listdir(input))
		with alive_bar(n_targs, title="Alignment completion status") as bar:
			for file in os.listdir(f"{self.output}/targets"):
				(base,ext)=os.path.splitext(file)
				for seqrecord in SeqIO.parse(os.path.join(f"{self.output}/targets",file), "fasta"):
					alignment_path = os.path.join(f"{self.output}/alignments", f"{base}.aln")
					if not os.path.exists(alignment_path):
						with open(f"{self.output}/alignments/{base}.aln","w") as out_h:
							subprocess.run(['mafft','--adjustdirection','--thread',threads,os.path.join(f"{self.output}/targets",file)], stdout=out_h,stderr=subprocess.DEVNULL)
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
		threshold = len(os.listdir(f"{self.output}/alignments")) * float(self.exclude_threshold)
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

	def alignment_to_xmfa(self):
		print("Generating XMFA from alignments and building concatenated alignment with 1000-gap padding per block")
		os.makedirs(f"{self.output}/supermatrix", exist_ok=True)
		alignment_dir = f"{self.output}/alignments"
		xmfa_path = f"{self.output}/supermatrix/core_aln.xmfa"
		fasta_path = f"{self.output}/supermatrix/concatenated_alignment.fasta"

		aln_files = sorted(f for f in os.listdir(alignment_dir) if f.endswith('.aln'))

		all_isolates = set()
		for aln_file in aln_files:
			alignment = AlignIO.read(os.path.join(alignment_dir, aln_file), "fasta")
			for record in alignment:
				isolate_name = record.id.split("_busco_out")[0]
				all_isolates.add(isolate_name)
		all_isolates = sorted(all_isolates) 

		concatenated = {iso: "" for iso in all_isolates}
		gap_padding = '-' * 1000

		with open(xmfa_path, 'w') as xmfa_out:
			for aln_file in aln_files:
				aln_path = os.path.join(alignment_dir, aln_file)
				alignment = AlignIO.read(aln_path, "fasta")
				block_len = alignment.get_alignment_length()

				block_records = {rec.id.split("_busco_out")[0]: str(rec.seq) for rec in alignment}

				for iso in all_isolates:
					seq = block_records.get(iso, "-" * block_len)
					padded_seq = seq + gap_padding
					xmfa_out.write(f">{iso}\n{padded_seq}\n")
					concatenated[iso] += padded_seq
				xmfa_out.write("=\n")

		with open(fasta_path, 'w') as fasta_out:
			for iso in all_isolates:
				fasta_out.write(f">{iso}\n{concatenated[iso]}\n")

		print(f"XMFA file written to {xmfa_path}")
		print(f"Concatenated alignment written to {fasta_path}")

	
	def alignment_concatenate_bak(self):
		if not os.path.exists(f"{self.output}/supermatrix"):
			os.mkdir(f"{self.output}/supermatrix")

		all_isolates = set()
		for file in os.listdir(f"{self.output}/alignments"):
			alignment = AlignIO.read(os.path.join(f"{self.output}/alignments",file), "fasta")
			for seq in alignment:
				seq.id=seq.id.split("_busco_out")[0]
				all_isolates.add(seq.id)
		
		missing_counts=self.identify_incomplete_isolates(all_isolates,f"{self.output}/alignments")
		to_remove=self.remove_threshold(missing_counts)
		

		print("Adding missing orthologs as gaps in alignments")
		n_to_cat=len(os.listdir(f"{self.output}/alignments"))
		with alive_bar(n_to_cat, title="Concatenation completion status") as bar:
			for file in os.listdir(f"{self.output}/alignments"):
				alignment = AlignIO.read(os.path.join(f"{self.output}/alignments",file), "fasta")
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
						
				if len(to_remove) > 0:
					alignment = MultipleSeqAlignment([seq for seq in alignment if seq.id not in to_remove])
					alignment.sort()

				alignment.sort(key=lambda x: x.id)

				if file == os.listdir(f"{self.output}/alignments")[0]:
					cat_algn = alignment
				else:
					cat_algn += alignment
				bar()

		print("Concatenated!")
		print("Concatenated alignment of length %i" % cat_algn.get_alignment_length())

		outfh = open(f"{self.output}/supermatrix/superaln.fasta", "w")
		AlignIO.write(cat_algn, outfh, "fasta")
		outfh.close()

	def alignment_concatenate(self):
		if not os.path.exists(f"{self.output}/supermatrix"):
			os.mkdir(f"{self.output}/supermatrix")

		all_isolates = set()
		alignment_dir = f"{self.output}/alignments"
		alignment_files = sorted(os.listdir(alignment_dir))

		for file in alignment_files:
			alignment = AlignIO.read(os.path.join(alignment_dir, file), "fasta")
			for seq in alignment:
				seq.id = seq.id.split("_busco_out")[0]
				seq.description = ""
				all_isolates.add(seq.id)
		all_isolates = sorted(all_isolates)

		missing_counts = self.identify_incomplete_isolates(all_isolates, alignment_dir)
		to_remove = self.remove_threshold(missing_counts)

		print("Adding missing orthologs as gaps in alignments")
		with alive_bar(len(alignment_files), title="Concatenation completion status") as bar:
			for idx, file in enumerate(alignment_files):
				alignment = AlignIO.read(os.path.join(alignment_dir, file), "fasta")
				alignment_dict = {}

				aln_len = alignment.get_alignment_length()
				for seq in alignment:
					seq.id = seq.id.split("_busco_out")[0]
					seq.description = ""
					alignment_dict[seq.id] = seq

				seqs_for_block = []
				for iso in all_isolates:
					if iso in alignment_dict:
						seqs_for_block.append(alignment_dict[iso])
					else:
						gap_seq = SeqRecord(Seq("-" * aln_len), id=iso, description="")
						seqs_for_block.append(gap_seq)

				if to_remove:
					seqs_for_block = [s for s in seqs_for_block if s.id not in to_remove]

				seqs_for_block.sort(key=lambda x: x.id)
				complete_alignment = MultipleSeqAlignment(seqs_for_block)

				if idx == 0:
					cat_algn = complete_alignment
				else:
					cat_algn += complete_alignment

				bar()

		print("Concatenated!")
		print("Concatenated alignment of length %i" % cat_algn.get_alignment_length())

		outfh = open(f"{self.output}/supermatrix/superaln.fasta", "w")
		AlignIO.write(cat_algn, outfh, "fasta")
		outfh.close()


	def run_clipkit(self,rc):
		print("Running clipkit")
		if not rc:
			subprocess.run(['clipkit', f"{self.output}/supermatrix/superaln.fasta", '-m', 'gappy', '-g', self.gap_threshold, '-o', f"{self.output}/supermatrix/superaln.degapped.fasta"])
		else:
			print("Doing this step now")
			subprocess.run(['clipkit', f"{self.output}/ClonalFrameML/BUSCO_alignment.rc_masked.mafft.aln", '-m', 'gappy', '-g', self.gap_threshold, '-o', f"{self.output}/ClonalFrameML/BUSCO_alignment.rc_masked.mafft.degapped.aln"])
		print("clipkit finished")

	def run_clonalframeml(self):
		if not os.path.exists(f"{self.output}/ClonalFrameML"):
			os.mkdir(f"{self.output}/ClonalFrameML")
		print("running ClonalFrameML")
		input_xmfa=f"{self.output}/supermatrix/core_aln.xmfa"
		input_tree=f"{self.output}/iqtree/iqtree.rc.initial.treefile"
		output_em=f"{self.output}/ClonalFrameML/BUSCO_alignment.cfml.em"
		command=["ClonalFrameML",input_tree,input_xmfa,output_em,"-em","true","-emsim","100","-xmfa_file","true"]
		if not os.path.exists(f"{output_em}.em.txt"):
			subprocess.run(command)
		else:
			print("Recombination filtering already exists. Skipping..")
	
	def run_mask_RC(self):
		aln_fasta=f"{self.output}/supermatrix/concatenated_alignment.fasta"
		input_em=f"{self.output}/ClonalFrameML/BUSCO_alignment.cfml.em"
		output_aln=f"{self.output}/ClonalFrameML/BUSCO_alignment.rc_masked.aln"
		output_svg=f"{self.output}/ClonalFrameML/BUSCO_alignment.rc_masked.svg"
		command=["maskrc-svg.py", "--aln", aln_fasta, "--out", output_aln, "--svg", output_svg, input_em]
		if not os.path.exists(output_svg):
			subprocess.run(command)
		else:
			print("Recombination masking already complete. Skipping..")
	
	def realign_supermatrix(self):
		input_alignment_file=f"{self.output}/ClonalFrameML/BUSCO_alignment.rc_masked.aln"
		mafft_output_file=f"{self.output}/ClonalFrameML/BUSCO_alignment.rc_masked.mafft.aln"
		mafft_command=["mafft", "--auto", input_alignment_file]
		if not os.path.exists(mafft_output_file):
			with open(mafft_output_file, "w") as out_f:
				result = subprocess.run(mafft_command, stdout=out_f, stderr=subprocess.PIPE, text=True, )
		else:
			print(f"{mafft_output_file} already exists. Skipping..")

	def ML(self,threads,recom,flag):
		if not os.path.exists(f"{self.output}/iqtree/"):
			os.mkdir(f"{self.output}/iqtree/")

		if recom:
			if os.path.exists(f"{self.output}/iqtree/iqtree.rc.initial.treefile") and flag=="final":
				subprocess.run(['iqtree', '-s', f"{self.output}/ClonalFrameML/BUSCO_alignment.rc_masked.mafft.degapped.aln", "--keep-ident","-T", threads,'-B', '1000', '--prefix', f'{self.output}/iqtree/iqtree.rc.final'])
			else:
				subprocess.run(['iqtree', '-s', f"{self.output}/supermatrix/concatenated_alignment.fasta", "--keep-ident","-T", threads, '--prefix', f'{self.output}/iqtree/iqtree.rc.initial'])
		else:
			subprocess.run(['iqtree', '-s', f"{self.output}/supermatrix/superaln.degapped.fasta", "--keep-ident","-T", threads, '-B', '1000', '--prefix', f'{self.output}/iqtree/iqtree'])

	def count_ortho_in_alignments(self):
		ortho_counts = {}
		alignments_directory=f"{self.output}/alignments"
		for alignment_file in os.listdir(alignments_directory):
			if alignment_file.endswith(".aln"): 
				alignment = AlignIO.read(os.path.join(alignments_directory, alignment_file), "fasta")
				for record in alignment:
					genome_id = record.id.split("_busco_out")[0]
					ortho_counts[genome_id] = ortho_counts.get(genome_id, 0) + 1
		
		with open(f"{self.output}/ortholog_counts.tsv", "w") as out_file:
			for genome_id, count in ortho_counts.items():
				out_file.write(f"{genome_id}\t{count}\n")
		return ortho_counts

	def plot_ortho_counts(self,ortho_counts):
		names = list(ortho_counts.keys())
		values = list(ortho_counts.values())

		output_file=f"{self.output}/Genome_ortholog_counts.png"
		fig_size = max(10, len(names)* 0.35)
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

		output_file = f"{self.output}/Genome_ortholog_counts.html"

		fig = go.Figure(data=[go.Bar(x=names, y=values)])
		fig.update_layout(
			title='Ortholog Counts in Each Genome',
			xaxis=dict(title='Genomes'),
			yaxis=dict(title='Ortholog Counts'),
		)


		pio.write_html(fig, file=output_file, auto_open=False)

	def run(self):
		if not os.path.exists(f"{self.output}"):
			os.mkdir(f"{self.output}")
		self.BUSCO(self.input,self.odb,self.threads)
		targets=self.database_parsing(self.odb)
		self.busco_targs(targets, self.odb)
		self.target_alignments(self.input,self.threads)
		ortho_counts=self.count_ortho_in_alignments()
		self.plot_ortho_counts(ortho_counts)
		self.plotly_ortho_counts(ortho_counts)
		if self.rc_filt:
			self.alignment_to_xmfa()
		else:
			self.alignment_concatenate()
			self.run_clipkit(self.rc_filt)
		self.ML(self.threads,self.rc_filt, "initial")
		if self.rc_filt:
			self.run_clonalframeml()
			self.run_mask_RC()
			self.realign_supermatrix()
			self.run_clipkit(self.rc_filt)
			self.ML(self.threads,self.rc_filt, "final")


def main():
	description=f"""
	Create Phylogenies from BUSCO output. \n
	Version: {__version__}"""

	parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter, add_help=False)

	required = parser.add_argument_group('Required Arguments')
	required.add_argument('-i', '--input', type=str, required=True, help="Input folder of Genomes")
	required.add_argument('-d', '--db', type=str, required=True, help="Location of odb database")
	required.add_argument('-o', '--output', type=str, required=True, help="Location of output [prefix]_Buscogeny_out")


	optional = parser.add_argument_group('Optional Arguments')
	optional.add_argument("-h", "--help", action="help", help="show this help message and exit")
	optional.add_argument("-t", "--threads", type=str, required=False, default='8', help="Number of threads to use. Default 8")
	optional.add_argument("-g", "--gappy_threshold", type=str, required=False, default='0.05', help="Specifies gaps threshold used by clipkit. Default 0.05")
	optional.add_argument("-e", "--exclude_threshold", type=str, required=False, default='0.2', help="Specifies proportion of alignments an isolate is allowed to be missing from. Default 0.2")
	optional.add_argument("-s", "--seq_type", type=str, choices=['prot', 'nucl'], default='prot', help="Alignment using: 'prot' for protein, 'nucl' for nucleotide. Default is 'prot'.")
	optional.add_argument("-r", "--rc_filt", action='store_true', help="Enable recombination filtering")

	args=parser.parse_args()
	input=args.input
	prefix = os.path.normpath(args.output)
	db=args.db
	threads=args.threads
	gap_threshold=args.gappy_threshold
	exclude_threshold=args.exclude_threshold
	seq_type=args.seq_type
	rc_filt=args.rc_filt

	if rc_filt and seq_type=="prot":
		print("Warning: Recombination filtering is only available on nucleotide alignments. Please use -s nucl instead. Exiting.")
		sys.exit()
	job = Buscogeny(input, prefix, db, threads, gap_threshold, exclude_threshold, seq_type, rc_filt)
	job.run()


if __name__ == '__main__':
	main()