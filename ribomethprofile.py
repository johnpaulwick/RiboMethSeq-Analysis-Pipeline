#!/usr/bin/env python3

# Author: John Paul Wick
# 16 September 2024


# Modules required for pipeline.py:
import argparse
import os
import subprocess
import pysam
import matplotlib.pyplot as plt


# Parse commandline argument (the fastq of ribomethseq reads):
parser = argparse.ArgumentParser()
parser.add_argument('fastq', type=str)
args = parser.parse_args()


# Defining filename variables for alignment pipeline:
sample = args.fastq # assign command line argument of input fastq to variable
samplename = sample[:-6] # strip fastq extension to get sample name to use in creating new file types
samplesam = str(f"{samplename}.sam")
samplebam = str(f"{samplename}.bam")
samplesortedbam = str(f"{samplename}.sorted.bam")


# Make directory to contain the temporary files
try:
    os.mkdir('temp')
except:
    True



############################  Alignment Pipeline ###################################


# Copy the files into temp (the human 5.8S RNA reference and the ribomethseq fastq must be present in the folder for the script to work):
subprocess.run(['cp', 'reference.fasta', 'temp/reference.fasta'])
subprocess.run(['cp', str(sample), f'temp/{sample}'])


# Index the reference gene inside temp:
subprocess.run(['bowtie2-build', 'temp/reference.fasta', 'temp/reference'])


# Do the alignment with the input fastq:
subprocess.run(['bowtie2', '-x', 'temp/reference', '-U', str(sample), '-S', f"temp/{str(samplesam)}"])


# Convert sam to bam:
os.system(f"samtools view -bS temp/{samplesam} > temp/{samplebam}")


# Sort bam:
os.system(f"samtools sort -m 100M -o temp/{samplesortedbam} temp/{samplebam}")


# Index sorted bam:
os.system(f"samtools index temp/{samplesortedbam}")



######################### Alignment Analysis with Pysam ############################


# Creating the pysam alignment file object:
samfile = pysam.AlignmentFile(f"temp/{str(samplesortedbam)}","rb")


# Obtaining the length of the reference sequence and creating an appropriately sized list of positions to iterate through:
with open("reference.fasta", 'r') as f:
    lines = f.readlines()
referencelength = len(lines[1].removesuffix('\n'))+1


################# First, I collect positions for 5' starts of the reads only ###########################


# This dictionary will store the counts at each position at the values, with the keys being each position
positions = list(range(0,referencelength-1)) # use 0 - length-1 numbering for starts (shift to left by one base)
countdictstarts = dict.fromkeys(positions, 0)


# Count fragment starts for each position:
for read in samfile.fetch():
    countdictstarts[read.reference_start] += 1


# Saving an image of the counts:
plt.figure(figsize=(13, 4))
plt.bar(list(countdictstarts.keys()), list(countdictstarts.values()), color="b")
plt.title(f"Counts of Start Positions - {samplename} Reads")
plt.xlabel("Position in Human 5.8S RNA Sequence")
plt.ylabel("Number of Reads Starting on the Position")
plt.savefig(f"{samplename}_starts.png")
plt.clf()


################# Next, I collect positions for 3' ends of the reads only ###########################


# This dictionary will store the counts at each position at the values, with the keys being each position
positions = list(range(1,referencelength)) # use 1 - length numbering for ends (no adjustment)
countdictends = dict.fromkeys(positions, 0)


# Count fragment ends for each position:
for read in samfile.fetch():
    countdictends[read.reference_end] += 1


# Saving an image of the counts:
plt.figure(figsize=(13, 4))
plt.bar(list(countdictends.keys()), list(countdictends.values()), color="r")
plt.title(f"Counts of End Positions - {samplename} Reads")
plt.xlabel("Position in Human 5.8S RNA Sequence")
plt.ylabel("Number of Reads Ending on the Position")
plt.savefig(f"{samplename}_ends.png")
plt.clf()


################# Finally, I combine total counts for 5' starts and 3' ends of the reads at each position ###########################


# This dictionary will store the counts at each position at the values, with the keys being each position
positions = list(range(0,referencelength)) # use 0 - length numbering (1 longer than the actual length this time) to include added counts of both ends
countdicttotal = dict.fromkeys(positions, 0)


# Add counts from the starts:
for key, value in countdictstarts.items():
    countdicttotal[key] += value


# Add counts from the ends:
for key, value in countdictends.items():
    countdicttotal[key] += value


# Saving an image of the counts:
plt.figure(figsize=(13, 4))
plt.bar(list(countdicttotal.keys()), list(countdicttotal.values()), color="g")
plt.title(f"Counts of Start and End Positions - {samplename} Reads")
plt.xlabel("Position in Human 5.8S RNA Sequence")
plt.ylabel("Number of Reads Starting or Ending on the Position")
plt.savefig(f"{samplename}_totals.png")
plt.clf()


######################### RELATIVE METHYLATION SCORES ##################################


# Getting total count of reads:
reads = 0
for read in samfile.fetch():
    reads += 1


# Calculating 5'/3' end representation fractions for each position:
scores = []
positions = list(range(1,referencelength)) # excluding the back-shifted position 0
for position in positions:
    scores.append(countdicttotal[position]/reads)


# Normalizing and processing the scores:
scores = [(score - min(scores)) / (max(scores) - min(scores)) for score in scores] # max min normalization on the scores to make them more intelligible
scores = [1-score for score in scores] # inverting it so that higher-methylated positions score higher


# Adding the scores to a new dictionary by position:
scoresdict = dict.fromkeys(positions, 0)
n = 0
for position in positions:
    scoresdict[position] = scores[n]
    n += 1


# Saving an image of the relative methylation scores:
plt.figure(figsize=(13, 4))
plt.bar(list(scoresdict.keys()), list(scoresdict.values()), color="y")
plt.title(f"Normalized Relative Methylation Level Scores - {samplename} Reads")
plt.xlabel("Position in Human 5.8S RNA Sequence")
plt.ylabel("Score")
plt.savefig(f"{samplename}_scores.png")
plt.clf()


#########################################################################################################################


# Getting rid of temp directory as it is no longer needed:
removetempdirectory = "rm -r temp"
os.system(removetempdirectory)


# Completion message:
print("\n\n\n\n\n################################\n\nRIBOMETHSEQ ANALYSIS COMPLETE\n\n################################\n\n")


#########################################################################################################################