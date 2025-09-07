```shell
# mafft

# MAFFT is suitable for aligning large numbers of reads and automatically selects appropriate algorithms. 
# For large sequence sets, MAFFT is faster than MUSCLE.

```shell
mamba install -c bioconda mafft
mafft --thread 120 sequences.fasta > aligned.fasta
```

## LTR

# Recommended MAFFT parameters for multiple sequence alignment of LTR or similar sequences:

# --auto: Lets MAFFT automatically select the best strategy based on sequence data size and complexity
# --genafpair and --maxiterate 1000: Used for complex alignment tasks, improves accuracy
# --adjustdirection: Automatically detects reverse complement sequences (useful for LTRs)
# --localpair and --maxiterate 1000: Provides high-precision alignment for smaller datasets
# --ep 0: Reduces gap opening penalty to handle indels more effectively
# --thread number: Enables multi-core processing for faster alignment

```shell
mafft --thread 120 --auto --genafpair --localpair --maxiterate 1000 --adjustdirection sequences.fasta > aligned.fasta
```

# Faster method (completed in 2-3 days, used in Arabidopsis)
```shell
cd-hit-est -i all_chr1_3.fa -o all_chr1_3.cluster.fa

mafft --retree 1 --maxiterate 0 --thread 190 fl_LTR.sample6k.dedup.fa > fl_LTR.sample6k.dedup.aligned.fa
sed -i 's/:/_/g' fl_LTR.sample6k.dedup.aligned.fa
FastTree -nt -gtr -log fl_LTR.sample6k.dedup.aligned.fa.log fl_LTR.sample6k.dedup.aligned.fa > fl_LTR.sample6k.dedup.aligned.fa.nwk &
```

## Tandem repeat

# Recommended MAFFT parameters for tandem repeat sequences with consistent length and potential reverse complements:

# --adjustdirectionaccurately: Enhanced version for handling reverse complement sequences
# --localpair: Uses local alignment algorithm for precise alignment of complex structures
# --maxiterate 1000: Increases iterations to optimize alignment results
# --retree 2: Improves tree reconstruction process for better alignment quality
# --thread number: Enables multi-core processing
# --globalpair: Uses global alignment mode for consistent alignment across sequences
# --op 400: Adjusts gap opening penalty (higher value reduces gap openings)

```shell
mafft --thread 120 --adjustdirectionaccurately --localpair --globalpair --retree 2 --op 400 sequences.fasta > aligned.fasta
```

# EMBOSS cons

# cons generates consensus sequences from multiple sequence alignments

```shell
mamba install -c bioconda emboss

cons -sequence aligned.fasta -outseq consensus.fasta
```

# IQ-TREE

# Builds phylogenetic trees and visualizations

## Installation
```shell
conda install -c bioconda iqtree
```

## Usage
```shell
iqtree -s example.phy -bb 1000 -nt AUTO -m MFP -pre output

# Alternative clustering step (commented out):
# cd-hit-est -i ../WEW.peri20M.flTE.fa -o ../WEW.peri20M.flTE.cluster.fa -c 0.9 -n 5 -T 0 -M 16000

# FastTree alternatives (commented out):
# FastTree -i example.phy > example.phy.nwk
# FastTreeMP example.phy > example.phy.nwk
```

# Parameter explanations:
# -s: Input file
# -bb: Bootstrap value (minimum 1000 recommended)
# -nt: Thread count (AUTO automatically selects optimal number)
# -m MFP: Finds best tree-building model (default in iqtree2 versions >1.5.4)
# -pre: Output file prefix

## Log file notes

# If sequence character composition significantly deviates from the MSA file average, 
# the sequence is marked as "failed" in the log file
# This test helps identify problematic sequences in the dataset
# Failed sequences are not automatically removed but can help explain unexpected tree topologies
```
