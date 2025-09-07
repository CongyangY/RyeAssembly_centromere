```shell
# Quality control
fastqc -t 10 -o qc_output *.fastq.gz
# -t: specify thread count
# -o: specify output directory

# Adapter trimming
######################################### Recommended parameters by Yang's group
cutadapt -j 24 -a AGATCGGAAGAGCACACGTCT -A AGATCGGAAGAGCGTCGTGTAG -o test_1_trim1.fq.gz -p test_2_trim1.fq.gz test_1.fastq.gz test_2.fastq.gz
# -j: thread count
# -a: adapter sequence to remove from R1; -A: adapter sequence to remove from R2 (3' adapter)
# According to the schematic, the first adapter is added to the 3' end of reads. For R2, the adapter on the left side (second sequencing) is removed.
# -o: output filtered R1 reads; -p: output filtered R2 reads

######################################### Recommended by Sun's group for ssDRIP-seq: use trimmomatic/fastqmcf for adapter removal
# Recommended to use trimmomatic, which can remove low-quality reads (or use after cutadapt to remove adapter and tail)
java -jar trimmomatic-0.36.jar PE -threads 12 -phred33 test_1_trim2.fastq.gz test_2_trim2.fastq.gz test_1_trim3.fastq.gz test_1_U.fastq.gz test_2_trim3.fastq.gz test_2_U.fastq.gz ILLUMINACLIP:adapter:2:30:10 HEADCROP:4 LEADING:1 TRAILING:1 SLIDINGWINDOW:4:1 MINLEN:50
# -threads: thread count
# -phred33: use phred33 format for quality scores
# Output files: paired R1, unpaired R1; paired R2, unpaired R2. Can use -baseout to specify prefix
# ILLUMINACLIP: start adapter removal (skip if already removed)
# adapter: file containing adapter sequences; 2: maximum allowed mismatches
# 30: minimum alignment score for palindrome clip mode in PE to trigger adapter removal
# 10: minimum alignment score for adapter removal, usually 7-15
# HEADCROP: remove bases from beginning of reads (not recommended here)
# LEADING: remove bases from start if quality below threshold
# TRAILING: remove bases from end if quality below threshold
# SLIDINGWINDOW: use sliding window, discard entire read if average quality in window drops below threshold
# MINLEN: minimum read length after trimming (not set here)

# Remove tail
######################################### Recommended parameters by Yang's group
cutadapt -j 24 -u 10 -u -10 -U 10 -U -10 -o test_1_trim1.fq.gz -p test_2_trim1.fq.gz test_1_trim2.fq.gz test_2_trim2.fq.gz
# -j: thread count
# -u: remove bases from beginning of R1; -U: remove bases from beginning of R2
# -o: output filtered R1 reads; -p: output filtered R2 reads

######################################### Recommended method by Sun's group for ssDRIP-seq
trim_galore --phred33 --fastqc --stringency 10 --gzip --length 50 --max_n 10 --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10 --paired test_1.fastq.gz test_2.fastq.gz
# --fastqc: perform QC after filtering
# --fastqc_args "<ARGS>": specify arguments for fastqc
# --stringency: stringency level for adapter removal, default 1 (very strict)
# --length: remove reads shorter than this after trimming
# --max_n: remove reads with more than this number of N bases
# --clip_R1: remove bases from 5' end of R1; --clip_R2: remove bases from 5' end of R2
# --three_prime_clip_R1: remove bases from 3' end of R1; --three_prime_clip_R2: remove bases from 3' end of R2
# --paired: perform operations in paired-end mode

# Alignment & sorting
bowtie2 --phred33 -p 24 -x ref.genome -1 test_1_trim3 -2 test_2_trim3 2 > test_align.info | samtools view -bS |samtools sort -@ 10 -m 5G -l 9 -o test.sorted.bam

# bowtie2 comments
# --phred33: input fastq quality format is phred33
# -x: reference genome
# -1: R1 of paired-end; -2: R2 of paired-end
# -U: input file for single-end, multiple files separated by commas
# -t: print alignment time and reference loading time (default to stderr)
# -S: specify output file name (default stdout)

# samtools comments
# -b: output in BAM format
# -S: input is SAM format (auto-detected in newer versions)
# -1: use fastest compression (number 1)
# -l: compression level (letter l), 0=uncompressed, 1=fastest, 9=best compression
# -@: number of threads
# -m: memory per thread
# -o: output file (BAM format)

# Remove low mapping quality reads
samtools view -q 20

# PCR duplicate removal
java -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=true METRICS_FILE=test.matrix INPUT=test.sorted.bam OUTPUT=test.sort.paird_dup.bam

# MarkDuplicates from Picard for duplicate removal
# REMOVE_DUPLICATES: true to remove duplicates, false to only mark
# METRICS_FILE: file to write duplicate metrics
# INPUT: input file
# OUTPUT: output file

# Separate strands
samtools view -b -f 128 -F 16 test.sort.paird_dup.bam > fwd1.bam
samtools view -b -f 80 test.sort.paird_dup.bam > fwd2.bam
samtools merge -f test_fwd.bam fwd1.bam fwd2.bam

samtools view -b -f 144 test.sort.paird_dup.bam > rev1.bam
samtools view -b -f 64 -F 16 test.sort.paird_dup.bam > rev2.bam
samtools merge -f test_rev.bam rev1.bam rev2.bam

samtools index test_fwd.bam # genomeCoverage requires index for generating bw file
samtools index test_rev.bam

# 128: R2; 64: R1; 16: reverse strand
# -f: include; -F: exclude
# -b: output in BAM format

# merge: merge multiple sam/bam files
# -f: overwrite if output exists

# Peak calling
macs2 callpeak -t test_rev.bam -f BAMPE -g 119300826 --keep-dup all -n test_rev -q 0.01
macs2 callpeak -t test_fwd.bam -f BAMPE -g 119300826 --keep-dup all -n test_fwd -q 0.01

# -t: treatment file
# -c: control file
# -f: file format, BAMPE for paired-end
# -g: effective genome size
# --keep-dup: keep all duplicates (already removed previously)
# -n: output file prefix
# --outdir: output directory
# -q: q-value cutoff

# Convert BAM to BigWig
######################################### Recommended by Yang's group
genomeCoverageBed -bga -ibam tes.bam -g genome.sizes -split -pc | sort - -k1,1 -k2,2n > test.bg
genomeCoverageBed -bga -ibam test_fwd.bam -g genome.sizes -split -pc | sort - -k1,1 -k2,2n > forward/test_fwd.bg
genomeCoverageBed -bga -ibam test_rev.bam -g genome.sizes -split -pc | sort - -k1,1 -k2,2n > reverse/test_rev.bg

# bga: output intervals even with score 0
# -ibam: input BAM file
# -g: genome sizes file
# -split: treat N or D in CIGAR as two intervals
# -pc: calculate coverage from left to right end of paired reads

awk -v OFS="\t" 'NR==FNR{sum+=($3-$2)*$4;len+=$3-$2}NR>FNR{print $1,$2,$3,$4*len/sum}' test.bg forward/test_fwd.bg > forward/test.norm_fwd.bg
awk -v OFS="\t" 'NR==FNR{sum+=($3-$2)*$4;len+=$3-$2}NR>FNR{print $1,$2,$3,$4*len/sum}' test.bg reverse/test_rev.bg > reverse/test.norm_rev.bg

# Calculate total sequencing amount (sum) and total length (len) from test.bg
# Normalize fwd/rev.bg by sequencing depth

######################################### Recommended by Sun's group
bamCoverage -v -p 10 -b test_fwd.bam -o test_fwd.bw --binSize 100 --effectiveGenomeSize 119300826 --normalizeUsing RPGC --ignoreForNormalization ChrM ChrC
bamCoverage -v -p 10 -b test_rev.bam -o test_rev.bw --binSize 100 --effectiveGenomeSize 119300826 --normalizeUsing RPGC --ignoreForNormalization ChrM ChrC

bamCoverage -v -p 10 -b test_fwd.bam -o test_fwd.bg --binSize 100 --effectiveGenomeSize 119300826 --normalizeUsing RPGC --ignoreForNormalization ChrM ChrC --outFileFormat bedgraph
bamCoverage -v -p 10 -b test_rev.bam -o test_rev.bg --binSize 100 --effectiveGenomeSize 119300826 --normalizeUsing RPGC --ignoreForNormalization ChrM ChrC --outFileFormat bedgraph

# -v: verbose mode
# -p: number of processes; -b: input BAM file; -o: output file
# --binSize: set bin size; --effectiveGenomeSize: set effective genome size
# --normalizeUsing: normalization method
# RPGC: number of reads per bin divided by scaling factor
# Scaling factor: (total mapped reads * fragment length) / effective genome size

# --scaleFactor: multiply final score by this factor
# Test if --scaleFactor conflicts with --normalizeUsing RPGC
```
