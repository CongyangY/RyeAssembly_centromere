# 01 3d-DNA

nohup ~/juicer/scripts/juicer.sh \
-z ~/juicer/references/genome.fa \
-p ~/juicer/restriction_sites/genome.chrom.sizes \
-y ~/juicer/restriction_sites/genome_DpnII.txt \
-s DpnII \
-d ~/juicer/work/ \
-D ~/juicer \
-t 40 > log.txt


~/sofw/3d-dna/run-asm-pipeline.sh -r 0 reference/genome.fa aligned/merged_nodups.txt > 3d.log 2>&1 &

/home/ubuntu/3d-dna/run-asm-pipeline-post-review.sh \
-r genome.review.assembly \
genome.fa aligned/merged_nodups.txt



# 02 yahs + chromap

# Set up directory structure and link input files
mkdir ../4_chromap; cd ../4_chromap
ln -s ../3_removeOther/result.fa raw.fa 

# Activate software environment
mamba activate chromap

# Build sequence index
samtools faidx raw.fa
chromap -i -r raw.fa -o raw.index
# Note: For large genomes, may need additional parameters: -k 17 -w 13

# Align Hi-C reads
chromap --preset hic -r raw.fa -x raw.index --remove-pcr-duplicates \
        -1 hiC.R1.fq.gz -2 hiC.R2.fq.gz --SAM -o aligned.sam -t 50
# Recommendation:
# 1. Keep threads ≤100 to avoid errors
# 2. Split large datasets (e.g., aligned_1.sam, aligned_2.sam) then merge:
#    samtools merge -@ 60 --write-index aligned.bam aligned_*.bam &

# Run YAHS scaffolding
yahs raw.fa aligned.bam
# Alternative: Disable contig error correction
# yahs -o no --no-contig-ec raw.fa aligned.bam

# Generate Hi-C contact maps with Juicer
~/sofw/yahs-main/juicer pre -a -o out_JBAT yahs.out.bin yahs.out_scaffolds_final.agp raw.fa.fai > out_JBAT.log 2>&1
cat out_JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}' > asm.size
java -Xmx200G -jar ~/juicer/juicer_tools_1.19.02.jar pre out_JBAT.txt out_JBAT.hic asm.size

# Post-processing
juicer post -o out_JBAT out_JBAT.review.assembly out_JBAT.liftoveroscopy raw.fa

############
# No contig error correction version
#############
~/sofw/yahs-main/juicer pre -a -o no no.bin no_scaffolds_final.agp raw.fa.fai > no.log 2>&1
cat no.log | grep PRE_C_SIZE | awk '{print $2" "$3}' > no_asm.size
java -Xmx200G -jar ~/juicer/juicer_tools_1.19.02.jar pre no.txt no.hic no_asm.size
juicer post -o no no.review.assembly no.liftover.agp raw.fa

############
# Chromosome renaming pipeline using GMAP
############
mkdir gmap_rename; cd gmap_rename

# Build GMAP database and align reference CDS
gmap_build -D . -d DB ../out_JBAT.FINAL.fa
gmap -D . -d DB -t 12 -f 2 -n 2 reference.cds.fasta > gmap.gff3 2>gmap.gff3.log
# Where reference.cds.fasta contains published CDS sequences

# Generate chromosome rename list
cut -f1 -d ';' gmap.gff3 | sed 's/ID=//g' | awk -v OFS="\t" '{print $1,$7}' | sort | uniq -c | sort -k1,1rn | head -n 7 > rename.list
# Adjust head -n 7 based on actual chromosome count

# Apply standardized naming (chr1R, chr2R etc.)
awk '{print "s/>"$2"$/>chr"$3"R/g"}' rename.list > sed.cmd  # Changed to chrNR format
sed -f sed.cmd ../out_JBAT.FINAL.fa | head -n7 | seqkit sort - > tmp.chr.fa

# Process unplaced contigs
awk '{print $2}' rename.list | seqkit grep -v -f - ../out_JBAT.FINAL.fa | sed 's/>.*/>contig/g' | seqkit rename - > tmp.contigs.fa

# Combine and create final assembly
cat tmp.chr.fa tmp.contigs.fa > result.fa

# Optional: Reverse complement problematic chromosomes (e.g., chr2R)
# seqkit grep -p 'chr2R' tmp.chr.fa | seqkit seq -r -p - > chr2R.fa
# seqkit grep -v -p 'chr2R' tmp.chr.fa | cat - chr2R.fa | seqkit sort - >> result.fa



# Create directory and link input files
mkdir ../5_ragtag; cd ../5_ragtag
ln -s ../3_removeOther/result.fa raw.fa
ln -s /path/to/reference_genome.fa ref.fa  # Replace with actual reference genome path



# 03 ragtag scffold
# Create directory and link input files
mkdir ../5_ragtag; cd ../5_ragtag
ln -s ../3_removeOther/result.fa raw.fa
ln -s /path/to/reference_genome.fa ref.fa  # Replace with actual reference genome path

# Activate software environment
mamba activate ragtag

# Run RagTag for scaffolding
ragtag.py scaffold ref.fa raw.fa -o ragtag_output -t 50 --aligner minimap2

# Check output results
ls ragtag_output/
# Main output file: ragtag_output/ragtag.scaffold.fasta

# 04 ragtag+yahs
# Create directory and link input files
mkdir ../6_ragtag_yahs; cd ../6_ragtag_yahs
ln -s ../5_ragtag/ragtag_output/ragtag.scaffold.fasta ragtag_scaffold.fa

# Activate software environment
mamba activate chromap

# Build sequence index
samtools faidx ragtag_scaffold.fa
chromap -i -r ragtag_scaffold.fa -o ragtag.index
# Note: For large genomes, may need to adjust index parameters: -k 17 -w 13

# Align Hi-C reads
chromap --preset hic -r ragtag_scaffold.fa -x ragtag.index --remove-pcr-duplicates \
        -1 hiC.R1.fq.gz -2 hiC.R2.fq.gz --SAM -o ragtag_aligned.sam -t 50

# Convert SAM to BAM and sort
samtools view -@ 10 -bS ragtag_aligned.sam | samtools sort -@ 10 -o ragtag_aligned.bam
samtools index ragtag_aligned.bam

# Run YAHS for Hi-C scaffolding
yahs ragtag_scaffold.fa ragtag_aligned.bam

# Generate Hi-C contact maps
~/sofw/yahs-main/juicer pre -a -o ragtag_yahs yahs.out.bin yahs.out_scaffolds_final.agp ragtag_scaffold.fa.fai > ragtag_yahs.log 2>&1
cat ragtag_yahs.log | grep PRE_C_SIZE | awk '{print $2" "$3}' > ragtag_asm.size
java -Xmx200G -jar ~/juicer/juicer_tools_1.19.02.jar pre ragtag_yahs.txt ragtag_yahs.hic ragtag_asm.size

# Post-processing
juicer post -o ragtag_yahs ragtag_yahs.review.assembly ragtag_yahs.liftover.agp ragtag_scaffold.fa

# Note: Manual adjustments may be needed based on CENH3 ChIP-seq signals and Hi-C interaction heatmaps in JuicerBox
# These adjustments can help correct potential misassemblies and improve chromosome-scale scaffolding

# Chromosome renaming (optional)
mkdir gmap_rename; cd gmap_rename
gmap_build -D . -d DB ../ragtag_yahs.FINAL.fa
gmap -D . -d DB -t 12 -f 2 -n 2 reference.cds.fasta > gmap.gff3 2>gmap.gff3.log

# Generate rename list
cut -f1 -d ';' gmap.gff3 | sed 's/ID=//g' | awk -v OFS="\t" '{print $1,$7}' | sort | uniq -c | sort -k1,1rn | head -n 7 > rename.list
awk '{print "s/>"$2"$/>chr"$3"R/g"}' rename.list > sed.cmd
sed -f sed.cmd ../ragtag_yahs.FINAL.fa | head -n7 | seqkit sort - > tmp.chr.fa

# Process unplaced contigs
awk '{print $2}' rename.list | seqkit grep -v -f - ../ragtag_yahs.FINAL.fa | sed 's/>.*/>contig/g' | seqkit rename - > tmp.contigs.fa

# Combine to create final assembly
cat tmp.chr.fa tmp.contigs.fa > result.fa
