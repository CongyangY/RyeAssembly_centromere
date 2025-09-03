# convert bam to fastq
ls HIFI_*.bam | while read id 
do
bam2fastq -o ${id/.bam/} $id
done

# assembly
hifiasm-o WN -t 190 -l 2 --ul ./Ont.fastq.gz *.fastq.gz 
gfatools  gfa2fa WN.bp.p_ctg.gfa > WN.bp.p_ctg.fasta
seqkit stat -a WN.bp.p_ctg.fasta > WN.bp.p_ctg.fasta.stat

# remove Organelles 
mkdir 1_Organelles 
cd 1_Organelles
ln -s  ../WN.bp.p_ctg.fasta raw.fa
minimap2 -x asm5 -t 100 ~/Public/Secale_cereale_mitochondrion.fasta raw.fa > mitochondrion.paf
minimap2 -x asm5 -t 100 ~/Public/Secale_cereale_chloroplast.fasta raw.fa > chloroplast.paf

ls *paf | while read id 
do
python ~/scripts/filter_pollute_contigs.py $id > ${id/paf/txt}
done

cat *txt | seqkit grep -v -f - raw.fa > result.fa

# remove otherPollute
mkdir ../2_OtherPollute
mmseqs createdb nt  mmseq_targetDB
mmseqs createindex mmseq_targetDB  tmp
ln -s ../1_
# 比对
mmseqs easy-search --format-mode 0 --search-type 3 query.fasta mmseq_targetDB queryIn_nt.mmseq2  tmp_dir
