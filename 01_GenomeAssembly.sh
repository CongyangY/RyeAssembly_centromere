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
ln -s ../1_organelles ./raw.fa 
# alignment 
mmseqs easy-search --format-mode 0 --search-type 3 query.fasta mmseq_targetDB queryIn_nt.mmseq2  tmp_dir
cut -f2 queryIn_nt.mmseq2 | sed 's/\..*//g' | awk -v OFS="\t" 'NR==FNR{a[$1]=$3} NR>FNR&&($1 in a){print $1,a[$1]}' ~/Public/nt/nucl_gb.accession2taxid - | cut -f 2 | taxonkit lineage | taxonkit reformat -F | cut -f 3 | paste queryIn_nt.mmseq2 - | awk -v OFS="\t" '$7<=$8{print $1,$7,$8,$13,$4} $7>$8{print $1,$8,$7,$13,$4}' > queryIn_nt.taxonomy.bed 

# Extract only Kingdom, Phylum, Class (first 3 taxonomic levels)
cut -f 1-3 -d ';' queryIn_nt.taxonomy.bed > queryIn_nt.taxonomy.bed2

# Split data by taxonomic group, calculate alignment length per contig per taxon
cut -f4 queryIn_nt.taxonomy.bed2 | sort -u | while read taxon; do 
  awk -v tax="$taxon" '$4==tax {print $0}' queryIn_nt.taxonomy.bed2 \
  | sortBed -i - \
  | mergeBed -i - \
  | awk -v tax="$taxon" -v OFS="\t" '{contig_len[$1]+=$3-$2} 
       END {for (contig in contig_len) print contig,tax,contig_len[contig]}' \
  >> queryIn_nt.taxonomy.pie &  
done

# Find the most abundant taxon per contig (by alignment length)
sort -k1,1 -k3rn,3 queryIn_nt.taxonomy.pie \
| awk -v OFS="\t" '($3>max_len[$1]){max_len[$1]=$3;top_taxon[$1]=$2} 
     END{for (contig in top_taxon) print contig,top_taxon[contig],max_len[contig]}' \
| sort -k1,1 > queryIn_nt.taxonomy.abundant.list 

# Identify potential contaminants (excluding Streptophyta)
grep -v 'Streptophyta' queryIn_nt.taxonomy.abundant.list \
| cut -f1 | sort -u > possible_contaminants.list

# Generate detailed reports for each potential contaminant
mkdir tmp
cat possible_contaminants.list | while read contig; do
  grep "$contig" queryIn_nt.taxonomy.pie > tmp/${contig}.list
done

# Generate list of contigs lacking taxonomic assignment
awk '!$4' queryIn_nt.taxonomy.bed | cut -f1 | sort -u > unassigned_contigs.list

# Select contigs to keep based on primary assignment being Insecta
grep 'Insecta' query_nt.short.taxonomy.abundant.list \
| awk 'NR==FNR{a[$1]} NR>FNR&&!($1 in a){print $1}' \
  - query_nt.short.saved_contig > query_nt.remove_contig.list

# Add prokaryotic contigs to removal list
egrep 'Bacteria|Archaea' query_nt.short.taxonomy.abundant.list \
| cut -f1 >> query_nt.remove_contig.list

# Generate final cleaned fasta file
seqkit grep -v -f query_nt.remove_contig.list contigs.noOrganelles.fasta \
> query_nt.clean.fa



verkko -d asm \
  --hifi hifi/*.fastq.gz \
  --nano ont/*.fastq.gz \
  --hic1 hic/*R1*fastq.gz  \
  --hic2 hic/*R2*fastq.gz
