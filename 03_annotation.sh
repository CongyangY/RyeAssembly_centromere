######################### genes
```shell
# RNA-seq based prediction

## mikado (highly recommended, integrates multiple assembly methods)
### Installation
#### Latest method
```shell
mamba create -n mikado2 -c bioconda sqlalchemy==1.4.46 mikado python=3.10 pandas=1.4
mamba activate mikado2
mamba install -c bioconda prodigal blast transdecoder diamond stringtie samtools hisat2 star scallop 

mamba create -n mikado -c bioconda cufflinks portcullis mikado
export PATH=$PATH:~/mambaforge/envs/mikado/bin/
```

```shell
git clone https://github.com/EI-CoreBioinformatics/mikado.git
cd mikado

# 
sed -i 's/sqlalchemy==.*/sqlalchemy==1.4.46/g' environment.yml
mamba env create -f environment.yml

conda activate mikado2
pip3 install wheel

sed -i 's/sqlalchemy==.*/sqlalchemy==1.4.46/g' requirements.txt
pip3 install -r requirements.txt

pip3 install Cython
python3 setup.py bdist_wheel
pip3 install dist/*.whl
mikado -h

mamba install -c bioconda prodigal blast transdecoder diamond stringtie samtools hisat2 star scallop  
```
#### Previous installation method

```shell
mamba create -n mikado -c bioconda prodigal blast transdecoder diamond stringtie cufflinks 
mamba activate mikado
pip3 install mikado
mamba install -c bioconda samtools hisat2 star scallop 

cpan URI::Escape

# Create separate environment for portcullis to avoid version conflicts
mamba create -n mikado2 -c bioconda cufflinks portcullis
mamba activate mikado2
pip install mikado

export PATH=$PATH:~/mambaforge/envs/mikado2/bin/

# Alternative installation with specific Python version
mamba create -n mikado -c bioconda prodigal blast transdecoder diamond stringtie mikado samtools hisat2 star scallop gmap python=3.9
# Note: cufflinks not included, gmap can be used for third-generation sequencing alignment
```

- Need to modify code at line 692 in ~/mambaforge/envs/mikado/lib/python3.6/site-packages/Mikado/loci/superlocus.py
- Python 3.6 doesn't support `asyncio.run()` function (introduced in Python 3.7)
- Modify code as follows:
```shell
# Original code
# data_dict = asyncio.run(data_dict)
# asyncio.run(intron_loader)

# Modified code
loop = asyncio.get_event_loop()
data_dict = loop.run_until_complete(data_dict)
intron_loader_result = loop.run_until_complete(intron_loader)
```

### Execution
1. Generate configuration file
```shell
daijin configure --scheduler "local" -o daijin.yaml -al hisat star -as stringtie scallop -t 60 --genome /path/to/genome.fa --use-transdecoder --sample-sheet sample.sheet --scoring plant.yaml -m permissive --copy-scoring plant.yaml --prot-db protein.aa

# Not recommended to use mikado with third-generation data due to potential errors
# For third-generation sequencing use lal, but second-generation settings cannot be removed
daijin configure --scheduler "local" -o daijin.yaml -al star hisat -as stringtie scallop -lal star gmap -t 60 --genome /path/to/genome.fa --use-transdecoder --sample-sheet sample.sheet --scoring plant.yaml -m permissive --copy-scoring plant.yaml --prot-db /path/to/uniprot_sprot.fasta
```
- `--scheduler ""`: Specify local execution instead of job submission system
- `-al hisat star`: Specify hisat2 and star as alignment software
- `-as cufflinks stringtie scallop`: Specify cufflinks, stringtie and scallop as assembly software
- `-t 60`: Specify thread count
- `--genome /path/to/genome.fa`: Specify genome to annotate
- `--use-transdecoder`: Use transdecoder for ORF identification
- `--sample-sheet sample-sheet`: Sequencing file configuration information
- `--scoring plants.yaml --copy-scoring plants.yaml`: Specify scoring configuration file and copy to current directory
- `-m permissive`: mikado running mode, options: nosplit, split, permissive, stringent, lenient
- `--prot-db protein.aa`: Specify homologous protein database

`sample.sheet` format:
```yaml
/path/to/2B-leaf-1RNA_trim2_1P.fq.gz	/path/to/2B-leaf-1RNA_trim2_2P.fq.gz	2B_rep1	fr-unstranded	False
```
- Columns: reads1 filename, reads2 filename, sample name, sequencing type, whether third-generation sequencing
- File uses tab separation, first column must be relative path

2. Assembly
```
daijin assemble -nd --cores 100 daijin.yaml
# -nd: don't use DRMAA mode (not available on most servers)
```
- May fail during plotting step - comment out plot functions in snakemake files
- Snakemake files location: `/path/to/mambaforge/envs/mikado/lib/python3.6/site-packages/Mikado/daijin/assemble.smk` and `/path/to/mambaforge/envs/mikado/lib/python3.6/site-packages/Mikado/daijin/mikado.smk`
- If software like hisat2 or star cannot be found, consider:
  - 1. Write commands to script file and run script (include environment activation)
  - 2. Use conda run command: `conda run daijin assemble -nd --cores 100 daijin.yaml`
- **For large genomes, need to modify multiple places: samtools index needs -c parameter, portcullis prep and junc parts also need -c parameter**
- If error occurs: `AttributeError: module 'numpy' has no attribute 'warnings'. Did you mean: 'hanning'?`, modify file "/path/to/mambaforge/envs/mikado2/lib/python3.10/site-packages/Mikado/scales/calculator.py"
- Change `np.warnings.filterwarnings("ignore")` to `warnings.filterwarnings("ignore")`

3. mikado further identification and filtering
```
daijin mikado -nd --threads 100 Daijin/mikado.yaml
```
Final file: `/5-mikado/pick/permissive/mikado-permissive.loci.gff3`

**Notes**
- Errors like `Error in library(seqLogo)` in transdecoder can be ignored (plotting component, doesn't affect results)
- If final mikado picker step fails, modify script at `~/sofw/mambaforge-pypy3/envs/mikado/lib/python3.5/site-packages/Mikado/loci/abstractlocus.py`
- Replace `graph.add_path(segments)` with `networkx.add_path(graph, segments)`

## stringtie
```
mamba activate rna_anno
mamba install -c bioconda stringtie -y
stringtie -i rna.bam -o stringtie.gtf 

mamba install bioconda::agat

# If error occurs: MSG: Each line of the file must be less than 65,536 characters.
# Fold fasta file: seqkit sort input.fa > output.fa 
gffread ./assembled.gtf -T -F -o ./assembled.gff
agat_convert_sp_gxf2gxf.pl --gff ./assembled.gff --output ./assembled.cleaned.gff

agat_sp_filter_gene_by_length.pl --size 200 --test ">=" --gff ./assembled.cleaned.gff --output ./assembled.filtered_size.gff3

agat_sp_filter_by_ORF_size.pl --size 100 --gff ./assembled.filtered_size.gff3 --output ./assembled.filtered_orf.gff3

agat_sp_fix_overlaping_genes.pl --gff assembled.filtered_orf3_sup100.gff --output fixed.gff3

agat_sp_fix_features_locations_duplicated.pl --gff fixed.gff3 --output ont_ready.gff3
```

## cufflinks 
```shell
mamba create -n rna_anno_py2 cufflinks 
mamba activate rna_anno_py2

cufflinks -i rna.bam -o cufflinks
```

## scallop
```
mamba activate rna_anno
mamba install -c bioconda scallop

scallop -i ../rna.bam -o scallop.gtf > scallop.lop 2>&1 &
```

## Pasapipeline
```shell
# Installation
mamba create -n pasa -c bioconda pasa
mamba activate pasa
mamba install -c bioconda trinity
```

```shell
# Usage
samtools merge -@ 40 --write-index all.bam /path/to/Daijin/2-alignments/output/*.sorted.bam

Trinity --genome_guided_bam all.bam \
        --genome_guided_max_intron 10000 \
        --max_memory 50G \
        --CPU 10 \
        --output Trinity_GG_out

# Maximum 12 threads to avoid errors
seqclean Trinity_GG_out/Trinity-GG.fasta -c 10

ln -s Trinity_GG_out/Trinity-GG.fasta ./transcripts.fa
ln -s Trinity_GG_out/Trinity-GG.fasta.clean ./transcripts.fa.clean

Launch_PASA_pipeline.pl -c alignAssembly.config -C -R --ALIGNERS blat,gmap,minimap2 -g genome.fa -t transcripts.fa.clean -T -u transcripts.fa --CPU 180 \
-A \# Key: add -A parameter to enable annotation 
--annots_gff3 existing_annotation.gff3  # Introduce reference annotation

# For large genomes, not recommended to use blat
# minimap2 can be used, but for large genomes, index building needs -c parameter

############### alignAssembly.config content

cat alignAssembly.config
## templated variables to be replaced exist as <__var_name__>

# Pathname of an SQLite database
# If the environment variable DSN_DRIVER=mysql then it is the name of a MySQL database
DATABASE=/tmp/mydb.sqlite

#######################################################
# Parameters to specify to specific scripts in pipeline
# create a key = "script_name" + ":" + "parameter"
# assign a value as done above.

#script validate_alignments_in_db.dbi
validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=80
validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=80
validate_alignments_in_db.dbi:--NUM_BP_PERFECT_SPLICE_BOUNDARY=0

#script subcluster_builder.dbi
subcluster_builder.dbi:-m=50
```

# Protein-based prediction
## exonerate
```shell
for i in $(seq 1 80); do exonerate --model protein2genome --percent 95 --showvulgar no --showalignment no --showquerygff yes --querychunkid ${i} --querychunktotal 80 /path/to/protein.faa /path/to/genome.fa > exonerate_chunk${i}.gff & done
```

## genomethreader
```shell
conda install bioconda::genomethreader

gth -genomic target.fa -cdna ref.cds -protein ref.aa -gff3out -skipalignmentout -o target.gff &
```

## genoma
```
mamba activate mikado2
mamba install -c bioconda gemoma

GeMoMa GeMoMaPipeline \
threads=60 AnnotationFinalizer.r=NO p=false o=true tblastn=false \
t=genome.fa outdir=gemoma \
s=own i=ref1 a=ref1.gff g=ref1.fasta \
s=own i=ref2 a=ref2.gff g=ref2.fasta &

# Example with multiple references
GeMoMa GeMoMaPipeline threads=190 AnnotationFinalizer.r=NO p=false o=true tblastn=false t=/path/to/genome.fa s=own i=ref1 a=Arabidopsis_thaliana.TAIR10.60.gff3 g=Arabidopsis_thaliana.TAIR10.dna.toplevel.fa s=own i=ref2 a=B73_v4.gff3 g=B73_v4.fa s=own i=ref3 a=CSGL.gene.gff3 g=CSGL.fa s=own i=ref4 a=GCA_900231445.1_Svevo.v1_genomic.gff g=GCA_900231445.1_Svevo.v1_genomic.fna s=own i=ref5 a=rice_v7.gff3 g=rice_v7.fa s=own i=ref6 a=wheat_WEW_AB.HighConf.gff3 g=wheat_WEW_AB.fa debug=true r=MAPPED ERE.s=FR_UNSTRANDED ERE.m=rna.bam ERE.mil=30 ERE.mmq=20 ERE.c=true d=DENOISE DenoiseIntrons.m=25000 DenoiseIntrons.me=0.03 DenoiseIntrons.c=10 threads=150 outdir=gemoma_results 

# Easily exceeds memory limits, recommend using which GeMoMa to find script location and modify Xmx parameters

GeMoMa_gff_to_gff3.pl GeMoMa/final_annotation.gff > gemoma.gene.gff
```

## Gmap 
```
mamba activate rna_anno
mamba install -c bioconda gmap
gmap_build -D . -d DB target.genome
# target.genome is the contig level of polyploid genome assembly

gmap -D . -d DB -t 12 -f 2 -n $N reference.cds.fasta > gmap.gff3
# target.genome is the contig level of polyploid genome assembly
# reference.cds.fasta is the coding sequences of diploid genome
# $N could be the ploidy of your target genome (e.g., $N=4 for tetraploid)
# For large genomes, use gmapl
```

# De novo prediction
## braker
```shell
# Stop docker
sudo systemctl stop docker
# Create new directory
mkdir /path/to/docker_data/

# Change Docker data directory in /etc/docker/daemon.json
{
    "registry-mirrors": ["https://docker.mirror.url"],
    "data-root": "/path/to/docker_data/"
}

# Restart docker
sudo systemctl start docker

# Verify
docker info | grep "Docker Root Dir"
```

```shell
# Set permissions
sudo chown -R 1000:100 /path/to/braker3
sudo chmod -R u+rwX /path/to/braker3

# Run docker container
docker run --user 1000:100 --user root --rm -it --privileged \
-v /path/to/data/:/home/jovyan/data/ --name braker3-container teambraker/braker3:latest bash

# Fix script issue
cp /opt/TSEBRA/bin/best_by_compleasm.py /opt/TSEBRA/bin/best_by_compleasm.py.bc
wget https://raw.githubusercontent.com/Gaius-Augustus/TSEBRA/refs/heads/fix_braker_issue_882/bin/best_by_compleasm.py
mv best_by_compleasm.py /opt/TSEBRA/bin/best_by_compleasm.py

# Copy and modify test script
cp /opt/BRAKER/example/docker-tests/test3.sh work.sh
# Modify thread count (don't exceed 40 to avoid errors)

# --addUTR=on recommended, requires:
# git clone https://github.com/Gaius-Augustus/GUSHR.git
# mamba install -c bioconda java-jdk=8 -y
# Add --GUSHR_PATH=/home/jovyan/GUSHR/ to braker command

# --skipAllTraining optional
# Example command:
# ( time braker.pl --genome=/home/jovyan/data/genome.fa --prot_seq=/home/jovyan/data/homolog/orthodb.fa --bam=/home/jovyan/data/Ont/all_add_Ont.bam --workingdir=$wd --threads 40 --addUTR=on --busco_lineage eukaryota_odb10 --GUSHR_PATH=/home/jovyan/GUSHR/ ) &> work.log

bash work.sh

# After docker exit, restore permissions
# sudo chown -R user:group /path/to/data/
```

## augustus
```shell
mamba activate rna_anno
mamba install -c bioconda augustus=3.2.3
# Incompatible with genemark_es, need separate environments

# Recommended to install via apt-get
sudo apt install augustus augustus-data augustus-doc

# Check available datasets
augustus --species=help

# Run with existing model
augustus --species=maize --gff3=on genome.masked.fa > augustus.gff & 

# Without existing model, train with related species annotation
autoAugTrain.pl --genome=species_ref.fa --trainingset=species_ref.gff --species=species_ref
augustus --species=species_ref --gff3=on genome.fasta > augustus.gff 
augustus_GFF3_to_EVM_GFF3.pl augustus.gff > augustus.gene.gff
# Important step: convert to EVM format
```

## GeneMark_es
```shell
# Download genemark_es and key from http://topaz.gatech.edu/GeneMark/license_download.cgi
# Follow INSTALL instructions to copy key

# Conda installation
mamba create -n genemark -c thiesgehrmann genemark_es -y

export PERL5LIB=/path/to/mambaforge/envs/rna_anno/lib/site_perl/5.26.2/:$PERL5LIB 

gmes_petap.pl --ES --cores 30 --sequence genome.fa &
# Thread count should not exceed 60

/path/to/EvmUtils/misc/GeneMarkHMM_GTF_to_EVM_GFF3.pl genemark.gtf > genemark.gene.gff3
```

# Integration
```shell
mamba activate rna_anno 
mamba install bioconda::evidencemodeler
mamba install bioconda::perl-db-file

export PERL5LIB=/path/to/mambaforge/envs/rna_anno/bin/PerlLib/:$PERL5LIB 

EVidenceModeler --genome masked.genome.fa \
--gene_predictions gene_predictions.gff3 \
--protein_alignments protein_alignments.gff3 \
--transcript_alignments transcript_alignments.gff3 \
--CPU 100 --segmentSize 100000 --overlapSize 10000 \
--weights weights.txt --sample_id output_prefix
```

Weight file format:
```shell
PROTEIN	protein_source1	1
PROTEIN	protein_source2	5
TRANSCRIPT	transcript_source1	1
TRANSCRIPT	transcript_source2	10
ABINITIO_PREDICTION	fgenesh	1
ABINITIO_PREDICTION	genemark	1
ABINITIO_PREDICTION	glimmerHMM	1
```
- Tab-separated, first column: PROTEIN, TRANSCRIPT, or ABINITIO_PREDICTION
- Second column: source identifier from annotation files
- Third column: weight value
- For multiple annotation files of same type, merge files and use second column to distinguish

# TE-related genes
If over 30% of a gene's amino acid sequence is annotated as transposon elements, it's considered a TE-related gene
```shell
# Software download
mkdir my_interproscan
cd my_interproscan
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.68-100.0/interproscan-5.68-100.0-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.68-100.0/interproscan-5.68-100.0-64-bit.tar.gz.md5

md5sum -c interproscan-5.68-100.0-64-bit.tar.gz.md5

tar -pxvzf interproscan-5.68-100.0-*-bit.tar.gz
python3 setup.py -f interproscan.properties

# Annotation
./interproscan.sh -i test_all_appl.fasta -f tsv -dp
./interproscan.sh -i test_all_appl.fasta -f tsv
# interproscan uses MD5 to check for duplicate sequences
# -dp parameter disables this check
```

# Functional annotation
- Installation (recommend direct download instead of conda)
```shell
mamba create -n interproscan -c bioconda interproscan
mamba activate interproscan

# Set versions
version_major=5.59
version_minor=91.0
CONDA_PREFIX=/path/to/mambaforge/envs/interproscan

# Get database MD5
wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${version_major}-${version_minor}/interproscan-${version_major}-${version_minor}-64-bit.tar.gz.md5
# Get databases (with core for faster download)
wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${version_major}-${version_minor}/interproscan-${version_major}-${version_minor}-64-bit.tar.gz
# Checksum
md5sum -c interproscan-${version_major}-${version_minor}-64-bit.tar.gz.md5
# Extract
tar xvzf interproscan-${version_major}-${version_minor}-64-bit.tar.gz
# Remove default sample DB
rm -rf $CONDA_PREFIX/share/InterProScan/data/
# Copy new DB
cp -r interproscan-${version_major}-${version_minor}/data $CONDA_PREFIX/share/InterProScan/
```

- Usage
```shell
mkdir function
interproscan.sh -i test_proteins.fasta -f tsv \
	-d ./function -cpu 100 -goterms -pa -hm -dp

/path/to/interproscan-5.59-91.0/interproscan.sh \
-i output_prefix.EVM.clean.pep \
-f tsv \
-d ./function \
-cpu 32 \
-goterms
```


########### TE

nohup ~/sofw/EDTA/EDTA.pl --genome genome.fa --curatedlib  ~/sofw/CLARI-TE/clariTeRep.convert.fa --sensitive 1 --anno 1 --threads 160 &
nohup python ./main.py --genome genome.fa --outdir ./TE/ --thread 120 --plant 1 --curated_lib ~/sofw/CLARI-TE/clariTeRep.convert.fa --recover 1 --annotate 1 --intact_anno 1 >genome.HiTE.log 2>&1 &
