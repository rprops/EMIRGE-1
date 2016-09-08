# EMIRGE
Set of scripts/guidelines for EMIRGE analysis of metagenomic data

#############################################################################
### How to use EMIRGE for reconstructing 16S sequences frome metagenomic data
### Publications of interest:
### 1. https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-5-r44
### 2. http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0056018
### 3. http://aem.asm.org/content/82/5/1423.abstract
##############################################################################

##################################################################################
### If you wish to run these processes in the background then add an & at the back
##################################################################################

### Copy metaG data from nfs to scratch
### This code copies only the fastq.gz in from the respective folders to the new location
### For some samples there are fastq of different shotgun sequencing methods
### I'll run EMIRGE for both of these.
rsync -a --include '*/' --include '*.fastq.gz' --exclude '*' /nfs/vdenef-lab/Shared/Sequence_data/CSP_LM13/LM13_JGI_MetaG /scratch/vdenef_fluxm/rprops/metaG --progress

### These folders didn't have a fast.gz, only the fastaq
cp /nfs/vdenef-lab/Shared/Sequence_data/CSP_LM13/LM13_JGI_MetaG/Fa13.BD.MM110.DN/8202.1.92752.GGACTCC-TATCCTC.fastq /scratch/vdenef_fluxm/rprops/metaG/LM13_JGI_MetaG/Fa13.BD.MM110.DN/
cp /nfs/vdenef-lab/Shared/Sequence_data/CSP_LM13/LM13_JGI_MetaG/Sp13.BD.MLB.SN/8199.1.94074.TGTGAA.anqdp.fastq /scratch/vdenef_fluxm/rprops/metaG/LM13_JGI_MetaG/Sp13.BD.MLB.SN

### Unzip all .gz fastq in the folders
### Also make a READme file for this (and store copy on nfs drive)
### -r stands for recursive, i.e. enter all directories in the provided path and unzip stuff
### Best to run this separately in a pbs script because it seems to automatically parallelize this over multiple CPUs
gunzip -r /scratch/vdenef_fluxm/rprops/metaG/

### Extract forward and reverse fastq from original fasta (script comb_to_rever_forw_fastq)
qsub comb_to_rever_forw_fastq.pbs

## Copy the fastq of your samples to the correct directory (example for one sample)
cp /nfs/vdenef-lab/Shared/Sequence_data/CSP_LM13/LM13_JGI_MetaG/Fa13.BD.MM110.DN/Fa13.BD.MM110.DN.10Mpairs.PE.2.fastq /scratch/vdenef_fluxm/rprops/metaG
cp /nfs/vdenef-lab/Shared/Sequence_data/CSP_LM13/LM13_JGI_MetaG/Fa13.BD.MM110.DN/Fa13.BD.MM110.DN.10Mpairs.PE.1.fastq /scratch/vdenef_fluxm/rprops/metaG

### Cluster silva_v123 NR99 database to 97 % with usearch
### In the original article they also prune the database from sequences with a length less than 1200 and more than 1900
usearch -cluster_fast SILVA_123_SSURef_Nr99_tax_silva2.fasta -id 0.97 -centroids SSURef_NR97_123_for_emirge.fasta -uc SSURef_NR97_123_for_emirge.clusters.uc

### Replace non-standard base characters in reference database (script in same directory as database)
### This was run with python-anaconda2/latest - could also work with the 201607 but not sure yet
python fix_nonstandard_chars.py < SSURef_NR97_123_for_emirge.fasta > /scratch/vdenef_fluxm/rprops/databases/SSURef_NR97_123_for_emirge2.fasta

### Create bowtie-index for your reference database
bowtie-build SSURef_NR97_123_for_emirge2.fasta SSU_candidate_db_btindex

### Make sure your these scripts are in the directory where all your sample folders are
### Run EMIRGE (see batch scripts for parallelizing)
### Make sure that --Phred33 is present (this is mandatory for new fastq files from Illumina)
### Takes approx. 20-40h on 10 cores per sample for 65 iterations
### Running with -j 1.0 (thus 100 % identity) which is different from the 0.97 that others use (unique seqs)
### There is also little point in parallelizing beyond 10 - 12 cores because of the dependency on usearch (uses 8 cores)
bash -x Batch_01.sh
bash -x Batch_02.sh
bash -x Batch_03.sh
bash -x Batch_04.sh
bash -x Batch_05.sh
bash -x Batch_06.sh

### Reformat fasta header to sort sequences according to abundances
### Usage: emirge_rename_fasta.py [options] <iter.DIR> > renamed.fasta
### Make sure biopython/1.60 is loaded !!!!
### Make sure that the files mapped.reads.txt, total.emirge.cons.fasta and reads.info.txt don't exist yet.
bash -x rename.sh

### Count number of mapping reads from bowtie bam files in final iteration folders
### This is needed to estimate the number of reads for each sequence
bash -x extract_mappings.sh

### Concatenate all the fasta files in one fasta for taxonomic classification
bash -x concat.emirge.sh

### Classify the sequences according to the TaxASS pipeline
# Step 1: All redundant information from the first line in the fasta file
# Retain sample names only
# Fasta file cannot be aligned
awk '{print $1}' total.emirge.renamed.fasta > total.emirge.renamed2.fasta
mv total.emirge.renamed2.fasta total.emirge.renamed.fasta

# Step 2: make BLAST database file (blast)
# makeblastdb -dbtype nucl -in FreshTrain18Aug2016.fasta -input_type fasta -parse_seqids -out FWonly_18Aug2016custom.db

### Taxonomy folder
# Step 3: Run blast and reformat output blast file
blastn -query total.emirge.renamed.fasta -task megablast -db /scratch/vdenef_fluxm/rprops/Emirge/BLAST/FWonly_18Aug2016custom.db -out custom.blast -outfmt 11 -max_target_seqs 5 system(blastn -query total.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta -task megablast -db /scratch/vdenef_fluxm/rprops/Emirge/BLAST/FWonly_18Aug2016custom.db -out custom.blast -outfmt 11 -max_target_seqs 5 -num_threads 40
blast_formatter -archive custom.blast -outfmt "6 qseqid pident length qlen qstart qend" -out otus.custom.blast.table

# Step 4: Correct BLAST pident using custom script
### This accounts for sequence length differences
Rscript calc_full_length_pident.R otus.custom.blast.table otus.custom.blast.table.modified

# Step 5: Filter BLAST results
Rscript filter_seqIDs_by_pident.R otus.custom.blast.table.modified ids.above.97 97 TRUE 
Rscript filter_seqIDs_by_pident.R otus.custom.blast.table.modified ids.below.97 97 FALSE

# Step 6: 
mkdir plots
Rscript plot_blast_hit_stats.R otus.custom.blast.table.modified 97 plots

# Step 7: recover sequence IDs left out of blast (python, bash)
python find_seqIDs_blast_removed.py total.emirge.renamed.fasta otus.custom.blast.table.modified ids.missing
cat ids.below.97 ids.missing > ids.below.97.all

# Step 8: create fasta files of desired sequence IDs (python)
python create_fastas_given_seqIDs.py ids.above.97 total.emirge.renamed.fasta otus.above.97.fasta
python create_fastas_given_seqIDs.py ids.below.97.all total.emirge.renamed.fasta otus.below.97.fasta

# Step 9: Screen for minimum length and max number of amibguous bases
# Assign taxonomy to each fasta file
### Only interested in unique seqs (otherwise you create non-existent sequence variants)
mothur "#unique.seqs(fasta=otus.below.97.fasta)"
mothur "#unique.seqs(fasta=otus.above.97.fasta)"
### Classify
mothur "#classify.seqs(fasta=otus.below.97.unique.fasta, template=/scratch/vdenef_fluxm/rprops/databases/silva.nr_v123.align, taxonomy=/scratch/vdenef_fluxm/rprops/databases/silva.nr_v123.tax, method=wang, probs=T, processors=10, cutoff=80)"
mothur "#classify.seqs(fasta=otus.above.97.unique.fasta, template=/scratch/vdenef_fluxm/rprops/databases/FreshTrain18Aug2016.fasta,  taxonomy=/scratch/vdenef_fluxm/rprops/databases/FreshTrain18Aug2016.taxonomy, method=wang, probs=T, processors=10, cutoff=0)"

# Step 10: combine taxonomy files and names files
cat otus.above.97.unique.FreshTrain18Aug2016.wang.taxonomy otus.below.97.unique.nr_v123.wang.taxonomy > otus.97.taxonomy
cat otus.above.97.names otus.below.97.names > otus.97.names

### Run R script to create Sequence table (same format as OTU table)
Rscript emirge_format.R mapped.reads.txt otus.97.taxonomy read.info.txt otus.97.names
