# EMIRGE
Set of scripts/guidelines for EMIRGE analysis of metagenomic data.

## The script has been tested with the following dependencies

```
bowtie v1.0.0
usearch v8.1
python-anaconda2/201607
emirge v0.60.3
mothur v1.38.1 
biopython v1.60
ncbi-blast v2.2.29
```

## Workflow
Copy metaG data from nfs to scratch. This code copies only the fastq.gz in from the respective folders to the new location. 
```
rsync -a --include '*/' --include '*.fastq.gz' --exclude '*' /nfs/vdenef-lab/Shared/Sequence_data/CSP_LM13/LM13_JGI_MetaG /scratch/vdenef_fluxm/rprops/metaG --progress
```

Unzip all .gz fastq in the folders. -r stands for recursive, i.e. enter all directories in the provided path and unzip stuff. Best to run this separately in a pbs script if possible, because it automatically parallelizes this over multiple cores, or run it as a background process with <code>nohup</code>.
```R
gunzip -r /scratch/vdenef_fluxm/rprops/metaG/
nohup gunzip -r /scratch/vdenef_fluxm/rprops/metaG/ &
```
Extract forward and reverse fastq from original fastq (script comb_to_rever_forw_fastq.pbs).
```
qsub comb_to_rever_forw_fastq.pbs
```

Cluster silva_v123 NR99 database to 97 % with usearch. In the original article they also prune the database from sequences with a length less than 1200 and more than 1900. 

```
usearch -cluster_fast SILVA_123_SSURef_Nr99_tax_silva2.fasta -id 0.97 -centroids SSURef_NR97_123_for_emirge.fasta -uc SSURef_NR97_123_for_emirge.clusters.uc
```

Replace non-standard base characters in reference database (script in same directory as database). This requires `biopython`.

```
python fix_nonstandard_chars.py < SSURef_NR97_123_for_emirge.fasta > SSURef_NR97_123_for_emirge2.fasta
```

Create bowtie-index for your reference database
```
bowtie-build SSURef_NR97_123_for_emirge2.fasta SSU_candidate_db_btindex
```
Make sure these scripts are in the directory where all your sample folders are. Run EMIRGE (see batch scripts for parallelizing). Make sure that --Phred33 is present (this is mandatory for new fastq files from Illumina). Takes approx. 20-40h on 10 cores per sample for 65 iterations. Running with -j 1.0 (thus 100 % identity) which is different from the 0.97 that others use (unique seqs). There is also little point in parallelizing beyond 10 - 12 cores because of the dependency on usearch (uses 8 cores). Standard settings are mean insert size of 200 bp (-i option) with standard deviation of 50 bp (-s option) and max. read length of 150 (-l option). Adjust these settings in emirge.pbs if so desired/required.
You can make several batch scripts (such as below) to run multiple emirge runs in parallel. Replace the directory names with yours.
```
set -e
for i in Fa13.BD.MLB.DN Fa13.BD.MM110.DN Fa13.BD.MM110.SD Fa13.BD.MM110.SN Fa13.BD.MM15.DN ; do
        cd $i
		# print working directory to the screen
        echo $i
		# copy mapping.pbs from above directory
        cp ../emirge.pbs .
		# Submit job
        qsub emirge.pbs
		# go to previous work directory
        cd -
done
```
## Process EMIRGE output
Reformat fasta header to sort sequences according to abundances. Usage: emirge_rename_fasta.py [options] <iter.DIR> > renamed.fasta
Make sure biopython/1.60 is loaded. Make sure that the files mapped.reads.txt, total.emirge.cons.fasta and reads.info.txt don't exist yet.
```
bash -x rename.sh
```

Count number of mapping reads from bowtie bam files in final iteration folders. This is needed to estimate the number of reads for each sequence

```
bash -x extract_mappings.sh
```

Concatenate all the fasta files in one fasta for taxonomic classification

```
bash -x concat.emirge.sh
```

## Classify the sequences according to the TaxASS pipeline for freshwater communities (<a href="https://github.com/McMahonLab/TaxAss">https://github.com/McMahonLab/TaxAss</a> )

### Step 1: Remove all redundant information from the first line in the fasta file. Retain sample names only. Fasta file must not be aligned.  

```
awk '{print $1}' total.emirge.renamed.fasta > total.emirge.renamed2.fasta
mv total.emirge.renamed2.fasta total.emirge.renamed.fasta
```
### Step 2: make BLAST database file (blast)

```
makeblastdb -dbtype nucl -in FreshTrain18Aug2016.fasta -input_type fasta -parse_seqids -out FWonly_18Aug2016custom.db
```

### Step 3: Run blast and reformat output blast file
```
blastn -query total.emirge.renamed.fasta -task megablast -db FWonly_18Aug2016custom.db -out custom.blast -outfmt 11 -max_target_seqs 5

blast_formatter -archive custom.blast -outfmt "6 qseqid pident length qlen qstart qend" -out otus.custom.blast.table
```
### Step 4: Correct BLAST pident using custom script
This accounts for sequence length differences

```
Rscript calc_full_length_pident.R otus.custom.blast.table otus.custom.blast.table.modified
```

### Step 5: Filter BLAST results

```
Rscript filter_seqIDs_by_pident.R otus.custom.blast.table.modified ids.above.97 97 TRUE

Rscript filter_seqIDs_by_pident.R otus.custom.blast.table.modified ids.below.97 97 FALSE
```
### Step 6: Make plots to evaluate blast run
```
mkdir plots

Rscript plot_blast_hit_stats.R otus.custom.blast.table.modified 97 plots
```
Plot is generated displaying the effect of correcting the Pid for length variations.
![blast_hits_used_for_pidents_0-100](https://cloud.githubusercontent.com/assets/19682548/18751230/3dc50a56-80ac-11e6-9cee-d68a58e69096.png)

### Step 7: recover sequence IDs left out of blast (python, bash)

```
python find_seqIDs_blast_removed.py total.emirge.renamed.fasta otus.custom.blast.table.modified ids.missing
cat ids.below.97 ids.missing > ids.below.97.all
```
### Step 8: create fasta files of desired sequence IDs (python)

```
python create_fastas_given_seqIDs.py ids.above.97 total.emirge.renamed.fasta otus.above.97.fasta

python create_fastas_given_seqIDs.py ids.below.97.all total.emirge.renamed.fasta otus.below.97.fasta
```
### Step 9: remove short, long sequences and those with ambiguous bases
```
mothur "#screen.seqs(fasta=otus.below.97.fasta,maxambig=0,minlength=1000,maxlength=1700)"
mothur "#screen.seqs(fasta=otus.above.97.fasta,maxambig=0,minlength=1000,maxlength=1700)"
```

### Step 10: extract unique sequences
```
mothur "#unique.seqs(fasta=otus.below.97.good.fasta)"
mothur "#unique.seqs(fasta=otus.above.97.good.fasta)"
```
### Step 11: classify sequences
```
mothur "#classify.seqs(fasta=otus.below.97.good.unique.fasta, template=/nfs/vdenef-lab/Shared/Ruben/databases_taxass/silva.nr_v123.align, taxonomy=/nfs/vdenef-lab/Shared/Ruben/databases_taxass/silva.nr_v123.tax, method=wang, probs=T, processors=10, cutoff=80)"

mothur "#classify.seqs(fasta=otus.above.97.good.unique.fasta, template=/nfs/vdenef-lab/Shared/Ruben/databases_taxass/FreshTrain18Aug2016.fasta,  taxonomy=/nfs/vdenef-lab/Shared/Ruben/databases_taxass/FreshTrain18Aug2016.taxonomy, method=wang, probs=T, processors=10, cutoff=80)"
```
### Step 12a: combine taxonomy files, names and fasta files
```
cat otus.above.97.good.unique.FreshTrain18Aug2016.wang.taxonomy otus.below.97.good.unique.nr_v123.wang.taxonomy > otus.97.taxonomy
cat otus.above.97.good.names otus.below.97.good.names > otus.97.names
cat otus.above.97.good.unique.fasta otus.below.97.good.unique.fasta > final.otu.emirge.fasta
```
### Step 12b: Run R script to create Sequence table (same format as OTU table)
```
Rscript emirge_format.R mapped.reads.txt otus.97.taxonomy read.info.txt otus.97.names
```
A histogram of the length distribution of your reconstructed sequencs is also automatically generated.
![length_distribution](https://cloud.githubusercontent.com/assets/19682548/18751134/ee91d996-80ab-11e6-8ca5-ecad0f6c7520.png)


### Extract taxonomy of interest
```
cat otus.97.taxonomy| sed 's/betI-A/betI_A/g' > otus.97.adjusted.taxonomy
```
#### in mothur:
```
get.lineage(taxonomy=otus.97.adjusted.taxonomy, taxon='Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;betI;betI_A-Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Comamonadaceae;Limnohabitans;',name=otus.97.names)
list.seqs(taxonomy=current)
get.seqs(accnos=current,fasta=final.otu.emirge.fasta)
```

### Add extra sequences from Czechs
```
cat reference.fasta final.otu.emirge.pick.fasta > final.otu.emirge.pick.merged.fasta
```
### Only retain first word on each sample header
```
awk '{print $1}' final.otu.emirge.pick.merged.fasta > total.emirge.renamed2.fasta
mv total.emirge.renamed2.fasta final.otu.emirge.pick.merged.fasta
```
### Remove non standard characters from sequence names
```
sed -e '/>/s/|/_/g' <final.otu.emirge.pick.merged.fasta >sequence.all.pcr.align.fna2
sed -e '/>/s/\.//g' <sequence.all.pcr.align.fna2 >final.otu.emirge.pick.merged.fasta
```
### Align sequences
```
align.seqs(seed=777,fasta=final.otu.emirge.pick.merged.fasta, reference=/scratch/vdenef_fluxm/rprops/databases/silva.seed_v123.align, processors=10,flip=T)
```
### Optional: cluster sequences 
```
cluster.split(fasta=final.otu.emirge.pick.align,taxonomy=otus.97.adjusted.pick.taxonomy, splitmethod=classify, taxlevel=4, cutoff=0.15, processors=20)
get.oturep(list=final.otu.emirge.pick.an.list,method=distance,fasta=final.otu.emirge.pick.align,name=otus.97.pick.names, label=0.03)
classify.otu(list=final.otu.emirge.pick.an.list, taxonomy=otus.97.adjusted.taxonomy,name=otus.97.pick.names, label=0.03)
```
### Make phylogenetic tree
```
FastTree -gtr -nt < final.otu.emirge.pick.merged.align > tree_file_NC_afterQC
```

### Use phyloassigner to place OTUs/oligotypes
```
module load gcc/4.8.5 openmpi/1.10.2/gcc/4.8.5 pplacer hmmer/3.0 bioperl/1.6.1 mothur phyloassigner/6.166
setupdb.pl tree_file_NC_afterQC final.otu.emirge.pick.merged.align "gi_33328318_gb_AF5304761_"
phyloassigner.pl final.otu.emirge.pick.merged.phyloassignerdb.tar.gz otu_oligo_total.fasta
```
