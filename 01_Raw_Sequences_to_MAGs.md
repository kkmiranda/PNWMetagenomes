# From Raw Sequences to MAGs

This is the workflow we used for the coastal metagenomes using Anvi'o 7.1. For the purpose of this tutorial, we are working just with the metagenome of Laminaria Setchellii (MG4). The same workflow was used for the other metagenomes. Only once all analyses were complete did we swap Reference ID for the more human-friendly Final ID.

### Metagenome Reference numbers:
| Reference ID  | Species                   | Tissue    | Final ID  |
|---------------|---------------------------|-----------|-----------|
| MG1           | Phyllospadix scouleri     | Sediment  |PSC_SED    |
| MG2           | Phyllospadix scouleri     | Blade     |PSC_BLD    |
| MG4           | Phyllospadix scouleri     | Rhizome   |PSC_RHZ    |
| MG4           | Laminaria setchellii      | Frond     |LSE_BLD    |
| MG5           | Nereocystis luetkeana     | Inner Bulb|NLE_BLB    |
| MG6           | Phyllospadix serrulatus   | Rhizome   |PSE_RHZ    |
| MG7           | Zostera marina            | Rhizome   |ZMA_RHZ    |
| MG8           | Zostera marina            | Sediment  |ZMA_SED    |

## 1. Quality Filtering
Begin with a tab-delimited file `lamsetch.txt` that points to `R1` and `R2` fasta files.

|sample     |r1                     |r2                 |
|-----------|-----------------------|-------------------|
|MG4        |MG4_S4_R1_001.fastq    |MG4_S4_R2_001.fastq|

Create a new directory for the quality-filtered `R1` and `R2` files, and then use `iu-gen-configs` to create config files for illumina-utils in it.
```
mkdir 01_QC
iu-gen-configs lamsetch.txt -o 01_QC
```
This above code created `MG4.ini` file. Forward and reverse raw sequences are filtered using `iu-filter-quality-minoche` (Minoche et al. 2011)

```
iu-filter-quality-minoche 01_QC/MG4.ini

# Stats Output should look like this
number of pairs analyzed      : 60666031
total pairs passed            : 48579975 (%80.08 of all pairs)
  total pair_1 trimmed        : 0 (%0.00 of all passed pairs)
  total pair_2 trimmed        : 925 (%0.00 of all passed pairs)
total pairs failed            : 12086056 (%19.92 of all pairs)
  pairs failed due to pair_1  : 910385 (%7.53 of all failed pairs)
  pairs failed due to pair_2  : 8620138 (%71.32 of all failed pairs)
  pairs failed due to both    : 2555533 (%21.14 of all failed pairs)
```

## 2. Assembly: From Quality Sequences to Contigs
Metagenomes were assembled into contigs using a snakemake pipeline in anvio, with the following settings for IDBA-UD (Peng et al. 2012):
```
 "idba_ud": {
        "--min_contig": 1000,
        "threads": 35,
        "run": true,
        "--mink": "",
        "--maxk": "",
        "--step": "",
        "--inner_mink": "",
        "--inner_step": "",
        "--prefix": "",
        "--min_count": "",
        "--min_support": "",
        "--seed_kmer": "",
        "--similar": "",
        "--max_mismatch": "",
        "--min_pairs": "",
        "--no_bubble": "",
        "--no_local": "",
        "--no_coverage": "",
        "--no_correct": "",
        "--pre_correction": "",
        "use_scaffolds": ""}
```
This generates the output file MG4-contigs.fa that we stored in a new directory `02_FASTA`. 

## 3. Cleaning contig fasta
To clean the contig names, we used `anvi-script-reformat-fasta` that also filtered only contigs with lengths > 2.5kbp. This is stored in the contigs directory `03_CONTIGS`
```
anvi-script-reformat-fasta 02_FASTA/MG4-contigs.fasta \
                            -o 03_CONTIGS/MG4-contigs.fa \
                            --min-len 2500 \
                            --simplify-names \
                            --report name_conversions.txt
```
This gave us an output like this
```
Input ........................................: 02_FASTA/MG4-contigs.fasta 
Output .......................................: 03_CONTIGS/MG4-contigs.fa
Minimum length ...............................: 2,500
Total num contigs ............................: 2,943,189
Total num nucleotides ........................: 1,089,641,695
Contigs removed ..............................: 2919142 (99.18% of all)
Nucleotides removed ..........................: 934366010 (85.75% of all)
Deflines simplified ..........................: True
```

## 4. Generate Anvi'o contigs database
Using the cleaned contigs fasta, we then used `anvi-gen-contigs-database` which computes k-mer frequency and identifies open-reading frames using Prodigal (Hyatt et al. 2011).
```
anvi-gen-contigs-database -f 03_CONTIGS/MG4-contigs.fa \
                            -o 03_CONTIGS/MG4-contigs.db
```

If you wish to view stats for the assembled contigs:
```
anvi-display-contigs-stats 03_CONTIGS/MG4-contigs.db
```

## 5. Mapping
Using bowtie, we map our short read sequences from our metagenome samples to our contigs assemblies from the previous step. First create a new directory for mapping
```
mkdir 04_MAPPING
```
Build an index for your contigs

```
bowtie2-build 03_CONTIGS/MG4-contigs.fa 04_MAPPING/lamsetch_contigs
```
Then run these commands to get an indexed BAM file for MG4
```
bowtie2 --threads 10 -x 04_MAPPING/lamsetch_contigs \
                        -1 01_QC/MG4-QUALITY_PASSED_R1.fastq \
                        -2 01_QC/MG4-QUALITY_PASSED_R2.fastq \
                        -S 04_MAPPING/lamsetch.sam
samtools view -F 4 -bS 04_MAPPING/lamsetch.sam > 04_MAPPING/lamsetch-RAW.bam
anvi-init-bam 04_MAPPING/lamsetch-RAW.bam -o 04_MAPPING/lamsetch.bam
```

## 6. Running HMMs
To decorate your contigs databases with HMM that uses multiple bacterial single-copy core gene collections to identify hits among your genes in `MG4-contigs.db` to these collections.
```
anvi-run-hmms -c 03_CONTIGS/MG4-contigs.db
```
The output should look like this~
```
Target found .................................: RNA:CONTIG                                                                   
Target found .................................: AA:GENE                                                                      
                                                                                                                             
HMM Profiling for Ribosomal_RNAs
===============================================
Reference ....................................: Seemann T, https://github.com/tseemann/barrnap
Kind .........................................: Ribosomal_RNAs
Alphabet .....................................: RNA
Context ......................................: CONTIG
Domain .......................................: N\A
HMM model path ...............................: /software/Anaconda3-2019.03-el7-x86_64/envs/anvio-6.1/lib/python3.6/site-packages/anvio/data/hmm/Ribosomal_RNAs/genes.hmm.gz
Number of genes ..............................: 12
Noise cutoff term(s) .........................: --cut_ga
Number of CPUs will be used for search .......: 1
Temporary work dir ...........................: /tmp/tmpl8tkgr8y
Log file .....................................: /tmp/tmpl8tkgr8y/00_log.txt
Number of raw hits ...........................: 11                                                                           
Pruned .......................................: 3 out of 11 hits were removed due to redundancy
Gene calls added to db .......................: 8 (from source "Ribosomal_RNAs")                                             

HMM Profiling for Protista_83
===============================================
Reference ....................................: Delmont, http://merenlab.org/delmont-euk-scgs
Kind .........................................: singlecopy
Alphabet .....................................: AA
Context ......................................: GENE
Domain .......................................: eukarya
HMM model path ...............................: /software/Anaconda3-2019.03-el7-x86_64/envs/anvio-6.1/lib/python3.6/site-packages/anvio/data/hmm/Protista_83/genes.hmm.gz
Number of genes ..............................: 83
Noise cutoff term(s) .........................: -E 1e-25
Number of CPUs will be used for search .......: 1
Temporary work dir ...........................: /tmp/tmpg78iqwey
Log file .....................................: /tmp/tmpg78iqwey/00_log.txt
Number of raw hits ...........................: 139                                                                          

HMM Profiling for Bacteria_71
===============================================
Reference ....................................: Lee modified, https://doi.org/10.1093/bioinformatics/btz188
Kind .........................................: singlecopy
Alphabet .....................................: AA
Context ......................................: GENE
Domain .......................................: bacteria
HMM model path ...............................: /software/Anaconda3-2019.03-el7-x86_64/envs/anvio-6.1/lib/python3.6/site-packages/anvio/data/hmm/Bacteria_71/genes.hmm.gz
Number of genes ..............................: 71
Noise cutoff term(s) .........................: --cut_ga
Number of CPUs will be used for search .......: 1
Temporary work dir ...........................: /tmp/tmpg78iqwey
Log file .....................................: /tmp/tmpg78iqwey/00_log.txt
Number of raw hits ...........................: 2,672                                                                        

HMM Profiling for Archaea_76
===============================================
Reference ....................................: Lee, https://doi.org/10.1093/bioinformatics/btz188
Kind .........................................: singlecopy
Alphabet .....................................: AA
Context ......................................: GENE
Domain .......................................: archaea
HMM model path ...............................: /software/Anaconda3-2019.03-el7-x86_64/envs/anvio-6.1/lib/python3.6/site-packages/anvio/data/hmm/Archaea_76/genes.hmm.gz
Number of genes ..............................: 76
Noise cutoff term(s) .........................: --cut_ga
Number of CPUs will be used for search .......: 1
Temporary work dir ...........................: /tmp/tmpg78iqwey
Log file .....................................: /tmp/tmpg78iqwey/00_log.txt
Number of raw hits ...........................: 1,371           
```

## Generating Anvio Profile Databases
Anvi-Profile databases store sample specific information about the corresponding contig. This takes generates a single profile off the BAM file you generated above that reports properties for each contigs based on your mapping results. We also set `min-contig-length` to 1 kbp which is lower than the default 2.5 kbp. This casts a more generous filter by processing contig sequences > 1 kbp. As we were using single profiles downstream, we use the flag `--cluster-contigs`
```
mkdir 05_ANVIO_PROFILE/MG4
anvi-profile -i 04_MAPPING/lamsetch.bam \
                -c 03_CONTIGS/MG4-contigs.db \
                -o 05_ANVIO_PROFILE/MG4/ \
                --cluster-contigs \
                --min-contig-length 1000 \
```

## Binning
Binning is an important step where we group sequences into clusters that correspond to their taxonomic unit. We employed a hybrid binning strategy where we used CONCOCT (Alneburg et al. 2014), an automatic binning software, to cluster our sequences into a maximum of 100 bins. 
### Concoct automatic binning
Refer to <a href='https://concoct.readthedocs.io/en/latest/installation.html'> CONCOCT installation documentation</a> to first download the software.

Then we enter the concoct environment and generate a new directory for concoct output files. 
```
source activate concoct
mkdir 06_CONCOCT
```
We then cut up the fasta into chunks of 10 kbp using the `-c` flag and generate a BEDfile. This gives more weight to larger contigs which mitigates the effect of assembly errors.
```
cut_up_fasta.py 03_CONTIGS/MG4-contigs.fa \
                -c 10000 \
                -o 0 \
                --merge_last \ 
                -b 06_CONCOCT/MG4_contigs_10K.bed > 06_CONCOCT/MG4_contigs_10K.fa
```
Generate Table with coverage depth information per sample and subcontig. This step takes the BAMfiles generated in the Mapping step and the BEDfile generated above.
```
concoct_coverage_table.py 06_CONCOCT/MG4_contigs_10K.bed 04_MAPPING/MG4.bam > MG4_coverage_table.tsv
```
Run CONCOCT with max number of clusters `-c` set to 100. `merge_cutup_clustering.py` merges subcontig clustering into original contig clustering.
```
concoct --coverage_file 06_CONCOCT/MG4_coverage_table.tsv \
        --composition_file 06_CONCOCT/MG4_contigs_10K.fa \
        -c 100 \
        -l 1000 \
        -b 06_CONCOCT/

merge_cutup_clustering.py 06_CONCOCT/clustering_gt1000.csv > 06_CONCOCT/MG4_concoct.csv
```
Reformat the csv to a tab delimited text file that can be imported into an anvio-profile.
```
06_CONCOCT/MG4_concoct.csv > 06_CONCOCT/MG4_concoct.txt  
sed -i 1d 06_CONCOCT/MG4_concoct.txt 
awk '$2="Bin_"$2' 06_CONCOCT/MG4_concoct.txt > 06_CONCOCT/MG4_concoct_wBin.txt
awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' 06_CONCOCT/MG4_concoct_wbin.txt > MG4_concoct.txt
```
Finally, import your collection into your anvio-profile under the collection name `MG4_CONCOCT`
```
anvi-import-collection 06_CONCOCT/MG4_concoct.txt \
                       -c 03_CONTIGS/MG4-contigs.db \
                       -p 05_ANVIO_PROFILE/MG4/PROFILE.db \
                       -C MG4_CONCOCT \
                       --contigs-mode
```
To see all the collections and bins, use:
```
anvi-show-collections-and-bins -p 05_ANVIO_PROFILE/MG4/PROFILE.db
```
To enter the anvio-interface with this collection pulled up:
```
anvi-interactive -c 03_CONTIGS/MG4-contigs.db \
                    -p 05_ANVIO_PROFILE/MG4/PROFILE.db \
                    -C MG4_CONCOCT --server-only -P 8080
```

### Manual Binning
Now we had 100 bins generated by CONCOCT that we could go into and manually refine one at a time based on sequence composition. This is the only step in the workflow that is not automated and relies on individual judgement. However, to ensure optimal reproducibility, we followed these resources:
- https://merenlab.org/tutorials/infant-gut/#manual-identification-of-genomes-in-the-infant-gut-dataset
- https://merenlab.org/data/tara-oceans-mags/
- https://merenlab.org/2017/01/03/loki-the-link-archaea-eukaryota/ 
- https://merenlab.org/2016/06/09/assessing-completion-and-contamination-of-MAGs/

Additional binning trivia:
- To aid in the human-guided binning, sometimes, the taxonomic makeup of each contig helps. For this we used `anvi-run-scg-taxonomy` on our `MG4-contigs.db` which annotates each gene using the Genome Taxonomy Database (GTDB) (Parks et al. 2018). 
- We used `anvi-interactive` as the interface for binning and `anvi-refine` to go into a single bin and refine it.
- The CONCOCT algorithm attempted to find at max 100 clusters but the output shows clusters completely empty. This is normal and these bins were immediately discarded.
- We classified MAGs as High Quality (<b>Completion > 90% and Redundancy < 10%</b>) or Low Quality (everything else).
- Binning could be an infinitely long procedure, so once any effort made didn't translate to an improved Completion:Redundancy ratio, binning on that bin was stopped.

## References
Minoche et al. 2011
Peng et al. 2012
Hyatt et al. 2011
Alneburg et al. 2014
Parks et al. 2018
