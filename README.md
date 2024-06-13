# RNAseq
The shell and R scripts used to RNAseq. 

## Software installation
The scripts are running in Ubuntu20.04.5LTS and conda environments.

The scripts aligned sequencing data to the genome using **HISAT2**, calculated expression counts using **featureCounts**, performed differential expression analysis using **DESeq2** and **edgeR**, and used **clusterProfiler** for enrichment analysis.

For the specific software installation process, please refer to `INSTALL.txt`.

## Files preparation
The mouse reference genome (Mu_musculus.GRCm39.dna.primary_assembly.fa.gz) and annotation file (Mus_musculus.GRCm39.108.gtf.gz) were downloaded from [ensembl database](https://asia.ensembl.org/info/data/ftp/index.html), and they need to be in the `./genome/` folder.

The annotation information `./genome/geneAnnotation.csv` was merged from ensembl's BioMart and [MGI](https://www.informatics.jax.org/batch) on 2022.11.12

The files in  `./genome/Pathway` were selected from GO (download [mgi.gaf.gz](http://current.geneontology.org/products/pages/downloads.html), gaf-version 2.2, on 2023.02.06) and KEGG (download by `clusterProfiler::download_KEGG` on 2023.02.17), as were the gene sets in `./genome/HeatmapGene`.

Before using HISAT2 alignment, you need to build an index as follows:
```
samtools faidx ./genome/Mus_musculus.GRCm39.dna.primary_assembly.fa
hisat2_extract_exons.py ./genome/Mus_musculus.GRCm39.108.gtf > ./genome/genome/genome.exon
hisat2_extract_splice_sites.py ./genome/Mus_musculus.GRCm39.108.gtf > ./genome/genome/genome.ss
hisat2-build -p 48 --exon ./genome/genome/genome.exon --ss ./genome/genome/genome.ss ./genome/Mus_musculus.GRCm39.dna.primary_assembly.fa ./genome/genome/genome > run_build_tran.log 2>&1
```
## Data preparation
Place the sequencing data in a position like `./.../groupName`, and add a condition file `colData.tsv` to it, the condition file refers to `colData.tsv`, please put the control group in front of the conditions.

## Start running
Modify the `run.sh` as needed, and then run (for example, the data is in `./anyone/bone`):
```
bash run.sh ./anyone/bone
```
After that, you will get the following result：
```
 ./
   genome/                  reference files
   anyone/
     bone/                  data group
       colData.tsv          the condition file
       KO1_1.clean.fq.gz    sequencing data
       KO2_2.clean.fq.gz    sequencing data
       ...
   result/                  results folder
     anyone/
       bone/                group
         bam/               the mapping files (*.bam)
         Mapping/           the statistics for the bam files
         Quantification/    gene expression files
         Differential/      the differential expression analysis results
         Enrichment/        the enrichment analysis results
```
