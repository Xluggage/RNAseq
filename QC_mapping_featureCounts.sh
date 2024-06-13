genome_hisat2=./genome/genome/genome
# Index of reference genome, where "./genome/genome/" indicates the directory location and "genome" is a prefix for the index file.
gtf=./genome/Mus_musculus.GRCm39.108.gtf
genome=./genome/Mus_musculus.GRCm39.dna.primary_assembly.fa

dataDir=$1
# The received parameter is the fastq data address (relative position), and replace the parameter's "./" with an empty string.

QCStatDir="./result/${dataDir}/QC/"
# the output folder of the statistics for quality control

cleanDataDir="./result/${dataDir}/clean/"
# the output folder of clean data

mappingStatDir="./result/${dataDir}/Mapping/"
# the output folder of the statistics for the bam files

mappingDir="./result/${dataDir}/bam/"
# the output folder of bam files

quantifyDir="./result/${dataDir}/Quantification/"
# the output folder of gene expression level files

counts="`basename ${dataDir}`_counts_gene.txt"
# the gene expression level files names

[ ! -d ${QCStatDir} ] && mkdir -p ${QCStatDir}
[ ! -d ${cleanDataDir} ] && mkdir -p ${cleanDataDir}
[ ! -d ${mappingStatDir} ] && mkdir -p ${mappingStatDir}
[ ! -d ${mappingDir} ] && mkdir -p ${mappingDir}
[ ! -d ${quantifyDir} ] && mkdir -p ${quantifyDir}

for i in `ls ${dataDir}/*_1.fq*`
do
    simpleName=`basename $i`
    simpleName=${simpleName%_1.fq*}
    fastp -w 12 -n 0 -q 5 -u 50 -x A --detect_adapter_for_pe -c -p -R "${simpleName}" --dont_overwrite -i ${dataDir}/${simpleName}_1.fq.gz -I ${dataDir}/${simpleName}_2.fq.gz -o ${cleanDataDir}/${simpleName}_1.clean.fq.gz -O ${cleanDataDir}/${simpleName}_2.clean.fq.gz -h ${QCStatDir}/${simpleName}_fastp.html -j ${QCStatDir}/${simpleName}_fastp.json 2> ${QCStatDir}/${simpleName}_fastp.txt
    hisat2 -p 12 --time --summary-file ${mappingStatDir}${simpleName}.log --new-summary --novel-splicesite-outfile ${mappingStatDir}${simpleName}.splice -x ${genome_hisat2} -1 ${cleanDataDir}/${simpleName}_1.clean.fq.gz -2 ${cleanDataDir}/${simpleName}_2.clean.fq.gz | samtools view -@ 8 -bh | samtools sort -@ 8 -O bam > ${mappingDir}${simpleName}.bam
    samtools index ${mappingDir}${simpleName}.bam
done

featureCounts -T 12 -p --extraAttributes gene_name,gene_biotype -J -G ${genome} -a ${gtf} -o ${quantifyDir}${counts} `ls ${mappingDir}*.bam`
# After version 2.0.2, the -p parameter only indicates that the data is paired-end, and --countReadPairs must also be specified
