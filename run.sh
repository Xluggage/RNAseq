# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/bioinfo/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/bioinfo/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/bioinfo/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/bioinfo/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
# This is used to make conda run in the script without errors. 
# Please replace "/home/bioinfo/miniconda3/" with the installation path of conda.

dataDir=${1/.\//}
# The received parameter is the fastq data address (relative position), and replace the parameter's "./" with an empty string.
# e.g. the received parameter is "./test/bone", dataDir is "test/bone"
# The "dataDir" will be part of the folder structure for the output results ("./result/${dataDir}/"), with the last level serving as the prefix for all output file names.

conda activate Expression
# source QC_mapping_featureCounts.sh $dataDir
# If you are using raw data instead of clean data, remove the "#" from the previous line and comment out the next line.
source mapping_featureCounts.sh $dataDir
# alignment and quantification
conda deactivate

conda activate Differential
Rscript mappingPlot.R ${dataDir}
# visualize the statistics for the bam files
Rscript DESeq2.R $dataDir
# differential expression analysis. if the sample number in any condition is less than 3, comment out the previous line and remove the "#" from the next line.
# Rscript edgeR.R $dataDir
# please put the group file "colData.tsv" in the fastq data folder.
conda deactivate

conda activate Enrichment
Rscript enrichment.R $dataDir
# enrichment analysis. If you are using edgeR instead of DESeq2, comment out the previous line and remove the "#" from the next line.
# Rscript enrichment_edgeR.R $dataDir
conda deactivate
