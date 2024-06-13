dataDir <- commandArgs(trailingOnly = TRUE)

sample <- basename(dataDir)
# the prefix for all output file names

savePath <- paste0("result/", dataDir, "/Differential/")
# the output folder

if (!dir.exists(savePath)) {
  dir.create(savePath, recursive = TRUE)
}

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(apeglm))
suppressPackageStartupMessages(library(IHW))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))

# read files ------------------------------------------------------------------

counts <-
  read.table(
    paste0(
      "result/",
      dataDir,
      "/Quantification/",
      basename(dataDir),
      "_counts_gene.txt"
    ),
    sep = "\t",
    header = T,
    row.names = 1
  )
# read the gene expression level files

colnames(counts) %<>% gsub("\\.[sb]am$", "", .) %>% gsub("^.*\\.", "", .)
# remove the prefix (path) and suffix (.bam) in the sample names

countsAnnotation <- counts[1:7]
# extract the annotation information

colData <-
  read.table(
    file.path(dataDir, "colData.tsv"),
    sep = "\t",
    header = T,
    row.names = 1
  )
# read the group information file

ensToGene <-
  read.delim("./genome/geneAnnotation.tsv", check.names = FALSE)
# read another annotation information

listHeatmapGene_GO <- list.files(
  path = paste0("./genome/HeatmapGene/GO/"),
  pattern = ".tsv$",
  full.names = TRUE
) %>% magrittr::set_names(., stringr::str_replace(basename(.), "\\.tsv$", "")) %>% lapply(., read.delim, head = FALSE)

listHeatmapGene_KEGG <- list.files(
  path = paste0("./genome/HeatmapGene/KEGG/"),
  pattern = ".tsv$",
  full.names = TRUE
) %>% magrittr::set_names(., stringr::str_replace(basename(.), "\\.tsv$", "")) %>% lapply(., read.delim, head = FALSE)
# read the function tables to draw heatmaps

if (!all(rownames(colData) %in% colnames(counts)))
{
  cat("The group is not in the gene expression files!\n")
  quit(save = "yes", status = 1)
} else  if (!all(rownames(colData) == colnames(counts)))
{
  cat(
    "The order of group is different from the expression files, now reorder the expression files. \n"
  )
  counts <- counts[, rownames(colData)]
}
if (!all(rownames(colData) == colnames(counts)))
{
  cat("The group is different from the expression files, quit!\n")
  quit(save = "yes", status = 1)
} else {
  cat("Start analysis...", sep = "\n")
}


# Function definition ---------------------------------------------------------

draw_PCA <- function(rlogData,
                     sampleName = sample,
                     saveDir = "./") {
  if (!dir.exists(saveDir)) {
    dir.create(saveDir, recursive = TRUE)
  }
  
  suppressPackageStartupMessages(library("DESeq2"))
  suppressPackageStartupMessages(library("ggplot2"))
  suppressPackageStartupMessages(library("ggrepel"))
  suppressPackageStartupMessages(library("RColorBrewer"))
  
  pcaData <- plotPCA(rlogData, ntop = 1000, returnData = TRUE)
  
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  p <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3, stroke = 0) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    geom_label_repel(
      aes(label = name),
      force = 10,
      min.segment.length = 0,
      max.overlaps = 30,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c("#4575B4", "#D73027"),
                       breaks = c(levels(pcaData$condition)[1], levels(pcaData$condition)[2])) +
    theme_classic() +
    expand_limits(x = c(
      min(pcaData$PC1, pcaData$PC2),
      max(pcaData$PC1, pcaData$PC2)
    ), y = c(
      min(pcaData$PC1, pcaData$PC2),
      max(pcaData$PC1, pcaData$PC2)
    )) +
    coord_fixed()
  
  ggsave(
    file = paste0(sampleName, "_PCA.png"),
    path = saveDir,
    plot = p,
    width = 8,
    height = 6,
    units = "in",
    dpi = 600
  )
  ggsave(
    file = paste0(sampleName, "_PCA.tiff"),
    path = saveDir,
    plot = p,
    device = "tiff",
    compression = "lzw",
    width = 8,
    height = 6,
    units = "in",
    dpi = 600
  )
  ggsave(
    file = paste0(sampleName, "_PCA.svg"),
    path = saveDir,
    plot = p,
    width = 8,
    height = 6,
    units = "in",
  )
  ggsave(
    file = paste0(sampleName, "_PCA.eps"),
    path = saveDir,
    plot = p,
    width = 8,
    height = 6,
    units = "in",
  )
  
}

draw_volcano <-
  function(result,
           cutoff = 1,
           labelNum = 10,
           adjustP = "padj",
           sampleName = sample,
           saveDir = "./") {
    # cutoff:     the range of difference values for the color tags
    # labelNum:   number of labels (Up/Down are counted separately, and the first n log2FoldChange labels are taken)
    # adjustP:    the column name of adjust P value (ihwP or padj)
    
    if (!dir.exists(saveDir)) {
      dir.create(saveDir, recursive = TRUE)
    }
    
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("ggrepel"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("magrittr"))
    suppressPackageStartupMessages(library("RColorBrewer"))
    
    if (adjustP == "ihwP") {
      result <- result[!is.na(result$ihwP) & result$ihwP != 1,]
    } else{
      result <- result[!is.na(result$padj), ]
    }
    
    result$change = as.factor(ifelse(
      result[, adjustP] < 0.05 &
        abs(result$log2FoldChange) > cutoff,
      ifelse(result$log2FoldChange > cutoff , "UP", "DOWN"),
      "STABLE"
    ))
    # set markers of up and down
    
    # · no label --------------------------------------------------------------
    
    p <-
      ggplot(data = result, aes(
        x = log2FoldChange,
        y = -log10(get(adjustP)),
        color = change
      )) +
      geom_point(size = 2, stroke = 0) +
      geom_vline(
        xintercept = c(cutoff, -cutoff),
        colour = "gray20",
        linetype = 3
      ) +
      geom_hline(
        yintercept = -log10(0.05),
        colour = "gray20",
        linetype = 3
      ) +
      theme_classic() +
      theme(legend.justification = c("right", "top"),
            legend.title = element_blank()) +
      xlab(expression(paste(log[2], "(fold change)"))) +
      ylab(expression(paste(-log[10], "(adjusted p-value)"))) +
      scale_colour_manual(values = c(
        UP = "#D7302780",
        DOWN = "#4575B480",
        STABLE = "#0000004C"
      ))
    ggsave(
      file = paste0(sampleName,
                    "_volcano.png"),
      path = saveDir,
      plot = p,
      width = 8,
      height = 6,
      units = "in",
      dpi = 600
    )
    ggsave(
      file = paste0(sampleName,
                    "_volcano.tiff"),
      path = saveDir,
      plot = p,
      device = "tiff",
      compression = "lzw",
      width = 8,
      height = 6,
      units = "in",
      dpi = 600
    )
    ggsave(
      file = paste0(sampleName,
                    "_volcano.svg"),
      path = saveDir,
      plot = p,
      width = 8,
      height = 6,
      units = "in",
    )
    ggsave(
      file = paste0(sampleName,
                    "_volcano.eps"),
      path = saveDir,
      plot = p,
      device = cairo_ps,
      fallback_resolution = 600,
      width = 8,
      height = 6,
      units = "in",
    )
    
    # · have label ------------------------------------------------------------
    
    if ("Gene name" %in% colnames(result)) {
      result %<>% mutate(Symbol = if_else(is.na(Symbol) |
                                            Symbol == "", `Gene name`, Symbol))
    }
    result %<>% mutate(Symbol = if_else(is.na(Symbol) |
                                          Symbol == "", row.names(.), Symbol))
    
    this_tile <- paste0(
      "Cutoff for logFC is ",
      cutoff,
      "\nThe number of up gene is ",
      nrow(result[result$change == "UP",]) ,
      "\nThe number of down gene is ",
      nrow(result[result$change == "DOWN",])
    )
    
    p2 <- p %+%
      ggtitle(this_tile) +
      geom_label_repel(
        data = {
          result[which(result$change != "STABLE"), ] %>% group_by(change) %>% slice_max(order_by = abs(log2FoldChange), n = labelNum)
        },
        aes(label = Symbol),
        force = 10,
        min.segment.length = 0,
        max.overlaps = 30,
        show.legend = FALSE
      )
    ggsave(
      file = paste0(sampleName,
                    "_volcano_label.png"),
      path = saveDir,
      plot = p2,
      width = 8,
      height = 6,
      units = "in",
      dpi = 600
    )
    ggsave(
      file = paste0(sampleName,
                    "_volcano_label.tiff"),
      path = saveDir,
      plot = p2,
      device = "tiff",
      compression = "lzw",
      width = 8,
      height = 6,
      units = "in",
      dpi = 600
    )
    ggsave(
      file = paste0(sampleName,
                    "_volcano_label.svg"),
      path = saveDir,
      plot = p2,
      width = 8,
      height = 6,
      units = "in",
    )
    ggsave(
      file = paste0(sampleName,
                    "_volcano_label.eps"),
      path = saveDir,
      plot = p2,
      device = cairo_ps,
      fallback_resolution = 600,
      width = 8,
      height = 6,
      units = "in",
    )
  }

draw_heatmap <-
  function(rlogData,
           diffGene,
           suffix = NULL,
           coldata = colData,
           sampleName = sample,
           saveDir = "./") {
    # suffix:     used to distinguish between result files to avoid duplicate names
    
    saveDir = paste0(saveDir, "/")
    
    if (!dir.exists(saveDir)) {
      dir.create(saveDir, recursive = TRUE)
    }
    
    if (!is.null(suffix)) {
      suffix = paste0("_", suffix)
    }
    
    suppressPackageStartupMessages(library("DESeq2"))
    suppressPackageStartupMessages(library("pheatmap"))
    suppressPackageStartupMessages(library("RColorBrewer"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("magrittr"))
    
    if ("Gene name" %in% colnames(diffGene)) {
      diffGene %<>% mutate(Symbol = if_else(is.na(Symbol) |
                                              Symbol == "", `Gene name`, Symbol))
    }
    diffGene %<>% mutate(Symbol = if_else(is.na(Symbol) |
                                            Symbol == "", row.names(.), Symbol))
    
    countsMatrix <- assay(rlogData) %>% .[row.names(diffGene),]
    
    p <- pheatmap::pheatmap(
      countsMatrix,
      scale = "row",
      color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(1000),
      border_color = NA,
      cellwidth = 30,
      cellheight = 10,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      treeheight_row = 30,
      treeheight_col = 20,
      show_rownames = TRUE,
      show_colnames = TRUE,
      labels_row = diffGene$Symbol,
      angle_col = 0,
      annotation_col = coldata,
      annotation_names_col = FALSE,
      legend = TRUE,
      legend_breaks = c(-2, 0, 2),
      filename = paste0(saveDir,
                        sampleName,
                        "_heatmap",
                        suffix,
                        ".pdf"),
      width = 8,
      # height = 6,
      units = "in"
    )
    
    ht_width <-
      grid::convertWidth(gtable::gtable_width(p$gtable), "inch")
    ht_height <-
      grid::convertHeight(gtable::gtable_height(p$gtable), "inch")
    # extract the width and height of heatmap object
    
    ggplot2::ggsave(
      file = paste0(sampleName,
                    "_heatmap",
                    suffix,
                    ".svg"),
      path = saveDir,
      plot = p,
      width = ht_width,
      height = ht_height,
      units = "in",
      limitsize = FALSE
    )
    
    ggplot2::ggsave(
      file = paste0(sampleName,
                    "_heatmap",
                    suffix,
                    ".eps"),
      path = saveDir,
      plot = p,
      width = ht_width,
      height = ht_height,
      units = "in",
      limitsize = FALSE
    )
    
    suppressPackageStartupMessages(library("pdftools"))
    
    bitmap <- pdf_render_page(
      paste0(saveDir,
             sampleName,
             "_heatmap",
             suffix,
             ".pdf"),
      page = 1,
      dpi = 300,
      numeric = TRUE
    )
    
    png::writePNG(bitmap,
                  paste0(saveDir,
                         sampleName,
                         "_heatmap",
                         suffix,
                         ".png"))
    
    tiff::writeTIFF(bitmap,
                    paste0(saveDir,
                           sampleName,
                           "_heatmap",
                           suffix,
                           ".tiff"),
                    compression = "LZW")
  }

draw_heatmap2 <-
  function(rlogData,
           diffGene,
           suffix = NULL,
           coldata = colData,
           sampleName = sample,
           saveDir = "./") {
    # suffix:     used to distinguish between result files to avoid duplicate names
    
    saveDir = paste0(saveDir, "/")
    
    if (!dir.exists(saveDir)) {
      dir.create(saveDir, recursive = TRUE)
    }
    
    if (!is.null(suffix)) {
      suffix = paste0("_", suffix)
    }
    
    suppressPackageStartupMessages(library("DESeq2"))
    suppressPackageStartupMessages(library("ComplexHeatmap"))
    suppressPackageStartupMessages(library("circlize"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("magrittr"))
    
    if ("Gene name" %in% colnames(diffGene)) {
      diffGene %<>% dplyr::mutate(Symbol = if_else(is.na(Symbol) |
                                                     Symbol == "", `Gene name`, Symbol))
    }
    diffGene %<>% dplyr::mutate(Symbol = if_else(is.na(Symbol) |
                                                   Symbol == "", row.names(.), Symbol))
    
    diffGene$symP <- symnum(
      diffGene$ihwP,
      corr = FALSE,
      cutpoints = c(0,  .001, .01, .05, 1),
      symbols = c("***", "**", "*", " ")
    )
    # Statistically significant extent of the P-value
    
    countsMatrix <- assay(rlogData) %>% .[row.names(diffGene),]
    
    countsMatrix <-
      countsMatrix[, row.names(coldata)] %>% t(.) %>% scale(.) %>% t(.)
    
    ha <- HeatmapAnnotation(
      df = coldata,
      col = list(condition = c("royalblue", "pink") %>% setNames(c(
        levels(coldata$condition)[1], levels(coldata$condition)[2]
      ))),
      show_annotation_name = FALSE
    )
    # the annotation of group
    
    if (sign(max(diffGene$log2FoldChange)) == sign(min(diffGene$log2FoldChange))) {
      Log2FCWidth = unit(1, "inch")
    } else{
      Log2FCWidth = unit(max(0.375 * (
        max(diffGene$log2FoldChange) - min(diffGene$log2FoldChange)
      ) /
        min(c(
          max(diffGene$log2FoldChange),
          -min(diffGene$log2FoldChange)
        )) + 0.125, 1), "inch")
    }
    # calculate the width of the log2FoldChange bar
    
    ha2 <- HeatmapAnnotation(
      pvalue = anno_simple(
        -log10(diffGene$ihwP),
        col = colorRamp2(
          breaks = c(0, -log10(0.05), 3),
          colors = c("#4575B4", "white", "#D73027")
        ),
        pch = diffGene$symP
      ),
      Log2FC = anno_numeric(
        diffGene$log2FoldChange,
        bg_gp = gpar(
          fill = if (sign(max(diffGene$log2FoldChange)) != sign(min(diffGene$log2FoldChange)))
            c("#4575B4", "#D73027")
          else
            ifelse(sign(max(
              diffGene$log2FoldChange
            )) == -1, "#4575B4", "#D73027")
        ),
        labels_format =  function(x)
          sprintf("%.2f", x),
        align_to = ifelse(
          sign(max(diffGene$log2FoldChange)) == sign(min(diffGene$log2FoldChange)),
          ifelse(sign(max(
            diffGene$log2FoldChange
          )) == -1, "right", "left"),
          0
        ),
        width = Log2FCWidth
      ),
      annotation_name_rot = 0,
      which = "row"
    )
    
    lgd_pvalue <- Legend(
      title = "p-value",
      col_fun = colorRamp2(
        breaks = c(0, -log10(0.05), 3),
        colors = c("#4575B4", "white", "#D73027")
      ),
      at = c(0, -log10(0.05), 3),
      labels = c("1", "0.05", "0.001")
    )
    # the annotations of P-value significance and difference value bar
    
    ht <- Heatmap(
      countsMatrix,
      name = "Expression",
      col = colorRamp2(c(-2, 0, 2), c("#4575B4", "white", "#D73027")),
      cluster_columns = FALSE,
      show_row_dend = FALSE,
      column_names_rot = 0,
      column_names_centered = TRUE,
      row_names_side = "left",
      row_labels = diffGene$Symbol,
      top_annotation = ha,
      right_annotation = ha2,
      height = unit(15, "pt") * nrow(countsMatrix),
      width = unit(45, "pt") * nrow(coldata)
    )
    # heatmap object
    
    dht <-
      invisible(draw(ht, annotation_legend_list = list(lgd_pvalue)))
    ht_height = sum(component_height(dht)) + unit(4, "mm")
    ht_width = sum(component_width(dht)) + unit(4, "mm")
    ht_height = max(convertHeight(ht_height, "inch", valueOnly = TRUE),
                    unit(2, "inch"))
    ht_width = convertWidth(ht_width, "inch", valueOnly = TRUE)
    dev.off()
    # extract the width and height of heatmap object
    
    pdf(
      file = paste0(saveDir,
                    sampleName,
                    "_heatmap",
                    suffix,
                    ".pdf"),
      width = ht_width,
      height = ht_height,
      bg = "white"
    )
    draw(ht, annotation_legend_list = list(lgd_pvalue))
    dev.off()
    
    setEPS()
    postscript(
      paste0(saveDir,
             sampleName,
             "_heatmap",
             suffix,
             ".eps"),
      width = ht_width,
      height = ht_height
    )
    draw(ht, annotation_legend_list = list(lgd_pvalue))
    dev.off()
    
    suppressPackageStartupMessages(library("pdftools"))
    bitmap <- pdf_render_page(
      paste0(saveDir,
             sampleName,
             "_heatmap",
             suffix,
             ".pdf"),
      page = 1,
      dpi = 300,
      numeric = TRUE
    )
    png::writePNG(bitmap,
                  paste0(saveDir,
                         sampleName,
                         "_heatmap",
                         suffix,
                         ".png"))
    tiff::writeTIFF(bitmap,
                    paste0(saveDir,
                           sampleName,
                           "_heatmap",
                           suffix,
                           ".tiff"),
                    compression = "LZW")
  }

draw_distance <-
  function(rlogData,
           sampleName = sample,
           saveDir = "./") {
    saveDir = paste0(saveDir, "/")
    if (!dir.exists(saveDir)) {
      dir.create(saveDir, recursive = TRUE)
    }
    
    suppressPackageStartupMessages(library("pheatmap"))
    suppressPackageStartupMessages(library("RColorBrewer"))
    
    sampleDists <- dist(t(assay(rlogData)))
    # calculate the distance between samples
    
    sampleDistMatrix <- as.matrix(sampleDists)
    # convert the dist object to a matrix
    
    rownames(sampleDistMatrix) <- rlogData$condition
    
    colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
    
    p <- pheatmap::pheatmap(
      sampleDistMatrix,
      clustering_distance_rows = sampleDists,
      clustering_distance_cols = sampleDists,
      color = colors,
      angle_col = 0,
      filename = paste0(saveDir, sampleName, "_sampleDistances.png")
    )
    
    tiff(
      paste0(saveDir, sampleName, "_sampleDistances.tiff"),
      width = 7,
      height = 7,
      units = "in",
      res = 300,
      compression = "lzw"
    )
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
    dev.off()
    
    setEPS()
    postscript(paste0(saveDir, sampleName, "_sampleDistances.eps"))
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
    dev.off()
  }

draw_cook <- function(x = dds,
                      sampleName = sample,
                      saveDir = "./") {
  saveDir = paste0(saveDir, "/")
  if (!dir.exists(saveDir)) {
    dir.create(saveDir, recursive = TRUE)
  }
  
  png(
    filename = paste0(saveDir, sampleName, "_cook.png"),
    width = 8,
    height = 6,
    units = "in",
    bg = "white",
    res = 300
  )
  par(mar = c(8, 5, 2, 2))
  boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2)
  dev.off()
}

draw_heatmapCounts <-
  function(rlogData,
           x = dds,
           sampleName = sample,
           saveDir = "./") {
    saveDir = paste0(saveDir, "/")
    if (!dir.exists(saveDir)) {
      dir.create(saveDir, recursive = TRUE)
    }
    
    suppressPackageStartupMessages(library("pheatmap"))
    suppressPackageStartupMessages(library("RColorBrewer"))
    
    select <-
      order(rowMeans(counts(x, normalized = TRUE)), decreasing = TRUE)[1:100]
    df <- as.data.frame(colData(x)[, "condition"])
    rownames(df) <- rownames(colData(x))
    colnames(df) <- "condition"
    
    pheatmap(
      assay(rlogData)[select,],
      color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(1000),
      cellwidth = 30,
      cellheight = 10,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = FALSE,
      annotation_col = df,
      annotation_names_col = FALSE,
      filename = paste0(saveDir, sampleName, "_heatmapCounts.png"),
      width = 8,
      units = "in",
    )
  }


# differential expression analysis --------------------------------------------

colData$condition %<>% factor %>% relevel(ref = colData[1, 1])
# set the control group to the group where the first sample is located (so please put the control group first)

dds <-
  DESeqDataSetFromMatrix(round(counts),
                         colData = colData,
                         design = ~ condition)
# construct objects for differential analysis

filterCounts <- rowSums(counts) >= ncol(dds)
dds <- dds[filterCounts,]
# filter out low-expression genes

dds <- DESeq(dds)
# differential analysis

res <- results(dds)
# extract result

res <-  lfcShrink(dds, coef = 2, type = "apeglm")
# shrinking log2FoldChange to reduce the log2FoldChange value of the low-expression gene

resIhwP <-
  res %>% as.data.frame %>% mutate(ihwP = ihw(pvalue ~ baseMean, data = ., alpha = 0.1) %>% adj_pvalues())
# adjusting the P value using the IHW package, which is more stringent than the BH method, amplifies the P value of the gene with low expression levels

# save files ------------------------------------------------------------------

counts[, levels(colData$condition)[1]] <-
  rowMeans(counts[, row.names(colData[which(colData$condition == levels(colData$condition)[1]), , drop = FALSE])])
counts[, levels(colData$condition)[2]] <-
  rowMeans(counts[, row.names(colData[which(colData$condition == levels(colData$condition)[2]), , drop = FALSE])])
# the average expression amount was calculated by group

# · Total table ---------------------------------------------------------------

resAll <- merge(
  resIhwP,
  countsAnnotation,
  by.x = "row.names",
  by.y = "row.names",
  all.x = TRUE,
  sort = FALSE
)  %>% merge(
  .,
  ensToGene[c(1, 6:10)],
  by.x = "Row.names",
  by.y = "Gene stable ID",
  all.x = TRUE,
  sort = FALSE
) %>%
  merge(
    .,
    counts,
    by.x = "Row.names",
    by.y = "row.names",
    all.x = TRUE,
    sort = FALSE
  ) %>% rename(
    "Gene ID" = "Row.names",
    "Gene name" = "gene_name",
    "Chromosome/scaffold name" = "Chr",
    "Gene type" = "gene_biotype"
  )

write.csv(resAll,
          file = paste0(savePath, sample, "_DEG_All.csv"),
          row.names = FALSE)

# · differential genes table --------------------------------------------------

resIhwP %<>% merge(
  .,
  ensToGene,
  by.x = "row.names",
  by.y = "Gene stable ID",
  all.x = TRUE,
  sort = FALSE
) %>% merge(
  .,
  counts,
  by.x = "Row.names",
  by.y = "row.names",
  all.x = TRUE,
  sort = FALSE
) %>% rename("Gene ID" = "Row.names")

diffIhwP <- subset(resIhwP, ihwP < 0.05 &
                     abs(log2FoldChange) > 1) %$% .[order(log2FoldChange, decreasing = TRUE),-c(2, 4)]
# differential genes were screened according to ihwP < 0.05 and log2FoldChange > ±1

write.csv(
  diffIhwP,
  file = paste0(savePath, sample, "_DEG_logFC1.0_ihwP.csv"),
  row.names = FALSE
)

diffIhwP_filtered <-
  subset(diffIhwP, !(get(levels(colData$condition)[1]) < 50 &
                       get(levels(colData$condition)[2]) < 50))
# differential genes were screened according to ihwP < 0.05 and log2FoldChange > ±1, and genes with a mean read of less than 50 in each group were filtered out

write.csv(
  diffIhwP_filtered,
  file = paste0(savePath, sample, "_DEG_logFC1.0_ihwP_filtered.csv"),
  row.names = FALSE
)

row.names(resIhwP) <- resIhwP[, 1]
resIhwP <- resIhwP[,-1]

row.names(diffIhwP) <- diffIhwP[, 1]
diffIhwP <- diffIhwP[,-1]

row.names(diffIhwP_filtered) <- diffIhwP_filtered[, 1]
diffIhwP_filtered <- diffIhwP_filtered[,-1]
# convert the Gene ID in the first column to the row name

rld <- rlog(dds, blind = FALSE)
# regularized log transformation is used to draw plot

save.image(file = paste0(savePath, "DEG.RData"))
# save temp file

# draw ------------------------------------------------------------------------

if (grepl("cell", sample, ignore.case = TRUE)) {
  diffGene <- diffIhwP
} else {
  diffGene <- diffIhwP_filtered
}
# If the sample contains the "cell" affix (i.e., the samples is a cell sample), use the result file of ihwP < 0.05 & |log2FoldChange| > 1 for differential analysis, 
# and use the file that filters out the lowest group average expression of 50 if it is not cell (i.e., the sample is a tissue sample)

draw_cook(dds, saveDir = savePath)
draw_heatmapCounts(rld, dds, saveDir = savePath)
# Check for abnormal samples

draw_PCA(rld, saveDir = savePath)
draw_distance(rld, saveDir = savePath)
# PCA and sample clustering

draw_volcano(
  resIhwP,
  cutoff = 1,
  adjustP = "ihwP",
  saveDir = paste0(savePath, "LogFC1_ihwP/")
)
draw_heatmap(rld, diffGene, saveDir = paste0(savePath, "LogFC1_ihwP/"))
# draw volcano plot and heatmap by differential gene

lapply(names(listHeatmapGene_GO), function(x) {
  i <-
    diffGene[which(diffGene$Symbol %in% listHeatmapGene_GO[[x]]$V1),]
  if (nrow(i) > 1)
    draw_heatmap(rld,
                 i,
                 suffix = x,
                 saveDir = paste0(savePath, "subHeatmapGO"))
}) %>% invisible()
lapply(names(listHeatmapGene_KEGG), function(x) {
  i <-
    diffGene[which(diffGene$`Entrez Gene ID` %in% listHeatmapGene_KEGG[[x]]$V1),]
  if (nrow(i) > 1)
    draw_heatmap(rld,
                 i,
                 suffix = x,
                 saveDir = paste0(savePath, "subHeatmapKEGG"))
}) %>% invisible()
# draw heatmaps, use the the pheatmap package, and enter genes are taken from the intersection of differential genes and different function tables

lapply(names(listHeatmapGene_GO), function(x) {
  i <-
    diffGene[which(diffGene$Symbol %in% listHeatmapGene_GO[[x]]$V1),]
  if (nrow(i) > 1)
    draw_heatmap2(rld,
                  i,
                  suffix = x,
                  saveDir = paste0(savePath, "/subHeatmapGO2"))
}) %>% invisible()
lapply(names(listHeatmapGene_KEGG), function(x) {
  i <-
    diffGene[which(diffGene$`Entrez Gene ID` %in% listHeatmapGene_KEGG[[x]]$V1),]
  if (nrow(i) > 1)
    draw_heatmap2(rld,
                  i,
                  suffix = x,
                  saveDir = paste0(savePath, "/subHeatmapKEGG2"))
}) %>% invisible()
# draw heatmaps, use the the ComplexHeatmap package, and enter genes are taken from the intersection of differential genes and different function tables

lapply(names(listHeatmapGene_GO), function(x) {
  i <-
    resIhwP[which(resIhwP$Symbol %in% listHeatmapGene_GO[[x]]$V1),]
  if (nrow(i) > 1)
    draw_heatmap2(rld,
                  i,
                  suffix = x,
                  saveDir = paste0(savePath, "HeatmapGO"))
}) %>% invisible()
lapply(names(listHeatmapGene_KEGG), function(x) {
  i <-
    resIhwP[which(resIhwP$`Entrez Gene ID` %in% listHeatmapGene_KEGG[[x]]$V1),]
  if (nrow(i) > 1)
    draw_heatmap2(rld,
                  i,
                  suffix = x,
                  saveDir = paste0(savePath, "HeatmapKEGG"))
}) %>% invisible()
# draw heatmaps, use the the ComplexHeatmap package, and enter genes are taken from the intersection of total genes and different function tables
