dataDir <- commandArgs(trailingOnly = TRUE)

sample <- basename(dataDir)
# the prefix for all output file names

savePath <- paste0("result/", dataDir, "/Differential_edgeR/")
# the output folder

if (!dir.exists(savePath)) {
  dir.create(savePath, recursive = TRUE)
}

suppressPackageStartupMessages(library(edgeR))
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

draw_volcano <-
  function(result,
           suffix = NULL,
           cutoff = 1,
           labelNum = 10,
           adjustP = "FDR",
           sampleName = sample,
           saveDir = "./") {
    # suffix:     used to distinguish between result files to avoid duplicate names
    # cutoff:     the range of difference values for the color tags
    # labelNum:   number of labels (Up/Down are counted separately, and the first n logFC labels are taken)
    # adjustP:    the column name of adjust P value
    
    if (!dir.exists(saveDir)) {
      dir.create(saveDir, recursive = TRUE)
    }
    
    if (!is.null(suffix)) {
      suffix = paste0("_", suffix)
    }
    
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("ggrepel"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("magrittr"))
    
    result$change = as.factor(ifelse(
      result[, adjustP] < 0.05 &
        abs(result$logFC) > cutoff,
      ifelse(result$logFC > cutoff , "UP", "DOWN"),
      "STABLE"
    ))
    # set markers of up and down
    
    # · no label --------------------------------------------------------------
    
    g <-
      ggplot(data = result, aes(
        x = logFC,
        y = -log10(get(adjustP)),
        color = change
      )) +
      geom_point(size = 2, stroke = 0) +
      geom_vline(
        xintercept = c(cutoff,-cutoff),
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
                    "_volcano",
                    suffix,
                    ".png"),
      path = saveDir,
      plot = g,
      width = 8,
      height = 6,
      units = "in",
      dpi = 600
    )
    
    ggsave(
      file = paste0(sampleName,
                    "_volcano",
                    suffix,
                    ".tiff"),
      path = saveDir,
      plot = g,
      device = "tiff",
      compression = "lzw",
      width = 8,
      height = 6,
      units = "in",
      dpi = 600
    )
    
    ggsave(
      file = paste0(sampleName,
                    "_volcano",
                    suffix,
                    ".svg"),
      path = saveDir,
      plot = g,
      width = 8,
      height = 6,
      units = "in"
    )
    
    ggsave(
      file = paste0(sampleName,
                    "_volcano",
                    suffix,
                    ".eps"),
      path = saveDir,
      plot = g,
      device = cairo_ps,
      fallback_resolution = 600,
      width = 8,
      height = 6,
      units = "in"
    )
    
    # · have label ------------------------------------------------------------
    
    if ("Gene name" %in% colnames(result)) {
      result %<>% mutate(Symbol = if_else(is.na(Symbol) |
                                            Symbol == "", `Gene name`, Symbol))
    }
    result %<>% mutate(Symbol = if_else(is.na(Symbol) |
                                          Symbol == "", row.names(.), Symbol))
    
    this_tile <- paste0(
      'Cutoff for logFC is ',
      cutoff,
      '\nThe number of up gene is ',
      nrow(result[result$change == 'UP',]) ,
      '\nThe number of down gene is ',
      nrow(result[result$change == 'DOWN',])
    )
    
    g2 <- g %+%
      ggtitle(this_tile) +
      geom_label_repel(
        data = {
          result[which(result$change != "STABLE"), ] %>% group_by(change) %>% slice_max(order_by = abs(logFC), n = labelNum)
        },
        aes(label = Symbol),
        force = 10,
        min.segment.length = 0,
        max.overlaps = 30,
        show.legend = FALSE
      )
    
    ggsave(
      file = paste0(sampleName,
                    "_volcano_label",
                    suffix,
                    ".png"),
      path = saveDir,
      plot = g2,
      width = 8,
      height = 6,
      units = "in",
      dpi = 600
    )
    
    ggsave(
      file = paste0(sampleName,
                    "_volcano_label",
                    suffix,
                    ".tiff"),
      path = saveDir,
      plot = g2,
      device = "tiff",
      compression = "lzw",
      width = 8,
      height = 6,
      units = "in",
      dpi = 600
    )
    
    ggsave(
      file = paste0(sampleName,
                    "_volcano_label",
                    suffix,
                    ".svg"),
      path = saveDir,
      plot = g2,
      width = 8,
      height = 6,
      units = "in"
    )
    
    ggsave(
      file = paste0(sampleName,
                    "_volcano_label",
                    suffix,
                    ".eps"),
      path = saveDir,
      plot = g2,
      device = cairo_ps,
      fallback_resolution = 600,
      width = 8,
      height = 6,
      units = "in"
    )
  }

draw_heatmap <-
  function(data,
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
    
    suppressPackageStartupMessages(library("edgeR"))
    suppressPackageStartupMessages(library("pheatmap"))
    suppressPackageStartupMessages(library("RColorBrewer"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("magrittr"))
    
    if ("Gene name" %in% colnames(diffGene)) {
      diffGene %<>% dplyr::mutate(Symbol = if_else(is.na(Symbol) |
                                                     Symbol == "", `Gene name`, Symbol))
    }
    diffGene %<>% dplyr::mutate(Symbol = if_else(is.na(Symbol) |
                                                   Symbol == "", row.names(.), Symbol))
    
    countsMatrix <-
      edgeR::cpm(data, log = TRUE) %>% .[row.names(diffGene),]
    
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
  function(data,
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
    
    suppressPackageStartupMessages(library("edgeR"))
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
      diffGene$FDR,
      corr = FALSE,
      cutpoints = c(0,  .001, .01, .05, 1),
      symbols = c("***", "**", "*", " ")
    )
    # Statistically significant extent of the P-value
    
    countsMatrix <-
      edgeR::cpm(data, log = TRUE) %>% .[row.names(diffGene),]
    
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
    
    if (sign(max(diffGene$logFC)) == sign(min(diffGene$logFC))) {
      Log2FCWidth = unit(1, "inch")
    } else{
      Log2FCWidth = unit(max(0.375 * (
        max(diffGene$logFC) - min(diffGene$logFC)
      ) /
        min(c(
          max(diffGene$logFC), -min(diffGene$logFC)
        )) + 0.125, 1), "inch")
    }
    # calculate the width of the LogFC bar
    
    ha2 <- HeatmapAnnotation(
      pvalue = anno_simple(
        -log10(diffGene$FDR),
        col = colorRamp2(
          breaks = c(0, -log10(0.05), 3),
          colors = c("#4575B4", "white", "#D73027")
        ),
        pch = diffGene$symP
      ),
      Log2FC = anno_numeric(
        diffGene$logFC,
        bg_gp = gpar(
          fill = if (sign(max(diffGene$logFC)) != sign(min(diffGene$logFC)))
            c("#4575B4", "#D73027")
          else
            ifelse(sign(max(diffGene$logFC)) == -1, "#4575B4", "#D73027")
        ),
        labels_format =  function(x)
          sprintf("%.2f", x),
        align_to = ifelse(
          sign(max(diffGene$logFC)) == sign(min(diffGene$logFC)),
          ifelse(sign(max(diffGene$logFC)) == -1, "right", "left"),
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

# differential expression analysis --------------------------------------------

colData$condition %<>% factor %>% relevel(ref = colData[1, 1])
# set the control group to the group where the first sample is located (so please put the control group first)

group <- colData[, 1]
design <- model.matrix(~ group)
# design of experiment

y <-
  DGEList(counts = counts,
          group = group,
          genes = countsAnnotation)
# construct objects for differential analysis

keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]
# filter out low-expression genes

y <- calcNormFactors(y)
# normalize the library sizes

y2 <- y
# copy object to another test

# · glm framework -------------------------------------------------------------

y <- estimateDisp(y, design, robust = TRUE)

fit <- glmQLFit(y, design, robust = TRUE)

qlf <- glmQLFTest(fit, coef = 2)
# differential analysis

res <- topTags(qlf,
               n = Inf,
               sort.by = "PValue",
               p.value = 1)
# extract result

summary(decideTestsDGE(qlf, lfc = 1))
# show the gene's number of up and down

# · classic pipeline ----------------------------------------------------------

y2 <- estimateDisp(y2, robust = TRUE)

et <- exactTest(y2)
# differential analysis

res_et <- topTags(et,
                  n = Inf,
                  sort.by = "PValue",
                  p.value = 1)
# extract result

summary(decideTestsDGE(et, lfc = 1))
# show the gene's number of up and down

# save files ------------------------------------------------------------------

counts[, levels(colData$condition)[1]] <-
  rowMeans(counts[, row.names(colData[which(colData$condition == levels(colData$condition)[1]), , drop = FALSE])])
counts[, levels(colData$condition)[2]] <-
  rowMeans(counts[, row.names(colData[which(colData$condition == levels(colData$condition)[2]), , drop = FALSE])])
# the average expression amount was calculated by group

# · glm -----------------------------------------------------------------------

resAll <- res[, c(8:12, 1:7)] %>%
  merge(
    .,
    ensToGene[c(1, 6:10)],
    by.x = "row.names",
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
# total table

write.csv(resAll,
          file = paste0(savePath, sample, "_DEG_All.csv"),
          row.names = FALSE)

diffGene <-
  subset(resAll, FDR < 0.05 &
           abs(logFC) > 1) %$% .[order(logFC, decreasing = TRUE),-c(3, 4)]
# differential genes were screened according to FDR < 0.05 and logFC > ±1

write.csv(
  diffGene,
  file = paste0(savePath, sample, "_DEG_logFC1.csv"),
  row.names = FALSE
)

diffGene_filtered <-
  subset(diffGene,!(get(levels(colData$condition)[1]) < 50 &
                      get(levels(colData$condition)[2]) < 50))
# differential genes were screened according to FDR < 0.05 and logFC > ±1, and genes with a mean read of less than 50 in each group were filtered out

write.csv(
  diffGene_filtered,
  file = paste0(savePath, sample, "_DEG_logFC1_filtered.csv"),
  row.names = FALSE
)

row.names(diffGene) <- diffGene[, 1]
diffGene <- diffGene[, -1]

row.names(diffGene_filtered) <- diffGene_filtered[, 1]
diffGene_filtered <- diffGene_filtered[, -1]

row.names(resAll) <- resAll[, 1]
resAll <- resAll[, -1]
# convert the Gene ID in the first column to the row name

# · classic -------------------------------------------------------------------

resAll_et <- res_et[, c(8:11, 1:7)] %>%
  merge(
    .,
    ensToGene[c(1, 6:10)],
    by.x = "row.names",
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
# total table

write.csv(
  resAll_et,
  file = paste0(savePath, sample, "_DEG_All_et.csv"),
  row.names = FALSE
)

diffGene_et <-
  subset(resAll_et, FDR < 0.05 &
           abs(logFC) > 1) %$% .[order(logFC, decreasing = TRUE),-3]
# differential genes were screened according to FDR < 0.05 and logFC > ±1

write.csv(
  diffGene_et,
  file = paste0(savePath, sample, "_DEG_logFC1_et.csv"),
  row.names = FALSE
)

diffGene_et_filtered <-
  subset(diffGene_et,!(get(levels(colData$condition)[1]) < 50 &
                         get(levels(colData$condition)[2]) < 50))
# differential genes were screened according to FDR < 0.05 and logFC > ±1, and genes with a mean read of less than 50 in each group were filtered out

write.csv(
  diffGene_et_filtered,
  file = paste0(savePath, sample, "_DEG_logFC1_et_filtered.csv"),
  row.names = FALSE
)

row.names(diffGene_et) <- diffGene_et[, 1]
diffGene_et <- diffGene_et[, -1]

row.names(diffGene_et_filtered) <- diffGene_et_filtered[, 1]
diffGene_et_filtered <- diffGene_et_filtered[, -1]

row.names(resAll_et) <- resAll_et[, 1]
resAll_et <- resAll_et[, -1]
# convert the Gene ID in the first column to the row name


save.image(file = paste0(savePath, "DEG.RData"))
# save temp file

# draw ------------------------------------------------------------------------

if (!grepl("cell", sample, ignore.case = TRUE)) {
  diffGene <- diffGene_filtered
  diffGene_et <- diffGene_et_filtered
}
# If the sample contains the "cell" affix (i.e., the samples is a cell sample), use the result file of FDR < 0.05 & |logFC| > 1 for differential analysis, 
# and use the file that filters out the lowest group average expression of 50 if it is not cell (i.e., the sample is a tissue sample)

png(
  filename = paste0(savePath, sample, "_Smear.png"),
  width = 8,
  height = 6,
  units = "in",
  bg = "white",
  res = 300
)
plotSmear(qlf)
dev.off()

png(
  filename = paste0(savePath, sample, "_QLDisp.png"),
  width = 8,
  height = 6,
  units = "in",
  bg = "white",
  res = 300
)
plotQLDisp(fit)
dev.off()

png(
  filename = paste0(savePath, sample, "_BCV.png"),
  width = 8,
  height = 6,
  units = "in",
  bg = "white",
  res = 300
)
plotBCV(y)
dev.off()

png(
  filename = paste0(savePath, sample, "_MDS.png"),
  width = 8,
  height = 6,
  units = "in",
  bg = "white",
  res = 300
)
plotMDS(cpm(y, log = TRUE), col = as.numeric(y$samples$group))
dev.off()

tiff(
  filename = paste0(savePath, sample, "_MDS.tiff"),
  width = 8,
  height = 6,
  units = "in",
  compression = "lzw",
  bg = "white",
  res = 300
)
plotMDS(cpm(y, log = TRUE), col = as.numeric(y$samples$group))
dev.off()

svg(
  filename = paste0(savePath, sample, "_MDS.svg"),
  width = 8,
  height = 6,
  bg = "white"
)
plotMDS(cpm(y, log = TRUE), col = as.numeric(y$samples$group))
dev.off()

setEPS()
postscript(
  file = paste0(savePath, sample, "_MDS.eps"),
  bg = "white",
  width = 8,
  height = 6
)
plotMDS(cpm(y, log = TRUE), col = as.numeric(y$samples$group))
dev.off()

draw_volcano(
  resAll,
  suffix = "glm",
  labelNum = 10,
  saveDir = paste0(savePath, "LogFC1/")
)
draw_volcano(
  resAll_et,
  suffix = "et",
  labelNum = 10,
  saveDir = paste0(savePath, "LogFC1/")
)
# draw volcano plot by differential gene

draw_heatmap(y,
             diffGene,
             suffix = "glm",
             saveDir = paste0(savePath, "LogFC1/"))
draw_heatmap(y,
             diffGene_et,
             suffix = "et",
             saveDir = paste0(savePath, "LogFC1/"))
# draw heatmap by differential gene

lapply(names(listHeatmapGene_GO), function(x) {
  i <-
    diffGene[which(diffGene$`Gene name` %in% listHeatmapGene_GO[[x]]$V1),]
  if (nrow(i) > 1)
    draw_heatmap(y,
                 i,
                 suffix = paste0("glm_", x),
                 saveDir = paste0(savePath, "subHeatmapGO/"))
}) %>% invisible()
lapply(names(listHeatmapGene_KEGG), function(x) {
  i <-
    diffGene[which(diffGene$`Entrez Gene ID` %in% listHeatmapGene_KEGG[[x]]$V1),]
  if (nrow(i) > 1)
    draw_heatmap(
      y,
      i,
      suffix = paste0("glm_", x),
      saveDir = paste0(savePath, "subHeatmapKEGG/")
    )
}) %>% invisible()
lapply(names(listHeatmapGene_GO), function(x) {
  i <-
    diffGene_et[which(diffGene_et$`Gene name` %in% listHeatmapGene_GO[[x]]$V1),]
  if (nrow(i) > 1)
    draw_heatmap(
      y,
      i,
      suffix = paste0("et_", x),
      saveDir = paste0(savePath, "subHeatmapGO_et/")
    )
}) %>% invisible()
lapply(names(listHeatmapGene_KEGG), function(x) {
  i <-
    diffGene_et[which(diffGene_et$`Entrez Gene ID` %in% listHeatmapGene_KEGG[[x]]$V1),]
  if (nrow(i) > 1)
    draw_heatmap(
      y,
      i,
      suffix = paste0("et_", x),
      saveDir = paste0(savePath, "subHeatmapKEGG_et/")
    )
}) %>% invisible()
# draw heatmaps, use the the pheatmap package, and enter genes are taken from the intersection of differential genes and different function tables

lapply(names(listHeatmapGene_GO), function(x) {
  i <-
    diffGene[which(diffGene$`Gene name` %in% listHeatmapGene_GO[[x]]$V1),]
  if (nrow(i) > 1)
    draw_heatmap2(y,
                  i,
                  suffix = paste0("glm_", x),
                  saveDir = paste0(savePath, "subHeatmapGO2/"))
}) %>% invisible()
lapply(names(listHeatmapGene_KEGG), function(x) {
  i <-
    diffGene[which(diffGene$`Entrez Gene ID` %in% listHeatmapGene_KEGG[[x]]$V1),]
  if (nrow(i) > 1)
    draw_heatmap2(
      y,
      i,
      suffix = paste0("glm_", x),
      saveDir = paste0(savePath, "subHeatmapKEGG2/")
    )
}) %>% invisible()
lapply(names(listHeatmapGene_GO), function(x) {
  i <-
    diffGene_et[which(diffGene_et$`Gene name` %in% listHeatmapGene_GO[[x]]$V1),]
  if (nrow(i) > 1)
    draw_heatmap2(
      y,
      i,
      suffix = paste0("et_", x),
      saveDir = paste0(savePath, "subHeatmapGO_et2/")
    )
}) %>% invisible()
lapply(names(listHeatmapGene_KEGG), function(x) {
  i <-
    diffGene_et[which(diffGene_et$`Entrez Gene ID` %in% listHeatmapGene_KEGG[[x]]$V1),]
  if (nrow(i) > 1)
    draw_heatmap2(
      y,
      i,
      suffix = paste0("et_", x),
      saveDir = paste0(savePath, "subHeatmapKEGG_et2/")
    )
}) %>% invisible()
# draw heatmaps, use the the ComplexHeatmap package, and enter genes are taken from the intersection of differential genes and different function tables

lapply(names(listHeatmapGene_GO), function(x) {
  i <-
    resAll[which(resAll$`Gene name` %in% listHeatmapGene_GO[[x]]$V1),]
  if (nrow(i) > 1)
    draw_heatmap2(y,
                  i,
                  suffix = paste0("glm_", x),
                  saveDir = paste0(savePath, "HeatmapGO/"))
}) %>% invisible()
lapply(names(listHeatmapGene_KEGG), function(x) {
  i <-
    resAll[which(resAll$`Entrez Gene ID` %in% listHeatmapGene_KEGG[[x]]$V1),]
  if (nrow(i) > 1)
    draw_heatmap2(y,
                  i,
                  suffix = paste0("glm_", x),
                  saveDir = paste0(savePath, "HeatmapKEGG/"))
}) %>% invisible()
lapply(names(listHeatmapGene_GO), function(x) {
  i <-
    resAll_et[which(resAll_et$`Gene name` %in% listHeatmapGene_GO[[x]]$V1),]
  if (nrow(i) > 1)
    draw_heatmap2(y,
                  i,
                  suffix = paste0("et_", x),
                  saveDir = paste0(savePath, "HeatmapGO_et/"))
}) %>% invisible()
lapply(names(listHeatmapGene_KEGG), function(x) {
  i <-
    resAll_et[which(resAll_et$`Entrez Gene ID` %in% listHeatmapGene_KEGG[[x]]$V1),]
  if (nrow(i) > 1)
    draw_heatmap2(
      y,
      i,
      suffix = paste0("et_", x),
      saveDir = paste0(savePath, "HeatmapKEGG_et/")
    )
}) %>% invisible()
# draw heatmaps, use the the ComplexHeatmap package, and enter genes are taken from the intersection of total genes and different function tables
