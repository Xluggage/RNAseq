dataDir <- commandArgs(trailingOnly = TRUE)

sample <- basename(dataDir)
# the prefix for all output file names

savePath <- paste0("result/", dataDir, "/Enrichment/")
# the output folder

if (!dir.exists(savePath)) {
  dir.create(savePath, recursive = TRUE)
}

suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(ReactomePA))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))

# read files ------------------------------------------------------------------

listPathway_GO_gene <- list.files(
  path = paste0("./genome/Pathway/GO/"),
  pattern = "gene.tsv$",
  full.names = TRUE
) %>% magrittr::set_names(., stringr::str_replace(basename(.), "\\.tsv$", "")) %>% lapply(., read.delim, head = FALSE)

listPathway_GO_term <- list.files(
  path = paste0("./genome/Pathway/GO/"),
  pattern = "term.tsv$",
  full.names = TRUE
) %>% magrittr::set_names(., stringr::str_replace(basename(.), "\\.tsv$", "")) %>% lapply(., read.delim, head = FALSE)

listPathway_KEGG_gene <- list.files(
  path = paste0("./genome/Pathway/KEGG/"),
  pattern = "gene.tsv$",
  full.names = TRUE
) %>% magrittr::set_names(., stringr::str_replace(basename(.), "\\.tsv$", "")) %>% lapply(., read.delim, head = FALSE)

listPathway_KEGG_name <- list.files(
  path = paste0("./genome/Pathway/KEGG/"),
  pattern = "name.tsv$",
  full.names = TRUE
) %>% magrittr::set_names(., stringr::str_replace(basename(.), "\\.tsv$", "")) %>% lapply(., read.delim, head = FALSE)
# read the gene set by different function

# · used for ORA --------------------------------------------------------------

if (grepl("cell", sample,ignore.case = TRUE)) {
  diffGene <- read.csv(
    paste0(
      "result/",
      dataDir,
      "/Differential_edgeR/",
      sample,
      "_DEG_logFC1.csv"
    ),
    header = TRUE,
    check.names = FALSE
  )
} else {
  diffGene <- read.csv(
    paste0(
      "result/",
      dataDir,
      "/Differential_edgeR/",
      sample,
      "_DEG_logFC1_filtered.csv"
    ),
    header = TRUE,
    check.names = FALSE
  )
}
# read differential genes file
# If the sample contains the "cell" affix (i.e., the samples is a cell sample), use the result file of FDR < 0.05 & |logFC| > 1 for differential analysis, 
# and use the file that filters out the lowest group average expression of 50 if it is not cell (i.e., the sample is a tissue sample)

diffGene %<>% .[which(diffGene$`Gene ID` != "ENSMUSG00000114515"), ]
diffGene %<>% .[which(diffGene$`Gene ID` != "ENSMUSG00000116275"), ]
# The Symbol of ENSMUSG00000114515 and ENSMUSG00000030695 is duplicated, both are Aldoa
# The Symbol of ENSMUSG00000116275 deleted by version 109 are duplicated with ENSMUSG00000102976, both are Zc3h11a  

symbol_all <- diffGene[, "Symbol"] %>% na.omit()
symbol_up <-
  diffGene[diffGene$logFC > 0, "Symbol"] %>% na.omit()
symbol_down <-
  diffGene[diffGene$logFC < 0, "Symbol"] %>% na.omit()
# extract differential gene Symbol

ensembl_all <- diffGene[, 1]
ensembl_up <-
  diffGene[diffGene$logFC > 0, 1]
ensembl_down <-
  diffGene[diffGene$logFC < 0, 1]
# extract differential gene ensembl ID

entrez_all <- diffGene[, "Entrez Gene ID"] %>% na.omit()
entrez_up <-
  diffGene[diffGene$logFC > 0, "Entrez Gene ID"] %>% na.omit()
entrez_down <-
  diffGene[diffGene$logFC < 0, "Entrez Gene ID"] %>% na.omit()
# extract differential gene Entrez ID


# · used for GSEA -------------------------------------------------------------

allGene <- read.csv(
  paste0("result/",
         dataDir,
         "/Differential_edgeR/",
         sample,
         "_DEG_All.csv"),
  header = TRUE,
  check.names = FALSE
)
# read differential analysis total genes file

allGene %<>% .[.$FDR < 0.05, ]
# extract genes by FDR < 0.05

allGene %<>% .[which(allGene$`Gene ID` != "ENSMUSG00000114515"), ]
allGene %<>% .[which(allGene$`Gene ID` != "ENSMUSG00000116275"), ]
# The Symbol of ENSMUSG00000114515 and ENSMUSG00000030695 is duplicated, both are Aldoa
# The Symbol of ENSMUSG00000116275 deleted by version 109 are duplicated with ENSMUSG00000102976, both are Zc3h11a  

symbolListTmp <-
  allGene[, c("Symbol", "logFC")] %>% .[order(.$logFC, decreasing = TRUE), ] %>% na.omit()
symbolList_up <-
  symbolListTmp[symbolListTmp$logFC > 0, "Symbol"]
symbolList_down <-
  symbolListTmp[symbolListTmp$logFC < 0, "Symbol"]
# extract gene Symbol which expression up or down

symbolList <- symbolListTmp[, 2]
names(symbolList) <- symbolListTmp[, 1]
rm(symbolListTmp)
# logFC value list, rownames is Symbol

ensemblListTmp <-
  allGene[, c("Gene ID", "logFC")] %>% .[order(.$logFC, decreasing = TRUE), ] %>% na.omit()
ensemblList <- ensemblListTmp[, 2]
names(ensemblList) <- ensemblListTmp[, 1]
rm(ensemblListTmp)
# logFC value list, rownames is Gene ID

entrezListTmp <-
  allGene[, c("Entrez Gene ID", "logFC")] %>% .[order(.$logFC, decreasing = TRUE), ] %>% na.omit()
entrezList <-  entrezListTmp[, 2]
names(entrezList) <- entrezListTmp[, 1]
rm(entrezListTmp)
# logFC value list, rownames is Entrez ID


# analysis --------------------------------------------------------------------

# · Function definition -------------------------------------------------------

# ·· analysis function --------------------------------------------------------

simply <- function(x,
                   saveDir = ".",
                   objectName = deparse1(substitute(x)),
                   sampleName = sample) {
  # Number of results for GO enrichment analysis to be refined, based on semantic similarity of different terms (parameter cutoff).https://yulab-smu.top/biomedical-knowledge-mining-book/semantic-similarity-overview.html
  # x:              enrichGO object or gseGO object
  # objectName:     used to distinguish between result files to avoid duplicate names
  # cutoff:         semantic similarity, default value of 0.7, i.e., merge gene sets with similarity greater than 0.7 in the enrichment result (use measure = "Wang")
  # by/select_fun:  specify how the similar gene set is processed, with "by" as the reference parameter and "select_fun" as the function to be used, "by = 'p.adjust', select_fun = min" is to select the term with the lowest p.adjust value in the similar gene set
  
  suppressPackageStartupMessages(require(clusterProfiler))
  
  if (is.null(x) || nrow(x) < 10) {
    return(x)
  } else{
    simply <- clusterProfiler::simplify(
      x,
      cutoff = 0.7,
      by = "p.adjust",
      select_fun = min,
      measure = "Wang"
    )
    
    saveDir <- paste0(saveDir, "/")
    if (!dir.exists(saveDir)) {
      dir.create(saveDir, recursive = TRUE)
    }
    write.csv(simply,
              file = paste0(saveDir,
                            sampleName,
                            "_",
                            objectName,
                            "_simply",
                            ".csv"))
    return(simply)
  }
}

enrichORA <-
  function(x,
           gene_up = symbol_up,
           gene_down = symbol_down,
           term2gene,
           term2name,
           objectName = "GO",
           saveDir = "./",
           sampleName = sample,
           dbType = "GO",
           keyType = "SYMBOL",
           pvalueCutoff = 1,
           qvalueCutoff = 1,
           min = 1,
           max = 500) {
    # the generic over representation enrichment analysis
    # x:              a list of differential genes in the form of a string vector (the gene name's type should be the same as term2gene, GO uses Symbol here) without sorting
    # gene_up/down:   a up/down subset of the differential gene list (Symbol), which is used to generate the up/down gene name and number after the results
    # term2gene:      GO ID (or pathway ID) and gene name mapping table
    # term2name:      GO ID (or pathway ID) and GO term name (or pathway name) mapping table
    # objectName:     used to distinguish between result files to avoid duplicate names
    # dbType:         the database type entered, if set to "GO", will add ONTOLOGY grouping information to the results
    # keyType:        the type of gene name entered, if not "SYMBOL", the enrichment result is converted to Symbol (KEGG type "ENTREZID", parameters: columns(org.Mm.eg.db))
    # pvalueCutoff:   the p-value cut-off value, default 1 (this is because the custom library is too small and the p-value is large, and it can be changed as appropriate)
    # qvalueCutoff:   the q-value cut-off value, default 1 (this is because the custom library is too small and the q-value is large, and it can be changed as appropriate)
    # min/max:        gene set size limits to screen for gene sets to participate in the analysis (screening out gene sets that are too small or too general)
    
    suppressPackageStartupMessages(require(clusterProfiler))
    suppressPackageStartupMessages(require(org.Mm.eg.db))
    suppressPackageStartupMessages(require(dplyr))
    suppressPackageStartupMessages(require(purrr))
    
    ORA <- try(enricher(
      gene = x,
      TERM2GENE = term2gene,
      TERM2NAME = term2name,
      universe = NULL,
      pAdjustMethod = "BH",
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      minGSSize = min,
      maxGSSize = max
    ))
    
    if ("try-error" %in% class(ORA)) {
      cat(
        "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
        paste("Sample:", sampleName),
        paste("Database: ", objectName, "diy"),
        ORA,
        "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
        sep = "\n",
        file = "error_enrichment.log",
        append = TRUE
      )
      return(NULL)
    } else if (is.null(ORA) || nrow(ORA) == 0) {
      cat(
        "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
        paste("Sample:", sampleName),
        paste("Database: ", objectName, "diy"),
        "The enrichment result is empty, nothing is enriched",
        "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
        sep = "\n",
        file = "error_enrichment.log",
        append = TRUE
      )
      return(NULL)
    } else{
      keyType <- toupper(keyType)
      if (keyType != "SYMBOL") {
        ORA <- setReadable(ORA, OrgDb = org.Mm.eg.db, keyType = keyType)
        # converts the GeneID in the enrichment result from the input gene name type (e.g., "ENTREZID") to Symbol
      }
      
      dbType <- toupper(dbType)
      if (dbType == "GO") {
        ORA@result <-
          go2ont(ORA@result$ID) %>% merge(
            ORA@result,
            .,
            by.x = "ID",
            by.y = "go_id",
            all.x = TRUE,
            sort = FALSE
          ) %>% relocate("ONTOLOGY" = "Ontology", .after = "ID") %>% {
            rownames(.) <- .$ID
            .
          }
        # the last "." is very important, don't delete
        # use the function "go2ont" to get ONTOLOGY, but it's possible to get NA (the NA means the term is not in org.Mm.eg.db)
      }
      
      if (!is.null(gene_up) & !is.null(gene_down)) {
        list <- strsplit(ORA@result$geneID, "/")
        ORA@result$geneUp <-
          purrr::map_chr(list,  ~ paste(intersect(.x, gene_up), collapse = "/"))
        ORA@result$geneDown <-
          purrr::map_chr(list,  ~ paste(intersect(.x, gene_down), collapse = "/"))
        ORA@result <-
          ORA@result %>% dplyr::relocate(., Count, .after = last_col())
        ORA@result$up <-
          unlist(purrr::map(list,  ~ length(intersect(.x, gene_up))))
        ORA@result$down <-
          unlist(purrr::map(list,  ~ length(intersect(.x, gene_down))))
        # add the Symbol and number of genes in UP and DOWN at the end of the enrichment result (the total Count is already there, so the order is shifted here)
      }
      
      saveDir <- paste0(saveDir, "/")
      if (!dir.exists(saveDir)) {
        dir.create(saveDir, recursive = TRUE)
      }
      write.csv(
        ORA,
        file = paste0(saveDir,
                      sampleName,
                      "_",
                      objectName,
                      ".csv"),
        row.names = FALSE
      )
      return(ORA)
    }
  }

enrichGSEA <-
  function(x,
           gene_up = symbolList_up,
           gene_down = symbolList_down,
           term2gene,
           term2name,
           objectName = "gseaGO",
           saveDir = ".",
           sampleName = sample,
           dbType = "GO",
           keyType = "SYMBOL",
           pvalueCutoff = 1,
           min = 1,
           max = 500) {
    # the generic gene set enrichment analysis
    # x:              a list of differential genes in the form of a string vector (the gene name's type should be the same as term2gene, GO uses Symbol here) without sorting
    # gene_up/down:   a up/down subset of the differential gene list (Symbol), which is used to generate the up/down gene name and number after the results
    # term2gene:      GO ID (or pathway ID) and gene name mapping table
    # term2name:      GO ID (or pathway ID) and GO term name (or pathway name) mapping table
    # objectName:     used to distinguish between result files to avoid duplicate names
    # dbType:         the database type entered, if set to "GO", will add ONTOLOGY grouping information to the results
    # keyType:        the type of gene name entered, if not "SYMBOL", the enrichment result is converted to Symbol (KEGG type "ENTREZID", parameters: columns(org.Mm.eg.db))
    # pvalueCutoff:   the p-value cut-off value, default 1 (this is because the custom library is too small and the p-value is large, and it can be changed as appropriate)
    # min/max:        gene set size limits to screen for gene sets to participate in the analysis (screening out gene sets that are too small or too general)
    
    suppressPackageStartupMessages(require(clusterProfiler))
    suppressPackageStartupMessages(require(org.Mm.eg.db))
    suppressPackageStartupMessages(require(dplyr))
    suppressPackageStartupMessages(require(purrr))
    
    set.seed(123)
    # set a random seed to make it easy to reproduce the results
    
    GSEA <- try(GSEA(
      geneList = x,
      TERM2GENE = term2gene,
      TERM2NAME = term2name,
      pvalueCutoff = pvalueCutoff,
      pAdjustMethod = "BH",
      exponent = 1,
      eps = 0,
      minGSSize = min,
      maxGSSize = max,
      verbose = FALSE,
      seed = TRUE,
      by = "fgsea"
    ))
    
    if ("try-error" %in% class(GSEA)) {
      cat(
        "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
        paste("Sample:", sampleName),
        paste("Database: ", objectName, "diy"),
        GSEA,
        "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
        sep = "\n",
        file = "error_enrichment.log",
        append = TRUE
      )
      return(NULL)
    } else if (is.null(GSEA) || nrow(GSEA) == 0) {
      cat(
        "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
        paste("Sample:", sampleName),
        paste("Database: ", objectName, "diy"),
        "The enrichment result is empty, nothing is enriched",
        "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
        sep = "\n",
        file = "error_enrichment.log",
        append = TRUE
      )
      return(NULL)
    } else{
      keyType <- toupper(keyType)
      if (keyType != "SYMBOL") {
        GSEA <-
          setReadable(GSEA, OrgDb = org.Mm.eg.db, keyType = keyType)
        # converts the GeneID in the enrichment result from the input gene name type (e.g., "ENTREZID") to Symbol
      }
      
      dbType <- toupper(dbType)
      if (dbType == "GO") {
        GSEA@result <-
          go2ont(GSEA@result$ID) %>% merge(
            GSEA@result,
            .,
            by.x = "ID",
            by.y = "go_id",
            all.x = TRUE,
            sort = FALSE
          ) %>% relocate("ONTOLOGY" = "Ontology", .after = "ID") %>% {
            rownames(.) <- .$ID
            .
          }
        # the last "." is very important, don't delete
        # use the function "go2ont" to get ONTOLOGY, but it's possible to get NA (the NA means the term is not in org.Mm.eg.db)
        }
      
      if (!is.null(gene_up) & !is.null(gene_down)) {
        list <- strsplit(GSEA@result$core_enrichment, "/")
        GSEA@result$geneUp <-
          purrr::map_chr(list,  ~ paste(intersect(.x, gene_up), collapse = "/"))
        GSEA@result$geneDown <-
          purrr::map_chr(list,  ~ paste(intersect(.x, gene_down), collapse = "/"))
        GSEA@result$Count <- unlist(purrr::map(list,  ~ length(.x)))
        GSEA@result$up <-
          unlist(purrr::map(list,  ~ length(intersect(.x, gene_up))))
        GSEA@result$down <-
          unlist(purrr::map(list,  ~ length(intersect(.x, gene_down))))
        # add the Symbol and number of core genes in UP and DOWN at the end of the enrichment result
      }
      
      saveDir <- paste0(saveDir, "/")
      if (!dir.exists(saveDir)) {
        dir.create(saveDir, recursive = TRUE)
      }
      write.csv(
        GSEA,
        file = paste0(saveDir,
                      sampleName,
                      "_",
                      objectName,
                      ".csv"),
        row.names = FALSE
      )
      return(GSEA)
    }
  }

oraGO <- function(x,
                  gene_up = symbol_up,
                  gene_down = symbol_down,
                  objectName = "GO",
                  saveDir = ".",
                  sampleName = sample,
                  ont = "ALL",
                  keyType = "ENSEMBL",
                  pvalue = 0.05,
                  qvalue = 0.2,
                  min = 10,
                  max = 500) {
  # over representation enrichment analysis using GO database
  # x:              the list of differential genes is in the form of a string of ensembl IDs string vectors, no need to sort, and for other gene IDs, see the parameter keyType
  # gene_up/down:   a up/down subset of the differential gene list (Symbol), which is used to generate the up/down gene name and number after the results
  # objectName:     used to distinguish between result files to avoid duplicate names
  # ont:            a subset of the GO database, divided into three types: "BP", "CC", "MF", or analyzed together using "ALL"
  # keyType:        the type of gene IDs entered, "ENSEMBL" is ensembl ID, "ENTREZID" is Entrez ID, parameters: columns(org.Mm.eg.db)
  # pvalue/qvalue:  the p-value and q-value cut-off values
  # min/max:        gene set size limit to screen for analysis (screening out overly general gene sets)
  # readable:       the results are output as gene names instead of gene IDs, and convertion is used by org.Mm.eg.db
  
  suppressPackageStartupMessages(require(clusterProfiler))
  suppressPackageStartupMessages(require(org.Mm.eg.db))
  suppressPackageStartupMessages(require(purrr))
  
  GO <-
    try(enrichGO(
      gene = x,
      OrgDb = org.Mm.eg.db,
      ont = ont,
      keyType = keyType,
      pvalueCutoff = pvalue,
      pAdjustMethod = "BH",
      qvalueCutoff = qvalue,
      minGSSize = min,
      maxGSSize = max,
      readable = TRUE
    ))
  if ("try-error" %in% class(GO)) {
    cat(
      "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
      paste("Sample:", sampleName),
      paste("Database:", objectName),
      GO,
      "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
      sep = "\n",
      file = "error_enrichment.log",
      append = TRUE
    )
    return(NULL)
  } else if (is.null(GO) || nrow(GO) == 0) {
    cat(
      "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
      paste("Sample:", sampleName),
      paste("Database:", objectName),
      "The enrichment result is empty, nothing is enriched",
      "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
      sep = "\n",
      file = "error_enrichment.log",
      append = TRUE
    )
    return(NULL)
  } else{
    if (!is.null(gene_up) & !is.null(gene_down)) {
      list <- strsplit(GO@result$geneID, "/")
      GO@result$geneUp <-
        purrr::map_chr(list,  ~ paste(intersect(.x, gene_up), collapse = "/"))
      GO@result$geneDown <-
        purrr::map_chr(list,  ~ paste(intersect(.x, gene_down), collapse = "/"))
      GO@result <-
        GO@result %>% dplyr::relocate(., Count, .after = last_col())
      GO@result$up <-
        unlist(purrr::map(list,  ~ length(intersect(.x, gene_up))))
      GO@result$down <-
        unlist(purrr::map(list,  ~ length(intersect(.x, gene_down))))
      # add the Symbol and number of genes in UP and DOWN at the end of the enrichment result (the total Count is already there, so the order is shifted here)
    }
    
    saveDir <- paste0(saveDir, "/")
    if (!dir.exists(saveDir)) {
      dir.create(saveDir, recursive = TRUE)
    }
    write.csv(GO,
              file = paste0(saveDir,
                            sampleName,
                            "_",
                            objectName,
                            ".csv"))
    return(GO)
  }
}

gseaGO <- function(x,
                   gene_up = symbolList_up,
                   gene_down = symbolList_down,
                   saveDir = ".",
                   sampleName = sample,
                   ont = "ALL",
                   keyType = "ENSEMBL",
                   pvalue = 0.05,
                   min = 10,
                   max = 500) {
  # gene set enrichment analysis using GO database
  # x:              the list of genes that is not filtered by difference value is in the form of a string of logFC value vectors, sorted in descending order, with the gene name (ensembl ID) as the name of each item, and the parameter keyType for other gene IDs
  # gene_up/down:   a up/down subset of the differential gene list (Symbol), which is used to generate the up/down gene name and number after the results
  # ont:            a subset of the GO database, divided into three types: "BP", "CC", "MF", or analyzed together using "ALL"
  # keyType:        the type of gene IDs entered, "ENSEMBL" is ensembl ID, "ENTREZID" is Entrez ID, parameters: columns(org.Mm.eg.db)
  # pvalue:         p valve cut-off value
  # min/max:        gene set size limit to screen for analysis (screening out overly general gene sets)
  
  suppressPackageStartupMessages(require(clusterProfiler))
  suppressPackageStartupMessages(require(org.Mm.eg.db))
  
  set.seed(123)
  gseaGO <- try(gseGO(
    geneList  = x,
    OrgDb  = org.Mm.eg.db,
    ont  = ont,
    keyType = keyType,
    pvalueCutoff = pvalue,
    pAdjustMethod = "BH",
    exponent = 1,
    minGSSize  = min,
    maxGSSize  = max,
    eps = 0,
    verbose = FALSE,
    seed = TRUE,
    by = "fgsea"
  ))
  if ("try-error" %in% class(gseaGO)) {
    cat(
      "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
      paste("Sample:", sampleName),
      paste("Database: GO\tGSEA\tGO-ont:",
            ont),
      gseaGO,
      "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
      sep = "\n",
      file = "error_enrichment.log",
      append = TRUE
    )
    return(NULL)
  } else if (is.null(gseaGO) || nrow(gseaGO) == 0) {
    cat(
      "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
      paste("Sample:", sampleName),
      paste("Database: GO\tGSEA\tGO-ont:",
            ont),
      "The enrichment result is empty, nothing is enriched",
      "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
      sep = "\n",
      file = "error_enrichment.log",
      append = TRUE
    )
    return(NULL)
  } else{
    if (keyType == "ENSEMBL") {
      gseaGO = setReadable(gseaGO, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL")
    } else if (keyType == "ENTREZID") {
      gseaGO = setReadable(gseaGO, OrgDb = "org.Mm.eg.db", keyType = "ENTREZID")
    }
    
    if (!is.null(gene_up) & !is.null(gene_down)) {
      list <- strsplit(gseaGO@result$core_enrichment, "/")
      gseaGO@result$geneUp <-
        purrr::map_chr(list,  ~ paste(intersect(.x, gene_up), collapse = "/"))
      gseaGO@result$geneDown <-
        purrr::map_chr(list,  ~ paste(intersect(.x, gene_down), collapse = "/"))
      gseaGO@result$Count <- unlist(purrr::map(list,  ~ length(.x)))
      gseaGO@result$up <-
        unlist(purrr::map(list,  ~ length(intersect(.x, gene_up))))
      gseaGO@result$down <-
        unlist(purrr::map(list,  ~ length(intersect(.x, gene_down))))
      # add the Symbol and number of core genes in UP and DOWN at the end of the enrichment result
    }
    
    saveDir <- paste0(saveDir, "/")
    if (!dir.exists(saveDir)) {
      dir.create(saveDir, recursive = TRUE)
    }
    write.csv(gseaGO,
              file = paste0(saveDir,
                            sampleName,
                            "_gseaGO.csv"))
    return(gseaGO)
  }
}

oraKEGG <- function(x,
                    gene_up = symbol_up,
                    gene_down = symbol_down,
                    objectName = "KEGG",
                    saveDir = ".",
                    sampleName = sample,
                    keyType = "kegg",
                    pvalue = 0.05,
                    qvalue = 0.2,
                    min = 10,
                    max = 500,
                    offline = TRUE) {
  # over representation enrichment analysis using KEGG pathway database
  # x:              the list of differential genes is in the form of a string of Entrez IDs string vectors, no need to sort, and for other gene IDs, see the parameter keyType
  # gene_up/down:   a up/down subset of the differential gene list (Symbol), which is used to generate the up/down gene name and number after the results
  # objectName:     used to distinguish between result files to avoid duplicate names
  # keyType:        the type of gene IDs, includes "kegg", "ncbi-geneid", "ncbi-proteinid" and "uniprot". the "kegg" is Entrez ID (eukaryotes) or Locus ID (prokaryotes)
  # pvalue/qvalue:  the p-value and q-value cut-off values
  # min/max:        gene set size limit to screen for analysis (screening out overly general gene sets)
  # organism:       organism species information, "ko" is KEGG Orthology, mouse is "mmu", other species check search_kegg_organism() or http://www.genome.jp/kegg/catalog/org_list.html
  # offline:        whether to use the data of the installed KEGG.db (modified version), the default is TRUE, and if you want the latest online data, changed to FALSE
  
  suppressPackageStartupMessages(require(clusterProfiler))
  suppressPackageStartupMessages(require(org.Mm.eg.db))
  
  KEGG <- try(enrichKEGG(
    gene = x,
    organism = "mmu",
    keyType = keyType,
    pvalueCutoff = pvalue,
    pAdjustMethod = "BH",
    qvalueCutoff = qvalue,
    minGSSize = min,
    maxGSSize = max,
    use_internal_data = offline
  ))
  if ("try-error" %in% class(KEGG)) {
    cat(
      "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
      paste("Sample:", sampleName),
      paste("Database", objectName),
      KEGG,
      "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
      sep = "\n",
      file = "error_enrichment.log",
      append = TRUE
    )
    return(NULL)
  } else if (is.null(KEGG) || nrow(KEGG) == 0) {
    cat(
      "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
      paste("Sample:", sampleName),
      paste("Database:", objectName),
      "The enrichment result is empty, nothing is enriched",
      "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
      sep = "\n",
      file = "error_enrichment.log",
      append = TRUE
    )
    return(NULL)
  } else{
    if (keyType == "kegg") {
      KEGG = setReadable(KEGG, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
      # converts the GeneID in the enrichment result from ENTREZID to Symbol
    }
    
    KEGG@result$Description <-
      sub("\\s-\\s[a-zA-Z ]+\\([a-zA-Z ]+\\)$",
          "",
          KEGG@result$Description)
    # Due to the update of the KEGG official website api, a species information will be added after the path name, which will be removed here(clusterProfiler updated 4.7.1.003 to add this function but it is wrong for mouse)
    
    if (!is.null(gene_up) & !is.null(gene_down)) {
      list <- strsplit(KEGG@result$geneID, "/")
      KEGG@result$geneUp <-
        purrr::map_chr(list,  ~ paste(intersect(.x, gene_up), collapse = "/"))
      KEGG@result$geneDown <-
        purrr::map_chr(list,  ~ paste(intersect(.x, gene_down), collapse = "/"))
      KEGG@result <-
        KEGG@result %>% dplyr::relocate(., Count, .after = last_col())
      KEGG@result$up <-
        unlist(purrr::map(list,  ~ length(intersect(.x, gene_up))))
      KEGG@result$down <-
        unlist(purrr::map(list,  ~ length(intersect(.x, gene_down))))
      # add the Symbol and number of genes in UP and DOWN at the end of the enrichment result (the total Count is already there, so the order is shifted here)
    }
    
    saveDir <- paste0(saveDir, "/")
    if (!dir.exists(saveDir)) {
      dir.create(saveDir, recursive = TRUE)
    }
    write.csv(KEGG,
              file = paste0(saveDir,
                            sampleName,
                            "_",
                            objectName,
                            ".csv"))
    return(KEGG)
  }
}

gseaKEGG <- function(x,
                     gene_up = symbolList_up,
                     gene_down = symbolList_down,
                     saveDir = ".",
                     sampleName = sample,
                     keyType = "kegg",
                     pvalue = 0.05,
                     min = 10,
                     max = 500,
                     offline = TRUE) {
  # gene set enrichment analysis using KEGG pathway database
  # x:              the list of genes that is not filtered by difference value is in the form of a string of logFC value vectors, sorted in descending order, with the gene name (Entrez ID) as the name of each item, and the parameter keyType for other gene IDs
  # gene_up/down:   a up/down subset of the differential gene list (Symbol), which is used to generate the up/down gene name and number after the results
  # keyType:        the type of gene IDs, includes "kegg", "ncbi-geneid", "ncbi-proteinid" and "uniprot". the "kegg" is Entrez ID (eukaryotes) or Locus ID (prokaryotes)
  # pvalue:         p valve cut-off value
  # min/max:        gene set size limit to screen for analysis (screening out overly general gene sets)
  # organism:       organism species information, "ko" is KEGG Orthology, mouse is "mmu", other species check search_kegg_organism() or http://www.genome.jp/kegg/catalog/org_list.html
  # offline:        whether to use the data of the installed KEGG.db (modified version), the default is TRUE, and if you want the latest online data, changed to FALSE
  
  suppressPackageStartupMessages(require(clusterProfiler))
  suppressPackageStartupMessages(require(org.Mm.eg.db))
  
  set.seed(123)
  gseaKEGG <- try(gseKEGG(
    geneList  = x,
    organism = "mmu",
    keyType = keyType,
    pvalueCutoff = pvalue,
    pAdjustMethod = "BH",
    exponent = 1,
    minGSSize  = min,
    maxGSSize  = max,
    eps = 0,
    verbose = FALSE,
    seed = TRUE,
    by = "fgsea",
    use_internal_data = offline
  ))
  if ("try-error" %in% class(gseaKEGG)) {
    cat(
      "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
      paste("Sample:", sampleName),
      "Database: KEGG\tGSEA",
      gseaKEGG,
      "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
      sep = "\n",
      file = "error_enrichment.log",
      append = TRUE
    )
    return(NULL)
  } else if (is.null(gseaKEGG) || nrow(gseaKEGG) == 0) {
    cat(
      "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
      paste("Sample:", sampleName),
      "Database: KEGG\tGSEA",
      "The enrichment result is empty, nothing is enriched",
      "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
      sep = "\n",
      file = "error_enrichment.log",
      append = TRUE
    )
    return(NULL)
  } else{
    if (keyType == "kegg") {
      gseaKEGG = setReadable(gseaKEGG, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
      # converts the GeneID in the enrichment result from ENTREZID to Symbol
    }
    
    gseaKEGG@result$Description <-
      sub("\\s-\\s[a-zA-Z ]+\\([a-zA-Z ]+\\)$",
          "",
          gseaKEGG@result$Description)
    # Due to the update of the KEGG official website api, a species information will be added after the path name, which will be removed here(clusterProfiler updated 4.7.1.003 to add this function but it is wrong for mouse)
    
    if (!is.null(gene_up) & !is.null(gene_down)) {
      list <- strsplit(gseaKEGG@result$core_enrichment, "/")
      gseaKEGG@result$geneUp <-
        purrr::map_chr(list,  ~ paste(intersect(.x, gene_up), collapse = "/"))
      gseaKEGG@result$geneDown <-
        purrr::map_chr(list,  ~ paste(intersect(.x, gene_down), collapse = "/"))
      gseaKEGG@result$Count <-
        unlist(purrr::map(list,  ~ length(.x)))
      gseaKEGG@result$up <-
        unlist(purrr::map(list,  ~ length(intersect(.x, gene_up))))
      gseaKEGG@result$down <-
        unlist(purrr::map(list,  ~ length(intersect(.x, gene_down))))
      # add the Symbol and number of core genes in UP and DOWN at the end of the enrichment result
    }
    
    saveDir <- paste0(saveDir, "/")
    if (!dir.exists(saveDir)) {
      dir.create(saveDir, recursive = TRUE)
    }
    write.csv(gseaKEGG,
              file = paste0(saveDir,
                            sampleName,
                            "_gseaKEGG.csv"))
    return(gseaKEGG)
  }
}


oraReactome <- function(x,
                        gene_up = symbol_up,
                        gene_down = symbol_down,
                        objectName = "Reactome",
                        saveDir = ".",
                        sampleName = sample,
                        pvalue = 0.05,
                        qvalue = 0.2,
                        min = 10,
                        max = 500) {
  # over representation enrichment analysis using Reactome database
  # x:              the list of differential genes is in the form of a string of Entrez IDs string vectors, no need to sort
  # gene_up/down:   a up/down subset of the differential gene list (Symbol), which is used to generate the up/down gene name and number after the results
  # objectName:     used to distinguish between result files to avoid duplicate names
  # pvalue/qvalue:  the p-value and q-value cut-off values
  # min/max:        gene set size limit to screen for analysis (screening out overly general gene sets)
  # organism:       organism species information, mouse is "mouse"
  # readable:       the results are output as gene names instead of gene IDs, and convertion is used by org.Mm.eg.db
  
  suppressPackageStartupMessages(require(clusterProfiler))
  suppressPackageStartupMessages(require(ReactomePA))
  
  Reactome <- try(enrichPathway(
    gene = x,
    organism = "mouse",
    pvalueCutoff = pvalue,
    pAdjustMethod = "BH",
    qvalueCutoff = qvalue,
    minGSSize = min,
    maxGSSize = max,
    readable = TRUE
  ))
  if ("try-error" %in% class(Reactome)) {
    cat(
      "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
      paste("Sample:", sampleName),
      paste("Database:", objectName),
      Reactome,
      "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
      sep = "\n",
      file = "error_enrichment.log",
      append = TRUE
    )
    return(NULL)
  } else if (is.null(Reactome) || nrow(Reactome) == 0) {
    cat(
      "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
      paste("Sample:", sampleName),
      paste("Database:", objectName),
      "The enrichment result is empty, nothing is enriched",
      "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
      sep = "\n",
      file = "error_enrichment.log",
      append = TRUE
    )
    return(NULL)
  } else{
    if (!is.null(gene_up) & !is.null(gene_down)) {
      list <- strsplit(Reactome@result$geneID, "/")
      Reactome@result$geneUp <-
        purrr::map_chr(list,  ~ paste(intersect(.x, gene_up), collapse = "/"))
      Reactome@result$geneDown <-
        purrr::map_chr(list,  ~ paste(intersect(.x, gene_down), collapse = "/"))
      Reactome@result <-
        Reactome@result %>% dplyr::relocate(., Count, .after = last_col())
      Reactome@result$up <-
        unlist(purrr::map(list,  ~ length(intersect(.x, gene_up))))
      Reactome@result$down <-
        unlist(purrr::map(list,  ~ length(intersect(.x, gene_down))))
      # add the Symbol and number of genes in UP and DOWN at the end of the enrichment result (the total Count is already there, so the order is shifted here)
    }
    
    saveDir <- paste0(saveDir, "/")
    if (!dir.exists(saveDir)) {
      dir.create(saveDir, recursive = TRUE)
    }
    write.csv(Reactome,
              file = paste0(saveDir,
                            sampleName,
                            "_",
                            objectName,
                            ".csv"))
    return(Reactome)
  }
}

gseaReactome <- function(x,
                         gene_up = symbolList_up,
                         gene_down = symbolList_down,
                         saveDir = ".",
                         sampleName = sample,
                         pvalue = 0.05,
                         min = 10,
                         max = 500) {
  # gene set enrichment analysis using Reactome database
  # x:              the list of genes that is not filtered by difference value is in the form of a string of logFC value vectors, sorted in descending order, with the gene name (Entrez ID) as the name of each item
  # gene_up/down:   a up/down subset of the differential gene list (Symbol), which is used to generate the up/down gene name and number after the results
  # pvalue:         p valve cut-off value
  # min/max:        gene set size limit to screen for analysis (screening out overly general gene sets)
  # organism:       organism species information, mouse is "mouse"
  
  suppressPackageStartupMessages(require(clusterProfiler))
  suppressPackageStartupMessages(require(ReactomePA))
  suppressPackageStartupMessages(require(org.Mm.eg.db))
  set.seed(123)
  gseaReactome <- try(gsePathway(
    geneList = x,
    organism = "mouse",
    pvalueCutoff = pvalue,
    pAdjustMethod = "BH",
    exponent = 1,
    minGSSize = min,
    maxGSSize = max,
    eps = 0,
    verbose = FALSE,
    seed = TRUE,
    by = "fgsea",
  ))
  if ("try-error" %in% class(gseaReactome)) {
    cat(
      "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
      paste("Sample:", sampleName),
      "Database: Reactome\tGSEA",
      gseaReactome,
      "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
      sep = "\n",
      file = "error_enrichment.log",
      append = TRUE
    )
    return(NULL)
  } else if (is.null(gseaReactome) || nrow(gseaReactome) == 0) {
    cat(
      "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
      paste("Sample:", sampleName),
      "Database: Reactome\tGSEA",
      "The enrichment result is empty, nothing is enriched",
      "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
      sep = "\n",
      file = "error_enrichment.log",
      append = TRUE
    )
    return(NULL)
  } else{
    gseaReactome = setReadable(gseaReactome, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    
    if (!is.null(gene_up) & !is.null(gene_down)) {
      list <- strsplit(gseaReactome@result$core_enrichment, "/")
      gseaReactome@result$geneUp <-
        purrr::map_chr(list,  ~ paste(intersect(.x, gene_up), collapse = "/"))
      gseaReactome@result$geneDown <-
        purrr::map_chr(list,  ~ paste(intersect(.x, gene_down), collapse = "/"))
      gseaReactome@result$Count <-
        unlist(purrr::map(list,  ~ length(.x)))
      gseaReactome@result$up <-
        unlist(purrr::map(list,  ~ length(intersect(.x, gene_up))))
      gseaReactome@result$down <-
        unlist(purrr::map(list,  ~ length(intersect(.x, gene_down))))
      # add the Symbol and number of core genes in UP and DOWN at the end of the enrichment result
    }
    
    saveDir <- paste0(saveDir, "/")
    if (!dir.exists(saveDir)) {
      dir.create(saveDir, recursive = TRUE)
    }
    write.csv(gseaReactome,
              file = paste0(saveDir,
                            sampleName,
                            "_gseaReactome.csv"))
    return(gseaReactome)
  }
}


# ·· draw function ------------------------------------------------------------

barPlot <-
  function(x,
           term = nrow(x),
           lateralAxis = "Count",
           labelLimit = 50,
           saveDir = ".",
           sampleName = sample,
           objectName = deparse1(substitute(x)),
           width = 8,
           height = 12,
           ...) {
    # x:              over representation enrichment analysis result object
    # term:           pick the terms to be used for drawing, which can be the first n lines or a specified term names string variable, and all lines are taken by default
    # lateralAxis:    y-axis, choose "Count" or "GeneRatio"
    # labelLimit:     set the length limit of the x-axis label, and if the label is too long, it will be folded
    # objectName:     used to distinguish between result files to avoid duplicate names
    # width/height:   the width and height of the output image
    
    if (is.null(x) || nrow(x) < 3) {
      cat(
        "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
        paste("Sample:", sampleName),
        paste("Warning:", objectName, "the number of rows is too low (<3), exit:barPlot"),
        "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
        sep = "\n",
        file = paste0("error_enrichment_plot.log"),
        append = TRUE
      )
      return(NULL)
    }
    
    suppressPackageStartupMessages(require(enrichplot))
    suppressPackageStartupMessages(require(ggplot2))
    
    p <- barplot(x,
                 showCategory = term,
                 x = lateralAxis,
                 label_format = labelLimit,
                 ...)
    # ... include font.size, title, color = "p.adjust"
    
    if (!dir.exists(saveDir)) {
      dir.create(saveDir, recursive = TRUE)
    }
    
    ggsave(
      filename = paste0(sampleName,
                        "_",
                        objectName,
                        "_bar.png"),
      path = saveDir,
      plot = p,
      bg = "white",
      width = width,
      height = height,
      units = "in",
      dpi = 300,
      limitsize = FALSE
    )
    
    ggsave(
      filename = paste0(sampleName,
                        "_",
                        objectName,
                        "_bar.tiff"),
      path = saveDir,
      plot = p,
      device = "tiff",
      compression = "lzw",
      bg = "white",
      width = width,
      height = height,
      units = "in",
      dpi = 300,
      limitsize = FALSE
    )
    
    ggsave(
      filename = paste0(sampleName,
                        "_",
                        objectName,
                        "_bar.svg"),
      path = saveDir,
      plot = p,
      bg = "white",
      width = width,
      height = height,
      units = "in",
      limitsize = FALSE
    )
    
    ggsave(
      filename = paste0(sampleName,
                        "_",
                        objectName,
                        "_bar.eps"),
      path = saveDir,
      plot = p,
      bg = "white",
      width = width,
      height = height,
      units = "in",
      limitsize = FALSE
    )
  }

barPlot1 <- function(x,
                     num = 10,
                     labelLimit = 50,
                     saveDir = ".",
                     sampleName = sample,
                     objectName = deparse1(substitute(x)),
                     width = 8,
                     height = 6) {
  # x:              over representation enrichment analysis result object
  # num:            pick the terms used for drawing, with the first n lines of each group
  # labelLimit:     set the length limit of the x-axis label, and if the label is too long, it will be folded
  # objectName:     used to distinguish between result files to avoid duplicate names
  # width/height:   the width and height of the output image
  
  if (is.null(x) || nrow(x) < 3) {
    cat(
      "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
      paste("Sample:", sampleName),
      paste("Warning:", objectName, "the number of rows is too low (<3), exit:barPlot1"),
      "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
      sep = "\n",
      file = paste0("error_enrichment_plot.log"),
      append = TRUE
    )
    return(NULL)
  }
  
  suppressPackageStartupMessages(require(dplyr))
  suppressPackageStartupMessages(require(tidyr))
  suppressPackageStartupMessages(require(stringr))
  suppressPackageStartupMessages(require(ggplot2))
  
  p <-
    x@result %>% drop_na(., ONTOLOGY) %>% group_by(ONTOLOGY) %>% slice_head(n = num) %>%  {
      ggplot(., aes(
        x = factor(Description, levels = Description),
        y = -log10(p.adjust),
        fill = ONTOLOGY
      )) + geom_col() +
        xlab("Description") + geom_text(aes(label = Count, vjust = -0.5)) +
        scale_x_discrete(
          labels = function(x)
            str_wrap(x, width = labelLimit)
        ) + scale_fill_manual(values = c("#D73027", "#00BA38", "#4575B4")) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        theme_classic() + theme(axis.text.x =
                                  element_text(
                                    face = "plain",
                                    color = "black",
                                    angle = 70,
                                    vjust = 1,
                                    hjust = 1
                                  )) + theme(plot.margin = unit(c(0.2, 0.1, 0.2, 0.1), "inches"))
      
    }
 
  if (!dir.exists(saveDir)) {
    dir.create(saveDir, recursive = TRUE)
  }
  
  ggsave(
    filename = paste0(sampleName,
                      "_",
                      objectName,
                      "_bar1.png"),
    path = saveDir,
    plot = p,
    bg = "white",
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )
  
  ggsave(
    filename = paste0(sampleName,
                      "_",
                      objectName,
                      "_bar1.tiff"),
    path = saveDir,
    plot = p,
    device = "tiff",
    compression = "lzw",
    bg = "white",
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )  
  ggsave(
    filename = paste0(sampleName,
                      "_",
                      objectName,
                      "_bar1.svg"),
    path = saveDir,
    plot = p,
    bg = "white",
    width = width,
    height = height,
    units = "in",
    limitsize = FALSE
  )
  
  ggsave(
    filename = paste0(sampleName,
                      "_",
                      objectName,
                      "_bar1.eps"),
    path = saveDir,
    plot = p,
    bg = "white",
    width = width,
    height = height,
    units = "in",
    limitsize = FALSE
  )
}

barPlot2 <- function(x,
                     term = 15,
                     labelLimit = 50,
                     saveDir = ".",
                     sampleName = sample,
                     objectName = deparse1(substitute(x)),
                     width = 8,
                     height = 6) {
  # x:              over representation enrichment analysis result object
  # term:           pick the terms to be used for drawing, which can be the first n lines or a specified term names string variable, and first 15 lines are taken by default
  # labelLimit:     set the length limit of the x-axis label, and if the label is too long, it will be folded
  # objectName:     used to distinguish between result files to avoid duplicate names
  # width/height:   the width and height of the output image
  
  if (is.null(x) || nrow(x) < 3) {
    cat(
      "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
      paste("Sample:", sampleName),
      paste("Warning:", objectName, "the number of rows is too low (<3), exit:barPlot2"),
      "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
      sep = "\n",
      file = paste0("error_enrichment_plot.log"),
      append = TRUE
    )
    return(NULL)
  }
  
  suppressPackageStartupMessages(require(tidyr))
  suppressPackageStartupMessages(require(stringr))
  
  
  if (is.numeric(term)) {
    if (nrow(x) < term) {
      term <- nrow(x)
    }
    data <-
      x@result[1:term, c("ID", "Description", "up", "down")] %>% pivot_longer(cols =
                                                                                c(3, 4),
                                                                              names_to = "group",
                                                                              values_to = "counts")
  } else{
    data <-
      x@result[which(x@result$Description %in% term), c("ID", "Description", "up", "down")] %>% pivot_longer(cols =
                                                                                                               c(3, 4),
                                                                                                             names_to = "group",
                                                                                                             values_to = "counts")
  }
  
  p <- ggplot(data, aes(
    x = factor(Description, levels = rev(unique(Description))),
    y = ifelse(group
               == "up", counts,-counts),
    fill = factor(group, levels = c("up", "down"))
  )) + geom_col() + coord_flip() + labs(
    fill = "",
    x = "",
    y = "Gene number",
    title = ""
  ) + geom_text(aes(label = counts, hjust = ifelse(group == "up",-0.5, 1.5))) + scale_x_discrete(
    labels = function(x)
      str_wrap(x, width = labelLimit),
    expand =
      expansion(mult = c(0.1, 0.1))
  ) + scale_fill_manual(values = c("up" = "#D73027", "down" = "#4575B4")) +
    scale_y_continuous(expand =
                         expansion(mult = c(0.1, 0.1)), labels = abs) + theme_classic() + theme(
                           axis.line.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           plot.margin = unit(c(0.2, 0.1, 0.2, 0.1), "inches")
                         )

  if (!dir.exists(saveDir)) {
    dir.create(saveDir, recursive = TRUE)
  }
  ggsave(
    filename = paste0(sampleName,
                      "_",
                      objectName,
                      "_bar2.png"),
    path = saveDir,
    plot = p,
    bg = "white",
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )
  
  ggsave(
    filename = paste0(sampleName,
                      "_",
                      objectName,
                      "_bar2.tiff"),
    path = saveDir,
    plot = p,
    device = "tiff",
    compression = "lzw",
    bg = "white",
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )
  
  ggsave(
    filename = paste0(sampleName,
                      "_",
                      objectName,
                      "_bar2.svg"),
    path = saveDir,
    plot = p,
    bg = "white",
    width = width,
    height = height,
    units = "in",
    limitsize = FALSE
  )
  
  ggsave(
    filename = paste0(sampleName,
                      "_",
                      objectName,
                      "_bar2.eps"),
    path = saveDir,
    plot = p,
    bg = "white",
    width = width,
    height = height,
    units = "in",
    limitsize = FALSE
  )
}

dotPlot <-
  function(x,
           term = nrow(x),
           lateralAxis = "GeneRatio",
           labelLimit = 50,
           saveDir = ".",
           sampleName = sample,
           objectName = deparse1(substitute(x)),
           width = 8,
           height = 12,
           ...) {
    # x:              enrichment analysis result object
    # term:           pick the terms to be used for drawing, which can be the first n lines or a specified term names string variable, and all lines are taken by default
    # lateralAxis:    y-axis, choose "Count" or "GeneRatio"
    # labelLimit:     set the length limit of the x-axis label, and if the label is too long, it will be folded
    # objectName:     used to distinguish between result files to avoid duplicate names
    # width/height:   the width and height of the output image
    
    if (is.null(x) || nrow(x) < 3) {
      cat(
        "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
        paste("Sample:", sampleName),
        paste("Warning:", objectName, "the number of rows is too low (<3), exit:dotPlot"),
        "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
        sep = "\n",
        file = paste0("error_enrichment_plot.log"),
        append = TRUE
      )
      return(NULL)
    }
    
    suppressPackageStartupMessages(require(enrichplot))
    suppressPackageStartupMessages(require(ggplot2))
    
    p <- dotplot(x,
                 showCategory = term,
                 x = lateralAxis,
                 label_format = labelLimit,
                 ...)
    # ... include font.size、title、color = "p.adjust"、size
    
    if (!dir.exists(saveDir)) {
      dir.create(saveDir, recursive = TRUE)
    }
    ggsave(
      filename = paste0(sampleName,
                        "_",
                        objectName,
                        "_dot.png"),
      path = saveDir,
      plot = p,
      bg = "white",
      width = width,
      height = height,
      units = "in",
      dpi = 300,
      limitsize = FALSE
    )
    
    ggsave(
      filename = paste0(sampleName,
                        "_",
                        objectName,
                        "_dot.tiff"),
      path = saveDir,
      plot = p,
      device = "tiff",
      compression = "lzw",
      bg = "white",
      width = width,
      height = height,
      units = "in",
      dpi = 300,
      limitsize = FALSE
    )
    
    ggsave(
      filename = paste0(sampleName,
                        "_",
                        objectName,
                        "_dot.svg"),
      path = saveDir,
      plot = p,
      bg = "white",
      width = width,
      height = height,
      units = "in",
      limitsize = FALSE
    )
    
    ggsave(
      filename = paste0(sampleName,
                        "_",
                        objectName,
                        "_dot.eps"),
      path = saveDir,
      plot = p,
      bg = "white",
      width = width,
      height = height,
      units = "in",
      limitsize = FALSE
    )
  }

cnetPlot <-
  function(x,
           term = nrow(x),
           colorEdge = TRUE,
           saveDir = ".",
           sampleName = sample,
           objectName = deparse1(substitute(x)),
           width = 16,
           height = 12,
           ...) {
    # x:              enrichment analysis result object
    # term:           pick the terms to be used for drawing, which can be the first n lines or a specified term names string variable, and all lines are taken by default
    # objectName:     used to distinguish between result files to avoid duplicate names
    # width/height:   the width and height of the output image
    
    if (is.null(x) || nrow(x) < 3) {
      cat(
        "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
        paste("Sample:", sampleName),
        paste("Warning:", objectName, "the number of rows is too low (<3), exit:centPlot"),
        "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
        sep = "\n",
        file = paste0("error_enrichment_plot.log"),
        append = TRUE
      )
      return(NULL)
    }
    
    if (is(x, "gseaResult")) {
      no_zero_index <- which(x@result$Count[1:term] != 0 & x@result$Count[1:term] != 1)
      if (length(no_zero_index) > 3) {
        term <- x@result$Description[no_zero_index]
      }else {
        return(NULL)
      }
    }
    # The number of core-enrichment genes may be 0 or 1, resulting in a drawing error, and the term that is not 0 or 1 is taken here
    
    suppressPackageStartupMessages(require(enrichplot))
    suppressPackageStartupMessages(require(ggplot2))
    
    p <- cnetplot(
      x,
      showCategory = term,
      color.params = list(edge = colorEdge),
      cex.params = list(category_label = 1.2, gene_label = 0.7),
      ...
    )
    # ... include color.params (foldChange, edge, category, gene), cex.params (foldChange, category_node, gene_node, category_label, gene_label), layout and node_label

    if (!dir.exists(saveDir)) {
      dir.create(saveDir, recursive = TRUE)
    }
    ggsave(
      filename = paste0(sampleName,
                        "_",
                        objectName,
                        "_cnet.png"),
      path = saveDir,
      plot = p,
      bg = "white",
      width = width,
      height = height,
      units = "in",
      dpi = 300,
      limitsize = FALSE
    )
    
    ggsave(
      filename = paste0(sampleName,
                        "_",
                        objectName,
                        "_cnet.tiff"),
      path = saveDir,
      plot = p,
      device = "tiff",
      compression = "lzw",
      bg = "white",
      width = width,
      height = height,
      units = "in",
      dpi = 300,
      limitsize = FALSE
    )
    
    ggsave(
      filename = paste0(sampleName,
                        "_",
                        objectName,
                        "_cnet.svg"),
      path = saveDir,
      plot = p,
      bg = "white",
      width = width,
      height = height,
      units = "in",
      limitsize = FALSE
    )
    
    ggsave(
      filename = paste0(sampleName,
                        "_",
                        objectName,
                        "_cnet.eps"),
      path = saveDir,
      plot = p,
      device = cairo_ps,
      fallback_resolution = 600,
      bg = "white",
      width = width,
      height = height,
      units = "in",
      limitsize = FALSE
    )
  }

testOnt <- function(x) {
  print(summary(factor(x@result$ONTOLOGY)))
}

dagPlot <-
  function(x,
           firstSigNodes = 10,
           saveDir = ".",
           useInfo = "all",
           sigForAll = TRUE,
           useFullNames = TRUE,
           sampleName = sample,
           objectName = deparse1(substitute(x)),
           width = 12,
           height = 8) {
    # This function is modified from clusterProfiler::plotGOgraph()
    # x:              enrichment analysis result object, which was filted by dplyr::filter(ONTOLOGY == "BP") (or "CC", "MF")
    # firstSigNodes:  specify the number of colored items in the plot (i.e., the term actually used for drawing) (marked is the enrichment result, unmarked is the parent related term)
    # useInfo:        whether to show p.adjust ("pval") and enrichment genes/term genes ("counts"), or both ("all")
    # objectName:     used to distinguish between result files to avoid duplicate names
    # width/height:   the width and height of the output image
    
    suppressPackageStartupMessages(require(topGO))
    suppressPackageStartupMessages(library(clusterProfiler))
    
    if (!class(x) %in% c("gseaResult", "enrichResult")) {
      stop("x should be output of gseGO or enrichGO...")
    }
    
    gs <- x@geneSets
    gs.df <- data.frame(gene = unlist(gs),
                        go   = rep(names(gs),
                                   times = sapply(gs, length)))
    gene2GO <- split(as.character(gs.df$go), gs.df$gene)
    
    if (is(x, "gseaResult")) {
      ont <- x@result$ONTOLOGY[1]
      allgenes <- x@geneList
      core_genes <- unique(unlist(geneInCategory(x)))
      allgenes[!names(allgenes) %in% core_genes] <- -1
      allgenes[core_genes] <- 1
    } else {
      ont <- x@result$ONTOLOGY[1]
      universe <- x@universe
      allgenes <- numeric(length(universe))
      names(allgenes) <- universe
      allgenes[x@gene] <- 1
    }
    
    selector <- function(scores)
      return(scores == 1)
    
    if (!ont %in% c("BP", "MF", "CC")) {
      stop("ONTOLOGY should be one of 'BP', 'MF' or 'CC'...")
    }
    
    
    pvalue <- x@result$p.adjust
    names(pvalue) <- x@result$ID
    
    groupGOTerms()
    
    GOdata <- new(
      "topGOdata",
      description = "clusterProfiler enrichment results",
      ontology = ont,
      allGenes = allgenes,
      geneSel = selector,
      annot = annFUN.gene2GO,
      gene2GO = gene2GO
    )
    
    firstSigNodes <- min(firstSigNodes, nrow(x))
    
    png(
      filename = paste0(saveDir,
                        "/",
                        sampleName,
                        "_",
                        objectName,
                        "_",
                        ont,
                        "_DAG.png"),
      bg = "white",
      width = width,
      height = height,
      units = "in",
      res = 300
    )
    showSigOfNodes(
      GOdata        = GOdata,
      termsP.value  = pvalue,
      firstSigNodes = firstSigNodes,
      useInfo       = useInfo,
      sigForAll     = sigForAll,
      useFullNames  = useFullNames
    )
    dev.off()
    
    tiff(
      filename = paste0(saveDir,
                        "/",
                        sampleName,
                        "_",
                        objectName,
                        "_",
                        ont,
                        "_DAG.tiff"),
      bg = "white",
      width = width,
      height = height,
      units = "in",
      res = 300,
      compression = "lzw"
    )
    showSigOfNodes(
      GOdata        = GOdata,
      termsP.value  = pvalue,
      firstSigNodes = firstSigNodes,
      useInfo       = useInfo,
      sigForAll     = sigForAll,
      useFullNames  = useFullNames
    )
    dev.off()
    
    svg(
      filename = paste0(saveDir,
                        "/",
                        sampleName,
                        "_",
                        objectName,
                        "_",
                        ont,
                        "_DAG.svg"),
      bg = "white",
      width = width,
      height = height,
    )
    showSigOfNodes(
      GOdata        = GOdata,
      termsP.value  = pvalue,
      firstSigNodes = firstSigNodes,
      useInfo       = useInfo,
      sigForAll     = sigForAll,
      useFullNames  = useFullNames
    )
    dev.off()
    
    pdf(
      file = paste0(saveDir,
                    "/",
                    sampleName,
                    "_",
                    objectName,
                    "_",
                    ont,
                    "_DAG.pdf"),
      bg = "white",
      width = width,
      height = height
    )
    showSigOfNodes(
      GOdata        = GOdata,
      termsP.value  = pvalue,
      firstSigNodes = firstSigNodes,
      useInfo       = useInfo,
      sigForAll     = sigForAll,
      useFullNames  = useFullNames
    )
    dev.off()
    
    setEPS()
    postscript(
      file = paste0(saveDir,
                    "/",
                    sampleName,
                    "_",
                    objectName,
                    "_",
                    ont,
                    "_DAG.eps"),
      bg = "white",
      width = width,
      height = height
    )
    showSigOfNodes(
      GOdata        = GOdata,
      termsP.value  = pvalue,
      firstSigNodes = firstSigNodes,
      useInfo       = useInfo,
      sigForAll     = sigForAll,
      useFullNames  = useFullNames
    )
    dev.off()
  }

ESPlot <-
  function(x,
           saveDir = ".",
           sampleName = sample,
           objectName = deparse1(substitute(x)),
           ...) {
    # x:              gene set enrichment analysis result object
    # objectName:     used to distinguish between result files to avoid duplicate names
    
    suppressPackageStartupMessages(require(enrichplot))
    suppressPackageStartupMessages(require(ggplot2))
    suppressPackageStartupMessages(require(stringr))
    
    if (!dir.exists(paste0(saveDir, "/ESplot/"))) {
      dir.create(paste0(saveDir, "/ESplot/"), recursive = TRUE)
    }
    # default create folder "ESplot" to save plots
    
    for (i in seq_along(x$ID)) {
      p <- gseaplot2(
        x,
        geneSetID = i,
        title = paste0(x@result$ID[i], " : ", x@result$Description[i]),
        pvalue_table = TRUE,
        ...
      )
      # ... include color, base_size, subplots, rel_heights, pvalue_table, ES_geom
     
      pdf(
        file = paste0(
          saveDir,
          "/ESplot/",
          sampleName,
          "_",
          objectName,
          "_",
          i,
          "_",
          str_replace(x$ID[i], ":", "-"),
          ".pdf"
        ),
        width = 12,
        height = 12,
        bg = "white"
      )
      print(p)
      dev.off()
      # because ":" cannot appear in the file name, replace the ":" in the term name with "-"
      
      png(
        filename = paste0(
          saveDir,
          "/ESplot/",
          sampleName,
          "_",
          objectName,
          "_",
          i,
          "_",
          str_replace(x$ID[i], ":", "-"),
          ".png"
        ),
        width = 12,
        height = 12,
        units = "in",
        res = 600,
        bg = "white"
      )
      print(p)
      dev.off()
      
      tiff(
        filename = paste0(
          saveDir,
          "/ESplot/",
          sampleName,
          "_",
          objectName,
          "_",
          i,
          "_",
          str_replace(x$ID[i], ":", "-"),
          ".tiff"
        ),
        width = 12,
        height = 12,
        units = "in",
        res = 600,
        bg = "white",
        compression = "lzw"
      )
      print(p)
      dev.off()
    }
  }

ridgePlot <-
  function(x,
           num = nrow(x),
           labelLimit = 50,
           decreasing = FALSE,
           saveDir = ".",
           sampleName = sample,
           objectName = deparse1(substitute(x)),
           width = 8,
           height = 12,
           ...) {
    # x:              over representation enrichment analysis result object
    # num:            pick the terms to be used for drawing, which can be the first n lines, and all lines are taken by default
    # labelLimit:     set the length limit of the x-axis label, and if the label is too long, it will be folded
    # decreasing:     terms are sorted in ascending or descending order (the default is ascending, i.e. positive from top to negative bottom)
    # objectName:     used to distinguish between result files to avoid duplicate names
    # width/height:   the width and height of the output image
    
    if (is.null(x) || nrow(x) < 3) {
      cat(
        "﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊﹊",
        paste("Sample:", sampleName),
        paste("Warning:", objectName, "the number of rows is too low (<3), exit:barPlot1"),
        "﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎﹎\n",
        sep = "\n",
        file = paste0("error_enrichment_plot.log"),
        append = TRUE
      )
      return(NULL)
    }
    
    suppressPackageStartupMessages(require(enrichplot))
    suppressPackageStartupMessages(require(ggplot2))
    
    p <- ridgeplot(
      x,
      showCategory = num,
      label_format = labelLimit,
      decreasing = decreasing,
      ...
    )
    # ... include fill = "p.adjust", core_enrichment = TRUE, orderBy = "NES"

    if (!dir.exists(saveDir)) {
      dir.create(saveDir, recursive = TRUE)
    }
    
    ggsave(
      filename = paste0(sampleName,
                        "_",
                        objectName,
                        "_ridge.png"),
      path = saveDir,
      plot = p,
      bg = "white",
      width = width,
      height = height,
      units = "in",
      dpi = 300,
      limitsize = FALSE
    )
    
    ggsave(
      filename = paste0(sampleName,
                        "_",
                        objectName,
                        "_ridge.tiff"),
      path = saveDir,
      plot = p,
      device = "tiff",
      compression = "lzw",
      bg = "white",
      width = width,
      height = height,
      units = "in",
      dpi = 300,
      limitsize = FALSE
    )
    
    ggsave(
      filename = paste0(sampleName,
                        "_",
                        objectName,
                        "_ridge.svg"),
      path = saveDir,
      plot = p,
      bg = "white",
      width = width,
      height = height,
      units = "in",
      limitsize = FALSE
    )
    
    ggsave(
      filename = paste0(sampleName,
                        "_",
                        objectName,
                        "_ridge.eps"),
      path = saveDir,
      plot = p,
      bg = "white",
      width = width,
      height = height,
      units = "in",
      limitsize = FALSE
    )
  }

pathPlot <- function(x,
                     diffGene,
                     class = "ORA",
                     limit = NULL,
                     tmpDir = "./genome/PathPlot",
                     saveDir = "PathPlot",
                     sampleName = sample,
                     objectName = deparse1(substitute(x)),
                     ...) {
  # KEGG pathway maps
  # x:              KEGG result (enrichKEGG object, gseKEGG object or dataframe)
  # diffGene:       the list of differential genes, with at least three column: "Symbol", "Entrez Gene ID" and "logFC"
  # class:          "ORA" or "GSEA", used to distinguish between the two enrichment methods to extract the enrichment genes column
  # limit:          limitation of numeric values in the legend, default is the maximum difference value (affects coloring)
  # tmpDir:         KEGG pathway maps cache directory for reuse
  # saveDir:        save address, default will create directory "PathPlot"
  # objectName:     used to distinguish between result files to avoid duplicate names
  
  suppressPackageStartupMessages(require(pathview))
  suppressPackageStartupMessages(require(purrr))
  suppressPackageStartupMessages(require(dplyr))
  suppressPackageStartupMessages(require(tibble))
  
  if (!dir.exists(saveDir)) {
    dir.create(saveDir, recursive = TRUE)
  }
  
  wd <- getwd()
  # save the working directory address for recovery
  on.exit(setwd(wd))
  # restore the original working directory when you exit
  
  tmpDir <- normalizePath(tmpDir)
  if (!dir.exists(tmpDir)) {
    dir.create(tmpDir, recursive = TRUE)
  }
  # KEGG pathway maps cache directory (not the same folder as saveDir's "PathPlot")
  
  setwd(saveDir)
  # reset the working directory to save files using the pathview package
  
  
  if (class == "ORA") {
    geneIDlist <- strsplit(x$geneID, "/")
  } else if (class == "GSEA") {
    geneIDlist <- strsplit(x$core_enrichment, "/")
  }
  # read enrichment genes list or core-enrichment genes list
  
  geneDataList <-
    purrr::map(
      geneIDlist,
      ~ filter(diffGene, Symbol %in% .) %>% .[, c("Entrez Gene ID", "logFC")] %>% remove_rownames(.) %>% column_to_rownames(., var =
                                                                                                                                       "Entrez Gene ID")
    )
  # the differential values corresponding to the enrichment genes were extracted from the list of differential genes and used for coloring the pathway map
  
  for (i in seq_along(x$ID)) {
    if (is.null(geneDataList) || nrow(geneDataList[[i]]) == 0) {
      next
    }
    if (is.null(limit)) {
      limitTmp <- list(gene = max(abs(geneDataList[[i]])), cpd = 1)
    } else{
      limitTmp <- list(gene = limit, cpd = 1)
    }
    # If the limit parameter is not provided, the numeric limit of the color legend uses the maximum difference value of the differential genes, and if it is provided, the provided value is used
    
    p <- pathview(
      gene.data  = geneDataList[[i]],
      pathway.id = x$ID[i],
      species    = "mmu",
      kegg.dir   = tmpDir,
      same.layer = FALSE,
      out.suffix = paste(sampleName, objectName, sep = "_"),
      limit      = limitTmp,
      ...
    )
    
    file.copy(from = paste0(tmpDir, "/", x$ID[i], ".png"),
              to = "./")
    # copy the original maps to the result folder
  }
}


# · GO ORA --------------------------------------------------------------------

saveDir <- paste0(savePath, "GO/ORA")

for (i in c("all", "up", "down")) {
  x <- paste0("ensembl_", i)
  assign(
    paste0("GO_", i),
    oraGO(
      get(x),
      objectName = paste0("GO_", i),
      gene_up = symbol_up,
      gene_down = symbol_down,
      saveDir = paste(saveDir, i, sep = "/")
    )
  )
  assign(paste0("GO_", i, "_simply"),
         simply(
           get(paste0("GO_", i)),
           saveDir = paste(saveDir, i, sep = "/"),
           objectName = paste0("GO_", i)
         ))
  
  barPlot(
    get(paste0("GO_", i, "_simply")),
    term = 30,
    objectName = paste0("GO_", i),
    saveDir = paste(saveDir, i, sep = "/")
  )
  barPlot1(
    get(paste0("GO_", i, "_simply")),
    objectName = paste0("GO_", i),
    saveDir = paste(saveDir, i, sep = "/")
  )
  if (i == "all") {
    barPlot2(
      get(paste0("GO_", i, "_simply")),
      term = 15,
      objectName = paste0("GO_", i),
      saveDir = paste(saveDir, i, sep = "/")
    )
  }
  dotPlot(
    get(paste0("GO_", i, "_simply")),
    term = 30,
    objectName = paste0("GO_", i),
    saveDir = paste(saveDir, i, sep = "/")
  )
  cnetPlot(
    get(paste0("GO_", i, "_simply")),
    term = 15,
    objectName = paste0("GO_", i),
    saveDir = paste(saveDir, i, sep = "/")
  )
  
  if (!is.null(get(paste0("GO_", i, "_simply")))) {
    test <- testOnt(get(paste0("GO_", i, "_simply")))
    # Test the distribution of ONTOLOGY
    
    if ("BP" %in% names(test)) {
      get(paste0("GO_", i, "_simply")) %>% filter(ONTOLOGY == "BP") %>% dagPlot(
        .,
        firstSigNodes = 10,
        objectName = paste0("GO_", i),
        saveDir = paste(saveDir, i, sep = "/")
      )
    }
    if ("CC" %in% names(test)) {
      get(paste0("GO_", i, "_simply")) %>% filter(ONTOLOGY == "CC") %>% dagPlot(
        .,
        firstSigNodes = 10,
        objectName = paste0("GO_", i),
        saveDir = paste(saveDir, i, sep = "/")
      )
    }
    if ("MF" %in% names(test)) {
      get(paste0("GO_", i, "_simply")) %>% filter(ONTOLOGY == "MF") %>% dagPlot(
        .,
        firstSigNodes = 10,
        objectName = paste0("GO_", i),
        saveDir = paste(saveDir, i, sep = "/")
      )
    }
    rm(test)
  }
}
rm(i)
rm(x)


# · GO GSEA -------------------------------------------------------------------

saveDir <- paste0(savePath, "GO/GSEA")

GO_gsea <- gseaGO(
  ensemblList,
  gene_up = symbolList_up,
  gene_down = symbolList_down,
  saveDir = saveDir
)
GO_gsea_simply <- simply(GO_gsea,
                         saveDir = saveDir, objectName = "gseaGO")
dotPlot(
  GO_gsea_simply,
  term = 30,
  saveDir = saveDir,
  objectName = "gseaGO"
)

cnetPlot(
  GO_gsea_simply,
  term = 15,
  saveDir = saveDir,
  objectName = "gseaGO"
)

if (!is.null(GO_gsea_simply)) {
  test <- testOnt(GO_gsea_simply)
  
  if ("BP" %in% names(test)) {
    GO_gsea_simply %>% filter(ONTOLOGY == "BP") %>% dagPlot(
      .,
      firstSigNodes = 10,
      objectName = "gseaGO",
      saveDir = saveDir
    )
  }
  if ("CC" %in% names(test)) {
    GO_gsea_simply %>% filter(ONTOLOGY == "CC") %>% dagPlot(
      .,
      firstSigNodes = 10,
      objectName = "gseaGO",
      saveDir = saveDir
    )
  }
  if ("MF" %in% names(test)) {
    GO_gsea_simply %>% filter(ONTOLOGY == "MF") %>% dagPlot(
      .,
      firstSigNodes = 10,
      objectName = "gseaGO",
      saveDir = saveDir
    )
  }
  rm(test)
}

ridgePlot(
  GO_gsea_simply,
  num = 30,
  saveDir = saveDir,
  objectName = "gseaGO"
)

if (!is.null(GO_gsea_simply)) {
  if (length(which(GO_gsea_simply@result$p.adjust < 0.05)) > 10) {
    GO_gsea_simply %>% filter(p.adjust < 0.05) %>% ESPlot(., objectName = "gseaGO", saveDir = saveDir)
  } else{
    GO_gsea_simply %>% slice(., 1:10) %>% ESPlot(., objectName = "gseaGO", saveDir = saveDir)
  }
}


# · KEGG ORA ------------------------------------------------------------------

saveDir <- paste0(savePath, "KEGG/ORA")
for (i in c("all", "up", "down")) {
  x <- paste0("entrez_", i)
  assign(
    paste0("KEGG_", i),
    oraKEGG(
      get(x),
      pvalue = 1,
      qvalue = 1,
      objectName = paste0("KEGG_", i),
      gene_up = symbol_up,
      gene_down = symbol_down,
      saveDir = paste(saveDir, i, sep = "/")
    )
  )
  
  barPlot(
    get(paste0("KEGG_", i)),
    term = 20,
    saveDir = paste(saveDir, i, sep = "/"),
    objectName = paste0("KEGG_", i)
  )
  
  if (i == "all") {
    barPlot2(
      get(paste0("KEGG_", i)),
      term = 15,
      objectName = paste0("KEGG_", i),
      saveDir = paste(saveDir, i, sep = "/")
    )
  }
  
  dotPlot(
    get(paste0("KEGG_", i)),
    term = 20,
    saveDir = paste(saveDir, i, sep = "/"),
    objectName = paste0("KEGG_", i)
  )
  
  cnetPlot(
    get(paste0("KEGG_", i)),
    term = 15,
    saveDir = paste(saveDir, i, sep = "/"),
    objectName = paste0("KEGG_", i)
  )
  
  if (!is.null(get(paste0("KEGG_", i)))) {
    diffGene_notNA <- switch(
      i,
      all = drop_na(diffGene, Symbol, `Entrez Gene ID`, logFC),
      up = drop_na(diffGene[diffGene$logFC > 0,], Symbol, `Entrez Gene ID`, logFC),
      down = drop_na(diffGene[diffGene$logFC < 0,], Symbol, `Entrez Gene ID`, logFC)
    )
    # extract data to add color for KEGG pathway plot
    
    pathPlot(
      get(paste0("KEGG_", i))[1:max(length(which(get(
        paste0("KEGG_", i)
      )@result$p.adjust < 0.05)), 10),],
      diffGene = diffGene_notNA,
      class = "ORA",
      saveDir = paste(saveDir, i, "PathPlot", sep = "/"),
      objectName = paste0("KEGG_", i)
    )
    # take p.adjust < 0.05 (or the first 10 if there are less than 10) to draw some pathway plots
    rm(diffGene_notNA)
  }
}
rm(i)
rm(x)


# · KEGG GSEA -----------------------------------------------------------------

saveDir <- paste0(savePath, "KEGG/GSEA")

KEGG_gsea <-
  gseaKEGG(
    entrezList,
    gene_up = symbolList_up,
    gene_down = symbolList_down,
    pvalue = 1,
    saveDir = saveDir
  )

dotPlot(KEGG_gsea,
        term = 20,
        saveDir = saveDir,
        objectName = "gseaKEGG")

cnetPlot(KEGG_gsea,
         term = 15,
         saveDir = saveDir,
         objectName = "gseaKEGG")

ridgePlot(KEGG_gsea,
          num = 20,
          saveDir = saveDir,
          objectName = "gseaKEGG")

if (!is.null(KEGG_gsea)) {
  if (length(which(KEGG_gsea@result$p.adjust < 0.05)) > 10) {
    KEGG_gsea %>% clusterProfiler::filter(p.adjust < 0.05) %>% ESPlot(., objectName = "gseaKEGG", saveDir = saveDir)
  } else{
    KEGG_gsea %>% clusterProfiler::slice(., 1:10) %>% ESPlot(., objectName = "gseaKEGG", saveDir = saveDir)
  }
}

if (!is.null(KEGG_gsea)) {
  diffGene_notNA <-
    drop_na(allGene, Symbol, `Entrez Gene ID`, logFC)
  # extract data to add color for KEGG pathway plot
  
  pathPlot(
    KEGG_gsea[1:max(length(which(KEGG_gsea@result$p.adjust < 0.05)), 10), ],
    diffGene = diffGene_notNA,
    class = "GSEA",
    saveDir = paste(saveDir, "PathPlot", sep = "/"),
    objectName = "gseaKEGG"
  )
  # take p.adjust < 0.05 (or the first 10 if there are less than 10) to draw some pathway plots
  rm(diffGene_notNA)
}


# · Reactome ORA --------------------------------------------------------------

saveDir <- paste0(savePath, "Reactome/ORA")

for (i in c("all", "up", "down")) {
  x <- paste0("entrez_", i)
  assign(
    paste0("Reactome_", i),
    oraReactome(
      get(x),
      objectName = paste0("Reactome_", i),
      gene_up = symbol_up,
      gene_down = symbol_down,
      saveDir = paste(saveDir, i, sep = "/"),
    )
  )
  
  barPlot(
    get(paste0("Reactome_", i)),
    term = 20,
    saveDir = paste(saveDir, i, sep = "/"),
    objectName = paste0("Reactome_", i)
  )
  
  if (i == "all") {
    barPlot2(
      get(paste0("Reactome_", i)),
      term = 15,
      objectName = paste0("Reactome_", i),
      saveDir = paste(saveDir, i, sep = "/")
    )
  }
  
  dotPlot(
    get(paste0("Reactome_", i)),
    term = 20,
    saveDir = paste(saveDir, i, sep = "/"),
    objectName = paste0("Reactome_", i)
  )
  
  cnetPlot(
    get(paste0("Reactome_", i)),
    term = 15,
    saveDir = paste(saveDir, i, sep = "/"),
    objectName = paste0("Reactome_", i)
  )
}
rm(i)
rm(x)


# · Reactome GSEA -------------------------------------------------------------

saveDir <- paste0(savePath, "Reactome/GSEA")

Reactome_gsea <-
  gseaReactome(
    entrezList,
    gene_up = symbolList_up,
    gene_down = symbolList_down,
    saveDir = saveDir
  )

dotPlot(
  Reactome_gsea,
  term = 20,
  saveDir = saveDir,
  objectName = "gseaReactome"
)

cnetPlot(
  Reactome_gsea,
  term = 15,
  saveDir = saveDir,
  objectName = "gseaReactome"
)

ridgePlot(
  Reactome_gsea,
  num = 20,
  saveDir = saveDir,
  objectName = "gseaReactome"
)

if (!is.null(Reactome_gsea)) {
  if (length(which(Reactome_gsea@result$p.adjust < 0.05)) > 10) {
    Reactome_gsea %>% filter(p.adjust < 0.05) %>% ESPlot(., objectName = "gseaReactome", saveDir = saveDir)
  } else{
    Reactome_gsea %>% slice(., 1:10) %>% ESPlot(., objectName = "gseaReactome", saveDir = saveDir)
  }
}


# · DIY dataset ---------------------------------------------------------------

# ·· GO -----------------------------------------------------------------------

lapply(names(listPathway_GO_gene), function(x) {
  # enrichment analysis is performed separately for each function
  class <- stringr::str_replace(x, "_gene$", "")
  GO_gene <- listPathway_GO_gene[[x]]
  GO_term <-
    listPathway_GO_term[[stringr::str_replace(x, "_gene$", "_term")]]
  
  # ··· GO ORA ----------------------------------------------------------------
  
  saveDir <- paste0(savePath, "DIY/", class, "/GO_ORA/")
  
  for (i in c("all", "up", "down")) {
    x <- paste0("GO_diy_", class, "_", i)
    assign(
      paste0("GO_diy_", class, "_", i),
      enrichORA(
        get(paste0("symbol_", i)),
        keyType = "Symbol",
        dbType = "GO",
        objectName = paste0("GO_", i),
        gene_up = symbol_up,
        gene_down = symbol_down,
        term2gene = GO_gene,
        term2name = GO_term,
        max = 2000,
        saveDir = paste(saveDir, i, sep = "/")
      ),
      envir = .GlobalEnv
    )
    
    barPlot(
      get(x),
      term = 30,
      objectName = paste0("GO_", i),
      saveDir = paste(saveDir, i, sep = "/")
    )
    barPlot1(
      get(x),
      objectName = paste0("GO_", i),
      saveDir = paste(saveDir, i, sep = "/")
    )
    if (i == "all") {
      barPlot2(
        get(x),
        term = 15,
        objectName = paste0("GO_", i),
        saveDir = paste(saveDir, i, sep = "/")
      )
    }
    dotPlot(
      get(x),
      term = 30,
      objectName = paste0("GO_", i),
      saveDir = paste(saveDir, i, sep = "/")
    )
    cnetPlot(
      get(x),
      term = 15,
      objectName = paste0("GO_", i),
      saveDir = paste(saveDir, i, sep = "/")
    )
    
    if (!is.null(get(x))) {
      test <- testOnt(get(x))
      # Test the distribution of ONTOLOGY
      
      if ("BP" %in% names(test)) {
        get(x) %>% filter(ONTOLOGY == "BP") %>% dagPlot(
          .,
          firstSigNodes = 10,
          objectName = paste0("GO_", i),
          saveDir = paste(saveDir, i, sep = "/")
        )
      }
      if ("CC" %in% names(test)) {
        get(x) %>% filter(ONTOLOGY == "CC") %>% dagPlot(
          .,
          firstSigNodes = 10,
          objectName = paste0("GO_", i),
          saveDir = paste(saveDir, i, sep = "/")
        )
      }
      if ("MF" %in% names(test)) {
        get(x) %>% filter(ONTOLOGY == "MF") %>% dagPlot(
          .,
          firstSigNodes = 10,
          objectName = paste0("GO_", i),
          saveDir = paste(saveDir, i, sep = "/")
        )
      }
      rm(test)
    }
  }
  rm(i)
  rm(x)
  
  # ··· GO GSEA ---------------------------------------------------------------
  
  saveDir <- paste0(savePath, "DIY/", class, "/GO_GSEA/")
  
  x <- paste0("GO_gsea_diy_", class)
  assign(
    paste0("GO_gsea_diy_", class),
    enrichGSEA(
      symbolList,
      gene_up = symbolList_up,
      gene_down = symbolList_down,
      keyType = "Symbol",
      dbType = "GO",
      term2gene = GO_gene,
      term2name = GO_term,
      max = 2000,
      saveDir = saveDir,
      objectName = "gseaGO"
    ),
    envir = .GlobalEnv
  )
  
  dotPlot(get(x),
          term = 30,
          saveDir = saveDir,
          objectName = "gseaGO")
  
  cnetPlot(get(x),
           term = 15,
           saveDir = saveDir,
           objectName = "gseaGO")
  
  if (!is.null(get(x))) {
    test <- testOnt(get(x))
    if ("BP" %in% names(test)) {
      get(x) %>% filter(ONTOLOGY == "BP") %>% dagPlot(
        .,
        firstSigNodes = 10,
        objectName = "gseaGO",
        saveDir = saveDir
      )
    }
    if ("CC" %in% names(test)) {
      get(x) %>% filter(ONTOLOGY == "CC") %>% dagPlot(
        .,
        firstSigNodes = 10,
        objectName = "gseaGO",
        saveDir = saveDir
      )
    }
    if ("MF" %in% names(test)) {
      get(x) %>% filter(ONTOLOGY == "MF") %>% dagPlot(
        .,
        firstSigNodes = 10,
        objectName = "gseaGO",
        saveDir = saveDir
      )
    }
    rm(test)
  }
  
  ridgePlot(get(x),
            num = 30,
            saveDir = saveDir,
            objectName = "gseaGO")
  
  if (!is.null(get(x))) {
    if (length(which(get(x)@result$p.adjust < 0.05)) > 10) {
      get(x) %>% filter(p.adjust < 0.05) %>% ESPlot(., objectName = "gseaGO", saveDir = saveDir)
    } else{
      get(x) %>% slice(., 1:10) %>% ESPlot(., objectName = "gseaGO", saveDir = saveDir)
    }
  }
  rm(x)
}) %>% invisible()


# ·· KEGG ---------------------------------------------------------------------

lapply(names(listPathway_KEGG_gene), function(x) {
  # enrichment analysis is performed separately for each function
  class <- stringr::str_replace(x, "_gene$", "")
  KEGG_gene <- listPathway_KEGG_gene[[x]]
  KEGG_name <-
    listPathway_KEGG_name[[stringr::str_replace(x, "_gene$", "_name")]]
  
  # ··· KEGG pathway ORA ------------------------------------------------------
  
  saveDir <- paste0(savePath, "DIY/", class, "/KEGG_ORA")
  
  for (i in c("all", "up", "down")) {
    x <- paste0("KEGG_diy_", class, "_", i)
    assign(
      paste0("KEGG_diy_", class, "_", i),
      enrichORA(
        get(paste0("entrez_", i)),
        keyType = "ENTREZID",
        dbType = "KEGG",
        objectName = paste0("KEGG_", i),
        gene_up = symbol_up,
        gene_down = symbol_down,
        term2gene = KEGG_gene,
        term2name = KEGG_name,
        saveDir = paste(saveDir, i, sep = "/"),
      ),
      envir = .GlobalEnv
    )
    
    barPlot(
      get(x),
      term = 20,
      saveDir = paste(saveDir, i, sep = "/"),
      objectName = paste0("KEGG_", i)
    )
    
    if (i == "all") {
      barPlot2(
        get(x),
        term = 15,
        objectName = paste0("KEGG_", i),
        saveDir = paste(saveDir, i, sep = "/")
      )
    }
    
    dotPlot(
      get(x),
      term = 20,
      saveDir = paste(saveDir, i, sep = "/"),
      objectName = paste0("KEGG_", i)
    )
    
    cnetPlot(
      get(x),
      term = 15,
      saveDir = paste(saveDir, i, sep = "/"),
      objectName = paste0("KEGG_", i)
    )
    
    if (!is.null(get(x))) {
      diffGene_notNA <- switch(
        i,
        all = drop_na(diffGene, Symbol, `Entrez Gene ID`, logFC),
        up = drop_na(diffGene[diffGene$logFC > 0, ], Symbol, `Entrez Gene ID`, logFC),
        down = drop_na(diffGene[diffGene$logFC < 0, ], Symbol, `Entrez Gene ID`, logFC)
      )
      # extract data to add color for KEGG pathway plot
      
      pathPlot(
        get(x)[1:max(length(which(get(x)@result$p.adjust < 0.05)), 10), ],
        diffGene = diffGene_notNA,
        class = "ORA",
        saveDir = paste(saveDir, i, "PathPlot", sep = "/"),
        objectName = paste0("KEGG_", i)
      )
      # take p.adjust < 0.05 (or the first 10 if there are less than 10) to draw some pathway plots
      rm(diffGene_notNA)
    }
  }
  rm(i)
  rm(x)
  
  # ··· KEGG pathway GSEA -----------------------------------------------------
  
  saveDir <- paste0(savePath, "DIY/", class, "/KEGG_GSEA")
  
  x <- paste0("KEGG_gsea_diy_", class)
  assign(
    paste0("KEGG_gsea_diy_", class),
    enrichGSEA(
      entrezList,
      gene_up = symbolList_up,
      gene_down = symbolList_down,
      keyType = "ENTREZID",
      dbType = "KEGG",
      objectName = "gseaKEGG",
      term2gene = KEGG_gene,
      term2name = KEGG_name,
      saveDir = saveDir
    ),
    envir = .GlobalEnv
  )
  
  dotPlot(get(x),
          term = 20,
          saveDir = saveDir,
          objectName = "gseaKEGG")
  
  cnetPlot(get(x),
           term = 15,
           saveDir = saveDir,
           objectName = "gseaKEGG")
  
  ridgePlot(get(x),
            num = 20,
            saveDir = saveDir,
            objectName = "gseaKEGG")
  
  if (!is.null(get(x))) {
    if (length(which(get(x)@result$p.adjust < 0.05)) > 10) {
      get(x) %>% filter(p.adjust < 0.05) %>% ESPlot(., objectName = "gseaKEGG", saveDir = saveDir)
    } else{
      get(x) %>% slice(., 1:10) %>% ESPlot(., objectName = "gseaKEGG", saveDir = saveDir)
    }
  }
  
  if (!is.null(get(x))) {
    diffGene_notNA <-
      drop_na(allGene, Symbol, `Entrez Gene ID`, logFC)
    # extract data to add color for KEGG pathway plot
    
    pathPlot(
      get(x)[1:max(length(which(get(x)@result$p.adjust < 0.05)), 10),],
      diffGene = diffGene_notNA,
      class = "GSEA",
      saveDir = paste(saveDir, "PathPlot", sep = "/"),
      objectName = "gseaKEGG"
    )
    # take p.adjust < 0.05 (or the first 10 if there are less than 10) to draw some pathway plots
    rm(diffGene_notNA)
  }
  rm(x)
}) %>% invisible()

# save temp file --------------------------------------------------------------

save.image(file = paste0(savePath, sample,
                         "_enrich.RData"))
