dataDir <- commandArgs(trailingOnly = TRUE)

sample <- basename(dataDir)
# the prefix for all output file names

savePath <- paste0("result/", dataDir, "/Mapping/")
# the output folder

library(magrittr)
library(stringr)
library(ggplot2)

colData <-
  read.table(
    file.path(dataDir, "colData.tsv"),
    sep = "\t",
    header = T,
    row.names = 1
  )
# read the group information file to sort the samples in the plot

list <-
  list.files(
    path = paste0("result/", dataDir, "/Mapping/"),
    pattern = ".*\\.log$",
    full.names = TRUE
  ) %>% magrittr::set_names(., stringr::str_replace(basename(.), "\\.log$", "")) %>% lapply(., readLines)
# read the the statistic files to a list

list2 <- lapply(list, function(x) {
  table <-
    x[c(3:6, 8:10)]  %>% str_trim(string = .) %>% str_replace_all(c(":.+\\(" = "\t", "%\\)" = "")) %>% read.table(
      text = .,
      sep = "\t",
      col.names = c("class", "proportion")
    )
  table[5:7, 1] %<>% str_c("Unpaired reads ", .) %>% str_replace("A", tolower)
  table[, 2] %<>% `/`(., 100)
  table[5:7, 2] %<>% `*`(table[1, 2], .)
  table <- table[-1,]
})
# extract the ratio information from text to table

df <- do.call(rbind, list2)
# merge all tables

df$sample <-
  unlist(lapply(seq_along(list2), function(i)
    rep(names(list2)[i], nrow(list2[[i]]))))
df$sample %<>%  factor(levels = row.names(colData))
# add a column "sample" to represent the groups

df$class %<>% factor(
  levels = c(
    "Aligned concordantly 1 time",
    "Aligned concordantly >1 times",
    "Aligned discordantly 1 time",
    "Unpaired reads aligned 1 time",
    "Unpaired reads aligned >1 times",
    "Unpaired reads aligned 0 time"
  )
)
# the order of class

p <-
  ggplot(df, aes(x = sample, y = proportion, fill = class)) + geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = rev(RColorBrewer::brewer.pal(6,"RdYlGn")))

ggsave(
  filename = paste0(sample, "_mappingProportion.png"),
  path = savePath,
  plot = p,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)

df$class %<>% factor(
  levels = c(
    "Unpaired reads aligned 0 time",
    "Unpaired reads aligned >1 times",
    "Unpaired reads aligned 1 time",
    "Aligned discordantly 1 time",
    "Aligned concordantly >1 times",
    "Aligned concordantly 1 time"
  )
)
# change the order of class

p <-
  ggplot(df, aes(x = sample, y = proportion, fill = class)) + geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(6,"RdYlGn"))

ggsave(
  filename = paste0(sample, "_mappingProportion2.png"),
  path = savePath,
  plot = p,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)
