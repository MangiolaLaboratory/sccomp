# library(tidyverse)

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

library(Seurat)
library(tidyseurat)
library(tidysc)
library(sccomp)
library(job)
library(patchwork)
library(future)

counts = readRDS("dev/data_integration/UVM_single_cell/counts.rds")

options(future.globals.maxSize= 100000*1024^2)
# rownames(counts) =  metadata = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
#                                                      keys = rownames(counts) ,
#                                                      keytype = "ENSEMBL",
#                                                      column = "SYMBOL",
#                                                      multiVals = "first"
# )
#
# counts = counts[!is.na(rownames(counts)),]

job({

  options(future.globals.maxSize= 100000*1024^2)
  plan(multisession, workers=8)

  counts_UVM =
    readRDS("dev/data_integration/UVM_single_cell/counts.rds") %>%
    drop_empty_and_dead() %>%
    adjust_abundance(
      ~integrate(sample) + mito_RPS + high_RPS,
      reference_samples =
        filter(., nCount_RNA == max(nCount_RNA)) %>%
        pull(sample) %>%
        unique()
    ) %>%
    RunPCA() %>%
    RunUMAP(dims=1:20)

})



%>%
  DimPlot(group.by = "cell_type")

