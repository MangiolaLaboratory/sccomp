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


options(future.globals.maxSize= 100000*1024^2)


job({

  oligo_adjusted_cell_type <- readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/oligo_breast/expanded_analyses_with_control/cancer_only_analyses/integrated_counts.rds")


  oligo_adjusted_cell_type %>%
    tidyseurat::mutate(cell_type = curated_cell_type) %>%
    saveRDS("dev/data_integration/UMAP_oligo.rds")
})


job({

  benign_adjusted_cell_type <- readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/oligo_breast/expanded_analyses_with_control/benign_adjusted_cell_type.rds")

  DefaultAssay(benign_adjusted_cell_type) = "integrated"

  benign_adjusted_cell_type %>%
    tidyseurat::mutate(cell_type = curated_cell_type_pretty) %>%
    tidyseurat::filter(!is.na(cell_type)) %>%
    RunPCA() %>%
    RunUMAP(dims=1:20) %>%
    saveRDS("dev/data_integration/UMAP_GSE115189_SCP345_SCP424_SCP591_SRR11038995_SRR7244582_10x6K_10x8K.rds")
})

job({
  options(future.globals.maxSize= 100000*1024^2)

    readRDS("dev/data_integration/UVM_single_cell/counts.rds")  |>
    tidyseurat::rename(type = `Sample Type`) %>%
      tidyseurat::filter(!is.na(cell_type)) %>%
      drop_dead() %>%
    adjust_abundance(
      ~integrate(sample) + mito_RPS + high_RPS,
      reference_samples =
        filter(., nCount_RNA == max(nCount_RNA)) %>%
        pull(sample) %>%
        unique()
    ) %>%
    RunPCA() %>%
    RunUMAP(dims=1:20) %>%
    saveRDS("dev/data_integration/UMAP_GSE139829_uveal_melanoma.rds")

})

job({
  options(future.globals.maxSize= 100000*1024^2)

    readRDS("dev/data_integration/SCP1288_renal_cell_carcinoma.rds")  |>
    tidyseurat::filter(!is.na(sample) & !is.na(cell_type) & !is.na(sex))  |>
      tidyseurat::filter(!is.na(cell_type)) %>%
      drop_dead() %>%
    adjust_abundance(
      ~integrate(sample) ,
      reference_samples =
        filter(., nFeature_originalexp == max(nFeature_originalexp)) %>%
        pull(sample) %>%
        unique()
    ) %>%
    RunPCA() %>%
    RunUMAP(dims=1:20) %>%
    saveRDS("dev/data_integration/UMAP_SCP1288_renal_cell_carcinoma.rds")
})

job({
  options(future.globals.maxSize= 100000*1024^2)

    readRDS("dev/data_integration/SCP1039_bc_cells.rds")  |>
    tidyseurat::select(-nCount_RNA, -nFeature_RNA) %>%
    tidyseurat::mutate(type = subtype=="TNBC") %>%
      tidyseurat::filter(!is.na(cell_type)) %>%
      drop_dead() %>%
    adjust_abundance(
      ~integrate(sample) + mito_RPS + high_RPS,
      reference_samples =
        filter(., nCount_originalexp == max(nCount_originalexp)) %>%
        pull(sample) %>%
        unique()
    ) %>%
    RunPCA() %>%
    RunUMAP(dims=1:20) %>%
    saveRDS("dev/data_integration/UMAP_SCP1039_bc_cells.rds")
})

job({
  options(future.globals.maxSize= 100000*1024^2)

    readRDS("dev/data_integration/s41587-020-0602-4_COVID_19.rds")  |>
      tidyseurat::mutate(is_critical = severity=="critical") %>%
      tidyseurat::filter(!is.na(cell_type)) %>%
    adjust_abundance(
      ~integrate(sample) ,
      reference_samples =
        filter(., nCount_RNA == max(nCount_RNA)) %>%
        pull(sample) %>%
        unique(),
      assay = "RNA"
    ) %>%
    RunPCA() %>%
    RunUMAP(dims=1:20) %>%
    saveRDS("dev/data_integration/UMAP_s41587-020-0602-4_COVID_19.rds")
})

job({
  options(future.globals.maxSize= 100000*1024^2)

    readRDS("dev/data_integration/GSE120575_melanoma.rds")  |>
      tidyseurat::filter(!is.na(cell_type)) %>%
      drop_dead() %>%
    adjust_abundance(
      ~integrate(sample) + mito_RPS + high_RPS,
      reference_samples =
        filter(., nCount_RNA == max(nCount_RNA)) %>%
        pull(sample) %>%
        unique()
    ) %>%
    RunPCA() %>%
    RunUMAP(dims=1:20) %>%
    saveRDS("dev/data_integration/UMAP_GSE120575_melanoma.rds")
})




