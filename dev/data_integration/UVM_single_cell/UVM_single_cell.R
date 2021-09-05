 
library(tidyverse)
library(glue)
library(Seurat)
library(tidyseurat)
# library(future)
# plan(multisession, workers = 2)


cell_type_df = 
	
	bind_rows(
		
		# One cell type one figure
		read_csv("dev/UVM_single_cell/scUM_NComms_Fig1B_metadata_MD.csv") %>%
			rename(cell = X1) %>%
			extract(cell, "cell", "(.+)-.+" ),
		
		# One cell type one figure
		read_csv("dev/UVM_single_cell/scUM_NComms_Fig4A_Immune.Subset_metadata_MD.csv") %>%
			rename(cell = X1) %>%
			extract(cell, "cell", "(.+)-.+" )
	) %>%
	rename(sample = orig.ident, cell_type = CellType) %>% 
	group_by(cell, sample) %>%
	slice(2) %>%
	ungroup() %>%
	
	# Rename cells
	mutate(cell = glue("{cell}-1")) %>%
	
	nest(cell_type = -sample)

clinical_information = read_csv("dev/UVM_single_cell/clinical_information.csv")

# raw_counts_nested  =
# 	dir("dev/UVM_single_cell/", pattern = ".gz", full.names = TRUE) %>%
# 	enframe(value = "path") %>%
# 	mutate(file_name = basename(path)) %>%
# 	extract( file_name, "sample", ".*_([A-Z0-9]+)_.*") %>%
# 	nest(data = -sample) %>%
# 	mutate(seurat = map(
# 		data,
# 		~ ReadMtx(mtx = .x$path[[3]], cells = .x$path[[1]], features = .x$path[[2]], feature.column = 1) %>%
# 			CreateSeuratObject()
# 	))
# 
# saveRDS(raw_counts_nested , "dev/UVM_single_cell/raw_counts_nested.rds")

raw_counts_nested = readRDS( "dev/UVM_single_cell/raw_counts_nested.rds")

counts = 
	raw_counts_nested %>%
	left_join(clinical_information) %>%
	left_join(cell_type_df) %>% 
	mutate(seurat = map2(
		seurat, cell_type, 
		~ left_join(.x, .y) %>% 
			filter(cell_type %>% is.na %>% `!`) 
	)) %>%
	select(-data, -cell_type) %>% 
	unnest_seurat(seurat)

counts %>% saveRDS("dev/UVM_single_cell/counts.rds")

library(sccomp)

sccomp_estimation = 
	counts %>%
	filter(`Sample Type` == "Primary") %>%
	sccomp_glm(  ~ Sex, sample, cell_type  )


data_for_plot = 
	sccomp_estimation %>% 
	tidyr::unnest(outliers) %>%
	group_by(sample) %>%
	mutate(proportion = (count+1)/sum(count+1)) %>%
	ungroup(sample) %>%
	left_join(counts %>% distinct(sample, Sex))

ggplot() +
	geom_boxplot(
		aes(Sex, proportion, fill=significant),
		outlier.shape = NA, 
		data = data_for_plot %>% filter(!outlier)
	) + 
	geom_jitter(aes(Sex, proportion, color=outlier), size = 1, data = data_for_plot) + 
	facet_wrap(~ interaction(cell_type), scale="free_y") +
	scale_y_continuous(trans="logit") +
	scale_color_manual(values = c("black", "#e11f28")) +
	scale_fill_manual(values = c("white", "#E2D379")) +
	xlab("Biological condition") + 
	ylab("Cell-group proportion") + 
	theme_bw() +cp-
	theme(strip.background =element_rect(fill="white"))

sccomp_estimation %>% 
	unnest(concentration) %>% 
	ggplot(aes(`.median_(Intercept)`, mean)) + 
	geom_errorbar(aes(ymin = `2.5%`, ymax=`97.5%`), color="#4DAF4A", alpha = 0.4) +
	geom_errorbar(aes(xmin = `.lower_(Intercept)`, xmax=`.upper_(Intercept)`), color="#4DAF4A", alpha = 0.4) +
	geom_point() +
	geom_abline(intercept = 5.7496330, slope = -0.9650953, linetype = "dashed", color="grey") +
	xlab("Category logit-proportion mean") +
	ylab("Category log-concentration") +
	theme_bw() 
