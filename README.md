# NPS_manuscript
Code used for 2026 Neuropsychiatric Symptoms in AD paper
NPS Paper - Figure Documentation
Cell-type Specific Differential Expression in AD and MDD
Project Overview
This document provides documentation for reproducing the figures in the Neuropsychiatric Symptoms (NPS) paper, focusing on cell-type and subtype-specific differential expression analysis in individuals with and without LOAD and co-morbid MDD.
Core Pipeline Scripts
Script	Description
12_AD_DEP_cell_type_DE_analysis_adjusted.R	NEBULA differential expression for AD cohort (+/- MDD). Uses Mayo cell type annotations.
MDD_cell_type_DE_analysis.R	NEBULA differential expression for young cohort (MDD vs Control).
14_AD_DEP_cell_type_subset_FDR.R	Two-step FDR correction: filter |log2FC|>0.2, then apply FDR to AD DEP results.
14_Normal_DEP_cell_type_subset_FDR.R	Two-step FDR correction for older cognitively normal cohort (+/- MDD).
MDD_cell_type_subset_FDR.R	Two-step FDR correction for young MDD cohort.
15_AD_DEP_vs_Normal_DEP_intersection_mayo.R	Venn diagram intersection: AD+MDD vs Older Normal+MDD DEGs by cell type.
15_AD_DEP_vs_MDD_new_intersection_mayo.R	Venn diagram intersection: AD+MDD vs Young MDD DEGs by cell type.
AD_MDD_new_intersection_mayo.R	Additional intersection analysis between AD and young MDD cohorts.
16_cellchat_AD_DEP.R	CellChat cell-cell communication analysis for AD cohort (+/- MDD).
17_cellchat_Normal_DEP.R	CellChat analysis for older cognitively normal cohort (+/- MDD).
cellchat_MDD_Normal.R	CellChat analysis for young MDD cohort.
18_cellchat_AD_DEP_fix_heatmaps.R	Regenerates CellChat heatmaps with improved formatting for AD cohort.
18_cellchat_Normal_DEP_fix_heatmaps.R	Regenerates CellChat heatmaps for Normal cohort.
 
Slides 2-4: DEGs by Cell Subtype Bar Charts
Description
Bar charts showing the number of differentially expressed genes (DEGs) by cell subtype for three cohorts: (1) AD cohort +/- MDD, (2) Older cognitively normal +/- MDD, (3) Young cohort +/- MDD.
Scripts Required
Slide 2 - AD cohort:
•	12_AD_DEP_cell_type_DE_analysis_adjusted.R - Runs NEBULA DE analysis
•	14_AD_DEP_cell_type_subset_FDR.R - Applies two-step FDR correction
Slide 3 - Older Normal cohort:
•	Similar DE analysis script for Normal_DEPR cohort (in nebula_results_Normal_DEPR/)
•	14_Normal_DEP_cell_type_subset_FDR.R - Applies two-step FDR correction
Slide 4 - Young MDD cohort:
•	MDD_cell_type_DE_analysis.R - Runs NEBULA DE for young cohort
•	MDD_cell_type_subset_FDR.R - Applies two-step FDR correction
Key Parameters
•	FDR threshold: < 0.1
•	Log2FC threshold: |log2FC| > 0.2
•	Cell type annotations: Mayo cell types
Slides 5-6: Glutamatergic Common Genes (Venn Diagrams)
Description
Venn diagrams and tables showing common and distinct glutamatergic DEGs between cohorts. Slide 5: Older CN vs AD (both +/- MDD). Slide 6: Young vs AD (both +/- MDD).
Scripts Required
Slide 5:
•	15_AD_DEP_vs_Normal_DEP_intersection_mayo.R - Generates Venn diagrams and intersection tables
Slide 6:
•	15_AD_DEP_vs_MDD_new_intersection_mayo.R - Compares AD DEGs with young MDD DEGs
Slide 7: Pathway Analysis - Glutamatergic AD+MDD
Description
Metascape pathway analysis heatmaps for up-regulated and down-regulated DEGs in glutamatergic neurons from the AD cohort (+/- MDD comparison).
Scripts Required
•	14_AD_DEP_cell_type_subset_FDR.R - Generates DEG lists for Metascape input
•	External tool: Metascape (https://metascape.org)
Analysis Steps
1.	Extract up-regulated genes (log2FC > 0.2, FDR < 0.1) from glutamatergic results
2.	Extract down-regulated genes (log2FC < -0.2, FDR < 0.1) from glutamatergic results
3.	Submit gene lists to Metascape for pathway enrichment
4.	Download HeatmapSelectedGO.png for visualization
 
Slide 8: Microglia Common Genes - Young vs Older
Description
Table showing 13 common microglial DEGs between young (+/- MDD) and older cognitively normal (+/- MDD) cohorts.
Scripts Required
•	Similar intersection script comparing Young MDD vs Older Normal DEP for microglia
•	MDD_cell_type_DE_analysis.R + 14_Normal_DEP_cell_type_subset_FDR.R for input data
Slides 9-10: GABAergic Common Genes (Venn Diagrams)
Description
Venn diagrams for GABAergic neuron DEGs. Slide 9: Young vs AD. Slide 10: Older CN vs AD.
Scripts Required
•	15_AD_DEP_vs_MDD_new_intersection_mayo.R - For Young vs AD comparison (Slide 9)
•	15_AD_DEP_vs_Normal_DEP_intersection_mayo.R - For Older CN vs AD comparison (Slide 10)
Slide 11: Pathway Analysis - Microglia AD+MDD
Description
Metascape pathway analysis for down-regulated genes in microglia from AD cohort (MDD vs no MDD).
Slides 12-13: Cell-Cell Communication Networks
Description
Chord diagrams showing differential cell-cell communication between MDD and non-MDD groups. Slide 12: Number of interactions. Slide 13: Strength of interactions.
Scripts Required
•	16_cellchat_AD_DEP.R - CellChat analysis for AD cohort
•	17_cellchat_Normal_DEP.R - CellChat analysis for Normal cohort
•	18_cellchat_AD_DEP_fix_heatmaps.R - Regenerates visualizations with better formatting
•	18_cellchat_Normal_DEP_fix_heatmaps.R - Regenerates Normal cohort visualizations
Figure-to-Script Summary Table
Slide	Content	Primary Scripts
2	DEGs by cell subtype: AD +/- MDD	12_AD_DEP_cell_type_DE_analysis_adjusted.R, 14_AD_DEP_cell_type_subset_FDR.R
3	DEGs by cell subtype: Older CN +/- MDD	DE analysis + 14_Normal_DEP_cell_type_subset_FDR.R
4	DEGs by cell subtype: Young +/- MDD	MDD_cell_type_DE_analysis.R, MDD_cell_type_subset_FDR.R
5	Glutamatergic: Older CN vs AD	15_AD_DEP_vs_Normal_DEP_intersection_mayo.R
6	Glutamatergic: Young vs AD	15_AD_DEP_vs_MDD_new_intersection_mayo.R
7	Pathway: Glutamatergic AD+MDD	14_AD_DEP_cell_type_subset_FDR.R + Metascape
8	Microglia: Young vs Older	Intersection scripts for microglia cell type
9	GABAergic: Young vs AD	15_AD_DEP_vs_MDD_new_intersection_mayo.R
10	GABAergic: Older CN vs AD	15_AD_DEP_vs_Normal_DEP_intersection_mayo.R
11	Pathway: Microglia AD+MDD	14_AD_DEP_cell_type_subset_FDR.R + Metascape
12	CellChat: Number of interactions	16_cellchat_AD_DEP.R, 17_cellchat_Normal_DEP.R, 18_*_fix_heatmaps.R
13	CellChat: Strength of interactions	16_cellchat_AD_DEP.R, 17_cellchat_Normal_DEP.R, 18_*_fix_heatmaps.R
Software Dependencies
R Packages
•	Seurat, Matrix, data.table, ggplot2, dplyr, plyr, stringr
•	nebula - For differential expression analysis
•	CellChat - For cell-cell communication analysis
•	ggvenn - For Venn diagram visualization
•	missForest - For missing data imputation
•	NMF, ggalluvial - For CellChat visualizations

