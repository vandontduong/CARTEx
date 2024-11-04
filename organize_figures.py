from pptx import Presentation
from pptx.util import Inches

experiment = "GSE136874"
description = "comparison of CD19 and GD2 CAR T cells on day 10 in culture"

prs = Presentation()
# set width and height to 16 and 9 inches.
prs.slide_width = Inches(10.5)
prs.slide_height = Inches(14)

blank_slide_layout = prs.slide_layouts[6]

def add_slide(prs, layout, subtitle):
	slide = prs.slides.add_slide(layout)
	txBox = slide.shapes.add_textbox(left = Inches(0), top = Inches(0), width = Inches(1), height = Inches(1))
	tf = txBox.text_frame
	tf.text = "CARTEx project"
	p = tf.add_paragraph()
	p.text = subtitle
	return slide





# new slide
slide = add_slide(prs, blank_slide_layout, "Figure 1. Bulk RNA sequencing of functional and exhausted CAR T cells across 21 days")

pic = slide.shapes.add_picture("./construction/plots/" + "plot_initial_pca.png", top = Inches(0.75), left = Inches(3.5), height = Inches(2))

pic = slide.shapes.add_picture("./construction/plots/" + "plot_canonical_exhaustion_markers.png", top = Inches(0.75), left = Inches(6.5), height = Inches(2))

# manually replace heatmap with gene expression; use placeholder for now

pic = slide.shapes.add_picture("./construction/plots/" + "cluster1.png", left = Inches(3.5), top = Inches(3), height = Inches(2))

pic = slide.shapes.add_picture("./construction/plots/" + "cluster2.png", left = Inches(3.5), top = Inches(4.65), height = Inches(2))

pic = slide.shapes.add_picture("./construction/plots/" + "cluster3.png", left = Inches(3.5), top = Inches(6.3), height = Inches(2))

pic = slide.shapes.add_picture("./construction/plots/" + "cluster4.png", left = Inches(3.5), top = Inches(7.95), height = Inches(2))

pic = slide.shapes.add_picture("./construction/plots/" + "cluster5.png", left = Inches(3.5), top = Inches(9.6), height = Inches(2))

pic = slide.shapes.add_picture("./construction/plots/" + "plot_weights.png", left = Inches(7), top = Inches(3), height = Inches(2.5))

pic = slide.shapes.add_picture("./construction/plots/" + "plt_CARTEx_200.png", left = Inches(6), top = Inches(6), height = Inches(3))




# new slide
slide = add_slide(prs, blank_slide_layout, "Figure 2. Biological context / Evaluation with other exhaustion and cell state signatures")



pic = slide.shapes.add_picture("./construction/plots/" + "plot_mean_zscore_by_signature.png", top = Inches(5), left = Inches(0.25), height = Inches(2))



# new slide
slide = add_slide(prs, blank_slide_layout, "Supplementary figures. ")

# Fraction of signature genes represented in experiments
pic = slide.shapes.add_picture("./construction/plots/" + "all_genes_repmap.png", top = Inches(0.75), left = Inches(0.25), height = Inches(2))

# Fraction of CARTEx (200) genes represented in each signature
pic = slide.shapes.add_picture("./construction/plots/" + "fraction_CARTEx_repchart.png", top = Inches(0.75), left = Inches(5), height = Inches(1.25))

# Fraction of signature genes represented in CARTEx (200)

pic = slide.shapes.add_picture("./construction/plots/" + "fraction_sig_repchart_CARTEx.png", top = Inches(2.25), left = Inches(5), height = Inches(1.25))





# new slide
slide = add_slide(prs, blank_slide_layout, "Figure 4. T cell exhaustion is distinct from aging, cell cycle, and differentiation")
# Naive and terminal effector T cells have age-invariant CARTEx scores
# Central memory and effector memory T cells have age-dependent CARTEx scores
# We derived relative controls for scRNAseq analyses

experiment = "GSE136184" # Cross-sectional aging experiment

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_cs_prepare_umap_age_group_2.png", top = Inches(0.75), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_cs_prepare_umap_CARTEx_200_group.png", top = Inches(0.75), left = Inches(2.5), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_cs_prepare_umap_phase.png", top = Inches(2.5), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_cs_prepare_barplot_phase_age_group_slim.png",  top = Inches(2.5), left = Inches(2.5), height = Inches(1.25))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_cs_prepare_umap_CARTEx_200_phase.png", top = Inches(2.5), left = Inches(4), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_cs_prepare_umap_predicted_monaco.png", top = Inches(4.25), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_cs_prepare_barplot_monaco_age_group_slim.png", top = Inches(4.25), left = Inches(2.5), height = Inches(1.25))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_cs_prepare_umap_CARTEx_200_monaco.png", top = Inches(4.25), left = Inches(4), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_cs_prepare_vlnplot_age_group_exhaustion_markers.png", top = Inches(6), left = Inches(0.25), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_aggplot_CARTEx_200_group_monaco_split_countsized.png", top = Inches(6), left = Inches(4), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_cs_extract_umap_highlight_combined.png", top = Inches(8.25), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_featplot_CARTEx_200_extract_ident.png", top = Inches(8.25), left = Inches(3.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_vlnplot_CARTEx_200_extract_ident.png", top = Inches(8.25), left = Inches(5.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_transition_volcano_OTvYN_CARTEx_630.png", top = Inches(8.25), left = Inches(6.75), height = Inches(2))



# new slide
slide = add_slide(prs, blank_slide_layout, "Supplementary Figures for aging dataset")

experiment = "GSE136184" # Cross-sectional aging experiment

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_cs_prepare_umap_CARTEx_200_group.png", top = Inches(0.75), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_cs_prepare_umap_CARTEx_200_phase.png", top = Inches(2.5), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_cs_prepare_umap_CARTEx_200_monaco.png", top = Inches(4.25), left = Inches(0.25), height = Inches(1.5))



pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_l_prepare_umap_age_group_2.png", top = Inches(6), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_l_prepare_umap_CARTEx_200_group.png", top = Inches(6), left = Inches(2.5), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_l_prepare_umap_phase.png", top = Inches(7.75), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_l_prepare_barplot_phase_age_group_slim.png",  top = Inches(7.75), left = Inches(2.5), height = Inches(1.25))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_l_prepare_umap_CARTEx_200_phase.png", top = Inches(7.75), left = Inches(4), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_l_prepare_umap_predicted_monaco.png", top = Inches(9.5), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_l_prepare_barplot_monaco_age_group_slim.png", top = Inches(9.5), left = Inches(2.5), height = Inches(1.25))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_l_prepare_umap_CARTEx_200_monaco.png", top = Inches(9.5), left = Inches(4), height = Inches(1.5))




# new slide
slide = add_slide(prs, blank_slide_layout, "Supplementary Figures for aging dataset")

experiment = "GSE136184" # Cross-sectional aging experiment

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_cs_prepare_umap_CARTEx_200.png", top = Inches(0.75), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_cs_prepare_umap_LCMV_Tex.png", top = Inches(0.75), left = Inches(2), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_cs_prepare_umap_NKlike_Tex.png", top = Inches(0.75), left = Inches(3.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_cs_prepare_umap_BBD_Tex.png", top = Inches(0.75), left = Inches(5.5), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_cs_prepare_umap_PD1_Tex.png", top = Inches(0.75), left = Inches(7.25), height = Inches(1.5))




# new slide
slide = add_slide(prs, blank_slide_layout, "Figure 5. Validation in exhaustion-prone models")

experiment = "GSE136874" # CD19 vs GD2 experiment

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CAR.png", left = Inches(0), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_phase.png", left = Inches(2.5), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_phase_CAR_slim.png", left = Inches(5), top = Inches(0.75), height = Inches(1.25))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_predicted_monaco.png", left = Inches(6.5), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_monaco_CAR_slim.png", left = Inches(9), top = Inches(0.75), height = Inches(1.25))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_vlnplot_CAR_exhaustion_markers.png", left = Inches(0), top = Inches(3), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_transition_volcano_GD2vCD19.png", left = Inches(2.75), top = Inches(3), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_featplot_CARTEx_200_CAR.png", left = Inches(5.25), top = Inches(3), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200.png", left = Inches(6.25), top = Inches(3), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_CARTEx_200.png", left = Inches(8.5), top = Inches(3), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_aggplot_CARTEx_200_CAR_monaco_split_countsized.png", left = Inches(0.25), top = Inches(5.25), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_dmap_CAR.png", left = Inches(2.5), top = Inches(5.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_dmap_CARTEx_200.png", left = Inches(4.75), top = Inches(5.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_transition_transitplot_CAR_CARTEx_200.png", left = Inches(7), top = Inches(5.25), height = Inches(2))

experiment = "GSE160160" # Continuous antigen exposure experiment

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_exposure.png", left = Inches(0), top = Inches(7.5), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_phase.png", left = Inches(2.5), top = Inches(7.5), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_phase_exposure_slim.png", left = Inches(5), top = Inches(7.5), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_predicted_monaco.png", left = Inches(6.5), top = Inches(7.5), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_monaco_exposure_slim.png", left = Inches(9), top = Inches(7.5), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_vlnplot_exposure_exhaustion_markers.png", left = Inches(0), top = Inches(9.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_explore_volcano_CAEvday0.png", left = Inches(2.75), top = Inches(9.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_featplot_CARTEx_200_exposure.png", left = Inches(5.25), top = Inches(9.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200.png", left = Inches(6.25), top = Inches(9.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_CARTEx_200.png", left = Inches(8.5), top = Inches(9.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_aggplot_CARTEx_200_exposure_monaco_split_countsized.png", left = Inches(0.25), top = Inches(12), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_dmap_exposure.png", left = Inches(2.5), top = Inches(12), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_dmap_CARTEx_200.png", left = Inches(4.75), top = Inches(12), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_transition_transitplot_exposure_CARTEx_200.png", left = Inches(7), top = Inches(12), height = Inches(2))


# new slide
slide = add_slide(prs, blank_slide_layout, "Supplementary Figures for validation in exhaustion-prone models")

experiment = "GSE136874" # CD19 vs GD2 experiment

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_transition_volcano_GD2vCD19_C2genes.png", top = Inches(0.75), left = Inches(0.25), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_CAR.png", top = Inches(3), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_phase.png", top = Inches(4.75), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_monaco.png", top = Inches(6.5), left = Inches(0.25), height = Inches(1.5))


experiment = "GSE160160" # Continuous antigen exposure experiment

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_explore_volcano_CAEvday0_CARTEx_C2genes.png", top = Inches(0.75), left = Inches(3.25), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_exposure.png", top = Inches(8.25), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_phase.png", top = Inches(10), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_monaco.png", top = Inches(11.75), left = Inches(0.25), height = Inches(1.5))




# new slide
slide = add_slide(prs, blank_slide_layout, "Supplementary Figures for validation in exhaustion-prone models")

experiment = "GSE136874" # CD19 vs GD2 experiment

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_LCMV_Tex.png", left = Inches(0), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_NKlike_Tex.png", left = Inches(2.5), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_BBD_Tex.png", left = Inches(5), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_PD1_Tex.png", left = Inches(7.5), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_scatter_CARTEx_200_LCMV.png", left = Inches(0), top = Inches(2.5), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_scatter_CARTEx_200_NKlike.png", left = Inches(2.5), top = Inches(2.5), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_scatter_CARTEx_200_BBD.png", left = Inches(5), top = Inches(2.5), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_scatter_CARTEx_200_PD1.png", left = Inches(7.5), top = Inches(2.5), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_LCMV_Tex.png", left = Inches(0), top = Inches(4.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_NKlike_Tex.png", left = Inches(2.5), top = Inches(4.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_BBD_Tex.png", left = Inches(5), top = Inches(4.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_PD1_Tex.png", left = Inches(7.5), top = Inches(4.75), height = Inches(2))


experiment = "GSE160160" # Continuous antigen exposure experiment

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_LCMV_Tex.png", left = Inches(0), top = Inches(7), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_NKlike_Tex.png", left = Inches(2.5), top = Inches(7), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_BBD_Tex.png", left = Inches(5), top = Inches(7), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_PD1_Tex.png", left = Inches(7.5), top = Inches(7), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_scatter_CARTEx_200_LCMV.png", left = Inches(0), top = Inches(8.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_scatter_CARTEx_200_NKlike.png", left = Inches(2.5), top = Inches(8.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_scatter_CARTEx_200_BBD.png", left = Inches(5), top = Inches(8.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_scatter_CARTEx_200_PD1.png", left = Inches(7.5), top = Inches(8.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_LCMV_Tex.png", left = Inches(0), top = Inches(11), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_NKlike_Tex.png", left = Inches(2.5), top = Inches(11), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_BBD_Tex.png", left = Inches(5), top = Inches(11), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_PD1_Tex.png", left = Inches(7.5), top = Inches(11), height = Inches(2))






# new slide
slide = add_slide(prs, blank_slide_layout, "Supplementary figures. Activation and resting T cells from healthy donors")
# important to assess baseline exhaustion

experiment = "GSE126030" # 

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_tissue.png", top = Inches(0.75), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_tissue.png", top = Inches(0.75), left = Inches(2.5), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_stimulation_status.png", top = Inches(2.5), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_stimulation_status.png", top = Inches(2.5), left = Inches(2.5), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_phase.png", top = Inches(4.25), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_phase_stimulation_status_slim.png", top = Inches(4.25), left = Inches(2.5), height = Inches(1.25))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_phase.png", top = Inches(4.25), left = Inches(4), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_predicted_monaco.png", left = Inches(0.25), top = Inches(6), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_monaco_stimulation_status_slim.png", left = Inches(2.5), top = Inches(6), height = Inches(1.25))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_monaco.png", top = Inches(6), left = Inches(4), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_featplot_CARTEx_200_stim.png", top = Inches(7.75), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_explore_vlnplot_stim_exhaustion_markers.png", top = Inches(7.75), left = Inches(2), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_aggplot_CARTEx_200_stimulation_status_monaco_split_countsized.png", top = Inches(7.75), left = Inches(4), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_CARTEx_200.png", top = Inches(7.75), left = Inches(7), height = Inches(2))


# new slide
slide = add_slide(prs, blank_slide_layout, "Supplementary figures. T cells from healthy donors and cancer patients (MD Anderson pan-cancer)")
# H, normal tissues from healthy donors; U, tumor-adjacent uninvolved tissues; P, primary tumor tissues; M, metastatic tumor tissues

experiment = "mdandersonTCM" # 

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_tissuetype.png", top = Inches(0.75), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_tissue_type.png", top = Inches(0.75), left = Inches(2.5), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_phase.png", top = Inches(2.5), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_phase_tissue_type_slim.png", top = Inches(2.5), left = Inches(2.5), height = Inches(1.25))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_phase.png", top = Inches(2.5), left = Inches(4), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_predicted_monaco.png", left = Inches(0.25), top = Inches(4.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_monaco_tissue_type_slim.png", left = Inches(2.5), top = Inches(4.25), height = Inches(1.25))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_monaco.png", top = Inches(4.25), left = Inches(4), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_featplot_CARTEx_200_tissue_type.png", top = Inches(6), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_vlnplot_tissue_type_exhaustion_markers.png", top = Inches(6), left = Inches(3), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_aggplot_CARTEx_200_TissueType_monaco_split_countsized.png", top = Inches(6), left = Inches(6), height = Inches(2))











# new slide
slide = add_slide(prs, blank_slide_layout, "Figure 6. CARTEx applied to clinical oncology")

experiment = "GSE120575" # Sade-Feldman

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_timepoint.png", left = Inches(0), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_response.png", left = Inches(2), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_phase.png", left = Inches(4.25), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_phase_response_slim.png", left = Inches(6.25), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_predicted_monaco.png", left = Inches(7.25), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_monaco_response_slim.png", left = Inches(9.25), top = Inches(0.75), height = Inches(1.5))

# pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_vlnplot_response_exhaustion_markers.png", left = Inches(0.25), top = Inches(3), height = Inches(2))

# pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_timepoint.png", left = Inches(0.25), top = Inches(3), height = Inches(2))

# pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_featplot_CARTEx_200_baseline_response.png", left = Inches(0.25), top = Inches(3), height = Inches(2))

# pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_baseline.png", left = Inches(1.25), top = Inches(3), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_featplot_CARTEx_200_response.png", left = Inches(0.25), top = Inches(3), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200.png", left = Inches(1.5), top = Inches(3), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_plot_volcano_baseline_response.png", left = Inches(3.5), top = Inches(3), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_CARTEx_200_baseline.png", left = Inches(6), top = Inches(3), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_aggplot_CARTEx_200_response_monaco_split_baseline_response_countsized.png", left = Inches(8.25), top = Inches(3), height = Inches(2))


experiment = "GSE151511" # CAR T cell infusion product

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_responder.png", left = Inches(0.25), top = Inches(5.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_phase.png", left = Inches(2.5), top = Inches(5.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_phase_responder_slim_noNA.png", left = Inches(4.75), top = Inches(5.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_predicted_monaco.png", left = Inches(5.5), top = Inches(5.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_monaco_responder_slim_noNA.png", left = Inches(7.75), top = Inches(5.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_featplot_CARTEx_200_responder.png", left = Inches(0.25), top = Inches(7.5), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200.png", left = Inches(1.5), top = Inches(7.5), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_plot_volcano_baseline_response.png", left = Inches(3.5), top = Inches(7.5), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_CARTEx_200.png", left = Inches(6), top = Inches(7.5), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_aggplot_CARTEx_200_response_monaco_split_countsized.png", left = Inches(8.25), top = Inches(7.5), height = Inches(2))


experiment = "GSE125881" # CAR T cell kinetics over 16 weeks

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_group2.png", left = Inches(0.25), top = Inches(9.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_phase.png", left = Inches(2.5), top = Inches(9.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_phase_group2_slim.png", left = Inches(4.75), top = Inches(9.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_predicted_monaco.png", left = Inches(5.5), top = Inches(9.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_monaco_group2_slim.png", left = Inches(7.75), top = Inches(9.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_featplot_CARTEx_200_group2.png", left = Inches(0.25), top = Inches(12), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200.png", left = Inches(1.5), top = Inches(12), height = Inches(1.5))

# pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_explore_plot_volcano_InVivo_IP.png", left = Inches(3.5), top = Inches(12), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_CARTEx_200.png", left = Inches(3.5), top = Inches(12), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_aggplot_CARTEx_200_group2_monaco_split_countsized.png", left = Inches(6.5), top = Inches(12), height = Inches(2))





# new slide
slide = add_slide(prs, blank_slide_layout, "Supplementary Figure for CARTEx applied to clinical oncology")

experiment = "GSE120575" # Sade-Feldman

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_vlnplot_response_exhaustion_markers.png", top = Inches(0.75), left = Inches(0.25), height = Inches(2))

experiment = "GSE151511" # CAR T cell infusion product

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_vlnplot_response_exhaustion_markers.png", top = Inches(0.75), left = Inches(3.25), height = Inches(2))

experiment = "GSE125881" # CAR T cell kinetics over 16 weeks


pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_group.png", top = Inches(3), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_phase.png", top = Inches(4.75), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_monaco.png", top = Inches(6.5), left = Inches(0.25), height = Inches(1.5))





# new slide
slide = add_slide(prs, blank_slide_layout, "Supplementary Figure for CARTEx applied to clinical oncology")

experiment = "GSE120575" # Sade-Feldman

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_LCMV_Tex.png", left = Inches(0), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_NKlike_Tex.png", left = Inches(2.5), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_BBD_Tex.png", left = Inches(5), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_PD1_Tex.png", left = Inches(7.5), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_scatter_CARTEx_200_LCMV_baseline_response.png", left = Inches(0), top = Inches(2.5), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_scatter_CARTEx_200_NKlike_baseline_response.png", left = Inches(2.5), top = Inches(2.5), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_scatter_CARTEx_200_BBD_baseline_response.png", left = Inches(5), top = Inches(2.5), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_scatter_CARTEx_200_PD1_baseline_response.png", left = Inches(7.5), top = Inches(2.5), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_LCMV_Tex_baseline.png", left = Inches(0), top = Inches(4.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_NKlike_Tex_baseline.png", left = Inches(2.5), top = Inches(4.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_BBD_Tex_baseline.png", left = Inches(5), top = Inches(4.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_PD1_Tex_baseline.png", left = Inches(7.5), top = Inches(4.75), height = Inches(2))


experiment = "GSE151511" # CAR T cell infusion product

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_LCMV_Tex.png", left = Inches(0), top = Inches(7), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_NKlike_Tex.png", left = Inches(2.5), top = Inches(7), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_BBD_Tex.png", left = Inches(5), top = Inches(7), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_PD1_Tex.png", left = Inches(7.5), top = Inches(7), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_scatter_CARTEx_200_LCMV.png", left = Inches(0), top = Inches(8.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_scatter_CARTEx_200_NKlike.png", left = Inches(2.5), top = Inches(8.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_scatter_CARTEx_200_BBD.png", left = Inches(5), top = Inches(8.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_scatter_CARTEx_200_PD1.png", left = Inches(7.5), top = Inches(8.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_LCMV_Tex.png", left = Inches(0), top = Inches(11), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_NKlike_Tex.png", left = Inches(2.5), top = Inches(11), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_BBD_Tex.png", left = Inches(5), top = Inches(11), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_PD1_Tex.png", left = Inches(7.5), top = Inches(11), height = Inches(2))











# new slide
slide = add_slide(prs, blank_slide_layout, "Supplementary Figure. COVID-19 severity and SARS-CoV-2-reactive T cell exhaustion")

experiment = "GSE153931" # COVID-19 severity and SARS-CoV-2-reactive T cell exhaustion

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_orig_virus.png", top = Inches(0.75), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_orig_virus.png", top = Inches(0.75), left = Inches(2.5), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_severity_mod.png", top = Inches(2.5), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_severity_mod.png", top = Inches(2.5), left = Inches(2.5), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_phase.png", top = Inches(4.25), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_phase_severity_mod_slim.png", top = Inches(4.25), left = Inches(2.5), height = Inches(1.25))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_phase.png", top = Inches(4.25), left = Inches(4), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_predicted_monaco.png", left = Inches(0.25), top = Inches(6), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_monaco_severity_mod_slim.png", left = Inches(2.5), top = Inches(6), height = Inches(1.25))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_monaco.png", top = Inches(6), left = Inches(4), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_featplot_CARTEx_200_severity_mod.png", top = Inches(7.75), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_aggplot_CARTEx_200_severity_mod_monaco_split_countsized.png", top = Inches(7.75), left = Inches(2), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_CARTEx_200.png", top = Inches(7.75), left = Inches(5), height = Inches(2))


# new slide
slide = add_slide(prs, blank_slide_layout, "Supplementary Figure for CARTEx applied to clinical oncology")

experiment = "GSE153931" # COVID-19 severity and SARS-CoV-2-reactive T cell exhaustion

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_LCMV_Tex.png", left = Inches(0), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_NKlike_Tex.png", left = Inches(2.5), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_BBD_Tex.png", left = Inches(5), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_PD1_Tex.png", left = Inches(7.5), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_scatter_CARTEx_200_LCMV.png", left = Inches(0), top = Inches(2.5), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_scatter_CARTEx_200_NKlike.png", left = Inches(2.5), top = Inches(2.5), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_scatter_CARTEx_200_BBD.png", left = Inches(5), top = Inches(2.5), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_scatter_CARTEx_200_PD1.png", left = Inches(7.5), top = Inches(2.5), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_LCMV_Tex.png", left = Inches(0), top = Inches(4.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_NKlike_Tex.png", left = Inches(2.5), top = Inches(4.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_BBD_Tex.png", left = Inches(5), top = Inches(4.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_PD1_Tex.png", left = Inches(7.5), top = Inches(4.75), height = Inches(2))





# new slide
slide = add_slide(prs, blank_slide_layout, "Figure 7. CARTEx applied to CD8 T cells in autoimmune disorders")

experiment = "GSE207935" # STAT3 GOF

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_affstatstim.png", left = Inches(0.25), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_phase.png", left = Inches(2.75), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_phase_affstatstim_slim.png", left = Inches(5), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_predicted_monaco.png", left = Inches(6.5), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_monaco_affstatstim_slim.png", left = Inches(8.75), top = Inches(0.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_featplot_CARTEx_200_affstatstim.png", left = Inches(0.25), top = Inches(3), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200.png", top = Inches(3), left = Inches(1.75), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_CARTEx_200.png", top = Inches(3), left = Inches(4.25), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_aggplot_CARTEx_200_affstatstim_monaco_split_countsized.png", top = Inches(3), left = Inches(7.25), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_plot_volcano_STAT3GOF_stim.png", top = Inches(5.25), left = Inches(0.25), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_plot_volcano_STAT3GOF_control.png", top = Inches(5.25), left = Inches(2.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_aggplot_heatmap_C5_C2.png", top = Inches(5.25), left = Inches(5.5), height = Inches(2))


# new slide
slide = add_slide(prs, blank_slide_layout, "Supplementary Figure. More STAT3 GOF")

experiment = "GSE207935" # STAT3 GOF


pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_vlnplot_affstatstim_exhaustion_markers.png", top = Inches(0.75), left = Inches(0.25), height = Inches(2))




# new slide
slide = add_slide(prs, blank_slide_layout, "Supplementary Figure. More STAT3 GOF")

experiment = "GSE207935" # STAT3 GOF

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_aggplot_vlnplot_all_C5_C2.png", left = Inches(0.25), top = Inches(0.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_sig_activation.png", left = Inches(0), top = Inches(3), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_sig_anergy.png", left = Inches(2), top = Inches(3), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_sig_senescence.png", left = Inches(4), top = Inches(3), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_sig_stemness.png", left = Inches(6), top = Inches(3), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_affstatstim.png", left = Inches(0), top = Inches(5.25), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_activation_affstatstim.png", left = Inches(0), top = Inches(7.5), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_anergy_affstatstim.png", left = Inches(0), top = Inches(9.75), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_senescence_affstatstim.png", left = Inches(0), top = Inches(12), height = Inches(2))






# new slide
slide = add_slide(prs, blank_slide_layout, "Supplementary Figure. Severe psoriasis is distinct from T cell exhaustion")

experiment = "GSE146264" # Severe psoriasis

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_Severity.png", top = Inches(0.75), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_Severity.png", top = Inches(0.75), left = Inches(2.5), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_phase.png", top = Inches(2.5), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_phase_severity_slim.png", top = Inches(2.5), left = Inches(2.5), height = Inches(1.25))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_phase.png", top = Inches(2.5), left = Inches(4), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_predicted_monaco.png", left = Inches(0.25), top = Inches(4.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_monaco_severity_slim.png", left = Inches(2.5), top = Inches(4.25), height = Inches(1.25))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_monaco.png", top = Inches(4.25), left = Inches(4), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_featplot_CARTEx_200_severity.png", top = Inches(6), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_aggplot_CARTEx_200_severity_monaco_split_countsized.png", top = Inches(6), left = Inches(2), height = Inches(2))

# pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_CARTEx_200.png", top = Inches(6), left = Inches(5), height = Inches(2))

# pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_vlnplot_severity_exhaustion_markers.png", top = Inches(3), left = Inches(0.25), height = Inches(2))





# new slide
slide = add_slide(prs, blank_slide_layout, "Supplementary Figure. BEACH domain deficiencies (autoimmunity caused by defective CTLA-4 trafficking)")

experiment = "GSE196606" # NBEAL2 and LRBA deficiencies

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_group.png", top = Inches(0.75), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_group.png", top = Inches(0.75), left = Inches(2.5), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_phase.png", top = Inches(2.5), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_phase_group_slim.png", top = Inches(2.5), left = Inches(2.5), height = Inches(1.25))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_phase.png", top = Inches(2.5), left = Inches(4), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_predicted_monaco.png", left = Inches(0.25), top = Inches(4.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_barplot_monaco_group_slim.png", left = Inches(2.5), top = Inches(4.25), height = Inches(1.25))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_umap_CARTEx_200_monaco.png", top = Inches(4.25), left = Inches(4), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_prepare_featplot_CARTEx_200_group.png", top = Inches(6), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_aggplot_CARTEx_200_group_monaco_split_countsized.png", top = Inches(6), left = Inches(2), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_query_agg_aggplot_CARTEx_200.png", top = Inches(6), left = Inches(5), height = Inches(2))









# new slide
slide = add_slide(prs, blank_slide_layout, "Supplementary Figure. T cell exhaustion in neurodegenerative diseases")

experiment = "Zenodo3993994" # 

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_blood_prepare_umap_disease.png", top = Inches(0.75), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_blood_prepare_umap_CARTEx_200_disease.png", top = Inches(0.75), left = Inches(2.5), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_blood_prepare_umap_phase.png", top = Inches(2.5), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_blood_prepare_barplot_phase_disease_slim.png",  top = Inches(2.5), left = Inches(2.5), height = Inches(1.25))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_blood_prepare_umap_CARTEx_200_phase.png", top = Inches(2.5), left = Inches(4), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_blood_prepare_umap_predicted_monaco.png", top = Inches(4.25), left = Inches(0.25), height = Inches(1.5))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_blood_prepare_barplot_monaco_disease_slim.png", top = Inches(4.25), left = Inches(2.5), height = Inches(1.25))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_blood_prepare_umap_CARTEx_200_monaco.png", top = Inches(4.25), left = Inches(4), height = Inches(1.5))

# pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_blood_prepare_vlnplot_disease_exhaustion_markers.png", top = Inches(6), left = Inches(0.25), height = Inches(2))

# pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_blood_aggplot_CARTEx_200_group_monaco_split_countsized.png", top = Inches(6), left = Inches(4), height = Inches(2))

pic = slide.shapes.add_picture("./experiments/" + experiment + "/plots/" + experiment + "_blood_explore_plot_volcano_Parkinson_blood.png", top = Inches(6), left = Inches(0.25), height = Inches(2))









prs.save('organized_figures.pptx')





