from pptx import Presentation
from pptx.util import Inches

experiment = "GSE136874"
description = "comparison of CD19 and GD2 CAR T cells on day 10 in culture"

prs = Presentation()
# set width and height to 16 and 9 inches.
prs.slide_width = Inches(16)
prs.slide_height = Inches(9)

blank_slide_layout = prs.slide_layouts[6]

def add_slide(prs, layout, subtitle):
	slide = prs.slides.add_slide(layout)
	txBox = slide.shapes.add_textbox(left = Inches(0), top = Inches(0), width = Inches(1), height = Inches(1))
	tf = txBox.text_frame
	tf.text = experiment + ": " + description
	p = tf.add_paragraph()
	p.text = subtitle
	return slide

# slide 1
slide = add_slide(prs, blank_slide_layout, "UMAPs")

pic = slide.shapes.add_picture("./plots/" + experiment + "_prepare_umap_CAR.jpeg", left = Inches(0.5), top = Inches(1), height = Inches(3))

pic = slide.shapes.add_picture("./plots/" + experiment + "_prepare_umap_phase.jpeg", left = Inches(5.5), top = Inches(1), height = Inches(3))

pic = slide.shapes.add_picture("./plots/" + experiment + "_prepare_umap_predicted_monaco.jpeg", left = Inches(10.5), top = Inches(1), height = Inches(3))

pic = slide.shapes.add_picture("./plots/" + experiment + "_prepare_barplot_phase_CAR.jpeg", left = Inches(5.5), top = Inches(4), height = Inches(3))

pic = slide.shapes.add_picture("./plots/" + experiment + "_prepare_barplot_monaco_CAR.jpeg", left = Inches(10.5), top = Inches(4), height = Inches(3))


# slide 2
slide = add_slide(prs, blank_slide_layout, "UMAPs")


prs.save(experiment + '_results.pptx')

