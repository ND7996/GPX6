LATEXMK = latexmk -pdf -interaction=nonstopmode -halt-on-error
SOURCE_DIR = Files
BUILD_DIR = build
FIGURE_DIR = Figures_resub
ORIGINAL_DIR = ../submission_JCIM/Files

EDITABLE_TEX = $(wildcard $(SOURCE_DIR)/*.tex)
REFERENCES = $(SOURCE_DIR)/references.bib
FIGURES = $(wildcard $(FIGURE_DIR)/*)
ORIGINAL_TEX = $(wildcard $(ORIGINAL_DIR)/*.tex)

FINAL_PDFS = \
	GPX6_manuscript_clean.pdf \
	GPX6_manuscript_marked.pdf \
	GPX6_suppinfo_clean.pdf \
	GPX6_suppinfo_marked.pdf \
	Cover_letter_and_response_JCIM.pdf

.PHONY: all manuscript si response clean distclean help

all: $(FINAL_PDFS)

manuscript: GPX6_manuscript_clean.pdf GPX6_manuscript_marked.pdf

si: GPX6_suppinfo_clean.pdf GPX6_suppinfo_marked.pdf

response: Cover_letter_and_response_JCIM.pdf

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/GPX6_suppinfo.pdf: $(EDITABLE_TEX) $(BUILD_DIR)/references.bib $(FIGURES) | $(BUILD_DIR)
	cd $(SOURCE_DIR) && $(LATEXMK) -outdir=../$(BUILD_DIR) GPX6_suppinfo.tex

$(BUILD_DIR)/GPX6_manuscript.pdf: $(EDITABLE_TEX) $(BUILD_DIR)/references.bib $(FIGURES) $(BUILD_DIR)/GPX6_suppinfo.pdf | $(BUILD_DIR)
	cd $(SOURCE_DIR) && $(LATEXMK) -outdir=../$(BUILD_DIR) GPX6_manuscript.tex

$(BUILD_DIR)/marked-sources.stamp: $(EDITABLE_TEX) $(REFERENCES) $(ORIGINAL_TEX) \
		../.codex/scripts/regenerate_resubmission_from_original.sh \
		../.codex/scripts/generate_conditional_si.sh | $(BUILD_DIR)
	../.codex/scripts/regenerate_resubmission_from_original.sh ..
	touch $@

$(BUILD_DIR)/references.bib: $(REFERENCES) | $(BUILD_DIR)
	cp $< $@

$(BUILD_DIR)/GPX6_suppinfo_marked.pdf: $(BUILD_DIR)/marked-sources.stamp $(BUILD_DIR)/references.bib $(FIGURES)
	cd $(BUILD_DIR) && $(LATEXMK) -jobname=GPX6_suppinfo_marked \
		-usepretex='\def\SHOWTRACKED{1}' GPX6_suppinfo_marked.tex

$(BUILD_DIR)/GPX6_manuscript_marked.pdf: $(BUILD_DIR)/marked-sources.stamp \
		$(BUILD_DIR)/references.bib $(BUILD_DIR)/GPX6_suppinfo_marked.pdf $(FIGURES)
	cd $(BUILD_DIR) && $(LATEXMK) -jobname=GPX6_manuscript_marked \
		-usepretex='\def\SHOWTRACKED{1}' GPX6_manuscript_marked.tex

$(BUILD_DIR)/Cover_letter_and_response_JCIM.pdf: \
		$(SOURCE_DIR)/Cover_letter_and_response_JCIM.tex | $(BUILD_DIR)
	cd $(SOURCE_DIR) && $(LATEXMK) -outdir=../$(BUILD_DIR) Cover_letter_and_response_JCIM.tex

GPX6_suppinfo_clean.pdf: $(BUILD_DIR)/GPX6_suppinfo.pdf
	cp $< $@

GPX6_manuscript_clean.pdf: $(BUILD_DIR)/GPX6_manuscript.pdf
	cp $< $@

GPX6_suppinfo_marked.pdf: $(BUILD_DIR)/GPX6_suppinfo_marked.pdf
	cp $< $@

GPX6_manuscript_marked.pdf: $(BUILD_DIR)/GPX6_manuscript_marked.pdf
	cp $< $@

Cover_letter_and_response_JCIM.pdf: $(BUILD_DIR)/Cover_letter_and_response_JCIM.pdf
	cp $< $@

clean:
	rm -rf $(BUILD_DIR)

distclean: clean
	rm -f $(FINAL_PDFS)

help:
	@echo "Edit only resubmission_JCIM/Files/*.tex"
	@echo "Run 'make' here to regenerate the five final PDFs."
