# Makefile for compiling the paper

# Variables
LATEX = pdflatex
BIBTEX = bibtex
MAIN = GPX6_manuscript
SI = GPX6_suppinfo
TEX_FILES = $(wildcard *.tex)
MAIN_FIGURES = \
	../Figures/Fig1A.png \
	../Figures/Fig1B.png \
	../Figures/Fig2.png \
	../Figures/Fig3A.png \
	../Figures/Fig3B.png \
	../Figures/Fig3C.png \
	../Figures/Fig3D.png \
	../Figures/Fig4A.png \
	../Figures/Fig4B.png \
	../Figures/Fig5A.png \
	../Figures/Fig5B.png \
	../Figures/Fig6.png \
	../Figures/Fig7.png \
	../Figures/finalTOC.png
SI_FIGURES = \
	../Figures/FigS1.png \
	../Figures/FigS2.png \
	../Figures/FigS3A.png \
	../Figures/FigS3B.png \
	../Figures/FigS3C.png \
	../Figures/FigS3D.png \
	../Figures/FigS4A.png \
	../Figures/FigS4B.png \
	../Figures/FigS5.png \
	../Figures/FigS6.png \
	../Figures/FigS7.png \
	../Figures/FigS8.png \
	../Figures/FigS9.png \
	../Figures/FigS10.png
LATEXMK = latexmk -pdf

.PRECIOUS: $(SI).pdf $(MAIN).pdf

# Rules
all: $(SI).pdf $(MAIN).pdf 

$(MAIN).pdf: $(SI).pdf $(TEX_FILES) $(MAIN_FIGURES) references.bib
	$(LATEXMK) $(MAIN).tex

$(SI).pdf: $(TEX_FILES) $(SI_FIGURES) references.bib
	$(LATEXMK) $(SI).tex

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *.pdf acs-* *.lot *.lof

distclean: clean
	rm -f $(MAIN).pdf
