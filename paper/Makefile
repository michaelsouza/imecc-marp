TEXFILE = proximal_bcd_graph_realization.tex
PDFFILE = proximal_bcd_graph_realization.pdf

FIGURES = $(wildcard figures/*)
TABLES = $(wildcard tables/*)

all: $(PDFFILE)

$(PDFFILE): $(TEXFILE) references.bib $(FIGURES) $(TABLES)
	latexmk -pdf -interaction=nonstopmode $(TEXFILE)

clean:
	latexmk -C $(TEXFILE)
	rm -f *.bbl *.run.xml
