R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

STEM = osga2021

#FIGS =

all: docs/$(STEM).pdf

docs/$(STEM).pdf: $(STEM).pdf
	cp $< $@

$(STEM).pdf: $(STEM).tex header.tex $(FIGS)
	xelatex $<
