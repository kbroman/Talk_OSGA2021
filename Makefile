R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

STEM = osga2021

FIGS = Figs/triple_asso.pdf \
	   Figs/b6btbr_plates.pdf \
	   Figs/b6btbr_expr_swaps.pdf

all: docs/$(STEM).pdf

docs/$(STEM).pdf: $(STEM).pdf
	cp $< $@

$(STEM).pdf: $(STEM).tex header.tex $(FIGS)
	xelatex $<

Figs/%.pdf: R/%.R
	cd R;R CMD BATCH $(R_OPTS) $(^F)
