R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

STEM = osga2021

FIGS = Figs/triple_asso.pdf \
	   Figs/b6btbr_plates.pdf \
	   Figs/b6btbr_expr_swaps.pdf \
	   Figs/missing_data.pdf \
	   Figs/xydosage.pdf \
	   Figs/hist_compare_geno.pdf \
	   Figs/protein_dups.png \
	   Figs/mrna_dups.png \
	   Figs/sex_mrna.pdf \
	   Figs/sex_protein.pdf \
	   Figs/evp_mixup_scheme.pdf \
	   Figs/gve_mixup_scheme.pdf \
	   Figs/gvp_mixup_scheme.pdf \
	   Figs/heatmap_mrna_v_protein.png \
	   Figs/heatmap_dna_v_mrna.png \
	   Figs/heatmap_dna_v_prot.png \
	   Figs/hs_light.pdf

all: docs/$(STEM).pdf docs/$(STEM)_notes.pdf

docs/%.pdf: %.pdf
	cp $< $@

$(STEM).pdf: $(STEM).tex header.tex $(FIGS)
	xelatex $<

$(STEM)_notes.pdf: $(STEM)_notes.tex header.tex $(FIGS)
	xelatex $<
	pdfnup $@ --nup 1x2 --no-landscape --paper letterpaper --frame true --scale 0.9
	mv $(STEM)_notes-nup.pdf $@

$(STEM)_notes.tex: $(STEM).tex Ruby/createVersionWithNotes.rb
	Ruby/createVersionWithNotes.rb $< $@

Figs/%.pdf: R/%.R
	cd R;R CMD BATCH $(R_OPTS) $(^F)

Figs/%.png: R/%.R
	cd R;R CMD BATCH $(R_OPTS) $(^F)

Figs/heatmap_mrna_v_protein.png: R/mrna_v_protein.R
	cd R;R CMD BATCH $(R_OPTS) $(^F)

Figs/heatmap_dna_v_mrna.png: R/dna_v_mrna.R
	cd R;R CMD BATCH $(R_OPTS) $(^F)

Figs/heatmap_dna_v_protein.png: R/dna_v_prot.R
	cd R;R CMD BATCH $(R_OPTS) $(^F)
