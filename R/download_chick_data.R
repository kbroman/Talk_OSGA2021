# The file contains a number of objects, including the following:
#
# - `dataset.mrna` - a list containing the mRNA expression data
#   (`r nrow(dataset.mrna$annot.samples)` mice)
# - `dataset.protein` - a list containing the proteomics data
#   (`r nrow(dataset.protein$annot.samples)` mice)
# - `genoprobs` - genotype probabilities across the genome from
#   `calc_genoprob()`, reduced to 8 allele dosages
#   (`r nrow(genoprobs[[1]])` mice)
#
# The mRNA and protein datasets are lists containing a number of
# objects. Here are the key pieces of `dataset.mrna`:
#
# - `annot.mrna` - metadata about the mRNA traits
#   (`r vec2string(colnames(dataset.mrna$annot.mrna))`)
# - `covar.matrix` -
#   `r nrow(dataset.mrna$covar.matrix)` &times; `r ncol(dataset.mrna$covar.matrix)` matrix
#   (including `sexM`, a 0/1 indicator with 1 = male)
# - `data` - a list with "norm", "raw", "rz", each a
#   `r nrow(dataset.mrna$data$norm)` &times; `r add_commas(ncol(dataset.mrna$data$norm))` matrix
# - `lod.peaks` - a list with "additive", "diet", "sex_int",
#   each a data frame with `r vec2string(colnames(dataset.mrna$lod.peaks$additive))`
#
# The `dataset.protein` object is similar. Here are the key pieces:
#
# - `annot.protein` - metadata about the proteins
#   (`r nrow(dataset.protein$annot.samples)` mice)
# - `covar.matrix` -
#   `r nrow(dataset.protein$covar.matrix)` &times; `r ncol(dataset.protein$covar.matrix)` matrix
#   (including `sexM`, a 0/1 indicator with 1 = male)
# - `data` - a matrix of size
#   `r nrow(dataset.protein$data)` &times; `r add_commas(ncol(dataset.protein$data))`
# - `lod.peaks` - a list with "additive", "diet", "sex_int",
#   each a data frame with `r vec2string(colnames(dataset.protein$lod.peaks$additive))`

library(here)

url <- "https://churchilllab.jax.org/qtlviewer/svenson/DOHFD/rdata"
file <- here("Data/Svenson_DO850_for_eQTL_viewer_v9.RData")
if(!file.exists(file)) {
    download.file(url, file)
}

load(file)
protein <- dataset.protein$data
saveRDS(protein, here("Data/chick_protein.rds"))
mrna <- dataset.mrna$data$norm
saveRDS(mrna, here("Data/chick_mrna.rds"))
saveRDS(map, here("Data/chick_map.rds"))
saveRDS(genoprobs, here("Data/chick_genoprobs.rds"))
