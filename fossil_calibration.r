library(phangorn)
library(ape)
library(phytools)

tree <- read.newick("~/Work/species_phylogeny_eucera/Whole_phylogeny.rooted.newick")
#calibrate <- makeChronosCalib(tree, node = "root", age.max = 179)

nodes <- c(
  getMRCA(tree, tip = c("Colletes_gigas", "Apis_melliFfera")),
  getMRCA(tree, tip = c("Bombus_impatiens", "Melipona_quadrifasciata")),
  getMRCA(tree, tip = c("Apis_melliFfera", "Euglossa_dilemma")),
  getMRCA(tree, tip = c("Bombus_impatiens", "Habropoda_laboriosa")),
  getMRCA(tree, tip = c("Megachile_rotundata", "Osmia_lignaria"))
)

age.min <- c(0, 53, 28, 60, 58)
age.max <- c(110, 110, 110, 110, 110)
soft.bounds <- c(FALSE, FALSE, FALSE, FALSE, FALSE)

mycalibration <- data.frame(nodes, age.min, age.max, soft.bounds)

timetree <- chronos(tree, lambda = 1, model = "relaxed", calibration = mycalibration, control = chronos.control() )
#timetree <- chronos(tree, lambda = 1, model = "correlated", calibration = mycalibration, control = chronos.control() )
#timetree <- chronos(tree, lambda = 1, model = "discrete", calibration = mycalibration, control = chronos.control() )

is.ultrametric(timetree)
write.tree(timetree, "~/Work/species_phylogeny_eucera/Whole_phylogeny.rooted.utrametric.newick")
