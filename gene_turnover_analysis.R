library(phytools)
library(ape)
library(geiger)
library(dplyr)
library(brms)
library(tidyverse)
library(magrittr)
library(nlme)
library(caper)
library(MASS)
library(GLMcat)
library(ggplot2)
library(gtools)
library(MuMIn)
library(phylolm)
library(glmmTMB)
library(multcomp)
library(ggtree)
library(emmeans)
library(RColorBrewer)
library(car)
library(stargazer)

################################################################################
## ORs
################################################################################################
setwd('~/Work/gene_turnover_analysis/')
# OR analysis:
# load tree and OR counts
tree <- read.tree("whole_phylogeny_rooted_ultrametric.newick")
prune <- c("Bombus_skorikovi", "Bombus_turneri", "Eufriesea_mexicana")

# filtered dataset:
#prune <- c("Bombus_skorikovi", "Bombus_turneri", "Eufriesea_mexicana",
#"Bombus_polaris", "Bombus_cullumanus", "Bombus_sibiricus",
#"Bombus_haemorrhoidalis", "Bombus_soroeensis", "Bombus_confusus",
#"Bombus_superbus", "Tetragonula_hockingsi","Tetragonula_mellipes",
#"Tetragonula_carbonaria","Tetragonula_clypearis", "Lepidotrigona_ventralis_hoosana",
#"Apis_dorsata", "Apis_florea", "Ceratina_australensis",
#"Osmia_bicornis_bicornis")

tree <- drop.tip(tree, prune)

g <- ggtree(tree)
align <- get_taxa_name(g)

OR_subfamily_counts <- read_tsv(file = "OR_subfamily_counts.txt")
OR_subfamily_counts <- subset(OR_subfamily_counts, select = -1)
OR_subfamily_counts <- t(OR_subfamily_counts)
OR_subfamily_counts <- OR_subfamily_counts[-1,]
OR_subfamily_counts <- OR_subfamily_counts[!(row.names(OR_subfamily_counts) %in% prune),]
Total <- rowSums(OR_subfamily_counts)

diet_breadth <- read.csv("diet_breadths_final.csv", header = T)
diet_breadth <- diet_breadth[match(align, diet_breadth$species),]

niimura_rates <- read.csv('OR_birth_death_calculations_input_only_tips.csv', header = T)
niimura_rates <- niimura_rates[match(align, niimura_rates$Species),]

counts <- read.csv("OR.consensus.aligned.counts", header = F)
colnames(counts) <- c("Counts", "Species")
counts <- counts[match(align, counts$Species),]

df <- cbind(diet_breadth, Total, niimura_rates$Mean.birth, niimura_rates$Mean.death, counts$Counts)
colnames(df) <- c("Species", "Plant_families", "Diet_breadth", "busco", "Total", "Niimura_birth_rate", "Niimura_death_rate", "Counts")
rownames(df) <- c(1:51)
df$Diet_breadth <- as.factor(df$Diet_breadth)

################################################################################
# Ancestral State reconstruction
db_asr <- diet_breadth[, c(1,3)]
row.names(db_asr) <- db_asr$species

db<-setNames(
  db_asr$diet_breadth,
  rownames(db_asr))
head(db)

db_er<-fitMk(tree,db,model="ER",pi="fitzjohn")
# db_sym<-fitMk(tree,db,model="SYM",pi="fitzjohn")
# db_ard<-fitMk(tree,db,model="ARD",pi="fitzjohn")
# db_aov <- anova(db_er, db_sym, db_ard)

db_ancr <- ancr(db_er, type = "marginal")
foo <- as.data.frame(db_ancr$ace)
foo$node <- rownames(foo)

write.csv(foo, file = "marginal_ancestral_state_probabilities.csv", row.names = FALSE)

################################################################################
# Phylogenetic signal diet breadth
cat_db <- df[, c(1,3)]
cat_db$diet_num <- as.numeric(as.factor(cat_db$Diet_breadth))
phylosig(tree, cat_db$diet_num, method = "K", test = T)

################################################################################
# PGLS
phylosig(tree, df$busco, method = "lambda", test = T)
phylosig(tree, df$Counts, method = "lambda", test = T)
phylosig(tree, df$Niimura_birth_rate, method = "lambda", test = T)
phylosig(tree, df$Niimura_death_rate, method = "lambda", test = T)

################################################################################
# PGLS with BUSCO
temptree <- tree
temptree$edge.length <- tree$edge.length * 100
rownames(df) <- align
name.check(tree, df)

pglsModel_busco_db <- gls(busco ~ Diet_breadth, correlation = corBrownian(1, phy = temptree, form = ~Species), data = df, method = "ML") 
#pglsModel_busco_db <- gls(busco ~ Diet_breadth, correlation = corPagel(1, phy = temptree, form = ~Species), data = df, method = "ML") 
#pglsModel_busco_db <- gls(busco ~ Diet_breadth, correlation = corBlomberg(1, phy = temptree, form = ~Species), data = df, method = "ML") 
summary(pglsModel_busco_db)
car::Anova(pglsModel_busco_db)

pglsModel_counts <- gls(Counts ~ Diet_breadth, correlation = corBrownian(1, phy = temptree, form = ~Species), data = df, method = "ML")
summary(pglsModel_counts)
car::Anova(pglsModel_counts)

pglsModel_busco <- gls(Counts ~ busco, correlation = corBrownian(1, phy = temptree, form = ~Species),data = df, method = "ML")
summary(pglsModel_busco)
car::Anova(pglsModel_busco)

pglsModel_counts_busco <- gls(Counts ~ Diet_breadth + busco, correlation = corBrownian(1, phy = temptree, form = ~Species), data = df, method = "ML")
summary(pglsModel_counts_busco)
car::Anova(pglsModel_counts_busco)
ph <- summary(glht(pglsModel_counts_busco, linfct = mcp(Diet_breadth = "Tukey")))

pglsModel_birth_busco <- gls(Niimura_birth_rate ~ busco, correlation = corBrownian(1, phy = temptree, form = ~Species), data = df, method = "ML") 
summary(pglsModel_birth_busco)
car::Anova(pglsModel_birth_busco)

pglsModel_death_busco <- gls(Niimura_death_rate ~ busco, correlation = corBrownian(1, phy = temptree, form = ~Species), data = df, method = "ML")
summary(pglsModel_death_busco)
car::Anova(pglsModel_death_busco)

pglsModel_birth <- gls(Niimura_birth_rate ~ Diet_breadth + busco, correlation = corBrownian(1, phy = temptree, form = ~Species), data = df, method = "ML")
summary(pglsModel_birth)
car::Anova(pglsModel_birth)

pglsModel_death <- gls(Niimura_death_rate ~ Diet_breadth + busco, correlation = corBrownian(1, phy = temptree, form = ~Species), data = df, method = "ML")
summary(pglsModel_death)
car::Anova(pglsModel_death)
ph <- summary(glht(pglsModel_death, linfct = mcp(Diet_breadth = "Tukey")))

################################################################################
## GRS
################################################################################
g <- ggtree(tree)
align <- get_taxa_name(g)
GR_subfamily_counts <- read_tsv(file = "GR_subfamily_counts.txt")
GR_subfamily_counts <- subset(GR_subfamily_counts, select = -1)
GR_subfamily_counts <- t(GR_subfamily_counts)
GR_subfamily_counts <- GR_subfamily_counts[-1,]
GR_subfamily_counts <- GR_subfamily_counts[!(row.names(GR_subfamily_counts) %in% prune),]
Total <- rowSums(GR_subfamily_counts)

diet_breadth <- read.csv("diet_breadths_final.csv", header = T)
diet_breadth <- diet_breadth[match(align, diet_breadth$species),]

niimura_rates <- read.csv('GR_birth_death_calculations_input_only_tips.csv', header = T)
niimura_rates <- niimura_rates[match(align, niimura_rates$Species),]

counts <- read.csv("GR.consensus.aligned.counts", header = F)
colnames(counts) <- c("Counts", "Species")
counts <- counts[match(align, counts$Species),]

df <- cbind(diet_breadth, Total, niimura_rates$Mean.birth, niimura_rates$Mean.death, counts$Counts)
colnames(df) <- c("Species", "Plant_families", "Diet_breadth", "busco", "Total", "Niimura_birth_rate", "Niimura_death_rate", "Counts")
rownames(df) <- c(1:51)
df$Diet_breadth <- as.factor(df$Diet_breadth)

################################################################################
# PGLS
phylosig(tree, df$Counts, method = "lambda", test = T)
phylosig(tree, df$Niimura_birth_rate, method = "lambda", test = T)
phylosig(tree, df$Niimura_death_rate, method = "lambda", test = T)

################################################################################
# PGLS with BUSCO
temptree <- tree
temptree$edge.length <- tree$edge.length * 100
rownames(df) <- align
name.check(tree, df)

pglsModel_busco <- gls(Counts ~ busco, correlation = corBrownian(1, phy = temptree, form = ~Species), data = df, method = "ML") 
summary(pglsModel_busco)
car::Anova(pglsModel_busco)

pglsModel_counts <- gls(Counts ~ Diet_breadth, correlation = corBrownian(1, phy = temptree, form = ~Species), data = df, method = "ML") 
summary(pglsModel_counts)
car::Anova(pglsModel_counts)

pglsModel_counts_busco <- gls(Counts ~ Diet_breadth + busco, correlation = corBrownian(1, phy = temptree, form = ~Species), data = df, method = "ML") 
car::Anova(pglsModel_counts_busco)

pglsModel_birth_busco <- gls(Niimura_birth_rate ~ busco, correlation = corBrownian(1, phy = temptree, form = ~Species), data = df, method = "ML") 
summary(pglsModel_birth_busco)
car::Anova(pglsModel_birth_busco)

pglsModel_death_busco <- gls(Niimura_death_rate ~ busco, correlation = corBrownian(1, phy = temptree, form = ~Species), data = df, method = "ML") 
summary(pglsModel_death_busco)
car::Anova(pglsModel_death_busco)

pglsModel_birth <- gls(Niimura_birth_rate ~ Diet_breadth + busco, correlation = corBrownian(1, phy = temptree, form = ~Species), data = df, method = "ML")
summary(pglsModel_birth)
car::Anova(pglsModel_birth)
ph <- summary(glht(pglsModel_birth, linfct = mcp(Diet_breadth = "Tukey")))

pglsModel_death <- gls(Niimura_death_rate ~ Diet_breadth + busco, correlation = corBrownian(1, phy = temptree, form = ~Species), data = df, method = "ML")
summary(pglsModel_death)
car::Anova(pglsModel_death)

################################################################################
## IRS
################################################################################
g <- ggtree(tree)
align <- get_taxa_name(g)
IR_subfamily_counts <- read_tsv(file = "IR_subfamily_counts.txt")
IR_subfamily_counts <- subset(IR_subfamily_counts, select = -1)
IR_subfamily_counts <- t(IR_subfamily_counts)
IR_subfamily_counts <- IR_subfamily_counts[-1,]
IR_subfamily_counts <- IR_subfamily_counts[!(row.names(IR_subfamily_counts) %in% prune),]
Total <- rowSums(IR_subfamily_counts)

diet_breadth <- read.csv("diet_breadths_final.csv", header = T)
diet_breadth <- diet_breadth[match(align, diet_breadth$species),]

niimura_rates <- read.csv('IR_birth_death_calculations_input_only_tips.csv', header = T)
niimura_rates <- niimura_rates[match(align, niimura_rates$Species),]

counts <- read.csv("IR.consensus.aligned.counts", header = F)
colnames(counts) <- c("Counts", "Species")
counts <- counts[match(align, counts$Species),]

df <- cbind(diet_breadth, Total, niimura_rates$Mean.birth, niimura_rates$Mean.death, counts$Counts)
colnames(df) <- c("Species", "Plant_families", "Diet_breadth", "busco", "Total", "Niimura_birth_rate", "Niimura_death_rate", "Counts")
rownames(df) <- c(1:51)
df$Diet_breadth <- as.factor(df$Diet_breadth)

################################################################################
# PGLS
phylosig(tree, df$Niimura_birth_rate, method = "lambda", test = T)
phylosig(tree, df$Niimura_death_rate, method = "lambda", test = T)
phylosig(tree, df$Counts, method = "lambda", test = T)

################################################################################
# PGLS with BUSCO
temptree <- tree
temptree$edge.length <- tree$edge.length * 100
rownames(df) <- align
name.check(tree, df)

pglsModel_busco <- gls(Counts ~ busco, correlation = corBrownian(1, phy = temptree, form = ~Species), data = df, method = "ML") 
summary(pglsModel_busco)
car::Anova(pglsModel_busco)

pglsModel_counts <- gls(Counts ~ Diet_breadth, correlation = corBrownian(1, phy = temptree, form = ~Species), data = df, method = "ML") 
summary(pglsModel_counts)
car::Anova(pglsModel_counts)

pglsModel_counts <- gls(Counts ~ Diet_breadth + busco, correlation = corBrownian(1, phy = temptree, form = ~Species), data = df, method = "ML") 
summary(pglsModel_counts)
car::Anova(pglsModel_counts)
ph <- summary(glht(pglsModel_counts, linfct = mcp(Diet_breadth = "Tukey")))

pglsModel_birth_busco <- gls(Niimura_birth_rate ~ busco, correlation = corBrownian(1, phy = temptree, form = ~Species), data = df, method = "ML") 
summary(pglsModel_birth_busco)
car::Anova(pglsModel_birth_busco)

pglsModel_death_busco <- gls(Niimura_death_rate ~ busco, correlation = corBrownian(1, phy = temptree, form = ~Species), data = df, method = "ML") 
summary(pglsModel_death_busco)
car::Anova(pglsModel_death_busco)

pglsModel_birth <- gls(Niimura_birth_rate ~ Diet_breadth + busco, correlation = corBrownian(1, phy = temptree, form = ~Species),
                       data = df, method = "ML")
summary(pglsModel_birth)
car::Anova(pglsModel_birth)

pglsModel_death <- gls(Niimura_death_rate ~ Diet_breadth + busco, correlation = corBrownian(1, phy = temptree, form = ~Species), data = df, method = "ML")
summary(pglsModel_death)
car::Anova(pglsModel_death)
ph <- summary(glht(pglsModel_death, linfct = mcp(Diet_breadth = "Tukey")))