# Code for variance composition analysis and betadiversity
# Knepp and Boothby 2022 FINAL DRAFT  13-09-2023

# ALPHA DIVERSITY METRICS
# -----------------------------------------------------------------------------#
#### SETUP: ####
# -----------------------------------------------------------------------------#
source("scripts/headers.R")

# Load data
knepp <- read.csv("raw-data/k-plants.csv")
boothby <- read.csv("raw-data/b-plants.csv")
tree <- read.tree("raw-data/Vascular_Plants_rooted.dated.tre")

# Initial format and get rid of  rows that contain *s
# KNEPP
knepp$Species <- stringi::stri_trans_general(knepp$Species, "latin-ascii")
knepp$Species <- gsub(" ", "_", knepp$Species, fixed=TRUE)
knepp[is.na(knepp)] <- 0
knepp <- knepp %>% 
  filter(! grepl('\\*', Plot))

# BOOTHBY
boothby$Species <- stringi::stri_trans_general(boothby$Species, "latin-ascii")
boothby$Species <-gsub(" ", "_", boothby$Species, fixed=TRUE)
boothby[is.na(boothby)] <- 0
boothby <- boothby %>% 
  filter(! grepl('\\*', Plot))

# "Build" phylogenies
ktree <- congeneric.merge(tree, unique(knepp$Species))
btree <- congeneric.merge(tree, unique(boothby$Species))

# Format for pez etc.
#KNEPP
k.comm <- with(knepp, tapply(Cover, list(paste(Site,Block,Group,Plot), Species), mean, rm.na=TRUE))
k.comm[is.na(k.comm)] <- 0
k.c.data <- comparative.comm(ktree, k.comm)

#BOOTHBY
b.comm <- with(boothby, tapply(Cover, list(paste(Site,Block,Plot), Species), mean, rm.na=TRUE))
b.comm[is.na(b.comm)] <- 0
b.c.data <- comparative.comm(btree, b.comm)

# Make hierarchical groupings for each site
#KNEPP
raw.groups <- strsplit(sites(k.c.data), " ")
site.no <- sapply(raw.groups, function(x) x[4])
env <- data.frame(
  block = sapply(raw.groups, function(x) x[2]),
  fractal = sapply(raw.groups, function(x) x[3]),
  major = substr(site.no, 0, 1),
  minor = substr(site.no, 2, 2)
)

env$fractal[env$fractal=="NA"] <- "A" # R is seeing "NA" and sometimes coercing to the concept of an NA
rownames(env) <- sites(k.c.data)
write.csv(env, file = "kneppfile.csv")

env$X <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0)
env$Y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1,1,0,0,0,0,1,0)
env$Z <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,1,1,1,1)

# Uncomment below to check the site names, species dropped, and species included
#sites(k.c.data) 
#k.c.data$dropped$comm.sp.lost
#species(k.c.data)

# Subsetting to update object 'knepp' - losing 'sticks', 'bare soil' etc
k_wanted_spp <- species(k.c.data)
knepp <- subset(knepp, Species %in% k_wanted_spp)
write.csv(knepp, file = "clean-data/knepp_clean.csv")
knepp <- read.csv("clean-data/knepp_clean.csv")

# Merge everything in a compative.comm dataset 
k.c.data <- comparative.comm(tree(k.c.data), comm(k.c.data), env=env)
saveRDS(k.c.data, "clean-data/k.c.data.RDS")

# BOOTHBY
raw.groups <- strsplit(sites(b.c.data), " ")
site.no <- sapply(raw.groups, function(x) x[3])
env <- data.frame(
  fractal = sapply(raw.groups, function(x) x[2]),
  major = substr(site.no, 0, 1),
  minor = substr(site.no, 2, 2)
)
env$fractal[env$fractal=="NA"] <- "A" 
rownames(env) <- sites(b.c.data)

#sites(b.c.data) 
#b.c.data$dropped$comm.sp.lost
#species(b.c.data)

b_wanted_spp <- species(b.c.data)
boothby <- subset(boothby, Species %in% b_wanted_spp)

write.csv(boothby, file = "clean-data/boothby_clean.csv")
boothby <- read.csv("clean-data/boothby_clean.csv")
# Merge everything in a compative.comm dataset 
b.c.data <- comparative.comm(tree(b.c.data), comm(b.c.data), env=env)
saveRDS(b.c.data, "clean-data/b.c.data.RDS")
