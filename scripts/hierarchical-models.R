source("src/headers.R")

# -----------------------------------------------------------------------------#
#### ALPHA DIVERSITY INDICES AND PLOTTING ####
# -----------------------------------------------------------------------------#

knepp <- read.csv("clean-data/knepp_clean.csv")
k.c.data <- readRDS("clean-data/k.c.data.RDS")

boothby <- read.csv("clean-data/boothby_clean.csv")
b.c.data <- readRDS("clean-data/b.c.data.RDS")

# ----------------------------------------------------------------------------#
#### Species richness ####
# ----------------------------------------------------------------------------#

barnamesknepp <- c("Block", "Fractal", "Major", "Minor", "Residual")
barnamesboothby <- c("Fractal", "Major", "Minor", "Residual")

par(mar=c(5.1,5.5,4.1,2.1))
k.rich <- .ses.mntd(k.c.data)$ntaxa

model <- stan_lmer(k.rich ~ (X)+(Y)+(Z)+ (1|block/fractal/major/minor), data=env(k.c.data)) 

getME(lmer(k.rich ~ (1|major/minor), data=env(k.c.data)), "mmList")

please <- with(k.c.data$env, cbind(X,Y,Z))
model <- stan_lmer(k.rich ~ (1|please), data=env(k.c.data)) 


#model <- stan_lmer(k.rich ~ (X)+(Y)+(Z)+ (1|block/fractal/major/minor), data=env(k.c.data)) 
#model <- stan_lmer(k.rich ~ (1|(X+Y+Z)/block/fractal/major/minor), data=env(k.c.data)) 
#model <- stan_lmer(k.rich ~ (1|(X+Y+Z)/block/fractal/minor), data=env(k.c.data)) # trying without 'major'
#model <- stan_lmer(k.rich ~ (1|(X+Y+Z)/block/major/minor), data=env(k.c.data)) # trying without 'fractal'
#model <- stan_lmer(k.rich ~ cbind(X+Y+Z) + (1|block/fractal/major/minor), data=env(k.c.data)) # with cbind

coefs <- as.data.frame(summary(model))
variances <- coefs[c("sigma","Sigma[minor:(major:(fractal:block)):(Intercept),(Intercept)]", "Sigma[major:(fractal:block):(Intercept),(Intercept)]","Sigma[fractal:block:(Intercept),(Intercept)]","Sigma[block:(Intercept),(Intercept)]"),"50%"]
names(variances) <- c("residual","minor","major","fractal","block")
variances[c("block","fractal","major","minor","residual")]

barplot(variances/sum(variances), ylim=c(0,1), ylab = "Proportion of variance", names.arg = barnamesknepp, col = "deeppink", density = 30, cex.axis = 1.5, cex.names=1.5, cex.lab = 2)
variances/sum(variances)


#pdf(file="kplot_SR.pdf")
#par(mar=c(5.1,5.5,4.1,2.1))
k.rich <- .ses.mntd(k.c.data)$ntaxa


#################################################
# WILL WAS ERE

# Wrapper function
calc.vars <- function(response){
    model <- stan_lmer(response ~ (1|block/fractal/major/minor), data=env(k.c.data)) # lose fractal / major - because captured in new variables, X, Y, Z
    #iter=10000, warmup=8000, adapt_delta=.999
    coefs <- as.data.frame(summary(model))
    variances <- coefs[c("sigma","Sigma[minor:(major:(fractal:block)):(Intercept),(Intercept)]", "Sigma[major:(fractal:block):(Intercept),(Intercept)]","Sigma[fractal:block:(Intercept),(Intercept)]","Sigma[block:(Intercept),(Intercept)]"),"50%"]
    names(variances) <- c("residual","minor","major","fractal","block")
    variances[c("block","fractal","major","minor","residual")]
    variances <- variances/sum(variances)
    return(variances)
}
k.rich.vars <- calc.vars(k.rich)


null.k.rich.vars <- matrix(NA, nrow=999, ncol=5)
for(i in 1:999)
    null.k.rich.vars[i,] <- calc.vars(sample(k.rich))
rank(c(k.rich.vars[1], null.k.rich.vars[,1]))[1]

hist(null.k.rich.vars[,1])
abline(v=k.rich.vars[1], col="red")


#################################################


barplot(variances/sum(variances), ylim=c(0,1), ylab = "Proportion of variance", names.arg = barnamesknepp, col = "deeppink", density = 30, cex.axis = 1.5, cex.names=1.5, cex.lab = 2)
variances/sum(variances)
#dev.off()

#Uncomment below to see stats on this metric:
#mean(k.rich)
#median(k.rich)
#Mode(k.rich)

#pdf(file="bplot_SR.pdf")
par(mar=c(5.1,5.5,4.1,2.1))
b.rich <- .ses.mntd(b.c.data)$ntaxa

model <- stan_lmer(b.rich ~ (1|fractal/major/minor), data=env(b.c.data))
#iter=10000, warmup=8000, adapt_delta=.999
coefs <- as.data.frame(summary(model))
variances <- coefs[c("sigma","Sigma[minor:(major:fractal):(Intercept),(Intercept)]","Sigma[major:fractal:(Intercept),(Intercept)]","Sigma[fractal:(Intercept),(Intercept)]"),"50%"]
names(variances) <- c("residual","minor","major","fractal")
variances <- variances[c("fractal","major","minor","residual")]

barplot(variances/sum(variances), ylim=c(0,1), ylab = "Proportion of variance", names.arg = barnamesboothby, col = "darkturquoise", density = 30, cex.axis = 1.5, cex.names=1.5, cex.lab = 2)
variances/sum(variances) 
#dev.off()


#mean(b.rich)
#median(b.rich)
#Mode(b.rich)

# ----------------------------------------------------------------------------#
#### Faith's PD ####
# ----------------------------------------------------------------------------#

#pdf(file="kplot_PD.pdf")
par(mar=c(5.1,5.5,4.1,2.1))
k.pd <- .pd(k.c.data)[,"pd"]
model <- stan_lmer(k.pd ~ (1|block/fractal/major/minor), data=env(k.c.data),
                   iter=10000, warmup=8000, adapt_delta=.999)

coefs <- as.data.frame(summary(model))
variances <- coefs[c("sigma","Sigma[minor:(major:(fractal:block)):(Intercept),(Intercept)]", "Sigma[major:(fractal:block):(Intercept),(Intercept)]","Sigma[fractal:block:(Intercept),(Intercept)]","Sigma[block:(Intercept),(Intercept)]"),"50%"]
names(variances) <- c("residual","minor","major","fractal","block")
variances <- variances[c("block","fractal","major","minor","residual")]

barplot(variances/sum(variances), ylim=c(0,1), ylab = "Proportion of variance", names.arg = barnamesknepp, col = "deeppink", density = 30, cex.axis = 1.5, cex.names=1.5, cex.lab = 2)
variances/sum(variances) 
#dev.off()

mean(k.pd)
#median(k.pd)
#Mode(k.pd)

#pdf(file="bplot_PD.pdf")
par(mar=c(5.1,5.5,4.1,2.1))
b.pd <- .pd(b.c.data)[,"pd"]
model <- stan_lmer(b.pd ~ (1|fractal/major/minor), data=env(b.c.data),
                   iter=10000, warmup=8000, adapt_delta=.999)
?.pd

coefs <- as.data.frame(summary(model))
variances <- coefs[c("sigma","Sigma[minor:(major:fractal):(Intercept),(Intercept)]","Sigma[major:fractal:(Intercept),(Intercept)]","Sigma[fractal:(Intercept),(Intercept)]"),"50%"]
names(variances) <- c("residual","minor","major","fractal")
variances <- variances[c("fractal","major","minor","residual")]

barplot(variances/sum(variances), ylim=c(0,1), ylab = "Proportion of variance", names.arg = barnamesboothby, col = "darkturquoise", density = 30, cex.axis = 1.5, cex.names=1.5, cex.lab = 2)
variances/sum(variances) 
#dev.off()

mean(b.pd)
#median(b.pd)
#Mode(b.pd)

# Uncomment below for correlation tests for SR and PD
# KNEPP
# Shapiro-Wilk normality test for SR
#shapiro.test(k.rich) # => p = 0.6784
# Shapiro-Wilk normality test for PD
#shapiro.test(k.pd) # => p = 0.1855
#hist(k.rich)
#hist(k.pd)

#ggqqplot(k.rich, ylab = "k.rich")
#ggqqplot(k.pd, ylab = "k.pd")

#cor.test(k.rich, k.pd, method = "pearson")
# cor 0.8591997
# t = 9.1979, df = 30, p-value = 3.09e-10

# BOOTHBY = not normal data distribution, both heavily right-skewed

# Shapiro-Wilk normality test for SR
#shapiro.test(b.rich) # => p = 0.000...
# Shapiro-Wilk normality test for PD
#shapiro.test(b.pd) # => p = 0.000....
#hist(b.rich)
#hist(b.pd)

#ggqqplot(b.rich, ylab = "b.rich")
#ggqqplot(b.pd, ylab = "b.pd")

# if the data are not normally distributed, it’s recommended to use the non-parametric correlation, 
# including Spearman and Kendall rank-based correlation tests.
#cor.test(b.rich, b.pd, method = "spearman")
# correlation coefficient = rho 0.8643651 
# S = 11144, p-value < 2.2e-16

# ----------------------------------------------------------------------------#
#### Simpson's ####
# ----------------------------------------------------------------------------#

#pdf(file="kplot_Simp.pdf")
par(mar=c(5.1,5.5,4.1,2.1))
k.simpson <- diversity(k.c.data$comm, "simpson") 

model <- stan_lmer(k.simpson ~ (1|block/fractal/major/minor), data=env(k.c.data),
                   iter=10000, warmup=8000, adapt_delta=.999)

coefs <- as.data.frame(summary(model))
variances <- coefs[c("sigma","Sigma[minor:(major:(fractal:block)):(Intercept),(Intercept)]", "Sigma[major:(fractal:block):(Intercept),(Intercept)]","Sigma[fractal:block:(Intercept),(Intercept)]","Sigma[block:(Intercept),(Intercept)]"),"50%"]
names(variances) <- c("residual","minor","major","fractal","block")
variances <- variances[c("block","fractal","major","minor","residual")]

barplot(variances/sum(variances), ylim=c(0,1), ylab = "Proportion of variance", names.arg = barnamesknepp, col = "deeppink", density = 30, cex.axis = 1.5, cex.names=1.5, cex.lab = 2)
variances/sum(variances) 
#dev.off()

#mean(k.simpson)
#median(k.simpson)

#pdf(file="bplot_Simp.pdf")
par(mar=c(5.1,5.5,4.1,2.1))
b.simpson <- diversity(b.c.data$comm, "simpson") 

b.simpson
model <- stan_lmer(b.simpson ~ (1|fractal/major/minor), data=env(b.c.data),
                   iter=10000, warmup=8000, adapt_delta=.999)

coefs <- as.data.frame(summary(model))
variances <- coefs[c("sigma","Sigma[minor:(major:fractal):(Intercept),(Intercept)]","Sigma[major:fractal:(Intercept),(Intercept)]","Sigma[fractal:(Intercept),(Intercept)]"),"50%"]
names(variances) <- c("residual","minor","major","fractal")
variances <- variances[c("fractal","major","minor","residual")]

barplot(variances/sum(variances), ylim=c(0,1), ylab = "Proportion of variance", names.arg = barnamesboothby, col = "darkturquoise", density = 30, cex.axis = 1.5, cex.names=1.5, cex.lab = 2)
variances/sum(variances)  
#dev.off()

#mean(b.simpson)
#median(b.simpson)
#Mode(b.simpson)

# ----------------------------------------------------------------------------#
#### SES MPD ####
# ----------------------------------------------------------------------------#

#pdf(file="kplot_SESmpd.pdf")
par(mar=c(5.1,5.5,4.1,2.1))
k.ses.mpd <- .ses.mpd(k.c.data)$mpd.obs.z
model <- stan_lmer(k.ses.mpd ~ (1|block/fractal/major/minor), data=env(k.c.data),
                   iter=10000, warmup=8000, adapt_delta=.999)

coefs <- as.data.frame(summary(model))
variances <- coefs[c("sigma","Sigma[minor:(major:(fractal:block)):(Intercept),(Intercept)]", "Sigma[major:(fractal:block):(Intercept),(Intercept)]","Sigma[fractal:block:(Intercept),(Intercept)]","Sigma[block:(Intercept),(Intercept)]"),"50%"]
names(variances) <- c("residual","minor","major","fractal","block")
variances <- variances[c("block","fractal","major","minor","residual")]

barplot(variances/sum(variances), ylim=c(0,1), ylab = "Proportion of variance", names.arg = barnamesknepp, col = "deeppink", density = 30, cex.axis = 1.5, cex.names=1.5, cex.lab = 2)
variances/sum(variances) 
#dev.off()

k.ses.mpd <- na.omit(k.ses.mpd)
mean(k.ses.mpd)
median(k.ses.mpd)
Mode(k.ses.mpd)

#pdf(file="bplot_SESmpd.pdf")
par(mar=c(5.1,5.5,4.1,2.1))
b.ses.mpd <- .ses.mpd(b.c.data)$mpd.obs.z
model <- stan_lmer(b.ses.mpd ~ (1|fractal/major/minor), data=env(b.c.data),
                   iter=10000, warmup=8000, adapt_delta=.999)

coefs <- as.data.frame(summary(model))
variances <- coefs[c("sigma","Sigma[minor:(major:fractal):(Intercept),(Intercept)]","Sigma[major:fractal:(Intercept),(Intercept)]","Sigma[fractal:(Intercept),(Intercept)]"),"50%"]
names(variances) <- c("residual","minor","major","fractal")
variances <- variances[c("fractal","major","minor","residual")]

barplot(variances/sum(variances), ylim=c(0,1), ylab = "Proportion of variance", names.arg = barnamesboothby, col = "darkturquoise", density = 30, cex.axis = 1.5, cex.names=1.5, cex.lab = 2)
variances/sum(variances) 
#dev.off()

b.ses.mpd <- na.omit(b.ses.mpd)
mean(b.ses.mpd)
median(b.ses.mpd)
Mode(b.ses.mpd)
# ----------------------------------------------------------------------------#
#### SESmntd ####
# ----------------------------------------------------------------------------#

#pdf(file="kplot_SESmntd.pdf")
par(mar=c(5.1,5.5,4.1,2.1))
k.ses.mntd <- .ses.mntd(k.c.data)$mntd.obs.z
model <- stan_lmer(k.ses.mntd ~ (1|block/fractal/major/minor), data=env(k.c.data),
                   iter=10000, warmup=8000, adapt_delta=.999)

coefs <- as.data.frame(summary(model))
variances <- coefs[c("sigma","Sigma[minor:(major:(fractal:block)):(Intercept),(Intercept)]", "Sigma[major:(fractal:block):(Intercept),(Intercept)]","Sigma[fractal:block:(Intercept),(Intercept)]","Sigma[block:(Intercept),(Intercept)]"),"50%"]
names(variances) <- c("residual","minor","major","fractal","block")
variances <- variances[c("block","fractal","major","minor","residual")]

barplot(variances/sum(variances), ylim=c(0,1), ylab = "Proportion of variance", names.arg = barnamesknepp, col = "deeppink", density = 30, cex.axis = 1.5, cex.names=1.5, cex.lab = 2)
variances/sum(variances) 
#dev.off()

k.ses.mntd <- na.omit(k.ses.mntd)
mean(k.ses.mntd)
median(k.ses.mntd)
Mode(k.ses.mntd)

#pdf(file="bplot_SESmntd.pdf")
par(mar=c(5.1,5.5,4.1,2.1))
b.ses.mntd <- .ses.mntd(b.c.data)$mntd.obs.z
model <- stan_lmer(b.ses.mntd ~ (1|fractal/major/minor), data=env(b.c.data),
                   iter=10000, warmup=8000, adapt_delta=.999)

coefs <- as.data.frame(summary(model))
variances <- coefs[c("sigma","Sigma[minor:(major:fractal):(Intercept),(Intercept)]","Sigma[major:fractal:(Intercept),(Intercept)]","Sigma[fractal:(Intercept),(Intercept)]"),"50%"]
names(variances) <- c("residual","minor","major","fractal")
variances <- variances[c("fractal","major","minor","residual")]

barplot(variances/sum(variances), ylim=c(0,1), ylab = "Proportion of variance", names.arg = barnamesboothby, col = "darkturquoise", density = 30, cex.axis = 1.5, cex.names=1.5, cex.lab = 2)
variances/sum(variances) 
#dev.off()

b.ses.mntd <- na.omit(b.ses.mntd)
mean(b.ses.mntd)
median(b.ses.mntd)
Mode(b.ses.mntd)

# Load data
knepp <- read.csv("clean-data/knepp_clean.csv")
boothby <- read.csv("clean-data/boothby_clean.csv")


# -----------------------------------------------------------------------------#
#### BETA DIVERSITY - TOTAL DIVERSITY ####
# -----------------------------------------------------------------------------#
#KNEPP
betaknepp <- with(knepp, tapply(Cover, list(paste(Block,Group,Plot, sep=""), Species), mean, rm.na=TRUE))
betaknepp[is.na(betaknepp)] <- 0
rownames(betaknepp) <- gsub("0", "", rownames(betaknepp), fixed = TRUE)
betaknepp <- as.data.frame(betaknepp)
saveRDS(betaknepp, "clean-data/betaknepp.RDS")

k.map <- read.csv("raw-data/k-map.csv")
k.map <- subset(k.map, k.map$uniqueID != "")
row.names(k.map) <- k.map$uniqueID
saveRDS(k.map, "clean-data/k.map.RDS")

# BOOTHBY
betabooth <- with(boothby, tapply(Cover, list(paste(Block,Plot, sep=""), Species), mean, rm.na=TRUE))
betabooth[is.na(betabooth)] <- 0
betabooth <- as.data.frame(betabooth)
saveRDS(betabooth, "clean-data/betabooth.RDS")

b.map <- read.csv("raw-data/b-map.csv")
b.map <- subset(b.map, b.map$name != "")
row.names(b.map) <- b.map$name
b.map <- subset(b.map, select = -c(layer,path,ID_SHAPE, ID_PART, ID_POINT, CLOCKWISE, LAKE, ID, fid, description))
saveRDS(b.map, "clean-data/b.map.RDS")

# KNEPP
betaknepp <- readRDS("clean-data/betaknepp.RDS")
k.map <- readRDS("clean-data/k.map.RDS")
# BOOTHBY
betabooth <- readRDS("clean-data/betabooth.RDS")
b.map <- readRDS("clean-data/b.map.RDS")

# Checking the only site differences between the map and the plant dataset are the * ones
setdiff(rownames(k.map), rownames(betaknepp))
setdiff(rownames(betaknepp), rownames(k.map))
betaknepp <- betaknepp[rownames(betaknepp) %in% rownames(k.map),]
k.map <- k.map[rownames(k.map) %in% rownames(betaknepp),]
betaknepp <- betaknepp[match(rownames(k.map), rownames(betaknepp)),]
identical(rownames(betaknepp), rownames(k.map))

saveRDS(k.map, "clean-data/k.map.RDS")
saveRDS(betaknepp, "clean-data/betaknepp.RDS")

setdiff(rownames(b.map), rownames(betabooth))
setdiff(rownames(betabooth), rownames(b.map))
betabooth <- betabooth[rownames(betabooth) %in% rownames(b.map),]
b.map <- b.map[rownames(b.map) %in% rownames(betabooth),]
betabooth <- betabooth[match(rownames(b.map), rownames(betabooth)),]
identical(rownames(betabooth), rownames(b.map))

saveRDS(b.map, "clean-data/b.map.RDS")
saveRDS(betabooth, "clean-data/betabooth.RDS")

# Getting BDTotal

# 1: convert to a matrix
k.comm <- as.matrix(betaknepp)
b.comm <- as.matrix(betabooth)
# 2: Turn the NAs into 0s (absences)
k.comm[is.na(k.comm)] <- 0
b.comm[is.na(b.comm)] <- 0

beta.div(k.comm, method = "sorensen", samp=FALSE, save.D=TRUE)
beta.div(b.comm, method = "sorensen", samp=FALSE, save.D=TRUE)

#### Plot of species: ####
# KNEPP
#kneppcount<- knepp%>%select(Site, Plot, Species)%>%
#  unique()%>%group_by(Species)%>%count()

#Uncomment pdf() and dev.off() lines to save plots to file
#pdf(file="knepp_species.pdf")
#ggplot(kneppcount, aes(x = reorder(Species,n, decreasing=TRUE), y = n,fill=n)) + 
#  geom_bar(stat="identity")+
#  xlab("Species") + 
#  ylab("Count")+ 
#  theme_classic()+
#  theme(axis.text.x = element_blank(), legend.position = "none")+
#  scale_fill_gradient(low="black", high="#FFB4DC")+
#  theme(text = element_text(size = 25))
#dev.off()

# BOOTHBY
#boothbycount<- boothby%>%select(Site, Plot, Species)%>%
#  unique()%>%group_by(Species)%>%count()

#pdf(file="bplot_species.pdf")
#ggplot(boothbycount, aes(x = reorder(Species,n, decreasing=TRUE), y = n,fill=n)) + 
#  geom_bar(stat="identity")+
#  xlab("Species") + 
#  ylab("Count")+ 
#  theme_classic()+
#  theme(axis.text.x = element_blank(), legend.position = "none")+
#  scale_fill_gradient(low="black", high="#A1ECEE")+
#  theme(text = element_text(size = 25))
#dev.off()
