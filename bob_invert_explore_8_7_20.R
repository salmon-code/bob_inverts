library(tidyverse)
library(vegan)
library(cowplot)
library(goeveg)

#you'll need to set the working directory for your own computer
setwd("/Users/StuartMunsch/Google Drive/Other/bob_inverts/")

#I've found it helpful to use the convention where all names are lowercase and separated by a _ for spaces
#...(then you don't have to remember each time you refer to a variable how it was capitalized)
inverts_sample_info <- read.csv("SW_18_Epi_groups_noXP.csv")
names(inverts_sample_info) <- c("label", "year", "date", "site", "strata", "surface", "strata_surface", "sample_number", "color")

inverts_sample_info$site_date <- paste(inverts_sample_info$site, inverts_sample_info$date, sep = "_")

inverts_sample_info$ob <- ifelse(inverts_sample_info$strata_surface == "OB", 1, 0)
inverts_sample_info$os <- ifelse(inverts_sample_info$strata_surface == "OS", 1, 0)
inverts_sample_info$ub <- ifelse(inverts_sample_info$strata_surface == "UB", 1, 0)
inverts_sample_info$us <- ifelse(inverts_sample_info$strata_surface == "US", 1, 0)

inverts_wide <- read.csv("SW_18_Epi_speciesabu_only_noXP_log.csv")

#first question is whether there are differences among 
  #1. outside bench, 2. outside seawall, 3. under bench, 4. under seawall

#appreciate sampling regime
table(inverts_sample_info$year) #all 2018

table(inverts_sample_info$date, 
      inverts_sample_info$strata_surface, 
      inverts_sample_info$site) #perfectly balanced among strata, sites, and dates

table(inverts_sample_info$site_date) #perfectly balanced among strata, sites, and dates

#step one:
  #get a common sense appreciation of factors driving patterns in your data

#remove taxa not present in at least 5% of samples

inverts_5 <- inverts_wide[,colSums(inverts_wide > 0) >= 16]
  #this says exclude rows with taxa that are in less than 16 samples. 16 is 5% of 320 (your total sample number)

inverts_pca <- prcomp(inverts_5)

inverts_pca_tibble <- tibble(pc1 = inverts_pca$x[,"PC1"],
                   pc2 = inverts_pca$x[,"PC2"],
                   site = inverts_sample_info$site,
                   date = inverts_sample_info$date,
                   strata_surface = inverts_sample_info$strata_surface,
                   strata = inverts_sample_info$strata,
                   surface = inverts_sample_info$surface)

ggplot(aes(x = pc1, y = pc2, col = strata_surface, shape = site), 
       data = inverts_pca_tibble) +
  geom_point()

plot_grid(
  ggplot(aes(x = pc1, y = pc2, col = site), 
         data = inverts_pca_tibble) +
    geom_point(),
  ggplot(aes(x = pc1, y = pc2, col = date), 
         data = inverts_pca_tibble) +
    geom_point(),
  ggplot(aes(x = pc1, y = pc2, col = strata), 
         data = inverts_pca_tibble) +
    geom_point(),
  ggplot(aes(x = pc1, y = pc2, col = surface), 
         data = inverts_pca_tibble) +
    geom_point(),
  ggplot(aes(x = pc1, y = pc2, col = strata_surface), 
         data = inverts_pca_tibble) +
    geom_point(), align = "vh", axis = "lr"
)

ggsave("bob inverts pca not prop.jpeg", height = 6, width = 12)

dim(inverts_5)
dim(inverts_10)


#excluding taxa based on increasingly strict requirements of presence in samples 
#(5, 10, 15 etc. percent)


inverts_30 <- inverts_wide[,colSums(inverts_wide > 0) >= 96]

which(rowSums(inverts_30) == 0)
inverts_30_del <- inverts_30[-58,]




#no covergence
inverts_nmds_5 <- metaMDS(inverts_5, k = 2, distance = "bray", autotransform = FALSE, trymax = 10000, previous.best = TRUE)

#no covergence
inverts_nmds_10 <- metaMDS(inverts_10, k = 2, distance = "bray", autotransform = FALSE, trymax = 1000, previous.best = TRUE)

#covergence in 772 tries
inverts_nmds_15 <- metaMDS(inverts_15, k = 2, distance = "bray", autotransform = FALSE, trymax = 1000, previous.best = TRUE)

#good, robust patterns
inverts_nmds_tibble <- tibble(nmds1 = inverts_nmds_15$points[,1], 
                       nmds2 =  inverts_nmds_15$points[,2],
                       site = inverts_sample_info$site,
                       date = inverts_sample_info$date,
                       strata_surface = inverts_sample_info$strata_surface,
                       strata = inverts_sample_info$strata,
                       surface = inverts_sample_info$surface)

plot_grid(
  ggplot(aes(x = nmds1, y = nmds2, col = site), 
         data = inverts_nmds_tibble) +
    geom_point(),
  ggplot(aes(x = nmds1, y = nmds2, col = date), 
         data = inverts_nmds_tibble) +
    geom_point(),
  ggplot(aes(x = nmds1, y = nmds2, col = strata), 
         data = inverts_nmds_tibble) +
    geom_point(),
  ggplot(aes(x = nmds1, y = nmds2, col = surface), 
         data = inverts_nmds_tibble) +
    geom_point(),
  ggplot(aes(x = nmds1, y = nmds2, col = strata_surface), 
         data = inverts_nmds_tibble) +
    geom_point(), align = "vh", axis = "lr"
)


install.packages("ARTool")
library(ARTool)
art(inverts_5)



## Bray-Curtis distances between samples
dis <- vegdist(inverts_5)

## Calculate multivariate dispersions
mod <- betadisper(dis, inverts_sample_info$strata_surface)
mod
anova(mod) #so you shouldn't run a permanova...




#run a permanova

#even when you tell it the model to account for four factors, the output is only three
#the issue is that it's comparing by factors by contrasts among them
#you need three contrasts to tell the difference between four levels of a factor
#can think of the model output sort of like ob vs. us, os vs. us, and ub vs. us where us is the baseline value.

#this permanova stratifies by a block of date x time, which should be appropriate in your study because it is perfectly balanced.
  #(at least this part, I think your other question might be dealing w/imbalanced data)
the_perm <- how(nperm = 999)
setBlocks(perm) <- with(inverts_sample_info, site_date)

the_permanova_blocked <- adonis2(inverts_5 ~ ob + os + ub + us, 
                          data = inverts_sample_info, 
                          permutations = the_perm)

#checking what happens when not stratified
the_permanova_not_blocked <- adonis(inverts_5 ~ ob + os + ub + us, 
                                     data = inverts_sample_info)

#does seem odd that nothing changes, but I think what's happening is that both outputs are quite significant, so you get the p value of 0.001 regardless of blocking

#here's the example that vegan provides. the output is exactly the same between blocked and unblocked permanovas except the pvalues, suggesting the output on our models is normal

### Example of use with strata, for nested (e.g., block) designs.
dat <- expand.grid(rep=gl(2,1), NO3=factor(c(0,10)),field=gl(3,1) )
dat
Agropyron <- with(dat, as.numeric(field) + as.numeric(NO3)+2) +rnorm(12)/2
Schizachyrium <- with(dat, as.numeric(field) - as.numeric(NO3)+2) +rnorm(12)/2
total <- Agropyron + Schizachyrium
dotplot(total ~ NO3, dat, jitter.x=TRUE, groups=field,
        type=c('p','a'), xlab="NO3", auto.key=list(columns=3, lines=TRUE) )

Y <- data.frame(Agropyron, Schizachyrium)
mod <- metaMDS(Y, trace = FALSE)
plot(mod)
### Ellipsoid hulls show treatment
with(dat, ordiellipse(mod, field, kind = "ehull", label = TRUE))
### Spider shows fields
with(dat, ordispider(mod, field, lty=3, col="red"))

### Incorrect (no strata)
perm <- how(nperm = 199)
adonis2 (Y ~ NO3, data = dat, permutations = perm)

## Correct with strata
setBlocks(perm) <- with(dat, field)
adonis2(Y ~ NO3, data = dat, permutations = perm)


