# MORPHOLOGY PCA -------------------------------------------------------------#

# read in phenotype-environment data
d <- read.csv("clean_pca_phenotype-enviro.csv", header = T, stringsAsFactors = F)
head(d)
rownames(d) <- d$island
d3 <- read.csv("clean_phenotype-enviro-morph.csv", header = T, stringsAsFactors = F)
head(d3)

# when you do a PCA, you need a matrix with ONLY quantitative trait data, and you cannot have any NA data within your matrix. I like to subset my data with all my traits of interest in columns and their identifier as the rownames. You can use the function complete.cases() to get rid of any NA data 
complete.cases(d)

# good to go! 
## since a lot of morphological traits can be correlated with body size, we can to 'body size correct' these data. You can do that by regressing each trait against your body size metric (in my case SVL) and store the residuals in their own dataframe or matrix

residuals <- d[-1]
for(i in 2:ncol(d)) {
  lm <- lm(d[,i] ~ d$svl, na.action = na.omit)
  print(summary(lm))
  residuals[, i-1] <- lm$residuals
}
# if you scroll quickly through this output you can see that all of our traits are significantly correlated with body size. yep!
head(residuals)

# now we can use this function to perfrom the PCA
# we can do the PCA on the raw trait data or the residuals. its kind of fun to toggle between them and see how much the PCA changes 

pca <- prcomp(d[, 2:ncol(d) ])
# the 2:ncol(morph2) just gets rid of the annoying X column 

library('factoextra', 'RColorBrewer')
# fviz_pca_ind(pca)
myplot <- fviz_pca_var(pca,
                       col.var = "contrib", # Color by contributions to the PC
                       gradient.cols = RColorBrewer::brewer.pal(9,'RdBu'),
                       repel = TRUE     # Avoid text overlapping
)
myplot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

## now we can really see that SVL explains a ton of the variation here, and which direction corresponds to increasing SVL in our PCA plots/eigenvectors


## now what I'm going to do is plot all of the values for each individual and color them by species. You could color them by morph. The way I do this is just by taking the vector of names/species/morphs, copying it into a new vector, and then replacing each element with the corresponding color. You could probably do it faster but this is what I do!

cols <- c("orange", "white", "yellow")
names(cols) <- c("o", "w", "y" )
PCcols <- paste0(d3$morph[d3$island %in% rownames(d) == TRUE],"_",d3$sex[d3$island %in% rownames(d) == TRUE])

PCcols[PCcols %in% c("o_m","o_f")] <- cols[1]
PCcols[PCcols %in% c("w_m","w_f")] <- cols[2]
PCcols[PCcols %in% c("y_m","y_f")] <- cols[3]
# check it
PCcols

# going to do the same thing but I'll use different point shapes to indicate different sexes 
pch <- paste0(d3$morph[d3$island %in% rownames(d) == TRUE],"_", d3$sex[d3$island %in% rownames(d) == TRUE])
pch[pch %in% c("o_m", "w_m", "y_m")] <- 21
pch[pch %in% c("o_f", "w_f", "y_f")] <- 24

## now we can plot our PCA with the colors and points corresponding to our species 
## PCA plot
plot(pca$x[,1], pca$x[,2], pch = as.numeric(pch), bg = PCcols, bty = "l", xlab = "PC1", ylab = "PC2")
legend("bottomright", legend = c("orange","white","yellow"), fill = cols, bty = "n")
legend("topright", legend = c("males","females"), pch = c(21, 24), bty = "n")

## Now, I like to draw my own convex polygons around species. There are methods to automate this in ggplot etc., but they often draw ellipses that I personally think are a little misleading (ie they eliminate some of the data in favor of showing non-overlapping shapes). This is personal preference but here is a sort of hack-y way to draw shapes around all of your points

col <- 1
for ( i in unique(d3$morph[d3$island %in% rownames(d) == TRUE]) ){
  positions <- which(d3$morph[d3$island %in% rownames(d) == TRUE] == i)	
  tmpDF <- data.frame(pca$x[,1][positions], pca$x[,2][positions])
  ch <- chull(tmpDF)
  coords <- tmpDF[c((ch), ch[1]), ]	
  polygon(coords, col = alpha(cols[which(names(cols) == i)], 0.1), border = cols[which(names(cols) == i)])
  col <- col + 2
}

pdf("naxosPCA.pdf")
plot(pca$x[,1], pca$x[,2], pch = as.numeric(pch), bg = PCcols, bty = "l", xlab = "PC1", ylab = "PC2")
legend("bottomright", legend = c("orange","white","yellow"), fill = cols, bty = "n")
legend("topright", legend = c("males","females"), pch = c(21, 24), bty = "n")
col <- 1
for ( i in unique(d3$morph[d3$island %in% rownames(d) == TRUE]) ){
  positions <- which(d3$morph[d3$island %in% rownames(d) == TRUE] == i)	
  tmpDF <- data.frame(pca$x[,1][positions], pca$x[,2][positions])
  ch <- chull(tmpDF)
  coords <- tmpDF[c((ch), ch[1]), ]	
  polygon(coords, col = alpha(cols[which(names(cols) == i)], 0.1), border = cols[which(names(cols) == i)])
  col <- col + 2
}
dev.off()