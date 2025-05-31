library("adegenet")

myCol <- c("darkblue","purple","green","orange","red","blue")
#T2T data
data_t2t <- read.csv('/Users/luciabazan/Documents/GitHub/cnv_sgdp/data/input_dapc_t2t_hprc.csv')
pop <- data_t2t$Superpopulation
data_t2t <- data_t2t[ -c(1,2) ]
dapc.data_t2t <-  dapc(data_t2t, grp=pop, var.contrib = TRUE, scale = FALSE, n.pca = 6, n.da = 100)

scatter(dapc.data_t2t, pch = 18:23, cstar = 0,  scree.da=FALSE,  leg=TRUE, clab=0, main = 'T2T')
# Centroids of every population
centroids_t2t = dapc.data_t2t[["grp.coord"]]
write.csv(centroids_t2t,"/Users/luciabazan/Documents/GitHub/cnv_sgdp/data/centroids_t2t.csv")
## 
#scatter(dapc.data_t2t, pch = 18:23, cstar = 0, mstree = TRUE, lwd = 2, lty = 2, scree.da=TRUE, posi.da="bottomleft", leg=TRUE, clab=0)

scatter(dapc.data_t2t, posi.da="bottomright",  bg="white",
        pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,
        posi.pca="bottomleft", main='T2T')

## STRUCTURE like plot
assignplot(dapc.data_t2t)

par(mar=c(5.1,4.1,1.1,1.1), xpd=TRUE)
compoplot(dapc.data_t2t, lab="", posi=list(x=12,y=-.01), cleg=.7, space=(0), main='T2T', xlab="individuals")


## Barplot
temp <- summary(dapc(data_t2t, grp=pop, var.contrib = TRUE, scale = FALSE, n.pca = 5, n.da = 6))$assign.per.pop*100
par(mar=c(4.5,7.5,1,1))
barplot(temp, xlab="% of reassignment to actual population",
        horiz=TRUE, las=1, main='T2T')


## Optimal PCA's 
xval <- xvalDapc(data_t2t, pop, n.pca.max = 11, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)

#############GRCh37 data
data_grch37 <- read.csv('/Users/luciabazan/Documents/GitHub/cnv_sgdp/data/input_dapc_grch38_hprc.csv')
pop <- data_grch37$Superpopulation
data_grch37 <- data_grch37[ -c(1,2)]
dapc.data_grch37 <-  dapc(data_grch37, grp=pop, var.contrib = TRUE, scale = FALSE, n.pca = 6, n.da = 100)
#centroids
centroids_grch38 = dapc.data_grch37[["grp.coord"]]
write.csv(centroids_grch38,"/Users/luciabazan/Documents/GitHub/cnv_sgdp/data/centroids_grch38.csv")

scatter(dapc.data_grch37, pch = 18:23, cstar = 0, scree.da=FALSE, leg=TRUE, clab=0)


scatter(dapc.data_grch37, cell = 0, pch = 18:23, cstar = 0, mstree = TRUE, lwd = 2, lty = 2, scree.da=TRUE, posi.da="bottomleft", leg=TRUE, clab=0)


## 
assignplot(dapc.data_grch37)

par(mar=c(5.1,4.1,1.1,1.1), xpd=TRUE)
compoplot(dapc.data_grch37, lab="", posi=list(x=12,y=-.01), cleg=.7, space=(0), main='GRCh38', xlab="individuals")


## BARPLOT
temp <- summary(dapc(data_grch37, grp=pop, var.contrib = TRUE, scale = FALSE, n.pca = 6, n.da = 100))$assign.per.pop*100
par(mar=c(4.5,7.5,1,1))
barplot(temp, xlab="% of reassignment to actual population",
        horiz=TRUE, las=1, main='GRCh38')


## Optimal PCA's 
xval <- xvalDapc(data_grch37, pop, n.pca.max = 11, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)



