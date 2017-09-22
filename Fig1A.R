# Script for figure 1A: sampling map

#### Panel A -- Sampling map ####

#### Data required
# The database
setwd("/Users/afr/Desktop/A/Phd_papers/QGD_submission/QGD_data/QGD_database")
QGD_DB <- read.delim("QGD_DB_clean.txt", sep="\t", h=T, stringsAsFactors = F)

#### Auxiliary variables
QGD_colors <- c("#CDAA7D","#7EC0EE")

#### Plot map (function)
# This function uses:
        # QGD_DB: the database
        # colp: color for Pleistocene records
        # colh: color for Holocene records
QGD.F1A.map <- function(QGD_DB, colp=QGD_colors[1], colh=QGD_colors[2]){
        library(rworldmap)
        library(mapproj)
        # remove the records without coordinates
        xQGD_DB <- QGD_DB[-which(is.na(QGD_DB$Longitude)),]
        # Create a table with factors for time bin and type of sequence
        xmap <- as.data.frame(matrix(numeric(), nrow = dim(xQGD_DB)[1], ncol = 5))
        colnames(xmap) <- c("Species", "Longitude", "Latitude", "Bin", "Type")
        xmap[, 1] <- xQGD_DB$Species
        xmap[, 2] <- as.numeric(xQGD_DB$Longitude)
        xmap[, 3] <- as.numeric(xQGD_DB$Latitude)
        xmap$Bin <- ifelse(xQGD_DB$Mean_Age >= 11700, colp, colh)
        # pch 24 (triangle) == DNA; pch 22 (squares) == Not DNA
        xmap$Type <- ifelse(!is.na(xQGD_DB$Sequence), 24, 22)
        # Plot
        map("world",proj="azequidist", resolution = 0.2, ylim = c(33,90),
            lwd=1.2, mar=c(0,2,2,2), interior = F)
        points(mapproject(xmap$Longitude, xmap$Latitude), bg=xmap$Bin, col="black",
               pch=xmap$Type, cex=0.7, lwd=0.7)
        # Legend
        legend(-1.25,-0.45, pch=c(24, 22), c("DNA","Not DNA"), bty="n", pt.cex = 1,
               cex=0.7, pt.lwd = 1.5, title="Record type", title.adj = c(0,0),
               y.intersp = 0.7, x.intersp = .9)
        legend(-1.25,-0.82, pch=15, c("Pleistocene", "Holocene"), bty = "n",
               col=c(colp, colh), pt.cex = 1, cex=0.7, pt.lwd = 3,title = "Time bin",
               title.adj = c(0,0), y.intersp = 0.7, x.intersp = .9)
}
#### To execute the script ####
pdf("Fig1_map.pdf", width = 4,height = 4)
QGD.F1A.map(QGD_DB)
dev.off()
