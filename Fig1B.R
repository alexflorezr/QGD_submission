# Script for figure 1B: climate change velocity per species
#### Panel B -- Climate change velocity per species ####

#### packages required
library(raster)

#### Data required
# The database
setwd("/Users/afr/Desktop/A/Phd_papers/QGD_submission/QGD_data/QGD_database")
QGD_DB <- read.delim("QGD_DB_clean.txt", sep="\t", h=T, stringsAsFactors = F)
# Remove marine species
bm <- which(QGD_DB$Species == "Balaena_mysticetus")
er <- which(QGD_DB$Species == "Eschrichtius_robustus")
QGD_DB_F1B <- QGD_DB[-c(bm, er),]
# The climate change velocity maps
setwd("/Users/afr/Desktop/A/Phd_papers/QGD_submission/QGD_data/QGD_climate/QGD_velocity_raw/")
CC_vel <- stack("climate_change_velocity_perYear_tmp_rescaled_truncated_0_000005.grd")

#### Auxiliary variables

#### Plot climate change velocity per species heat map ####
#### internal functions ####
# QGD.F1B.bin.layer gets bins and the number of the raster layers for each record
# This function uses:
        # QGD_DB: the database
        # layer(boolean): TRUE gets stack the layer number for each time period
QGD.F1B.bin.layer <- function(QGD_DB, layer=FALSE){
        xvctr <- vector(length = length(QGD_DB))
        # vector of time bins
        x50to0 <- c(seq(50000, 22000, by=-2000), seq(21000, 0, by=-1000))
        # vector of raster stack layers corresponding to the time bins vector
        xlayers <- seq(26, 61, by=1)
        for (bin in seq_along(x50to0)[-37]){
                # which records fall within an specific time bin
                xr <- which(QGD_DB < x50to0[bin] & QGD_DB >= x50to0[bin+1])
                if(layer == FALSE){
                        xvctr[xr] <- x50to0[bin+1]
                }else{
                        xvctr[xr] <- xlayers[bin]  
                }
        }
        xvctr
}

# QGD.F1B.extract extracts CC value for each record
# This function uses:
        # QGD_DB: the database
        # CC_vel: stack of climate change velocity maps
        # QGD.F1B.extract uses QGD.F1B.bin.layer
QGD.F1B.extract <- function(QGD_DB, CC_vel){
        xcn <- c("Species","Longitude", "Latitude", "Mean_date", "Time_bin",
                 "Layer", "Rec_type", "Cell", "Cli_vel_tmp")
        xpnts <- as.data.frame(matrix(nrow=nrow(QGD_DB), ncol=length(xcn)))
        colnames(xpnts) <- xcn
        xpnts$Species   <- QGD_DB$Species
        xpnts$Longitude <- as.numeric(QGD_DB$Longitude)
        xpnts$Latitude  <- as.numeric(QGD_DB$Latitude)
        xpnts$Mean_date <- as.numeric(QGD_DB$Mean_Age)
        xpnts$Rec_type  <- ifelse(nchar(QGD_DB$Sequence) > 1, "Seq","Fossil")
        # QGD.F1B.bin.layer extracts the time bin
        xpnts$Time_bin  <- QGD.F1B.bin.layer(xpnts$Mean_date)
        # QGD.F1B.bin.layer extracts the raster stack layer
        xpnts$Layer     <- QGD.F1B.bin.layer(xpnts$Mean_date, layer=TRUE)
        xcoor           <- cbind(xpnts$Longitude, xpnts$Latitude)
        xextract        <- extract(CC_vel, layer=xpnts$Layer, nl=1, y=xcoor, cellnumbers=T)
        xpnts$Cell      <- xextract[,1]
        xpnts$Cli_vel_tmp <- xextract[,2]
        xpnts
}

# QGD.F1B.median estimates median values per time bin per species
# This function uses:
        # CC_vel_sp == the output from QGD.F1B.extract
QGD.F1B.median <- function(CC_vel_sp){
        xsp <- unique(CC_vel_sp$Species)
        x50to0 <- c(seq(50000, 22000, by=-2000), seq(21000, 0, by=-1000))
        xmedian <- matrix(nrow=length(xsp), ncol=length(x50to0))
        rownames(xmedian) <- xsp
        colnames(xmedian) <- rev(x50to0)
        for(sp in seq_along(xsp)){
                xsp_DB    <- subset(CC_vel_sp, CC_vel_sp$Species == xsp[sp])
                # if all the CC vel values are NA then jump to the next sp
                if(sum(is.na(xsp_DB$Cli_vel_tmp)) == length(xsp_DB$Cli_vel_tmp)){
                        next()
                }else{
                        xsp_vctr  <- aggregate(Cli_vel_tmp ~ Time_bin, data=xsp_DB, median, na.rm=T)
                        xm        <- match(xsp_vctr$Time_bin, colnames(xmedian))
                        xmedian[sp,xm] <- xsp_vctr$Cli_vel_tmp
                }
        }
        xmedian
}

# QGD.F1B.quantile.sp transforms the median matrix to quantiles using values per species
# This function uses: 
        # CC_vel_median == the QGD.F1B_median output
        # sp_q == CC vel quantiles for each species 
        # vColor == a vector with 5 colors
QGD.F1B.quantile.sp <- function(CC_vel_median, sp_q, vColor){
        # Vector for the temporal axis and species
        x50to0 <- c(seq(50000, 22000, by=-2000), seq(21000, 0, by=-1000))
        xsp <- rownames(CC_vel_median)
        # the matrix to plot
        ximage <- matrix(nrow=dim(CC_vel_median)[1], ncol=length(x50to0))
        colnames(ximage) <- rev(x50to0)
        rownames(ximage) <- rownames(CC_vel_median)
        # loop through time bins and species
        for(sp in seq_along(xsp)){
                xb <- CC_vel_median[sp,]
                xq <- sp_q[sp,2]
                for(bin in seq_along(x50to0)){
                        xc <- which(colnames(ximage) == x50to0[bin])
                        if(is.na(xb[xc])){
                                ximage[sp, xc] <- vColor[1]
                        }else{
                                if(xb[xc] <= xq[2]){
                                        ximage[sp, xc] <- vColor[2]
                                }
                                if(xb[xc] > xq[2] & xb[xc] <= xq[3]){
                                        ximage[sp, xc] <- vColor[3]
                                }
                                if(xb[xc] > xq[3] & xb[xc] <= xq[4]){
                                        ximage[sp, xc] <- vColor[4]
                                }
                                if(xb[xc] > xq[4] & xb[xc] <= xq[5]){
                                        ximage[sp, xc] <- vColor[5]
                                } 
                        }
                }
        }
        ximage
}

# QGD.F1B.plot plots the image using per species, per time bin quartiles coloring
# This function uses: 
        # CC_vel_quant == output of QGD.F1B.quantile.sp
QGD.F1B.plot <- function(CC_vel_quant){
        #l <- layout(matrix(c(1,2),nrow=1,ncol=2), heights = c(1), width=c(80,20))
        #par(mar=c(2,2,2,0), oma=c(2,2,2,0))
        plot(1,1, type="n", frame.plot=F, axes=F, xlim=c(0,50000), ylim=c(0, 53125), 
             xlab="", ylab="", xaxs="i")
        sp <- nrow(CC_vel_quant):1
        xheight <- seq_along(CC_vel_quant[,1])*3125
        for (r in seq_along(rownames(CC_vel_quant))){
                for(c in 1:(ncol(CC_vel_quant))){
                        x0 <- as.numeric(colnames(CC_vel_quant)[c])
                        x1 <- as.numeric(colnames(CC_vel_quant)[c+1])
                        y0 <- xheight[r]
                        y1 <- xheight[r] + 3125
                        xc <- CC_vel_quant[r,c]
                        polygon(c(x1,x0,x0,x1), c(y1, y1, y0, y0), col = xc,
                                border = "white", lwd=0.2)
                }
                r <- r +1
        }
        axis(1, at=(colnames(CC_vel_quant)), labels =colnames(CC_vel_quant), las=2, pos = 1, cex.axis=0.8)
        xsp <-sapply(lapply(strsplit(rownames(CC_vel_quant), split = "_"), substr, 1, 1),
                     paste0, collapse="_") 
        axis(2, at = xheight+1.25, rev(xsp), las=2)
        mtext(side=1, "Years before present", line=2.5)
        #par(mar=c(4,1,2,2))
        #plot(1,1, type="n", xlim=c(0,1), ylim=c(0,8), frame=F, axes=F, xlab="")
        #lt <- c("NA", "25%", "50%", "75%", "100%")
        #for(lg in 3:7){
        #      lx0 <- 0.4
        #       lx1 <- 0.1
        #       ly0 <- lg
        #        ly1 <- lg + 1
        #        polygon(c(lx0,lx1,lx1,lx0), c(ly0,ly0,ly1,ly1), col=vC[lg-2])
        #        text(lt[lg-2], x = 0.75, y= ly0 +0.5 )
        
        #}
}
#### To execute the script ####
QGD_F1B_table <- QGD.F1B.extract(CC_vel = CC_vel, QGD_DB = QGD_DB_F1B)
QGD_F1B_median <- QGD.F1B.median(QGD_F1B_table)
fig3_col <- c("#FFFFFF","#FFDAB9","#FFA54F","#CD6600","#8B4500")
sp_q <- aggregate(Cli_vel_tmp ~ Species, data=QGD_F1B_table, quantile, na.rm=T)
QGD_F1B_image_sp <- QGD.F1B.quantile.sp(QGD_F1B_median, sp_q, fig3_col)
QGD.F1B.plot(QGD_F1B_image_sp)
