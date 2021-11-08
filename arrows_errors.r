



par(mfrow=c(1,1), mar=c(5,5,3,3))
year <- 2000:2020
setwd('/home/stevie/Change_detection/KIUIC/tifs/')
kiuic_rasters <- read.csv('total_change_kiuic.csv')
plot(year, kiuic_rasters$V1, type="l", col = "blue",lwd =3, ylab = "Forest cover loss (ha)", ylim =c(0,50000), xlab="Year", cex.lab=1.5, cex.axis = 1.2) # kiuic
setwd('/home/stevie/Change_detection/PALMAR/tifs/')

palmar_rasters <- read.csv('total_change_palmar.csv')

lines(year, palmar_rasters$V1, type="l", col = "red", lwd=3) # palmar

setwd('/home/stevie/Change_detection/FCP/tifs/')
fcp_rasters <- read.csv('total_change_fcp.csv')
lines(year, rasters, type="l", col = "green", lwd = 3) # fcp
legend("topleft", legend=c("Semi-deciduous", "Semi-evergreen", "Deciduous"), col=c("blue", "green", "red"), lty=1, lwd = 3)



setwd('/home/stevie/Change_detection/')

kiuic_errors2006 <- read.csv('KIUIC/2006_errors')
kiuic_errors2015 <- read.csv('KIUIC/2015_errors')


setwd('/home/stevie/Change_detection/')
palmar_errors2008 <- read.csv('PALMAR/2008_errors')
palmar_errors2015 <- read.csv('PALMAR/2015_errors')

setwd('/home/stevie/Change_detection/')
fcp_errors2010 <- read.csv('FCP/2010_errors')
fcp_errors2015 <- read.csv('FCP/2015_errors')


points(2006,kiuic_errors2006$defo_area, col = "blue", pch=20, cex=1)
arrows(2006,kiuic_errors2006$upperic,2006,kiuic_errors2006$loweric, angle = 90, code = 3, length = 0.1, cex=2, col = "blue")
points(2015,kiuic_errors2015$defo_area, col = "blue", pch=20, cex=1)
arrows(2015,kiuic_errors2015$loweric,2015,kiuic_errors2015$upperic, angle = 90, code = 3, length = 0.1, cex=2, col = "blue")



points(2008,palmar_errors2008$defo_area, col = "red", pch=20, cex=1)
arrows(2008,palmar_errors2008$loweric,2008,palmar_errors2008$upperic, angle = 90, code = 3, length = 0.1, cex=2, col = "red")

points(2015,palmar_errors2015$defo_area, col = "red", pch=20, cex=1)
arrows(2015,palmar_errors2015$loweric,2015,palmar_errors2015$upperic, angle = 90, code = 3, length = 0.1, cex=2, col = "red")



points(2010,fcp_errors2010$defo_area, col = "green", pch=20, cex=1)
arrows(2010,fcp_errors2010$loweric,2010,fcp_errors2010$upperic, angle = 90, code = 3, length = 0.1, cex=2, col = "green")

points(2015,fcp_errors2015$defo_area, col = "green", pch=20, cex=1)
arrows(2015,fcp_errors2015$loweric,2015,fcp_errors2015$upperic, angle = 90, code = 3, length = 0.1, cex=2, col = "green")

