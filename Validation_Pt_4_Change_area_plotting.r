### This script plots change quartiles from Changesize output against time
## 17/06/2020 (1/1)

## set working directory
setwd('/home/stevie/Change_detection/FCP/')

#   install.packages("viridis")
## Colour libraries
library(ggplot2)
library(viridis)

# load all change size stat tables
stats_2020 <- read.csv('changeSize_stats 20 .csv')
stats_2019 <- read.csv('changeSize_stats 19 .csv')
stats_2018 <- read.csv('changeSize_stats 18 .csv')
stats_2017 <- read.csv('changeSize_stats 17 .csv')
stats_2016 <- read.csv('changeSize_stats 16 .csv')
stats_2015 <- read.csv('changeSize_stats 15 .csv')
stats_2014 <- read.csv('changeSize_stats 14 .csv')
stats_2013 <- read.csv('changeSize_stats 13 .csv')
stats_2012 <- read.csv('changeSize_stats 12 .csv')
stats_2011 <- read.csv('changeSize_stats 11 .csv')
stats_2010 <- read.csv('changeSize_stats 10 .csv')
stats_2009 <- read.csv('changeSize_stats 09 .csv')
stats_2008 <- read.csv('changeSize_stats 08 .csv')
stats_2007 <- read.csv('changeSize_stats 07 .csv')
stats_2006 <- read.csv('changeSize_stats 06 .csv')
stats_2005 <- read.csv('changeSize_stats 05 .csv')
stats_2004 <- read.csv('changeSize_stats 04 .csv')
stats_2003 <- read.csv('changeSize_stats 03 .csv')
stats_2002 <- read.csv('changeSize_stats 02 .csv')
stats_2001 <- read.csv('changeSize_stats 01 .csv')
stats_2000 <- read.csv('changeSize_stats 00 .csv')

# 


time = c(2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010,
         2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020)

time2 = as.character(time)


df <- data.frame(stats_2000$clump.size, stats_2001$clump.size, stats_2002$clump.size, stats_2003$clump.size, 
                 stats_2004$clump.size, stats_2005$clump.size, stats_2006$clump.size, stats_2007$clump.size, stats_2008$clump.size,
                 stats_2009$clump.size, stats_2010$clump.size, stats_2011$clump.size, stats_2012$clump.size, stats_2013$clump.size,
                 stats_2014$clump.size, stats_2015$clump.size, stats_2016$clump.size, stats_2017$clump.size, stats_2018$clump.size,
                 stats_2019$clump.size, stats_2020$clump.size)

stats <- t(df)

mean_change <- stats[,1]
min_change <- stats[,2]
Q1_change <- stats[,3]
median_change <- stats[,4]
Q3_change <- stats[,5]
max_change <- stats[,6]

timedf <- data.frame(time, mean_change, min_change, Q1_change, median_change, Q3_change, max_change)

# ## all hisorical period (1995 - 2020)
timedf_FCP <- data.frame(time, mean_change, min_change, Q1_change, median_change, Q3_change, max_change,"3")
maxdf_FCP <- data.frame(max_change)
maxdf_FCP$x <- "3"
# ## 2000 -  2020
timedf_FCP <- timedf_FCP[0:21,]
# 
# 


## set working directory
setwd('/home/stevie/Change_detection/KIUIC/')

# load all change size stat tables
stats_2020 <- read.csv('changeSize_stats 20 .csv')
stats_2019 <- read.csv('changeSize_stats 19 .csv')
stats_2018 <- read.csv('changeSize_stats 18 .csv')
stats_2017 <- read.csv('changeSize_stats 17 .csv')
stats_2016 <- read.csv('changeSize_stats 16 .csv')
stats_2015 <- read.csv('changeSize_stats 15 .csv')
stats_2014 <- read.csv('changeSize_stats 14 .csv')
stats_2013 <- read.csv('changeSize_stats 13 .csv')
stats_2012 <- read.csv('changeSize_stats 12 .csv')
stats_2011 <- read.csv('changeSize_stats 11 .csv')
stats_2010 <- read.csv('changeSize_stats 10 .csv')
stats_2009 <- read.csv('changeSize_stats 09 .csv')
stats_2008 <- read.csv('changeSize_stats 08 .csv')
stats_2007 <- read.csv('changeSize_stats 07 .csv')
stats_2006 <- read.csv('changeSize_stats 06 .csv')
stats_2005 <- read.csv('changeSize_stats 05 .csv')
stats_2004 <- read.csv('changeSize_stats 04 .csv')
stats_2003 <- read.csv('changeSize_stats 03 .csv')
stats_2002 <- read.csv('changeSize_stats 02 .csv')
stats_2001 <- read.csv('changeSize_stats 01 .csv')
stats_2000 <- read.csv('changeSize_stats 00 .csv')

# 


time = c(2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010,
         2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020)

time2 = as.character(time)


df <- data.frame(stats_2000$clump.size, stats_2001$clump.size, stats_2002$clump.size, stats_2003$clump.size, 
                 stats_2004$clump.size, stats_2005$clump.size, stats_2006$clump.size, stats_2007$clump.size, stats_2008$clump.size,
                 stats_2009$clump.size, stats_2010$clump.size, stats_2011$clump.size, stats_2012$clump.size, stats_2013$clump.size,
                 stats_2014$clump.size, stats_2015$clump.size, stats_2016$clump.size, stats_2017$clump.size, stats_2018$clump.size,
                 stats_2019$clump.size, stats_2020$clump.size)

stats <- t(df)

mean_change <- stats[,1]
min_change <- stats[,2]
Q1_change <- stats[,3]
median_change <- stats[,4]
Q3_change <- stats[,5]
max_change <- stats[,6]

timedf <- data.frame(time, mean_change, min_change, Q1_change, median_change, Q3_change, max_change)




# 
timedf_KIUIC<- data.frame(time, mean_change, min_change, Q1_change, median_change, Q3_change, max_change, "2")
maxdf_KIUIC <- data.frame(time,max_change)
maxdf_KIUIC$x <- "2"
# ## 2000 -  2020
timedf_KIUIC <- timedf_KIUIC[0:21,]
# 
# 

## set working directory
setwd('/home/stevie/Change_detection/PALMAR/')

# load all change size stat tables
stats_2020 <- read.csv('changeSize_stats 20 .csv')
stats_2019 <- read.csv('changeSize_stats 19 .csv')
stats_2018 <- read.csv('changeSize_stats 18 .csv')
stats_2017 <- read.csv('changeSize_stats 17 .csv')
stats_2016 <- read.csv('changeSize_stats 16 .csv')
stats_2015 <- read.csv('changeSize_stats 15 .csv')
stats_2014 <- read.csv('changeSize_stats 14 .csv')
stats_2013 <- read.csv('changeSize_stats 13 .csv')
stats_2012 <- read.csv('changeSize_stats 12 .csv')
stats_2011 <- read.csv('changeSize_stats 11 .csv')
stats_2010 <- read.csv('changeSize_stats 10 .csv')
stats_2009 <- read.csv('changeSize_stats 09 .csv')
stats_2008 <- read.csv('changeSize_stats 08 .csv')
stats_2007 <- read.csv('changeSize_stats 07 .csv')
stats_2006 <- read.csv('changeSize_stats 06 .csv')
stats_2005 <- read.csv('changeSize_stats 05 .csv')
stats_2004 <- read.csv('changeSize_stats 04 .csv')
stats_2003 <- read.csv('changeSize_stats 03 .csv')
stats_2002 <- read.csv('changeSize_stats 02 .csv')
stats_2001 <- read.csv('changeSize_stats 01 .csv')
stats_2000 <- read.csv('changeSize_stats 00 .csv')

# 


time = c(2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010,
         2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020)

time2 = as.character(time)


df <- data.frame(stats_2000$clump.size, stats_2001$clump.size, stats_2002$clump.size, stats_2003$clump.size, 
                 stats_2004$clump.size, stats_2005$clump.size, stats_2006$clump.size, stats_2007$clump.size, stats_2008$clump.size,
                 stats_2009$clump.size, stats_2010$clump.size, stats_2011$clump.size, stats_2012$clump.size, stats_2013$clump.size,
                 stats_2014$clump.size, stats_2015$clump.size, stats_2016$clump.size, stats_2017$clump.size, stats_2018$clump.size,
                 stats_2019$clump.size, stats_2020$clump.size)

stats <- t(df)

mean_change <- stats[,1]
min_change <- stats[,2]
Q1_change <- stats[,3]
median_change <- stats[,4]
Q3_change <- stats[,5]
max_change <- stats[,6]

timedf <- data.frame(time, mean_change, min_change, Q1_change, median_change, Q3_change, max_change)



# 
timedf_PALMAR <- data.frame(time, mean_change, min_change, Q1_change, median_change, Q3_change, max_change,"1")
maxdf_PALMAR <- data.frame(max_change)
maxdf_PALMAR$x <- "1"
# ## 2000 -  2020
timedf_PALMAR <- timedf_PALMAR[0:21,]

par(mfrow=c(1,1))

 
# PLOT ALL THREE SITES
ggplot() + geom_line(data = timedf_KIUIC, aes(x = time, y = max_change, colour = "KK"))  +
   geom_line(data =timedf_PALMAR, aes(x =time, y = max_change, colour = "EP")) +
   geom_line(data =timedf_FCP, aes(x=time, y=max_change, colour = "FCP"), na.rm = T, stat="identity") +
   labs(colour = "",x ="Year of change", y ="Area of change (ha)") +
   theme_light()

par(mfrow=c(1,1), mar = c(5,5,1,1))
## option 2   
plot(timedf_KIUIC$time, timedf_KIUIC$max_change, type= "l", lwd =3, ylab = "Forest cover loss (ha)",xlim = c(2000, 2020),  xlab="Year", col ="blue", cex.lab=1.5, cex.axis = 1.2)
lines(timedf_FCP$time, timedf_FCP$max_change, type = "l", lwd = 3, col = "green")
lines(timedf_PALMAR$time, timedf_PALMAR$max_change, type = "l", lwd = 3, col = "red")
legend("topleft", legend=c("Semi-deciduous", "Semi-evergreen", "Deciduous"),
       col=c("blue", "green", "red"), lty=1, lwd = 3)

## plot change area barplot
ggplot(timedf, aes(x = time, y = mean_change)) +                            # bar plot
geom_col(size = 1, color = "darkblue", fill = "white")+ geom_line(col="red")

## plot maximum change
ggplot(timedf) +  geom_line(aes(x = time, y = max_change), size = 1.5, color="dark red", group = 1) + ggtitle("Maximum area of change")
    

## plot quartiles
  ggplot(timedf) + geom_line(aes(x = time, y = Q3_change, colour="3rd_Quartile"), size = 1, group = 1) + 
    labs(y = "Change area (ha)") + ggtitle("Area of change (Quartiles)") +
    geom_line(aes(x = time, y = median_change, color="Median"), size = 1, linetype = "dashed", group = 1) +
    geom_line(aes(x = time, y = mean_change, color="Mean"), size = 1, group = 1) +
    geom_line(aes(x = time, y = Q1_change, color="1st_Quartile"), size = 1, group = 1) +
    geom_line(aes(x = time, y = min_change, color="Minimum"), size = 1, group = 1) +
  scale_color_viridis(discrete = TRUE, option = "plasma")  + theme_light()

 
  
