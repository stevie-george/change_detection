### binomial logistic regression (bernoulli)
## script by Stephanie P. George 
# this script performs binomial logistic regression to  associate presence and absence of deforestation events
# and finds the infleccion point, or the magnitude at which there is 50% chance find a deforestation event

#install.packages("popbio")
library(popbio)

#### define your parameters for 
# 0 = forest (persistent)
# 1 = deforested area

## define useful paths
path <- file.path("/home/stevie/Change_detection")
## set working directory
setwd(path)

# read in data 
data = read.csv('deforestation_training_sites.csv')
head(data)
tail(data)

## general linear model family binomial
modelo.logistico <- glm(deforestation ~ magnitude, family="binomial", data=data)
plot(deforestation~magnitude, data=data, pch= 20, ylab = "Deforestation probability (p)", xlab = "Magnitude", cex.axis = 1.2, cex.lab = 1.5)

# show the curve of the fit 
curve(predict(modelo.logistico,data.frame(magnitude=x),type="resp"),add=TRUE, lwd = 3)
## logi.hist.plot(data$magitude_re,data$deforestation,boxp=FALSE,type="hist",col="gray")
# The binomial glm using the logit link uses a logistic regression equation of the form
#
#    ln(p/(1-p)) = a + bx
#
# where b and a are the slope and the intercept respectively. These can be returned from the 
# glm function using 

coef(modelo.logistico)

#then we rearrange the equation above to be in terms of x, and use p = 0.5. 
p <- 0.5
x <- (log(p/(1-p)) - coef(modelo.logistico)[1]) / coef(modelo.logistico)[2]
x # is our magnitude threshold !
###

abline(.5,0, col=c("blue"), lty=c(2), lwd=3)
abline (-.061, 0, col = c("gray"), lty = 2, lwd = 3)
text(-0.2, .55,"Magnitude Threshold (p 0.5) = -0.061", cex = 1.2)

