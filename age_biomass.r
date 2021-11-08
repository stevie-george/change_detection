








chrono_data <- read.csv("/home/stevie/CARDAMOM/chrono/Kiuic_Chrono_2010_DALEC.csv")



model <- lm(log(chrono_data$ALL) ~ log(chrono_data$Edad.Estimada))
summary(model)
plot(log(chrono_data$ALL) ~log(chrono_data$Edad.Estimada),pch=20 , xlab =" log(Age)", ylab="log(Aboveground biomass)")
abline(model)



# Call:
#   lm(formula = log(chrono_data$ALL) ~ log(chrono_data$Edad.Estimada))
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.30969 -0.22571  0.02199  0.29319  0.96731 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                     2.73827    0.10467   26.16   <2e-16 ***
#   log(chrono_data$Edad.Estimada)  0.61335    0.03551   17.27   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4276 on 274 degrees of freedom
# Multiple R-squared:  0.5213,	Adjusted R-squared:  0.5196 
# F-statistic: 298.4 on 1 and 274 DF,  p-value: < 2.2e-16
