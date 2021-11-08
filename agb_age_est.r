## read in chronosequence data for obtaining model
chrono_data<- read.csv("path/to/data.csv")

lm1 <- lm((chrono_data$Biomasa)~ log(chrono_data$Edad))
summary(lm1)
lm2 <- lm(log(chrono_data$Biomasa)~ log(chrono_data$Edad))
summary(lm2)
lm3 <- lm(log(chrono_data$Biomasa)~ (chrono_data$Edad))
summary(lm3)
lm4 <- lm((chrono_data$Biomasa)~ log(chrono_data$Edad))
summary(lm4)
lm5 <- lm(sqrt(chrono_data$Biomasa)~ sqrt(chrono_data$Edad))
summary(lm5)

##### Get AGB from age - biomass relationship obtained from chronosequence
### function to obtain AGB from age using a log-transformed agb-age relationship

get_agb <- function(age) {
  exp(lm2$coefficients[1] + log(age)*lm2$coefficients[2])
}

## apply the model
age_agb <- get_agb(age)