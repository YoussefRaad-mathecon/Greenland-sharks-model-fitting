library(imager)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(momentuHMM)
library(expm)
library(MASS)
library(clock)
library(qdapRegex)
library(mixR)
library(kableExtra)
library(pbs)


# Null 2 states
# distributions for observation processes
dist <- list(Depth = "gamma")
### Shark 30
dat30 <- prepData(df30, coordNames = NULL)
# initial parameters
Par0_m30 <- list(Depth=c(mean30_Depth[1],mean30_Depth[2],
                         sd30_Depth[1],sd30_Depth[2]))
# fit model
m30 <- fitHMM(data = dat30, nbStates = 2, dist = dist, Par0 = Par0_m30)





### Shark 17
dat17 <- prepData(df17, coordNames = NULL)

# initial parameters
Par0_m17 <- list(Depth=c(mean17_Depth[1],mean17_Depth[2],
                         sd17_Depth[1],sd17_Depth[2]))
# fit model
m17 <- fitHMM(data = dat17, nbStates = 2, dist = dist, Par0 = Par0_m17)




### Shark 6
dat6 <- prepData(df6, coordNames = NULL)

# initial parameters
Par0_m6 <- list(Depth=c(mean6_Depth[1],mean6_Depth[2],
                        sd6_Depth[1],sd6_Depth[2]))
# fit model
m6 <- fitHMM(data = dat6, nbStates = 2, dist = dist, Par0 = Par0_m6)



# Null 3 states


### Shark 30
# initial parameters
Par0_m30_3 <- list(Depth=c(mean30_Depth[1],mean30_Depth[2],mean30_Depth[3],
                           sd30_Depth[1],sd30_Depth[2],sd30_Depth[3]))
# fit model
m30_3 <- fitHMM(data = dat30, nbStates = 3, dist = dist, Par0 = Par0_m30_3)



### Shark 17
dat17 <- prepData(df17, coordNames = NULL)

# initial parameters
Par0_m17_3 <- list(Depth=c(mean17_Depth[1],mean17_Depth[2],mean17_Depth[3],
                           sd17_Depth[1],sd17_Depth[2],sd17_Depth[3]))
# fit model
m17_3 <- fitHMM(data = dat17, nbStates = 3, dist = dist, Par0 = Par0_m17_3)




### Shark 6
dat6 <- prepData(df6, coordNames = NULL)

# initial parameters
Par0_m6_3 <- list(Depth=c(mean6_Depth[1],mean6_Depth[2],mean6_Depth[3],
                          sd6_Depth[1],sd6_Depth[2],sd6_Depth[3]))
# fit model
m6_3 <- fitHMM(data = dat6, nbStates = 3, dist = dist, Par0 = Par0_m6_3)








# Covariate - Temperature - for transition probabilities - 2 states


### Shark 30
#remove NA's
dat30$Temperature <- imputeTS::na_locf(dat30$Temperature)
# initial parameters (obtained from nested model m1)
Par0_m30_Temp <- getPar0(model=m30, formula=~Temperature)

# fit model with formula for transition probabilities
m30_Temp <- fitHMM(data = dat30, nbStates = 2, dist = dist, Par0 = Par0_m30_Temp$Par,
                   beta0=Par0_m30_Temp$beta, formula=~Temperature)






### Shark 6
# remove NA
dat6$Temperature <- imputeTS::na_locf(dat6$Temperature)

# initial parameters (obtained from nested model m1)
Par0_m6_Temp <- getPar0(model=m6, formula=~Temperature)

# fit model with formula for transition probabilities
m6_Temp <- fitHMM(data = dat6, nbStates = 2, dist = dist, Par0 = Par0_m6_Temp$Par,
                  beta0=Par0_m6_Temp$beta, formula=~Temperature)




# Covariate - Temperature - for transition probabilities - 3 states



### Shark 30
# remove NA
dat30$Temperature <- imputeTS::na_locf(dat30$Temperature)

# initial parameters (obtained from nested model m1)
Par0_m30_Temp_3 <- getPar0(model=m30_3, formula=~Temperature)

# fit model with formula for transition probabilities
m30_3_Temp <- fitHMM(data = dat30, nbStates = 3, dist = dist, Par0 = Par0_m30_Temp_3$Par,
                     beta0=Par0_m30_Temp_3$beta, formula=~Temperature)





### Shark 6
# remove NA
dat6$Temperature <- imputeTS::na_locf(dat6$Temperature)

# initial parameters (obtained from nested model m1)
Par0_m6_Temp_3 <- getPar0(model=m6_3, formula=~Temperature)

# fit model with formula for transition probabilities
m6_3_Temp <- fitHMM(data = dat6, nbStates = 3, dist = dist, Par0 = Par0_m6_Temp_3$Par,
                    beta0=Par0_m6_Temp_3$beta, formula=~Temperature)



# Covariate - periodic spline of Time Of Day - for transition probabilities - 2 states 


### Shark 30
dat30 <- prepData(df30, coordNames = NULL)
levels(dat30$ID) <- "shark 30"
# initial parameters (obtained from nested model m1)
Par0_m30_TOD <- getPar0(model=m30, formula=~(pbs::pbs(TOD, df=3, periodic=TRUE, Boundary.knots=c(0,24))))

# fit model
m30_TOD <- fitHMM(data = dat30, nbStates = 2, dist = dist, Par0 = Par0_m30_TOD$Par,
                  beta0=Par0_m30_TOD$beta, formula=~(pbs::pbs(TOD, df=3, periodic=TRUE, Boundary.knots=c(0,24))))




### Shark 17
# initial parameters (obtained from nested model m1)
Par0_m17_TOD <- getPar0(model=m17, formula=~(pbs::pbs(TOD, df=3, periodic=TRUE, Boundary.knots=c(0,24))))

# fit model
m17_TOD <- fitHMM(data = dat17, nbStates = 2, dist = dist, Par0 = Par0_m17_TOD$Par,
                  beta0=Par0_m17_TOD$beta, formula=~(pbs::pbs(TOD, df=3, periodic=TRUE, Boundary.knots=c(0,24))))




### Shark 6 TOD 2 state
#getpar0 from m6 null model // intercept beta_0
Par0_m6_TOD <- getPar0(model=m6, 
                       formula=~(pbs::pbs(TOD, df=3, periodic=TRUE, Boundary.knots=c(0,24))))

#Fithmm 
m6_TOD <- fitHMM(data = dat6, nbStates = 2, dist = dist,
                 Par0 = Par0_m6_TOD$Par,
                 beta0=Par0_m6_TOD$beta, 
                 formula=~(pbs::pbs(TOD, df=3, periodic=TRUE, Boundary.knots=c(0,24))))




# Covariate - periodic spline of Time Of Day - for transition probabilities - 3 states 


### Shark 30
# initial parameters (obtained from nested model m1)
Par0_m30_3_TOD <- getPar0(model=m30_3, formula=~(splines2::mSpline(TOD, df=3, periodic=TRUE, Boundary.knots=c(0,24))))

# fit model
m30_3_TOD <- fitHMM(data = dat30, nbStates = 3, dist = dist, Par0 = Par0_m30_3_TOD$Par,
                    beta0=Par0_m30_3_TOD$beta, formula=~(splines2::mSpline(TOD, df=3, periodic=TRUE, Boundary.knots=c(0,24))))









### Shark 17
# initial parameters (obtained from nested model m1)
Par0_m17_3_TOD <- getPar0(model=m17_3, formula=~(pbs::pbs(TOD, df=3, periodic=TRUE, Boundary.knots=c(0,24))))

# fit model
m17_3_TOD <- fitHMM(data = dat17, nbStates = 3, dist = dist, Par0 = Par0_m17_3_TOD$Par,
                    beta0=Par0_m17_3_TOD$beta, formula=~(pbs::pbs(TOD, df=3, periodic=TRUE, Boundary.knots=c(0,24))))




### Shark 6
# initial parameters (obtained from nested model m1)
Par0_m6_3_TOD <- getPar0(model=m6_3, formula=~(pbs::pbs(TOD, df=3, periodic=TRUE, Boundary.knots=c(0,24))))

# fit model
m6_3_TOD <- fitHMM(data = dat6, nbStates = 3, dist = dist, Par0 = Par0_m6_3_TOD$Par,
                   beta0=Par0_m6_3_TOD$beta, formula=~(pbs::pbs(TOD, df=3, periodic=TRUE, Boundary.knots=c(0,24))))




# Covariate - periodic time of year - for transition probabilities - 2 states 


### Shark 30
dat30 <- prepData(df30, coordNames = NULL)
levels(dat30$ID) <- "shark 30"
# initial parameters (obtained from nested model m1)
Par0_m30_month <- getPar0(model=m30, formula=~(pbs::pbs(TOY, df=3, periodic=TRUE, Boundary.knots=c(1,13))))

# fit model with formula for transition probabilities
m30_month <- fitHMM(data = dat30, nbStates = 2, dist = dist, Par0 = Par0_m30_month$Par,
                    beta0=Par0_m30_month$beta, formula=~(pbs::pbs(TOY, df=3, periodic=TRUE, Boundary.knots=c(1,13))))





### Shark 17
#initial parameters (obtained from nested model m1)
Par0_m17_month <- getPar0(model=m17, formula=~pbs::pbs(TOY, df=3, periodic=TRUE, Boundary.knots=c(1,13)))

# fit model with formula for transition probabilities
m17_month <- fitHMM(data = dat17, nbStates = 2, dist = dist, Par0 = Par0_m17_month$Par,
                    beta0=Par0_m17_month$beta, formula=~pbs::pbs(TOY, df=3, periodic=TRUE, Boundary.knots=c(1,13)))




### Shark 6
# initial parameters (obtained from nested model m1)
Par0_m6_month <- getPar0(model=m6, formula=~pbs::pbs(TOY, df=3, periodic=TRUE, Boundary.knots=c(1,13)))

#fit model with formula for transition probabilities
m6_month <- fitHMM(data = dat6, nbStates = 2, dist = dist, Par0 = Par0_m6_month$Par,
                   beta0=Par0_m6_month$beta, formula=~pbs::pbs(TOY, df=3, periodic=TRUE, Boundary.knots=c(1,13)))



# Covariate - periodic time of year - for transition probabilities - 3 states 



### Shark 30
#initial parameters (obtained from nested model m1)
Par0_m30_3_month <- getPar0(model=m30_3, formula=~pbs::pbs(TOY, df=3, periodic=TRUE, Boundary.knots=c(1,13)))

# fit model with formula for transition probabilities
m30_3_month <- fitHMM(data = dat30, nbStates = 3, dist = dist, Par0 = Par0_m30_3_month$Par,
                      beta0=Par0_m30_3_month$beta, formula=~pbs::pbs(TOY, df=3, periodic=TRUE, Boundary.knots=c(1,13)))





### Shark 17
# initial parameters (obtained from nested model m1)
Par0_m17_3_month <- getPar0(model=m17_3, formula=~pbs::pbs(TOY, df=3, periodic=TRUE, Boundary.knots=c(1,13)))

# fit model with formula for transition probabilities
m17_3_month <- fitHMM(data = dat17, nbStates = 3, dist = dist, Par0 = Par0_m17_3_month$Par,
                      beta0=Par0_m17_3_month$beta, formula=~pbs::pbs(TOY, df=3, periodic=TRUE, Boundary.knots=c(1,13)))





### Shark 6
# initial parameters (obtained from nested model m1)
Par0_m6_3_month <- getPar0(model=m6_3, formula=~pbs::pbs(TOY, df=3, periodic=TRUE, Boundary.knots=c(1,13)))

# fit model with formula for transition probabilities
m6_3_month <- fitHMM(data = dat6, nbStates = 3, dist = dist, Par0 = Par0_m6_3_month$Par,
                     beta0=Par0_m6_3_month$beta, formula=~pbs::pbs(TOY, df=3, periodic=TRUE, Boundary.knots=c(1,13)))



# temperature + TOD - 2 states


### Shark 30
# initial parameters (obtained from nested model m1)
Par0_m30_TODTemp<- getPar0(model=m30, formula=~((pbs::pbs(TOD, df=3, periodic=TRUE, Boundary.knots=c(0,24)))+Temperature))

# fit model
m30_TODTemp <- fitHMM(data = dat30, nbStates = 2, dist = dist, Par0 = Par0_m30_TODTemp$Par,
                      beta0=Par0_m30_TODTemp$beta, formula=~((pbs::pbs(TOD, df=3, periodic=TRUE, Boundary.knots=c(0,24)))+Temperature))





### Shark 6
# initial parameters (obtained from nested model m1)
Par0_m6_TODTemp<- getPar0(model=m6, formula=~((pbs::pbs(TOD, df=3, periodic=TRUE, Boundary.knots=c(0,24)))+Temperature))

# fit model
m6_TODTemp <- fitHMM(data = dat6, nbStates = 2, dist = dist, Par0 = Par0_m6_TODTemp$Par,
                     beta0=Par0_m6_TODTemp$beta, formula=~((pbs::pbs(TOD, df=3, periodic=TRUE, Boundary.knots=c(0,24)))+Temperature))




# temperature + TOD - 3 states

### Shark 30
# initial parameters (obtained from nested model m1)
Par0_m30_3_TODTemp <- getPar0(model=m30_3, formula=~((pbs::pbs(TOD, df=3, periodic=TRUE, Boundary.knots=c(0,24)))+Temperature))

# fit model
m30_3_TODTemp <- fitHMM(data = dat30, nbStates = 3, dist = dist, Par0 = Par0_m30_3_TODTemp$Par,
                        beta0=Par0_m30_3_TODTemp$beta, formula=~((pbs::pbs(TOD, df=3, periodic=TRUE, Boundary.knots=c(0,24)))+Temperature))





### Shark 6
# initial parameters (obtained from nested model m1)
Par0_m6_3_TODTemp <- getPar0(model=m6_3, formula=~((pbs::pbs(TOD, df=3, periodic=TRUE, Boundary.knots=c(0,24)))+Temperature))

# fit model
m6_3_TODTemp <- fitHMM(data = dat6, nbStates = 3, dist = dist, Par0 = Par0_m6_3_TODTemp$Par,
                       beta0=Par0_m6_3_TODTemp$beta, formula=~((pbs::pbs(TOD, df=3, periodic=TRUE, Boundary.knots=c(0,24)))+Temperature))



