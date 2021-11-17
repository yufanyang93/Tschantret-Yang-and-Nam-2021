library(countrycode)
library(brglm2)
library(plm)
library(MASS)
library(caret)
library(lme4)
library(haven)
library(dplyr)
library(tidyr)
library(miceadds)
library(logistf)
library(margins)
library(ggplot2)
library(texreg)
library(stargazer)
library(mice)
library(lmtest)
library(multiwayvcov)
library(permute)
library(sandwich)
library(mitools)
install.packages("remotes")
remotes::install_github("LukasWallrich/rNuggets")
library(rNuggets)
library(splitstackshape)

rm(list = ls())

full_data<-load("final_data.Rda")

###              Multiple imputation with mice                               ###
###              The analysis below is the one in the paper                  ###

library(mice)
library(doRNG)

missing_mice<-subset(full_data, select = c(e_migdppc, e_migdpgro, e_mipopula, e_miurbani, e_polity2, e_regiongeo, e_regionpol, 
                                           v2csprtcpt, v2pepwrsoc, e_miinterc, e_miinteco, v2x_polyarchy, v2x_libdem, 
                                           v2x_partipdem, v2x_delibdem, v2x_egaldem, COWcode, year, nationalist_d, leftist_d, 
                                           rightist_d, religious_d, anarchist_d, other_d, time, time_s, time_c, ccode_un, y_un, 
                                           terr_sum, nationalist, leftist, rightist, religious, anarchist, other, terronset))%>%
  distinct(COWcode, year, .keep_all=T)%>%
  filter(complete.cases(COWcode))

missing_mice<-missing_mice%>%
  mutate(Europe=ifelse(e_regiongeo==1|e_regiongeo==2|e_regiongeo==3|e_regiongeo==4, 1, 0))%>%
  mutate(North_America=ifelse(e_regiongeo==16, 1, 0))%>%
  mutate(SSA=ifelse(e_regionpol==4, 1, 0))%>%
  mutate(Latin_America=ifelse(e_regionpol==2|e_regionpol==10, 1, 0))%>%
  mutate(Asia=ifelse(e_regiongeo==10|e_regiongeo==11|e_regiongeo==12|e_regiongeo==13|e_regiongeo==14, 1, 0))%>%
  mutate(Oceania=ifelse(e_regiongeo==15, 1, 0))%>%
  mutate(MENA=ifelse(e_regionpol==3, 1, 0))

allVars <- names(missing_mice)
missVars<-names(missing_mice)[colSums(is.na(missing_mice)) > 0]

predictorMatrix <- matrix(0, ncol = length(allVars), nrow = length(allVars))
rownames(predictorMatrix) <- allVars
colnames(predictorMatrix) <- allVars

imputerVars <- c("v2x_polyarchy", "v2x_libdem", "v2x_partipdem", "v2x_delibdem", "v2x_egaldem",
                 "COWcode", "year", "e_miinterc", "e_miinteco", "v2csprtcpt", "v2pepwrsoc")

imputerVars <- intersect(unique(imputerVars), allVars)
imputerMatrix <- predictorMatrix
imputerMatrix[,imputerVars] <- 1

imputedOnlyVars <- c("e_migdppc", "e_migdpgro", "e_mipopula", "e_miurbani", "e_polity2")
imputedVars <- intersect(unique(c(imputedOnlyVars, imputerVars)), missVars)
imputedMatrix <- predictorMatrix
imputedMatrix[imputedVars,] <- 1

predictorMatrix <- imputerMatrix * imputedMatrix
diag(predictorMatrix) <- 0

dryMice <- mice(data = missing_mice, m = 1, predictorMatrix = predictorMatrix, maxit = 0)
predictorMatrix <- dryMice$predictorMatrix
imputerVars <- colnames(predictorMatrix)[colSums(predictorMatrix) > 0]
imputedVars <- rownames(predictorMatrix)[rowSums(predictorMatrix) > 0]
setdiff(imputerVars, imputedVars)
intersect(imputerVars, imputedVars)
setdiff(imputedVars, imputerVars)
setdiff(missVars, imputedVars)
predictorMatrix[rowSums(predictorMatrix) > 0, colSums(predictorMatrix) > 0]

dryMice$method[setdiff(allVars, imputedVars)] <- ""
dryMice$method[sapply(dryMice$method, nchar) > 0]
M <- 4

set.seed(12345)
miceout <- foreach(i = seq_len(M), .combine = ibind) %dorng% {
  miceout <- mice(data = missing_mice, m = 5, pmaxit=5, rint = TRUE,
                  predictorMatrix = predictorMatrix, method = dryMice$method,
                  MaxNWts = 2000)
  miceout
}

actuallyImputedVars <-
  setdiff(names(missing_mice)[colSums(is.na(missing_mice)) > 0],
          names(complete(miceout, action = 1))[colSums(is.na(complete(miceout, action = 1))) > 0])
actuallyImputedVars

names(complete(miceout, action = 1))[colSums(is.na(complete(miceout, action = 1))) > 0]

###Split data by pre/post- WWII###
d_long<-mice::complete(miceout, "long", include=T)
d_long_pre1945<-d_long[which(d_long$year<=1945),]
d_long_post1945<-d_long[which(d_long$year>=1946),]
pre1945<-as.mids(d_long_pre1945)
post1945<-as.mids(d_long_post1945)

###Basic Model##
model1 <- with(miceout, {glm(terronset ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt + 
                               e_miinterc + factor(ccode_un) + factor(y_un),
                             family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})
model1_pool<-pool(model1)

model2 <- with(miceout, {glm(nationalist_d ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) + e_miurbani + v2csprtcpt + 
                               e_miinterc + factor(ccode_un) + factor(y_un),
                             family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})
model2_pool<-pool(model2)

model3 <- with(miceout, {glm(leftist_d ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt +  
                               e_miinterc + factor(ccode_un) + factor(y_un),
                             family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})
model3_pool<-pool(model3)

model4 <- with(miceout, {glm(anarchist_d ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt + 
                               e_miinterc + factor(ccode_un) + factor(y_un),
                             family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})
model4_pool<-pool(model4)

model5<-with(miceout, {glm(other_d ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt + 
                             e_miinterc + factor(ccode_un) + factor(y_un),
                           family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})

model5_pool<-pool(model5)

texreg(list(model1_pool, model2_pool, model3_pool, model4_pool, model5_pool), 
       custom.model.names = c("All Groups", "Nationalist Groups", "Leftist Groups", "Anarchist Groups", "Other Groups"),
       omit.coef = "factor",
       custom.coef.names = c("(Intercept)", "Polity", "Polity Squared", "GDP per capita (logged)", "GDP per capita Growth", 
                             "Population (logged)", "Urbanization", "Civil Society Participation", "Civil Conflict"),
       caption = "Penalized Regression Results", file="Panelized Regression with Imputed Data.txt",
       stars = c(0.001, 0.01, 0.05))

wordreg(list(model1_pool, model2_pool, model3_pool, model4_pool, model5_pool), 
        custom.model.names = c("All Groups", "Nationalist Groups", "Leftist Groups", "Anarchist Groups", "Other Groups"),
        custom.coef.names = c("(Intercept)", "Polity", "Polity Squared", "GDP per capita (logged)", "GDP per capita Growth", 
                              "Population (logged)", "Urbanization","Civil Society Participation", 
                              "Civil Conflict"),
        caption = "Penalized Regression Results", file="Panelized Regression with Imputed Data.doc",
        omit.coef = "factor", stars = c(0.001, 0.01, 0.05))

###Pre-1945 Analysis##
pre1 <- with(pre1945, {glm(terronset ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt + 
                             e_miinterc + factor(ccode_un) + factor(y_un),
                           family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})
pre1_pool<-pool(pre1)

pre2 <- with(pre1945, {glm(nationalist_d ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt + 
                             e_miinterc + factor(ccode_un) + factor(y_un),
                           family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})
pre2_pool<-pool(pre2)

pre3 <- with(pre1945, {glm(leftist_d ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt +  
                             e_miinterc + factor(ccode_un) + factor(y_un),
                           family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})
pre3_pool<-pool(pre3)

pre4 <- with(pre1945, {glm(anarchist_d ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt + 
                             e_miinterc + factor(ccode_un) + factor(y_un),
                           family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})
pre4_pool<-pool(pre4)

pre5<-with(pre1945, {glm(other_d ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt + 
                           e_miinterc + factor(ccode_un) + factor(y_un),
                         family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})
pre5_pool<-pool(pre5)

texreg(list(pre1_pool, pre2_pool, pre3_pool, pre4_pool, pre5_pool), 
       custom.model.names = c("All Groups", "Nationalist Groups", "Leftist Groups", "Anarchist Groups", "Other Groups"),
       custom.coef.names = c("(Intercept)", "Polity", "Polity Squared", "GDP per capita (logged)", "GDP per capita Growth", 
                             "Population (logged)", "Urbanization","Civil Society Participation", 
                             "Civil Conflict"),
       caption = "Penalized Regression Results (pre-1945)", file="Panelized Regression with Imputed Data (pre-1945).txt",
       omit.coef = "factor", stars = c(0.001, 0.01, 0.05))

wordreg(list(pre1_pool, pre2_pool, pre3_pool, pre4_pool, pre5_pool), 
        custom.model.names = c("All Groups", "Nationalist Groups", "Leftist Groups", "Anarchist Groups", "Other Groups"),
        custom.coef.names = c("(Intercept)", "Polity", "Polity Squared", "GDP per capita (logged)", "GDP per capita Growth", 
                              "Population (logged)", "Urbanization","Civil Society Participation", 
                              "Civil Conflict"),
        caption = "Penalized Regression Results (pre-1945)", file="Panelized Regression with Imputed Data (pre-1945).doc",
        omit.coef = "factor", stars = c(0.001, 0.01, 0.05))

pre6 <- with(pre1945, {glm(terronset ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt +  
                             e_miinterc +  North_America + Europe + SSA + Asia + MENA + Latin_America + time + time_s+ time_c,
                           family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})
pre6_pool<-pool(pre6)

pre7 <- with(pre1945, {glm(nationalist_d ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt +  
                             e_miinterc +  North_America + Europe + SSA + Asia + MENA + Latin_America + time + time_s+ time_c,
                           family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})
pre7_pool<-pool(pre7)

pre8 <- with(pre1945, {glm(leftist_d ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt + 
                             e_miinterc +  North_America + Europe + SSA + Asia + MENA + Latin_America + time + time_s+ time_c,
                           family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})
pre8_pool<-pool(pre8)

pre9 <- with(pre1945, {glm(anarchist_d ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt + 
                             e_miinterc +  North_America + Europe + SSA + Asia + MENA + Latin_America+ time + time_s+ time_c,
                           family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})
pre9_pool<-pool(pre9)

pre10<-with(pre1945, {glm(other_d ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt +  
                            e_miinterc + North_America + Europe + SSA + Asia + MENA + Latin_America + time + time_s+ time_c,
                          family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})

pre10_pool<-pool(pre10)

texreg(list(pre6_pool, pre7_pool, pre8_pool, pre9_pool, pre10_pool), 
       custom.model.names = c("All Groups", "Nationalist Groups", "Leftist Groups", "Anarchist Groups", "Other Groups"),
       custom.coef.names = c("(Intercept)", "Polity", "Polity Squared", "GDP per capita (logged)", "GDP per capita Growth", 
                             "Population (logged)", "Urbanization","Civil Society Participation", 
                             "Civil Conflict", "North America", "Europe", "Sub-Saharan Africa", "Asia", "MENA", 
                             "Latin America", "Time", "Time Squared", "Time Cubic"),
       caption = "Penalized Regression Results (pre-1945, no FE)", file="Panelized Regression with Imputed Data (pre-1945, no FE).txt",
       omit.coef = "factor", stars = c(0.001, 0.01, 0.05))

wordreg(list(pre6_pool, pre7_pool, pre8_pool, pre9_pool, pre10_pool), 
        custom.model.names = c("All Groups", "Nationalist Groups", "Leftist Groups", "Anarchist Groups", "Other Groups"),
        custom.coef.names = c("(Intercept)", "Polity", "Polity Squared", "GDP per capita (logged)", "GDP per capita Growth", 
                              "Population (logged)", "Urbanization","Civil Society Participation", 
                              "Civil Conflict", "North America", "Europe", "Sub-Saharan Africa", "Asia", "MENA", "
                              Latin America", "Time", "Time Squared", "Time Cubic"),
        caption = "Penalized Regression Results (pre-1945, no FE)", file="Panelized Regression with Imputed Data (pre-1945, no FE).doc",
        omit.coef = "factor", stars = c(0.001, 0.01, 0.05))

###Post-1945 Analysis###
post1 <- with(post1945, {glm(terronset ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt + 
                               e_miinterc + factor(ccode_un) + factor(y_un),
                             family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})
post1_pool<-pool(post1)

post2 <- with(post1945, {glm(nationalist_d ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt + 
                               e_miinterc + factor(ccode_un) + factor(y_un),
                             family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})
post2_pool<-pool(post2)

post3 <- with(post1945, {glm(leftist_d ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt +  
                               e_miinterc + factor(ccode_un) + factor(y_un),
                             family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})
post3_pool<-pool(post3)

post4 <- with(post1945, {glm(anarchist_d ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt + 
                               e_miinterc + factor(ccode_un) + factor(y_un),
                             family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})
post4_pool<-pool(post4)

post5<-with(post1945, {glm(other_d ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt + 
                             e_miinterc + factor(ccode_un) + factor(y_un),
                           family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})

post5_pool<-pool(post5)

texreg(list(post1_pool, post2_pool, post3_pool, post4_pool, post5_pool), 
       custom.model.names = c("All Groups", "Nationalist Groups", "Leftist Groups", "Anarchist Groups", "Other Groups"),
       custom.coef.names = c("(Intercept)", "Polity", "Polity Squared", "GDP per capita (logged)", "GDP per capita Growth", 
                             "Population (logged)", "Urbanization","Civil Society Participation", 
                             "Civil Conflict"),
       caption = "Penalized Regression Results (post-1945)", file="Panelized Regression with Imputed Data (post-1945).txt",
       omit.coef = "factor", stars = c(0.001, 0.01, 0.05))

wordreg(list(post1_pool, post2_pool, post3_pool, post4_pool, post5_pool), 
        custom.model.names = c("All Groups", "Nationalist Groups", "Leftist Groups", "Anarchist Groups", "Other Groups"),
        custom.coef.names = c("(Intercept)", "Polity", "Polity Squared", "GDP per capita (logged)", "GDP per capita Growth", 
                              "Population (logged)", "Urbanization","Civil Society Participation", 
                              "Civil Conflict"),
        caption = "Penalized Regression Results (post-1945)", file="Panelized Regression with Imputed Data (post-1945).doc",
        omit.coef = "factor", stars = c(0.001, 0.01, 0.05))

post6 <- with(post1945, {glm(terronset ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt +  
                               e_miinterc +  North_America + Europe + SSA + Asia + MENA + Latin_America + time + time_s+ time_c,
                             family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})
post6_pool<-pool(post6)

post7 <- with(post1945, {glm(nationalist_d ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt +  
                               e_miinterc +  North_America + Europe + SSA + Asia + MENA + Latin_America + time + time_s+ time_c,
                             family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})
post7_pool<-pool(post7)

post8 <- with(post1945, {glm(leftist_d ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt + 
                               e_miinterc +  North_America + Europe + SSA + Asia + MENA + Latin_America + time + time_s+ time_c,
                             family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})
post8_pool<-pool(post8)

post9 <- with(post1945, {glm(anarchist_d ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt + 
                               e_miinterc +  North_America + Europe + SSA + Asia + MENA + Latin_America + time + time_s+ time_c,
                             family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})
post9_pool<-pool(post9)

post10<-with(post1945, {glm(other_d ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt +  
                              e_miinterc + North_America + Europe + SSA + Asia + MENA + Latin_America + time + time_s+ time_c,
                            family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T)})

post10_pool<-pool(post10)

texreg(list(post6_pool, post7_pool, post8_pool, post9_pool, post10_pool), 
       custom.model.names = c("All Groups", "Nationalist Groups", "Leftist Groups", "Anarchist Groups", "Other Groups"),
       custom.coef.names = c("(Intercept)", "Polity", "Polity Squared", "GDP per capita (logged)", "GDP per capita Growth", 
                             "Population (logged)", "Urbanization","Civil Society Participation", 
                             "Civil Conflict", "North America", "Europe", "Sub-Saharan Africa", "Asia", "MENA", 
                             "Latin America", "Time", "Time Squared", "Time Cubic"),
       caption = "Penalized Regression Results (post-1945, no FE)", file="Panelized Regression with Imputed Data (post-1945, no FE).txt",
       omit.coef = "factor", stars = c(0.001, 0.01, 0.05))

wordreg(list(post6_pool, post7_pool, post8_pool, post9_pool, post10_pool), 
        custom.model.names = c("All Groups", "Nationalist Groups", "Leftist Groups", "Anarchist Groups", "Other Groups"),
        custom.coef.names = c("(Intercept)", "Polity", "Polity Squared", "GDP per capita (logged)", "GDP per capita Growth", 
                              "Population (logged)", "Urbanization","Civil Society Participation", 
                              "Civil Conflict", "North America", "Europe", "Sub-Saharan Africa", "Asia", "MENA", 
                              "Latin America", "Time", "Time Squared", "Time Cubic"),
        caption = "Penalized Regression Results (post-1945, no FE)", file="Panelized Regression with Imputed Data (post-1945, no FE).doc",
        omit.coef = "factor", stars = c(0.001, 0.01, 0.05))

###     Appendix: Survival Analysis                                          ###
###     Note: To run the code below, you will need vdem data downloaded      ###
library(survival)
library(SurvRegCensCov)
library(survminer)
library(flexsurv)
library(ggfortify)
library(eha)

terr_group<-read.csv("terrorist-groups.csv")
terr_group$EndYear<-ifelse(terr_group$EndYear==".", 2019, terr_group$EndYear)
terr_group<-terr_group%>%mutate(duration=as.numeric(EndYear)-StrYear+1)%>%
  mutate(terronset=1)%>%rename(year=StrYear)
terr_group$Ideology<-factor(terr_group$Ideology, levels=c("Nationalist", "Leftist", "Anarchist", "Rightist", "Religious", "Other"))

survival<-terr_group%>%complete(year=seq(min(1860), max(1969), by=1))
survival$duration[is.na(survival$duration)]<-0
survival<-expandRows(survival, "duration", drop = FALSE)
survival<-survival%>%rename(StrYear=year)
survival<-survival%>%
  group_by(Group)%>%
  mutate(count_year=row_number())%>%
  mutate(year=StrYear+count_year-1)%>%
  mutate(fail=ifelse(count_year==duration, 1, 0))%>%
  arrange(Group, year)

survival$COWcode<-countrycode(survival$Country, "country.name", "cown")
survival$COWcode<-ifelse(survival$Country=="Serbia", 345, survival$COWcode)
survival$COWcode<-ifelse(survival$Country=="Palestine", 666, survival$COWcode)

vdem<-readRDS("vdem.rds")
vdem<-vdem%>%filter(year<1971&year>1859)%>%
  select(e_migdppc, e_migdpgro, e_mipopula, e_miurbani, e_polity2, e_regiongeo, 
         e_regionpol, v2csprtcpt, v2pepwrsoc, e_miinterc, e_miinteco, v2x_polyarchy, 
         v2x_libdem, v2x_partipdem, v2x_delibdem, v2x_egaldem, COWcode, year)

survival<-merge(survival, vdem, by = c("COWcode", "year"))
survival$exit<-as.numeric(ifelse(survival$EndYear>1970, NA, survival$EndYear))

fit1<-coxph(Surv(count_year, exit, fail)~e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt +  
              e_miinterc + factor(Ideology), data = survival, ties = "breslow")

fit2<-coxph(Surv(count_year, exit, fail)~e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt +  
              e_miinterc + factor(Ideology) + cluster(COWcode), data = survival, ties = "breslow")

fit3<-weibreg(Surv(count_year, fail)~e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt +  
                e_miinterc + factor(Ideology), data = survival)

fit4<-weibreg(Surv(count_year, fail)~e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt +  
                e_miinterc + factor(Ideology) + cluster(COWcode), data = survival)

stargazer(fit3, fit4, fit1, fit2, type = "html", out = "survival_analysis.htm")

wordreg(list(fit1, fit2), 
        custom.model.names = c("Cox Regression", "Cox Regression with Clustered SE"),
        custom.coef.names = c("Polity", "Polity Squared", "GDP per capita (logged)", "GDP per capita Growth", 
                              "Population (logged)", "Urbanization", "Civil Society Participation", "Civil Conflict", 
                              "Leftist", "Anarchist", "Rightist", "Religious", "Other Ideology"),
        caption = "Terrorist Groups Survival Analysis", file="Cox Analysis",
        stars = c(0.001, 0.01, 0.05))


##Graph##
g1<-survfit(coxph(Surv(count_year, exit, fail)~e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt +  
                    e_miinterc + strata(Ideology), data = survival, ties = "breslow"))

g2<-survfit(coxph(Surv(count_year, exit, fail)~e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt +  
                    e_miinterc + strata(Ideology) + cluster(COWcode), data = survival, ties = "breslow"))

ggsurvplot(g1, data = survival, censor.size=0.1, size=0.6, xlim=c(1900, 1970), linetype = "dashed",
           palette = c("#361D32", "#309975", "#AB3D39", "#F55951", "#EDD2CB", "#FFADAA", coef.int=TRUE,
                       conf.int.fill="grey"))

autoplot(g2, xlim = c(1860, 1970), xlab = "Year", ylab = "Survival", surv.linetype = "dashed", conf.int=FALSE, 
         censor.colour = "#687980", censor.size = 3, surv.alpha = 0.8, surv.size = 0.7)+theme_bw()+
  scale_color_manual(values = c("#999B84", "#206A5D", "#CA8A8B", "#FFCC29", "#8AB6D6", "#FF8474"))

ggsave("survival_curve.png", width=12.5, height=8.25, dpi=300)

p = seq(0.99, 0.01, by=-.01)

###substantive effects##
library(ggstatsplot)
library(emmeans)
library(ggeffects)

pre1_grid<-ref_grid(pre1, transform = "response")
mar_pre1<-emmeans(pre1_grid, ~v2csprtcpt|e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + 
                    (log(e_mipopula+1)) +e_miurbani + e_miinterc)

pre2_grid<-ref_grid(pre2, transform = "response")
mar_pre2<-emmeans(pre2_grid, ~v2csprtcpt|e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + 
                    (log(e_mipopula+1)) +e_miurbani + e_miinterc)

pre3_grid<-ref_grid(pre3, transform = "response")
mar_pre3<-emmeans(pre3_grid, ~v2csprtcpt|e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + 
                    (log(e_mipopula+1)) +e_miurbani + e_miinterc)

pre4_grid<-ref_grid(pre4, transform = "response")
mar_pre4<-emmeans(pre4_grid, ~v2csprtcpt|e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + 
                    (log(e_mipopula+1)) +e_miurbani + e_miinterc)

pre5_grid<-ref_grid(pre5, transform = "response")
mar_pre5<-emmeans(pre5_grid, ~v2csprtcpt|e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + 
                    (log(e_mipopula+1)) +e_miurbani + e_miinterc)

ggcoefstats(mar_pre_all, statistic = "z", conf.int = TRUE, conf.level = .90,
            point.args = list(color = c("#999B84", "#206A5D", "#CA8A8B", "#FFCC29", "#8AB6D6"), shape=18, size = 3),
            xlab = "Marginal Effects of Civil Society Participation", ylab = "Dependent Variable",
            xlim = c(0.010, 0.020))+ theme_bw()+
  ggplot2::scale_y_discrete(labels = c("All Groups", "Nationalist Groups", "Leftist Groups", 
                                       "Anarchist Groups", "Other Groups"))+
  ggplot2::scale_x_continuous(limits = c(0.010, 0.015))

pre_1945_dat<-complete(pre1945)
pre_1945_dat<-pre_1945_dat%>%dplyr::select(terronset, e_polity2, e_migdppc, e_migdpgro, e_mipopula, e_miurbani, 
                                           v2csprtcpt, e_miinterc, ccode_un, y_un, nationalist_d, leftist_d, anarchist_d,
                                           other_d)

post_1945_dat<-complete(post1945)
post_1945_dat<-post_1945_dat%>%dplyr::select(terronset, e_polity2, e_migdppc, e_migdpgro, e_mipopula, e_miurbani, 
                                             v2csprtcpt, e_miinterc, ccode_un, y_un, nationalist_d, leftist_d, anarchist_d,
                                             other_d)

predict_1 <- glm(terronset ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt + 
                   e_miinterc + factor(ccode_un) + factor(y_un),
                 family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T, data = pre_1945_dat)

graph1<-cplot(predict_1, x = "v2csprtcpt", what = c("effect"), type = c("response"),
              n = 5, xlab = "Civil Society Participation", ylab = "Probability of Terrorist Group Onset")

predict_2<-glm(terronset ~ e_polity2 + I(e_polity2^2) + (log(e_migdppc+1)) + e_migdpgro + (log(e_mipopula+1)) +e_miurbani + v2csprtcpt + 
                 e_miinterc + factor(ccode_un) + factor(y_un),
               family=binomial(link="logit"), maxit = 100, method = "brglmFit", p1=T, data = post_1945_dat)

graph2<-cplot(predict_2, x = "v2csprtcpt", what = c("effect"), type = c("response"),
              n = 5, xlab = "Civil Society Participation", ylab = "Probability of Terrorist Group Onset")

###summary statistics###
library(pastecs)
all<-complete(miceout)
all<-all%>%dplyr::select(terronset, e_polity2, e_migdppc, e_migdpgro, e_mipopula, e_miurbani, 
                         v2csprtcpt, e_miinterc, nationalist_d, leftist_d, anarchist_d,
                         other_d)
sum_stat<-stat.desc(all)
sum_stat<-round(sum_stat, 2)

