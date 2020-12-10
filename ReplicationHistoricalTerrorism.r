#################################################################################
#
#
          # Main replication file Modernization and Terrorism
#
#################################################################################
library(plm)
library(MASS)
library(lme4)
library(caret)
library(brglm)
library(dummies)
library(ggplot2)
library(lme4)
library(haven)
library(dplyr)
library(logistf)
library(dotwhisker)
library(margins)
library(survival)

rm(list = ls())

#Set up your own working directory here#

#Read in the dataset#
hist_terr <- read_dta("Historical Terrorism_v2.dta")

#Generate the PML fixed effects
y_panel <- aggregate(terronset ~ ccode, hist_terr, mean)
y_panel$censor <- ifelse(y_panel$terronset <= 0, 1, 0)
y_panel$uncensor <- ifelse(y_panel$terronset > 0, 1, 0)
for (i in 1:nrow(y_panel)) {
    hist_terr$censor[hist_terr$ccode == y_panel$ccode[i]] <- y_panel$censor[i]
}

hist_terr$uncensor <- 1 - hist_terr$censor
hist_terr$ccode_un <- hist_terr$ccode*hist_terr$uncensor

y_panel <- aggregate(terronset ~ year, hist_terr, mean)
y_panel$ycensor <- ifelse(y_panel$terronset <= 0, 1, 0)
y_panel$yuncensor <- ifelse(y_panel$terronset > 0, 1, 0)
for (i in 1:nrow(y_panel)) {
    hist_terr$ycensor[hist_terr$year == y_panel$year[i]] <- y_panel$ycensor[i]
}

hist_terr$yuncensor <- 1 - hist_terr$ycensor
hist_terr$y_un <- hist_terr$year*hist_terr$yuncensor

#Transform variables
hist_terr$lngdp <- log(hist_terr$e_migdppc + 1)
hist_terr$lngdpg <- log(hist_terr$e_migdpgro + 1)
hist_terr$lnpop <- log(hist_terr$e_population + 1)
hist_terr$urban <- log(hist_terr$e_miurbpop + 1)
hist_terr$lnpopc <- hist_terr$lnpop
hist_terr$urbanc <- hist_terr$urban
hist_terr$growthdiff <- hist_terr$lngdpg

hist_terr <- hist_terr %>% group_by(ccode) %>%
             mutate_at(vars(urbanc:growthdiff), list(~ .x - lag(.x)))


#Model 1
model1 <- brglm(terronset ~ e_polity2 + I(e_polity2^2) + lngdp + lngdpg + growthdiff + lnpop + urban + e_miinterc + e_miinteco + factor(ccode_un) + factor(y_un),
                family=binomial(link="logit"), 
                control.brglm=brglm.control(br.maxit=1000),
                method = "brglm.fit", p1 = T, data = hist_terr
                )

model2 <- glm(terronset ~ e_polity2 + I(e_polity2^2) + lngdp + lngdpg + growthdiff + lnpop + urban + e_miinterc + e_miinteco,
                family=binomial(link="logit"),
                data = hist_terr
                )

#Coefficient plot
m1_df <- broom::tidy(model2)
m1_df <- m1_df[-seq(nrow(m1_df),nrow(m1_df)-65),]

dwplot(m1_df,
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) +
       theme_bw() +
       xlab("Coefficient Estimates") +
       theme(plot.title = element_text(face="bold"),
          legend.position = c(0.75, 0.01),
          legend.justification = c(0, 0),
          legend.background = element_rect(colour="grey80"),
          legend.title.align = .5) +
          scale_colour_grey(start = .2, end = .2)

#Margins plot
margins1 <- margins(model1)
cplot(model2, x = "urbanc", se.type = "shade")

cplot(model2, x = "lngdp", dx = "lngdpg", what = "effect", se.type = "shade", level = 0.90)


# Missing
library(Amelia)

m_df <-  subset(hist_terr, select = c(terronset, e_boix_regime, lngdp, lngdpg, e_peedgini, e_peaveduc, 
                               e_miinterc, e_polity2, ccode, ccode_un, year))

df.out <- amelia(m_df, m = 5, ts = "year", cs = "ccode")


## From HERE ##
hist_terr$neg_growth<-ifelse(hist_terr$e_migdpgro<0, hist_terr$e_migdpgro, 0)

hist_terr<-hist_terr%>%
    group_by(ccode)%>%
    mutate(start_time=0:(n()-1))

library(survival)
model3<-coxph(Surv(start_time, year, terronset)~e_polity2 + I(e_polity2^2) + lngdp + neg_growth + lnpop + urban + e_miinterc + e_miinteco, data = hist_terr)
summary(model3)

model4<-coxph(Surv(start_time, year, terronset)~e_polity2 + I(e_polity2^2) + lngdp + lngdpg + lnpop + urban + e_miinterc + e_miinteco, data = hist_terr)
summary(model4)


