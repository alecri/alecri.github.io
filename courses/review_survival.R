## ----setup, include = FALSE---------------------------------------------------------------------------------------------
library(knitr)
library(kfigr)
opts_chunk$set(comment = NA, fig.align = "center", warning = FALSE)
options(width = 120, show.signif.stars = F)


## ----extract code (ignore), eval = FALSE, echo = FALSE------------------------------------------------------------------
knitr::purl("review_survival.qmd")


## ----load packages, message=FALSE---------------------------------------------------------------------------------------
if (!require(pacman)) install.packages("pacman")
pacman::p_load(tidyverse, plotly, shiny, scales, patchwork, DiagrammeR,
               survival, KMsurv, flexsurv, eha, ctqr, cmprsk, mstate, ggsurvfit, survminer, 
               Epi, epitools, gnm, splines)
theme_set(theme_classic() + theme(legend.position = "bottom"))


## ----read data----------------------------------------------------------------------------------------------------------
orca <- read.table("http://www.stats4life.se/data/oralca.txt", stringsAsFactors = TRUE)


## ----describe data------------------------------------------------------------------------------------------------------
head(orca)
glimpse(orca)
summary(orca)


## ----boxDiagram, fig.cap = "Figure 1: Box diagram for transitions.", anchor = "Figure"----------------------------------
grViz("
digraph boxes_and_circles {

  # a 'graph' statement
  graph [overlap = true, fontsize = 10]

  # 'node' statements
  node [shape = box, fontname = Helvetica]
  'A) Alive'; Death; 'B) Alive'; 'Other death'; 'Oral ca. death'
  
  # 'edge' statements
  'A) Alive' -> Death 
  'B) Alive' -> 'Other death' 'B) Alive' ->'Oral ca. death'
}
")


## ----table events-------------------------------------------------------------------------------------------------------
table(orca$event)
orca <- mutate(orca, all = event != "Alive")
table(orca$all)


## ----follow-up, fig.cap="Figure 2: Possible representations of follow-up time.", anchor = "Figure"----------------------
ggplotly(
  orca %>%
    mutate(
      text = paste("Subject ID = ", id, "<br>", "Time = ", time, "<br>", "Event = ",  
                   event, "<br>", "Age = ", round(age, 2), "<br>", "Stage = ", stage)
    ) %>%
  ggplot(aes(x = id, y = time, text = text)) +
    geom_linerange(aes(ymin = 0, ymax = time)) +
    geom_point(aes(shape = event, color = event), stroke = 1, cex = 2) +
    scale_shape_manual(values = c(1, 3, 4)) +
    labs(y = "Time (years)", x = "Subject ID") + coord_flip(),
  tooltip = "text"
)


## ----time-scales, fig.width=10------------------------------------------------------------------------------------------
wrap_plots(
  ggplot(orca, aes(x = id, y = time)) +
    geom_linerange(aes(ymin = 0, ymax = time)) +
    geom_point(aes(shape = event, color = event), stroke = 1, cex = 2) +
    scale_shape_manual(values = c(1, 3, 4)) +
    labs(y = "Time (years)", x = "Subject ID") + coord_flip(),
  orca %>%
    mutate(age_orig = age,
           age_end = age + time) %>%
    ggplot(aes(x = id, y = age_end)) +
    geom_linerange(aes(ymin = age_orig, ymax = age_end)) +
    geom_point(aes(shape = event, color = event), stroke = 1, cex = 2) +
    scale_shape_manual(values = c(1, 3, 4)) + 
    labs(y = "Age (years)", x = "Subject ID") + coord_flip()
)


## ----survObj------------------------------------------------------------------------------------------------------------
su_obj <- Surv(orca$time, orca$all)
str(su_obj)


## ----KMmethod-----------------------------------------------------------------------------------------------------------
fit_km <- survfit(Surv(time, all) ~ 1, data = orca)
print(fit_km, print.rmean = TRUE)


## ----suvTable KM--------------------------------------------------------------------------------------------------------
dat_km <- tidy(fit_km)
head(dat_km)


## ----plot KM------------------------------------------------------------------------------------------------------------
ggsurvfit(fit_km) + 
  add_censor_mark() +
  add_confidence_interval() +
  add_risktable() +
  labs(x = "Time (years)")


## ----alternative plots KM, fig.width = 12, results='hide'---------------------------------------------------------------
wrap_plots(
  ggsurvfit(fit_km, type = "risk") + add_confidence_interval(),
  ggsurvfit(fit_km, type = "cumhaz") + add_confidence_interval()
)


## ----aggregating data---------------------------------------------------------------------------------------------------
cuts <- seq(0, 23, 1)
lifetab_dat <- orca %>%
  mutate(time_cat = cut(time, cuts)) %>%
  group_by(time_cat) %>%
  summarise(nlost = sum(all == 0),
            nevent = sum(all == 1))


## ----suvTable lifetable-------------------------------------------------------------------------------------------------
dat_lt <- with(lifetab_dat, lifetab(tis = cuts, ninit = nrow(orca), 
                                    nlost = nlost, nevent = nevent))
round(dat_lt, 4)


## ----suvTable FH--------------------------------------------------------------------------------------------------------
fit_fh <- survfit(Surv(time, all) ~ 1, data = orca, type = "fleming-harrington", conf.type = "log-log")
dat_fh <- tidy(fit_fh)
## for the Nelson-Aalen estimator of the cumulative hazard
#dat_fh <- tidy_survfit(fit_fh, type = "cumhaz")
head(dat_fh)


## ----comaring Surv estimates--------------------------------------------------------------------------------------------
ggplotly(
  ggplot() +
  geom_step(data = dat_km, aes(x = time, y = estimate, colour = "K-M")) +
  geom_step(data = dat_fh, aes(x = time, y = estimate, colour = "N-A")) +
  geom_step(data = dat_lt, aes(x = cuts[-length(cuts)], y = surv, colour = "LT")) +
  labs(x = "Time (years)", y = "Survival", colour = "Estimator")
)


## ----measure of central tendency----------------------------------------------------------------------------------------
(mc <- data.frame(q = c(.25, .5, .75),
                 km = quantile(fit_km),
                 fh = quantile(fit_fh)))


## ----surv quantiles plot------------------------------------------------------------------------------------------------
ggsurvfit(fit_km) +
  add_quantile(y_value = 0.75, linetype = "dotted") +
  add_quantile(y_value = 0.5, linetype = "dotted") +
  add_quantile(y_value = 0.25, linetype = "dotted") +
  scale_ggsurvfit() + labs(x = "Time (years)")


## ----parametric surv----------------------------------------------------------------------------------------------------
fit_exp <- flexsurvreg(Surv(time, all) ~ 1, data = orca, dist = "exponential")
fit_exp
fit_w <- flexsurvreg(Surv(time, all) ~ 1, data = orca, dist = "weibull")
fit_ll <- flexsurvreg(Surv(time, all) ~ 1, data = orca, dist = "llogis")
fit_sp <- flexsurvspline(Surv(time, all) ~ 1, data = orca, k = 1, scale = "odds")


## ----ggflexsurvplot-----------------------------------------------------------------------------------------------------
ggflexsurvplot(fit_exp, xlab = "Time (years)", censor = F)


## ----comparison parametric surv, fig.width=12---------------------------------------------------------------------------
wrap_plots(
  ggplot(data.frame(summary(fit_exp)), aes(x = time)) + 
    geom_line(aes(y = est, col = "Exponential")) +
    geom_line(data = data.frame(summary(fit_w)), aes(y = est, col = "Weibull")) +
    geom_line(data = data.frame(summary(fit_ll)), aes(y = est, col = "Log-logistic")) +
    geom_line(data = data.frame(summary(fit_sp)), aes(y = est, col = "Flex splines")) +
    labs(x = "Time (years)", y = "Survival", col = "Distributions"),
  ggplot(data.frame(summary(fit_exp, type = "hazard")), aes(x = time)) + 
    geom_line(aes(y = est, col = "Exponential")) +
    geom_line(data = data.frame(summary(fit_w, type = "hazard")), aes(y = est, col = "Weibull")) +
    geom_line(data = data.frame(summary(fit_ll, type = "hazard")), aes(y = est, col = "Log-logistic")) +
    geom_line(data = data.frame(summary(fit_sp, type = "hazard")), aes(y = est, col = "Flex splines")) +
    labs(x = "Time (years)", y = "Hazard", col = "Distributions"),
  ncol = 2, guides = "collect"
)


## ----incidence rates----------------------------------------------------------------------------------------------------
#ci.exp(glm(all ~ 0 + stage, data = orca, family = "poisson", offset = log(time)))
group_by(orca, stage) %>%
  summarise(
    D = sum(all),
    Y = sum(time)
  ) %>%
  cbind(pois.approx(x = .$D, pt = .$Y))


## ----comp surv stages---------------------------------------------------------------------------------------------------
su_stg  <- survfit(Surv(time, all) ~ stage, data = orca)
su_stg


## ----KM stages----------------------------------------------------------------------------------------------------------
ggsurvfit(su_stg) + 
  add_confidence_interval() +
  add_risktable() +
  labs(x = "Time (years)")


## ----surv table stages--------------------------------------------------------------------------------------------------
lifetab_stg <- tidy(su_stg)
lifetab_stg %>%
  group_by(strata) %>%
  do(head(., n = 3))


## ----cumhaz surv stages, fig.width = 12---------------------------------------------------------------------------------
wrap_plots(
  ggsurvfit(su_stg, type = "risk") + add_confidence_interval(),
  ggsurvfit(su_stg, type = "cumhaz") + add_confidence_interval()
)


## ----M-H logrank test---------------------------------------------------------------------------------------------------
survdiff(Surv(time, all) ~ stage, data = orca)


## ----Peto & Peto test---------------------------------------------------------------------------------------------------
survdiff(Surv(time, all) ~ stage, data = orca, rho = 1)


## ----cox model----------------------------------------------------------------------------------------------------------
m1 <- coxph(Surv(time, all) ~ sex + I((age-65)/10) + stage, data = orca)
summary(m1)


## ----test ph------------------------------------------------------------------------------------------------------------
cox.zph.m1 <- cox.zph(m1)
cox.zph.m1


## ----plot_ph, fig.height=12---------------------------------------------------------------------------------------------
ggcoxzph(cox.zph.m1)


## ----plot_diagnostics, fig.height=6, message=FALSE----------------------------------------------------------------------
ggcoxdiagnostics(m1, type = "dfbeta", linear.predictions = FALSE)

## ----plot_diagnostics2, fig.height=6, message=FALSE---------------------------------------------------------------------
ggcoxdiagnostics(m1, type = "deviance", linear.predictions = FALSE)


## ----plot_martingale----------------------------------------------------------------------------------------------------
cox_martingale <- coxph(Surv(time, all) ~ log(age) + sqrt(age), data = orca)
ggcoxfunctional(cox_martingale, data = orca)


## ----updated cox--------------------------------------------------------------------------------------------------------
orca2 <- orca %>%
  filter(stage != "unkn") %>%
  mutate(st3 = Relevel(droplevels(stage), list(1:2, 3, 4)))
m2 <- coxph(Surv(time, all) ~ sex + I((age-65)/10) + st3, data = orca2, ties = "breslow")
round(ci.exp(m2), 4)


## ----forest plot, fig.width=12------------------------------------------------------------------------------------------
wrap_plots(
  ggforest(m1, data = orca),
  ggforest(m2, data = orca2),
  ncol = 2
)


## ----newdata for predictions--------------------------------------------------------------------------------------------
newd <- expand.grid(sex = c("Male", "Female"), age = c(40, 80), st3 = levels(orca2$st3))
newd$id <- 1:12
newd


## ----survival from cox, fig.width=12------------------------------------------------------------------------------------
surv_summary(survfit(m2, newdata = newd)) %>%
  merge(newd, by.x = "strata", by.y = "id") %>%
  ggplot(aes(x = time, y = surv, col = sex, linetype = factor(age))) +
  geom_step() + facet_grid(. ~ st3) +
  labs(x = "Time (years)", y = "Survival probability", linetype = "Age (years")


## ----weibull model------------------------------------------------------------------------------------------------------
m2w <- flexsurvreg(Surv(time, all) ~ sex + I((age-65)/10) + st3, data = orca2, dist = "weibull")
m2w


## -----------------------------------------------------------------------------------------------------------------------
m2wph <- weibreg(Surv(time, all) ~ sex + I((age-65)/10) + st3, data = orca2)
summary(m2wph)


## ----median surv weibull------------------------------------------------------------------------------------------------
median.weibull <- function(shape, scale) qweibull(0.5, shape = shape, scale = scale)
set.seed(2153)
newd <- data.frame(sex = c("Male", "Female"), age = 65, st3 = "I+II")
summary(m2w, newdata = newd, fn = median.weibull, t = 1, B = 10000)


## ----median surv cox----------------------------------------------------------------------------------------------------
survfit(m2, newdata = newd)


## -----------------------------------------------------------------------------------------------------------------------
cuts <- sort(unique(orca2$time[orca2$all == 1]))
orca_splitted <- survSplit(Surv(time, all) ~ ., data = orca2, cut = cuts, episode = "tgroup")
head(orca_splitted, 15)


## ---- conditional poisson-----------------------------------------------------------------------------------------------
mod_poi <- gnm(all ~ sex + I((age-65)/10) + st3, data = orca_splitted, 
               family = poisson, eliminate = factor(time))
summary(mod_poi)


## ----comp poisson cox---------------------------------------------------------------------------------------------------
round(data.frame(cox = ci.exp(m2), poisson = ci.exp(mod_poi)), 4)


## ----poisson basehazard (takes time), echo=c(1:3), eval=-c(2, 3), eval=FALSE--------------------------------------------
orca_splitted$dur <- with(orca_splitted, time - tstart)
mod_poi2 <- glm(all ~ -1 + factor(time) + sex + I((age-65)/10) + st3, 
                data = orca_splitted, family = poisson, offset = log(dur))


## ----loading from url, echo=FALSE---------------------------------------------------------------------------------------
orca_splitted$dur <- with(orca_splitted, time - tstart)
#save(mod_poi2, file = "data/mod_poi2.Rdata")
load(url("http://www.stats4life.se/data/mod_poi2.Rdata"))


## ----step hazard, warning=FALSE-----------------------------------------------------------------------------------------
newd <- data.frame(time = cuts, dur = 1,
                   sex = "Female", age = 65, st3 = "I+II")
blhaz <- data.frame(ci.pred(mod_poi2, newdata = newd))
ggplot(blhaz, aes(x = c(0, cuts[-138]), y = Estimate, xend = cuts, yend = Estimate)) + geom_segment() +
  scale_y_continuous(trans = "log", limits = c(.05, 5), breaks = pretty_breaks()) +
  labs(x = "Time (years)", y = "Baseline hazard")


## ----flexible hazard----------------------------------------------------------------------------------------------------
k <- quantile(orca2$time, 1:5/6)
mod_poi2s <- glm(all ~ ns(time, knots = k) + sex + I((age-65)/10) + st3, 
                 data = orca_splitted, family = poisson, offset = log(dur))
round(ci.exp(mod_poi2s), 3)
blhazs <- data.frame(ci.pred(mod_poi2s, newdata = newd))
ggplot(blhazs, aes(x = newd$time, y = Estimate)) + geom_line() +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = .2) +
  scale_y_continuous(trans = "log", breaks = pretty_breaks()) +
  labs(x = "Time (years)", y = "Baseline hazard")


## ----predicted survivals------------------------------------------------------------------------------------------------
newd <- data.frame(sex = "Female", age = 65, st3 = "I+II")
surv_cox <- tidy(survfit(m2, newdata = newd))
surv_weibull <- summary(m2w, newdata = newd, tidy = TRUE)
## For the poisson model we need some extra steps
tbmid <- sort(unique(.5*(orca_splitted$tstart + orca_splitted$time)))
mat <- cbind(1, ns(tbmid, knots = k), 0, 0, 0, 0)
Lambda <- ci.cum(mod_poi2s, ctr.mat = mat, intl = diff(c(0, tbmid)))
surv_poisson <- data.frame(exp(-Lambda))


## ----plot predicted survivals-------------------------------------------------------------------------------------------
ggplot(surv_cox, aes(time, estimate)) + geom_step(aes(col = "Cox")) +
  geom_line(data = surv_weibull, aes(y = est, col = "Weibull")) +
  geom_line(data = surv_poisson, aes(x = c(0, tbmid[-1]), y = Estimate, col = "Poisson")) +
  labs(x = "Time (years)", y = "Survival", col = "Models")


## ---- eval = FALSE------------------------------------------------------------------------------------------------------
library(shiny)
library(tidyverse)
library(splines)
library(Epi)
library(survival)
library(flexsurv)

shinyApp(
  
  ui = fluidPage(
    h2("Choose covariate pattern:"),
    selectInput("sex", label = h3("Sex"), 
                choices = list("Female" = "Female" , "Male" = "Male")),
    sliderInput("age", label = h3("Age"),
                min = 20, max = 80, value = 65),
    selectInput("st3", label = h3("Stage (3 levels)"), 
                choices = list("Stage I and II" = "I+II", "Stage III" = "III",
                               "Stage IV" = "IV")),
    plotOutput("survPlot")
  ),
  
  server = function(input, output){
    
    orca <- read.table("http://www.stats4life.se/data/oralca.txt", header = T,
                       stringsAsFactors = TRUE)
    orca2 <- orca %>%
      filter(stage != "unkn") %>%
      mutate(st3 = Relevel(droplevels(stage), list(1:2, 3, 4)),
             all = 1*(event != "Alive"))
    m2 <- coxph(Surv(time, all) ~ sex + I((age-65)/10) + st3, data = orca2, ties = "breslow")
    m2w <- flexsurvreg(Surv(time, all) ~ sex + I((age-65)/10) + st3, data = orca2, dist = "weibull")
    cuts <- sort(unique(orca2$time[orca2$all == 1]))
    orca_splitted <- survSplit(Surv(time, all) ~ ., data = orca2, cut = cuts, episode = "tgroup")
    orca_splitted$dur <- with(orca_splitted, time - tstart)
    k <- quantile(orca2$time, 1:5/6)
    mod_poi2s <- glm(all ~ ns(time, knots = k) + sex + I((age-65)/10) + st3, 
                     data = orca_splitted, family = poisson, offset = log(dur))
    
    newd <- reactive({
      data.frame(sex = input$sex, age = input$age, st3 = input$st3)      
    })
    
    output$survPlot <- renderPlot({
      newd <- newd()
      surv_cox <- tidy(survfit(m2, newdata = newd))
      surv_weibull <- summary(m2w, newdata = newd, tidy = TRUE)
      tbmid <- sort(unique(.5*(orca_splitted$tstart + orca_splitted$time)))
      mat <- cbind(1, ns(tbmid, knots = k), 1*(input$sex == "Male"), (input$age - 65)/10, 
                   1*(input$st3 == "III"), 1*(input$st3 == "IV"))
      Lambda <- ci.cum(mod_poi2s, ctr.mat = mat, intl = diff(c(0, tbmid)))
      surv_poisson <- data.frame(exp(-Lambda))
      
      ggplot(surv_cox, aes(time, estimate)) + geom_step(aes(col = "Cox")) +
        geom_line(data = surv_weibull, aes(y = est, col = "Weibull")) +
        geom_line(data = surv_poisson, aes(x = c(0, tbmid[-1]), y = Estimate, col = "Poisson")) +
        labs(x = "Time (years)", y = "Survival", col = "Models")
    })
  }
)


## ----quadratic age------------------------------------------------------------------------------------------------------
m3 <- coxph(Surv(time, all) ~ sex + I(age-65) + I((age-65)^2) + st3, data = orca2)
summary(m3)


## ----plot non-linear HR-------------------------------------------------------------------------------------------------
age <- seq(20, 80, 1) - 65
hrtab <- ci.exp(m3, ctr.mat = cbind(0, age, age^2, 0, 0))
ggplot(data.frame(hrtab), aes(x = age+65, y = exp.Est.., ymin = X2.5., ymax = X97.5.)) +
  geom_line() + geom_ribbon(alpha = .1) +
  scale_y_continuous(trans = "log", breaks =  pretty_breaks()) +
  labs(x = "Age (years)", y = "Hazard ratio") +
  geom_vline(xintercept = 65, lty = 2) + geom_hline(yintercept = 1, lty = 2)


## ----plot zph for sex---------------------------------------------------------------------------------------------------
plot(cox.zph.m1[1])
abline(h = m1$coef[1], col = 2, lty = 2, lwd = 2)


## ----time-dependent coef------------------------------------------------------------------------------------------------
orca3 <- survSplit(Surv(time, all) ~ ., data = orca2, cut = c(5, 15), episode = "tgroup")
head(orca3)
m3 <- coxph(Surv(tstart, time, all) ~ relevel(sex, 2):strata(tgroup) + I((age-65)/10) + st3, data = orca3)
m3


## -----------------------------------------------------------------------------------------------------------------------
p <- .25
fit_qs <- ctqr(Surv(time, all) ~ st3, p = p,  data = orca2)
fit_qs


## -----------------------------------------------------------------------------------------------------------------------
fit_kms <- survfit(Surv(time, all) ~ st3, data = orca2)
pc_qs <- tibble(
  st3 = levels(orca2$st3),
  coef = coef(fit_qs),
  time_p = c(NA, coef[-1] + coef[1]),
  time_ref = coef[1],
  p = c(p, p - .005, p + .005)
)[-1, ]

surv_summary(fit_kms, data = orca2) %>%
  ggplot(aes(time, surv, col = st3)) +
  geom_step() +
  geom_segment(data = pc_qs, aes(x = time_p, y = 1 - p, xend = time_ref, 
                                 yend = 1 - p))


## -----------------------------------------------------------------------------------------------------------------------
fit_q <- ctqr(Surv(time, all) ~ sex + I((age-65)/10) + st3, p = seq(.1, .7, .1), 
              data = orca2)
fit_q


## ----plots percentiles--------------------------------------------------------------------------------------------------
coef_q <- data.frame(coef(fit_q)) %>%
  rownames_to_column(var = "coef") %>%
  gather(p, beta, -coef) %>%
  mutate(
    p = as.double(gsub("p...", "", p)),
    se = c(sapply(fit_q$covar, function(x) sqrt(diag(x)))),
    low = beta - 1.96 * se,
    up = beta + 1.96 * se
    )
ggplot(coef_q, aes(p, beta)) +
  geom_line() +
  geom_ribbon(aes(ymin = low, ymax = up), alpha = .1) +
  facet_wrap(~ coef, nrow = 2, scales = "free")


## -----------------------------------------------------------------------------------------------------------------------
expand.grid(sex = c("Male", "Female"), age = c(40, 80), st3 = levels(orca2$st3)) %>% 
  cbind(round(predict(fit_q, newdata = .), 3))


## ----CIF----------------------------------------------------------------------------------------------------------------
cif <- Cuminc(time = "time", status = "event", data = orca2)
head(cif)


## ----CIF and survival plot----------------------------------------------------------------------------------------------
ggplot(cif, aes(time)) +
  geom_step(aes(y = `CI.Oral ca. death`, colour = "Cancer death")) +
  geom_step(aes(y = `CI.Oral ca. death`+ `CI.Other death`, colour = "Total")) +
  geom_step(aes(y = Surv, colour = "Overall survival")) + 
  labs(x = "Time (years)", y = "Proportion", colour = "")


## ----CIF st3------------------------------------------------------------------------------------------------------------
cif_stage <- Cuminc(time = "time", status = "event", group = "st3", data = orca2)
cif_stage %>% 
  group_by(group) %>%
  do(head(., n = 3))


## ----plot CIF st3-------------------------------------------------------------------------------------------------------
wrap_plots(
  ggplot(cif_stage, aes(time)) +
    geom_step(aes(y = `CI.Oral ca. death`, colour = group)) +
    labs(x = "Time (years)", y = "Proportion", title = "Cancer death by stage"),
  ggplot(cif_stage, aes(time)) +
    geom_step(aes(y = `CI.Other death`, colour = group)) +
    labs(x = "Time (years)", y = "Proportion", title = "Other deaths by stage"),
  ncol = 2, guides = "collect"
)


## ----cox competing------------------------------------------------------------------------------------------------------
m2haz1 <- coxph(Surv(time, event == "Oral ca. death") ~ sex + I((age-65)/10) + st3, data = orca2)
round(ci.exp(m2haz1), 4)
m2haz2 <- coxph(Surv(time, event == "Other death") ~ sex + I((age-65)/10) + st3, data = orca2)
round(ci.exp(m2haz2), 4)


## ----Fine–Gray model----------------------------------------------------------------------------------------------------
m2fg1 <- with(orca2, crr(time, event, cov1 = model.matrix(m2), failcode = "Oral ca. death"))
summary(m2fg1, Exp = T)
m2fg2 <- with(orca2, crr(time, event, cov1 = model.matrix(m2), failcode = "Other death"))
summary(m2fg2, Exp = T)

