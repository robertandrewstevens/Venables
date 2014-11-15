# R Workshop Session 3:  Generalized Linear and Generalized Additive Models
## Bill Venables, CSIRO, Australia
## UseR! 2012, Nashville, 11 June, 2012

# Contents
## 1 An example from MASS: low birth weight
### 1.1 Automated screening of variables
### 1.2 An extended model with smooth terms
### 1.3 Looking at the terms
### 1.4 A helper function: the most frequent value
### 1.5 The main two-way interaction
## 2 Tiger prawn species split
### 2.1 An initial GLM
### 2.2 A long-term trend?
### 2.3 A working GAM with new technology
### 2.4 The spatio-temporal effect
## 3 Technical highlights
## References
## Session information

# 1 An example from MASS: low birth weight
## From Venables and Ripley (2002, Chap. 7).
### low:  indicator of birth weight less than 2.5 kg.
### age:  mother's age in years.
### lwt:  mother's weight in pounds at last menstrual period.
### race:  mother's race ('1' = white, '2' = black, '3' = other).
### smoke:  smoking status during pregnancy.
### ptl:  number of previous premature labours.
### ht:  history of hypertension.
### ui:  presence of uterine irritability.
### ftv:  number of physician visits during the rst trimester.
### bwt:  birth weight in grams.

### The original MASS code:

library(MASS)
attach(birthwt)
race <- factor(race, labels = c("white", "black", "other"))
table(ptl)
ptd <- factor(ptl > 0)
table(ftv)
ftv <- factor(ftv)
levels(ftv)[-(1:2)] <- "2+"
table(ftv) # as a check
bwt <- data.frame(low = factor(low), age, lwt, race,
                  smoke = (smoke > 0), ptd, ht = (ht > 0), ui = (ui > 0), ftv)
detach() 
rm(race, ptd, ftv)

### My preference now:

suppressPackageStartupMessages(library(SOAR)) # picky!
suppressPackageStartupMessages(library(MASS))
BirthWt <- within(birthwt, {
  race <- factor(race, labels = c("white", "black", "other"))
  ptl <- ptl > 0
  ftv <- factor(ftv)
  levels(ftv)[-(1:2)] <- "2+"
  low <- factor(low, labels = c("normal", "low"))
  smoke <- (smoke > 0)
  ht <- (ht > 0)
  ui <- (ui > 0)
  bwt <- NULL ## remove actual birth weight
})
Store(BirthWt)
head(BirthWt, 2)

### Advice from the fortunes package:
###   If I were to be treated by a cure created by stepwise regression,
###   I would prefer voodoo.
###   - Dieter Menne (in a thread about regressions with many
###     variables) R-help (October 2009)
### * Automated screening is more defensible in cases of pure prediction.
### * Automated screening is dangerous if used for inference.
### Caveat emptor!
### We reduce some of the clutter with:
options(show.signif.stars = FALSE)
stepAIC <- function(..., trace = FALSE) ## change default
  MASS::stepAIC(..., trace = trace)
dropterm <- function(..., sorted = TRUE) ## change default
  MASS::dropterm(..., sorted = sorted)

## 1.1 Automated screening of variables

### A starting point, main effects only:

BW0 <- glm(low ~ ., binomial, BirthWt)
dropterm(BW0, test = "Chisq")

### Screen for possible interactions:

sBW0 <- stepAIC(BW0, scope = list(lower = ~1, upper = ~.^2))
dropterm(sBW0, test = "Chisq")

## 1.2 An extended model with smooth terms

### We consider some flexibility in the age term and its interaction with ftv.

suppressPackageStartupMessages(require(mgcv))
BW1 <- gam(low ~ smoke*ui + ht + s(lwt) + ptl + s(age) +
           poly(age, 2)*ftv, family = binomial, data = BirthWt)
anova(BW1)

## 1.3 Looking at the terms

### First the“main effect”terms. This is a bit tricky...

(nam <- names(model.frame(BW1)))
layout(matrix(1:8, 2, 4, byrow = TRUE))
termplot(BW1, terms = nam[2:7], se=TRUE) ## fixed
plot(BW1)                                ## smooth

## The first nam is the response and the 8th and 9th refer to the smooth terms. termplot can handle non-smoothed terms, but smooth terms must be handled by the plot method for gam objects.

## 1.4 A helper function: the most frequent value

mostFreq <- function(x, ...) UseMethod("mostFreq")
mostFreq.numeric <- stats::median.default  ## check argument names
mostFreq.logical <- function(x, ...) {
  tx <- as.vector(table(x))
  tx[2] > tx[1]
}
mostFreq.character <- function(x, ...) {
  tx <- table(x)
  names(tx)[which.max(tx)]
}
mostFreq.factor <- function(x, ...)
  mostFreq.character(as.character(x))
Store(list = ls(pattern = "^mostFreq"))

## 1.5 The main two-way interaction

### Predict the probability of low birth weight with varying age and ftv, and other variables at or near their modal value.

all.vars(formula(BW1))
pBirthWt <- with(BirthWt,
                 expand.grid(smoke = mostFreq(smoke), ui = mostFreq(ui),
                             ht = mostFreq(ht), lwt = mostFreq(lwt),
                             ptl = mostFreq(ptl), age = min(age):max(age),
                             ftv = levels(ftv)))
pBirthWt$pBW1 <- predict(BW1, pBirthWt, type = "response")
library(lattice)
(ageFtv <- xyplot(pBW1 ~ age|ftv, pBirthWt, layout = c(3,1),
                  type = "l", ylab = "Pr(low birth wt)",
                  ylim = 0:1, aspect = 1))

### To bring the predictions closer to the actual data, confine the predictions to age ranges that apply within the levels of ftv

pBirthWt <- within(pBirthWt, {
  rngs <- do.call(cbind, with(BirthWt, tapply(age, ftv, range)))
  pBW1a <- pBW1
  is.na(pBW1a[age < rngs[1, ftv] | age > rngs[2, ftv]]) <- TRUE
  rm(rngs)
})
(ageFtvA <- xyplot(pBW1a ~ age|ftv, pBirthWt, layout = c(3,1),
                   type = "l", ylab = "Pr(low birth wt)",
                   ylim = 0:1, aspect = 1))

## 2 Tiger prawn species split

### The Northern Prawn Fishery: Tiger prawn effort and “the blue line”

### Background
### * Two species of Tiger prawns are caught together.
### * Both species require separate Stock Assessment.
### * The assessment model requires catches of Tiger prawns to be split (by weight)
### * Problem: Build a model for partitioning catches into the two component species.
### * Data: independent surveys (12 in all) where catches have been split into the two species, Penaeus semisulcatus (Grooved) and P. esculentus (Brown).
### * Both species have annual offshore migration patterns.

### Variables available:
### Response: Psem, Pesc, (Total = Psem + Pesc) in gms;
### Predictors:
### * Longitude, Latitude – of trawl shot
### * Coast, Sea – alternative spatial coordinates
### * Depth – of trawl shot
### * Mud – the % mud in the substrate
### * DayOfYear – to allow for annual migration periodicity
### * ElapsedDays – days since 1970-01-01, for long term trend
### * Survey – used for a random effect extension

### Strategy:
### * Build a simpler GLM using mainly splines, with a term in DayOfYear and Sea to allow for temporal (annual migration) effects
### * Develop a more sophisticated GAM to take advantage of more recent modelling technology
### * Look at a long-term trend term as a perturbation to the model
### * Consider GLMMs with random terms for Survey, eventually

### Model terms, GLM:
### * Spline in Coast surrogate for large-scale benthic changes
### * Splines in Sea, Depth and Mud – more local spatial effects
### * Periodic term in DayOfYear and its interaction with Sea – annual migration effects
### * Spline in ElapsedDays – testing for long-term stability

## 2.1 An initial GLM

### The GLM fitting process is slow to converge under the normal algorithm. Two possible alternatives:
### * Use the glm2 library, which has a modified convergence process, (Marschner, 2011).
### * In this case the problem is with the variable weights needed, so fit a model ignoring weights and use the linear predictor as a starting value for the weighted fit.

### The periodic terms will use Fourier polynomials:

Annual <- function (day, k = 4) {  ## day of the year, starting from 0
  theta <- 2*base::pi*day/364.25
  X <- matrix(0, length(theta), 2 * k)
  nam <- as.vector(outer(c("c", "s"), 1:k, paste, sep = ""))
  dimnames(X) <- list(names(day), nam)
  m <- 0
  for (j in 1:k) {
    X[, (m <- m + 1)] <- cos(j * theta)
    X[, (m <- m + 1)] <- sin(j * theta)
  }
  X
}

library(splines)

Tigers <- read.csv("Tigers.csv", header = TRUE, sep = ";")

temp <- glm(Psem/Total ~ ns(Coast, 10) + ns(Sea, 5) + ns(Depth, 5) +
            ns(Mud, 4) + Annual(DayOfYear, 4)*Sea, family = quasibinomial,
            data = Tigers, trace = TRUE)  ## unweighted
Tigers$eta <- predict(temp)
TModelGLM <- update(temp, etastart = eta, weights = Total)
rm(temp)
Tigers$eta <- predict(TModelGLM)
TModelGLM$call$trace <- NULL  ## for future updating
Store(TModelGLM)
(nam <- names(model.frame(TModelGLM)))  ## for term plotting
nam <- nam[2:7]  ## terms to plot

### Look at the shape of the main effect terms, to see implications:

layout(matrix(1:6, 2, 3, byrow = TRUE))   ## 2 x 3 array of plots
termplot(TModelGLM, terms = nam, se = TRUE, rug = TRUE)

## 2.2 A long-term trend?

### The stability of species ratios over time is important. We can check for this by including a spline term in ElapsedDays:

TM2 <- update(TModelGLM, . ~ . + ns(ElapsedDays, 7))
anova(TModelGLM, TM2, test = "F")

### Significant, but is it important?

layout(matrix(1, 1, 1, byrow = TRUE))   ## return to 1 plot per screen/page
termplot(TM2, terms = "ns(ElapsedDays, 7)", se = TRUE, rug = TRUE)

## 2.3 A working GAM with new technology

### The mgcv package represents a major advance in smooth model fitting technology in sevaral respects, including
### Smoothed terms in multiple predictors can now be handled
### * A wide variety of basis functions is available, including e.g. thin plate splines, cyclic spline bases, &c
### * A powerful visualisation tool in projections of predictor variable space is avialable in addition to tools for inspection of individual terms
### The price is:
### * The package is still under development and new versions are fairly common (though becoming less so)
### * The implementation is to some extent non-standard R

### The working model:

require(mgcv)
Attach()
TModelGAM <- gam(Psem/Total ~ s(Longitude, Latitude) +
                 te(DayOfYear, Sea, k = c(5, 5), bs = c("cc", "cs")) +
                 te(DayOfYear, Depth, k = c(5, 5), bs = c("cc", "cs")) +
                 te(Sea, Depth, k = c(5, 5), bs = "cs") +
                 s(Mud, k = 5),
                 family = quasibinomial, data = Tigers,
                 knots = list(DayOfYear = seq(0, 365.25, length=5)),
                 weights = Total, control = gam.control(trace = TRUE))
TModelGAM_NS <- update(TModelGAM, . ~ . + s(ElapsedDays))  ## non-stationary


# Session information

sessionInfo()