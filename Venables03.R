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

# Session information

sessionInfo()
