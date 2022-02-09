if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, ggplot2, dplyr, lubridate, stringr, readxl, data.table, 
               gdata, MatchIt, cobalt)

install.packages('Matching')

## Simulate data
set.seed(12345678)
n <- 10000
select.dat <- tibble(
  x = rnorm(n, .2, 1),
  z = rnorm(n, -.1, 1),
  d = (x + rnorm(n,0,1) >0),
  y0 = -2.5 + 1.5*x + rnorm(n, 0, 1),
  y1 = y0 + 4,
  y = y1*d + y0*(1-d),
  d_alt = ( x + z + rnorm(n,0,1) > 0),
  y0_alt = -2.5 + 1.5*x + 2.25*z + rnorm(n, 0, 1),
  y1_alt = y0_alt+4,
  y_alt = y1_alt*d_alt + y0_alt*(1-d_alt)
)
  
## nearest neighbor matching with euclidean distance weights
nn.est1 <- Matching::Match(Y=select.dat$y,
                 Tr=select.dat$d,
                 X=select.dat$x,
                 M=1,
                 Weight=1,
                 estimand="ATE")
summary(nn.est1)


## nearest neighbor matching with mahalanobis weights
nn.est2 <- Matching::Match(Y=select.dat$y,
                           Tr=select.dat$d,
                           X=select.dat$x,
                           M=1,
                           Weight=2,
                           estimand="ATE")
summary(nn.est2)

## regression
reg1.dat <- select.dat %>% filter(d==1)
reg1 <- lm(y ~ x, data=reg1.dat)

reg0.dat <- select.dat %>% filter(d==0)
reg0 <- lm(y ~ x, data=reg0.dat)
pred1 <- predict(reg1,new=select.dat)
pred0 <- predict(reg0,new=select.dat)
mean(pred1-pred0)


## violation of selection on observables
## nearest neighbor
nn.est3 <- Matching::Match(Y=select.dat$y_alt,
                           Tr=select.dat$d_alt,
                           X=select.dat$x,
                           M=1,
                           Weight=2,
                           estimand="ATE")
summary(nn.est3)

## regression
reg1.dat <- select.dat %>% filter(d_alt==1)
reg1 <- lm(y_alt ~ x, data=reg1.dat)

reg0.dat <- select.dat %>% filter(d_alt==0)
reg0 <- lm(y_alt ~ x, data=reg0.dat)
pred1_alt <- predict(reg1,new=select.dat)
pred0_alt <- predict(reg0,new=select.dat)
mean(pred1_alt-pred0_alt)


## restoring selection on observables
match.mat <- select.dat %>% select(x,z)
nn.est4 <- Matching::Match(Y=select.dat$y_alt,
                           Tr=select.dat$d_alt,
                           X=match.mat,
                           M=1,
                           Weight=2,
                           estimand="ATE")
summary(nn.est4)

## regression
reg1.dat <- select.dat %>% filter(d_alt==1)
reg1 <- lm(y_alt ~ x + z, data=reg1.dat)

reg0.dat <- select.dat %>% filter(d_alt==0)
reg0 <- lm(y_alt ~ x + z, data=reg0.dat)
pred1_alt <- predict(reg1,new=select.dat)
pred0_alt <- predict(reg0,new=select.dat)
mean(pred1_alt-pred0_alt)

## regression in a single step
select.dat <- select.dat %>%
  mutate(xbar=mean(x),
         zbar=mean(z))
reg <-  lm(y_alt ~ d_alt + x + z + d_alt*(x-xbar) + d_alt*(z-zbar), data=select.dat)
summary(reg)


##propensity score
logit.reg <- glm(d_alt ~ x+z,
         data = select.dat, family = binomial(link = 'logit'))
select.dat <- select.dat %>%
  mutate(ps = predict(logit.reg, type = 'response')) %>%
  filter(ps>0 & ps<1)


# Create IPW weights
select.dat <- select.dat %>%
  mutate(ipw = case_when(
    d_alt == 1 ~ 1/ps,
    d_alt == 0 ~ 1/(1-ps)
    ),
    out_weight=y_alt*ipw
  )

select.dat %>% group_by(d_alt) %>% summarize(out_weight)

## Assessing balance
baseline.match <- matchit(d_alt~x+z, data=select.dat, method=NULL, distance="mahalanobis")
plot(summary(baseline.match))

update.match <- matchit(d_alt~x+z, data=select.dat, method="nearest", distance="mahalanobis", replace=TRUE)
match.mod <- summary(update.match)
library(cobalt)
bal.tab(update.match)
love.plot(update.match)
bal.plot(update.match, var.name="x", which="both")
bal.plot(update.match, var.name="z", which="both")
