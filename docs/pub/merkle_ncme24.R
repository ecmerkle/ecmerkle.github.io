#################################################################################
## PACKAGES
#################################################################################
library(lavaan)

## install github version of blavaan:
## remotes::install_github("ecmerkle/blavaan", INSTALL_opts = "--no-multiarch")
library(blavaan)

library(ggplot2)

library(bayesplot)

options(width = 90)



#################################################################################
## CFA Section
#################################################################################
data("data.pisaMath", package = "sirt")

dat <- data.pisaMath$data

patts <- with(dat, paste0(M192Q01, M406Q01, M423Q01, M496Q01, M564Q01,
                          M571Q01, M603Q01))

summary( as.factor(patts) )

summary( rgamma(1e5, 1, 1) )

mypriors <- dpriors(nu = "normal(.5, .5)", lambda = "normal(.5, .25)",
                    theta = "gamma(1, 1)[sd]")

## specifying my model:
m1 <- ' f1 =~ M192Q01 + M406Q01 + M423Q01 + M496Q01 + M564Q01 + M571Q01 + M603Q01 '

## drawing prior samples (100 for each of three chains):
m1pri <- bcfa(m1, data = dat, burnin = 100, sample = 100, std.lv = TRUE, prisamp = TRUE,
              dp = mypriors)

pridat <- sampleData(m1pri, simplify = TRUE)

dataset <- pridat[[ 1 ]]
hist(dataset[, 1], main = "")

cors <- sapply(pridat, function(x) cor( x[,1], x[,2] ))
hist(cors, main = "")

m1est <- bcfa(m1, data = dat, burnin = 1000, sample = 1000, std.lv = TRUE,
              dp = mypriors)

summary(m1est)

res <- blavFitIndices(m1est, pD = "loo")

summary(res)

postdat <- sampleData(m1est, nrep = 1000, simplify = TRUE)

y1rep <- t( sapply(postdat, function(x) x[,1]) )

ppc_dens_overlay(y = blavInspect(m1est, 'data')[, 1], yrep = y1rep[1:200, ])

## Percent correct vs posterior predictions, all items
obsdat <- blavInspect(m1est, 'data')

yrep <- lapply(1:7, function(i) t( sapply(postdat, function(x) x[, i]) ))
yrep <- do.call("cbind", yrep)

ppc_stat_grouped(y = as.numeric(obsdat), yrep = yrep, group = rep(1:7, each = nobs(m1est)))

## Observed correlations vs posterior correlations, all pairs of items
cormat <- cor(obsdat)
obscor <- cormat[ lower.tri(cormat) ]

correp <- lapply(postdat, cor)
correp <- lapply(correp, function(x) x[ lower.tri(x) ])
correp <- do.call("rbind", correp)

yobs <- as.numeric(obscor)

grp <- as.factor( unlist( sapply(1:6, function(i) paste0("y", i, (i+1):7)) ) )
ppc_stat_grouped(y = yobs, yrep = correp, group = grp, facet_args = list(nrow = 4)) +
  theme(legend.position = "none")

ppps <- sapply(1:21, function(i) min(mean(correp[ ,i] < yobs[i]), mean(correp[, i] > yobs[i])))

plot(obscor, ppps, xlab = "Observed correlation", ylab = "ppp", pch = 20)


#################################################################################
## IRT Section
#################################################################################
mypriors <- dpriors(tau = "normal(0, 1.5)", lambda = "normal(1, .5)")

## the model is the same as before:
m2 <- ' f1 =~ M192Q01 + M406Q01 + M423Q01 + M496Q01 + M564Q01 + M571Q01 + M603Q01 '

## drawing prior samples (100 for each of three chains):
m2pri <- bcfa(m2, data = dat, burnin = 100, sample = 100, std.lv = TRUE, prisamp = TRUE,
              dp = mypriors, ordered = TRUE)

pridat <- sampleData(m2pri, simplify = TRUE, type = "link")

dataset <- pridat[[ 1 ]]
hist(dataset[, 1], main = "")

hist(pnorm( dataset[, 1] ), main = "")

hist(pnorm( do.call("rbind", pridat)[, 1] ), main = "")

mypriors <- dpriors(tau = "normal(0, .28)", lambda = "normal(.4, .2)")

m2fit <- bcfa(m2, data = dat, std.lv = TRUE, ordered = TRUE, dp = mypriors)

m2nifit <- bcfa(m2, data = dat, std.lv = TRUE, ordered = TRUE)

summary(m2fit)

fitMeasures(m2fit)

fitMeasures(m2nifit)

## related to code at https://github.com/stan-dev/bayesplot/issues/232
combined <- rbind(plot(m2fit, 1:7, 'intervals_data', showplot = FALSE),
                  plot(m2nifit, 1:7, 'intervals_data', showplot = FALSE))
combined$model <- rep(c("Informative", "Default"), each = nrow(combined)/2)

pos <- position_nudge(y = ifelse(combined$model == "Default", 0, 0.2))

ggplot(combined, aes(x = m, y = parameter, color = model)) + 
  geom_linerange(aes(xmin = l, xmax = h), position = pos, linewidth = 2)+
  geom_linerange(aes(xmin = ll, xmax = hh), position = pos)+
  geom_point(position = pos, color="black")

it_tot <- function(fit) {
  tmpdata <- fit@Data@X[[1]]
  sapply(1:ncol(tmpdata),
         function(i) cor(tmpdata[,i], rowSums(tmpdata[,-i])))
}

itt2 <- ppmc(m2fit, discFUN = it_tot)

summary(itt2)


#################################################################################
## Explanatory IRT Section
#################################################################################
m3 <- ' f1 =~ M192Q01 + M406Q01 + M423Q01 + M496Q01 + M564Q01 + M571Q01 + M603Q01
        f1 ~ female + hisei + female:hisei '

m3fit <- bsem(m3, data = dat, std.lv = TRUE, ordered = TRUE, dp = mypriors,
              fixed.x = TRUE, save.lvs = TRUE)

summary(m3fit)

lvmeans <- blavInspect(m3fit, 'lvmeans')
dat$f1pred <- lvmeans[, 1]

regwts <- coef(m3fit)[grep("^f1~", names(coef(m3fit)))]
regdf <- cbind.data.frame(female = c(0, 1), int = c(0, regwts[1]), slp = c(regwts[2], sum(regwts[2:3])))

ggplot(dat, aes(x = hisei, y = f1pred)) + geom_point() + geom_abline(data = regdf, aes(slope = slp, intercept = int)) +
  facet_wrap( ~ female, labeller = label_both) + xlab("SES") + ylab("Proficiency")

p <- ggplot(dat, aes(x = hisei, y = f1pred)) + geom_point() + facet_wrap( ~ female, labeller = label_both) +
  xlab("SES") + ylab("Proficiency")

samps <- do.call("rbind", blavInspect(m3fit, 'mcmc'))
ndraws <- 50
regdf <- cbind.data.frame(female = rep(c(0, 1), ndraws), int = rep(0, ndraws * 2), slp = rep(0, ndraws * 2))
draws <- sample(1:nrow(samps), ndraws)
for (i in 1:length(draws)) {
  regwts <- samps[draws[i], grep("^f1~", colnames(samps))]
  regdf$int[i * 2] <- regwts[1]
  regdf$slp[(i - 1)*2 + 1:2] <- c(regwts[2], sum(regwts[2:3]))
}
  
p + geom_abline(data = regdf, aes(slope = slp, intercept = int), alpha=.2)

dat$female <- dat$female - .5

m4fit <- bsem(m3, data = dat, std.lv = TRUE, ordered = TRUE, dp = mypriors,
              fixed.x = TRUE, save.lvs = TRUE)

summary(m4fit)

lvmeans <- blavInspect(m4fit, 'lvmeans')
dat$f1pred <- lvmeans[, 1]

p <- ggplot(dat, aes(x = hisei, y = f1pred)) + geom_jitter() + facet_wrap( ~ female, labeller = label_both) +
  xlab("SES") + ylab("Proficiency")

samps <- do.call("rbind", blavInspect(m4fit, 'mcmc'))
ndraws <- 100
regdf <- cbind.data.frame(female = rep(c(-.5, .5), ndraws), int = rep(0, ndraws * 2), slp = rep(0, ndraws * 2))
draws <- sample(1:nrow(samps), ndraws)
for (i in 1:length(draws)) {
  regwts <- samps[draws[i], grep("^f1~", colnames(samps))]
  regdf$int[(i - 1)*2 + 1:2] <- regwts[1] * c(-.5, .5)
  regdf$slp[(i - 1)*2 + 1:2] <- regwts[2] + regwts[3] * c(-.5, .5)
}
  
p + geom_abline(data = regdf, aes(slope = slp, intercept = int), alpha=.2)


## each panel represents one posterior sample
library("patchwork")

lvs <- do.call("rbind", blavInspect(m4fit, 'lvs'))
ndraws <- 6
draws <- sample(1:nrow(samps), ndraws)

ps <- vector("list", length(draws))

for (i in 1:ndraws) {
  dat$f1pred <- lvs[draws[i], ]

  regwts <- samps[draws[i], grep("^f1~", colnames(samps))]
  regdf <- cbind.data.frame(female = c(-.5, .5), int = regwts[1] * c(-.5, .5), slp = regwts[2] + regwts[3] * c(-.5, .5))

  ps[[i]] <- ggplot(dat, aes(x = hisei, y = f1pred)) + geom_jitter() +
    geom_abline(data = regdf, aes(slope = slp, intercept = int)) +
    facet_wrap( ~ female, labeller = label_both) +
    xlab("SES") + ylab("Proficiency")
}

Reduce("+", ps)


#################################################################################
## Multilevel SEM Section
#################################################################################
m5 <- '
    level: within
        fw =~ enjoy1 + enjoy2 + enjoy3 + enjoy4
        math ~ fw + support

    level: between
        fenth =~ enth1 + enth2 + enth3
        fenj =~ enjoy1 + enjoy2 + enjoy3 + enjoy4
        math ~ fenth + fenj + support
'

load("pisa_span_sub.rda")
dat <- pisa_span_sub

mypris <- dpriors(lambda = "normal(1, .5)", theta = "gamma(.5, .5)[sd]",
                  psi = "gamma(.5, .5)[sd]", beta = "normal(0, 1)")

m5fit <- bsem(m5, data = dat, dp = mypris, cluster = "school", save.lvs = TRUE)

summary(m5fit)

## Visualize person-level relationship between enjoyment and achievement:
lvmeans <- blavInspect(m5fit, "lvmeans")
dat$lvmeans <- lvmeans[, 1]

schmath <- with(dat, tapply(math, school, mean))
dat$schmath <- schmath[match( dat$school, as.numeric(names(schmath)) )]
dat$math_within <- dat$math - dat$schmath

ggplot(dat, aes(x = lvmeans, y = math_within)) + geom_point() + xlab("Enjoyment (person lv)") +
  ylab("Math Achievement (person)")

## Visualize school-level relationship betwen enjoyment and achievement:
lv2 <- blavInspect(m5fit, "lvmeans", level = 2)

## ensure the school means are ordered in the same way that blavaan orders clusters
schmath <- schmath[ match(blavInspect(m5fit, 'cluster.id'), as.numeric(names(schmath))) ]

l2dat <- cbind.data.frame(schmath, lv2)
ggplot(l2dat, aes(x = fenj, y = schmath)) + geom_point() + xlab("Enjoyment (school lv)") +
  ylab("Math Achievement (school)")
