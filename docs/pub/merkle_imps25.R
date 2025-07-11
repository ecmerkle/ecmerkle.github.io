#############################################################
## Packages + options
#############################################################
library("lavaan")

library("blavaan")

library("ggplot2")

library("bayesplot")

library("mvtnorm")

library("mnormt")

library("viridis")

options(width = 90)

options(future.globals.maxSize = 1000 * 1024^2)


#############################################################
## Part 1
#############################################################
## First figure
set.seed(1165)
mn <- c(45, 3)
sdv <- c(5, .4)
cor1 <- matrix(1, 2, 2)
cor1[1,2] <- cor1[2,1] <- -.5
V <- diag(sdv) %*% cor1 %*% diag(sdv)

samps <- as.data.frame(rmvnorm(1e5, mean = mn, sigma = V))
names(samps) <- c("x", "y")

mysamp <- cbind.data.frame(xtrav = c(38, 42, 47), ytrav = c(2.4, 2.7, 2.6), draw = c(1, 2, 3))

ggplot(samps, aes(x = x, y = y)) + stat_density_2d(aes(fill = ..level..), geom = "polygon", n = 70) + scale_fill_viridis(alpha = .2) + theme_bw() + theme(legend.position = "none") + geom_point(aes(x = mn[1], y = mn[2])) + xlab("Intercept") + ylab("Slope")


## Second figure
ggplot(samps, aes(x = x, y = y)) + stat_density_2d(aes(fill = ..level..), geom = "polygon", n = 70) + scale_fill_viridis(alpha = .2) + theme_bw() + theme(legend.position = "none") + geom_point(aes(x = mn[1], y = mn[2])) + xlab("Intercept") + ylab("Slope") + geom_line(data = mysamp, mapping = aes(x = xtrav, y = ytrav)) + geom_point(data = mysamp, mapping = aes(x = xtrav, y = ytrav)) + geom_text(data = mysamp, mapping = aes(x = xtrav, y = ytrav, label = draw), nudge_x = -.2, nudge_y = .05)


#############################################################
## Part 2
#############################################################
## Data + response patterns
data("data.pisaMath", package = "sirt")

dat <- data.pisaMath$data

patts <- with(dat, paste0(M192Q01, M406Q01, M423Q01, M496Q01, M564Q01,
                          M571Q01, M603Q01))

table(patts)


## Priors and prior checks
## blavaan defaults:
dpriors()

## replacing the defaults for thresholds and loadings:
mypriors <- dpriors(tau = "normal(0, 1.5)", lambda = "normal(1, .5)")


## specify the model:
m2 <- ' f1 =~ M192Q01 + M406Q01 + M423Q01 + M496Q01 + M564Q01 + M571Q01 + M603Q01 '

## drawing prior samples (100 for each of three chains):
m2pri <- bcfa(m2, data = dat, burnin = 100, sample = 100, std.lv = TRUE, prisamp = TRUE,
              dp = mypriors, ordered = TRUE)


pridat <- sampleData(m2pri, simplify = TRUE, type = "link")


dataset <- pridat[[ 1 ]]
hist(dataset[, 1], main = "")


hist(pnorm( dataset[, 1] ), main = "")


hist(pnorm( do.call("rbind", pridat)[, 1] ), main = "")


## Prior specification and model estimation
mypriors <- dpriors(tau = "normal(0, .28)", lambda = "normal(.4, .2)")


m2fit <- bcfa(m2, data = dat, std.lv = TRUE, ordered = TRUE, dp = mypriors)


m2nifit <- bcfa(m2, data = dat, std.lv = TRUE, ordered = TRUE)


summary(m2fit)


fitMeasures(m2fit)

fitMeasures(m2nifit)


## Model summaries
### related to code at https://github.com/stan-dev/bayesplot/issues/232
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


## Explanatory IRT estimation
m3 <- ' f1 =~ M192Q01 + M406Q01 + M423Q01 + M496Q01 + M564Q01 + M571Q01 + M603Q01
        f1 ~ female + hisei + female:hisei '


m3fit <- bsem(m3, data = dat, std.lv = TRUE, ordered = TRUE, dp = mypriors,
              fixed.x = TRUE, save.lvs = TRUE)


summary(m3fit)


## Explanatory IRT summary
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


## Centering female and re-estimating
dat$female <- dat$female - .5


m4fit <- bsem(m3, data = dat, std.lv = TRUE, ordered = TRUE, dp = mypriors,
              fixed.x = TRUE, save.lvs = TRUE)


## New model summaries and graphs
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



#############################################################
## Part 3
#############################################################
## The same graph from Part 1
set.seed(1165)
mn <- c(45, 3)
sdv <- c(5, .4)
cor1 <- matrix(1, 2, 2)
cor1[1,2] <- cor1[2,1] <- -.5
V <- diag(sdv) %*% cor1 %*% diag(sdv)

samps <- as.data.frame(rmvnorm(1e5, mean = mn, sigma = V))
names(samps) <- c("x", "y")

mysamp <- cbind.data.frame(xtrav = c(38, 42, 47), ytrav = c(2.4, 2.7, 2.6), draw = c(1, 2, 3))

ggplot(samps, aes(x = x, y = y)) + stat_density_2d(aes(fill = ..level..), geom = "polygon", n = 70) + scale_fill_viridis(alpha = .2) + theme_bw() + theme(legend.position = "none") + geom_point(aes(x = mn[1], y = mn[2])) + xlab("Intercept") + ylab("Slope") + geom_line(data = mysamp, mapping = aes(x = xtrav, y = ytrav)) + geom_point(data = mysamp, mapping = aes(x = xtrav, y = ytrav)) + geom_text(data = mysamp, mapping = aes(x = xtrav, y = ytrav, label = draw), nudge_x = -.2, nudge_y = .05)



#############################################################
## Part 4
#############################################################

## Leapfrog function
leapfrog <- function(theta, momentum, grfun, data, eps, Minv) {
  newmom <- momentum + (eps/2) * grfun(theta, data)
  newtheta <- theta + eps * Minv %*% newmom
  newmom <- newmom + (eps/2) * grfun(newtheta, data)

  list(theta = newtheta, momentum = newmom)
}
