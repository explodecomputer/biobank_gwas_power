library(pwr)
library(dplyr)
library(parallel)

n <- 50000 * 10
ndiscovery <- c(25000, 30000, 35000, 40000, 45000) * 10
nreplication <- n - ndiscovery
maf <- c(0.01, 0.05, 0.2, 0.5)

params <- expand.grid(
	maf=maf,
	power=c(0.25, 0.5, 0.75, 0.99),
	ndiscovery=ndiscovery,
	discovery_threshold=c(5e-5, 5e-6, 5e-7, 5e-8)
)

params$nreplication <- n - params$ndiscovery
params$nfdr <- pmax(1000000 * params$discovery_threshold, 1)
params$replication_threshold <- 0.05 / params$nfdr
params$discovery_power <- NA
params$replication_power <- NA

for(i in 1:nrow(params))
{
	message(i)
	params$fstat[i] <- pwr.f2.test(u=1, v=n, sig.level=5e-8, power=params$power[i])$f2
	params$rsq[i] <- params$fstat[i] / (1 + params$fstat[i])
	params$eff[i] <- params$rsq[i] / (2 * params$maf[i] * (1-params$maf[i]))
	params$discovery_power[i] <- pwr.f2.test(u=1, v=params$ndiscovery[i], sig.level=params$discovery_threshold[i], f2=params$fstat[i])$power
	params$replication_power[i] <- pwr.f2.test(u=1, v=params$nreplication[i], sig.level=params$replication_threshold[i], f2=params$rsq[i])$power
}

params$total_power <- params$discovery_power * params$replication_power



## Simulations
# For each row do 10 Simulations
# Keep the stats, including mean effect size

fastAssoc <- function(y, x)
{
	index <- is.finite(y) & is.finite(x)
	n <- sum(index)
	y <- y[index]
	x <- x[index]
	vx <- var(x)
	bhat <- cov(y, x) / vx
	ahat <- mean(y) - bhat * mean(x)
	# fitted <- ahat + x * bhat
	# residuals <- y - fitted
	# SSR <- sum((residuals - mean(residuals))^2)
	# SSF <- sum((fitted - mean(fitted))^2)

	rsq <- (bhat * vx)^2 / (vx * var(y))
	fval <- rsq * (n-2) / (1-rsq)
	tval <- sqrt(fval)
	se <- abs(bhat / tval)

	# Fval <- (SSF) / (SSR/(n-2))
	# pval <- pf(Fval, 1, n-2, lowe=F)
	p <- pf(fval, 1, n-2, lowe=F)
	return(list(
		ahat=ahat, bhat=bhat, se=se, fval=fval, pval=p
	))
}

do_sim <- function(params, i)
{
	g1 <- rbinom(params$ndiscovery[i], 2, params$maf[i])
	g2 <- rbinom(params$nreplication[i], 2, params$maf[i])

	phen1 <- g1 * params$eff[i]
	ve1 <- var(phen1) / (params$rsq[i]) - var(phen1)
	phen1 <- scale(phen1 + rnorm(params$ndiscovery[i], 0, sqrt(ve1)))

	phen2 <- g2 * params$eff[i]
	ve2 <- var(phen2) / (params$rsq[i]) - var(phen2)
	phen2 <- scale(phen2 + rnorm(params$nreplication[i], 0, sqrt(ve2)))

	# mod1 <- summary(lm(phen1 ~ g1))$coefficients
	# mod2 <- summary(lm(phen2 ~ g2))$coefficients
	# mod3 <- summary(lm(c(phen1,phen2) ~ c(g1,g2)))$coefficients
	mod1 <- fastAssoc(phen1, g1)
	mod2 <- fastAssoc(phen2, g2)
	mod3 <- fastAssoc(c(phen1,phen2), c(g1,g2))

	return(as.data.frame(list(
		eff_discovery = mod1$bhat,
		se_discovery = mod1$se,
		p_discovery = mod1$pval,
		eff_replication = mod2$bhat,
		se_replication = mod2$se,
		p_replication = mod2$pval,
		significant = mod1$pval < params$discovery_threshold[i] & mod2$pval < params$replication_threshold[i],
		eff_total = mod3$bhat,
		se_total = mod3$se,
		p_total = mod3$pval,
		significant_t = mod3$pval < 5e-8
	)))
}


do_sim_multi <- function(params, i, nsim)
{
	params$eff_discovery <- NA
	params$eff_replication <- NA
	params$simulation_power <- NA
	params$eff_total <- NA
	params$simulation_power_t <- NA
	l <- list()
	for(j in 1:nsim)
	{
		l[[j]] <- do_sim(params, i)
	}
	l <- bind_rows(l)
	params$eff_discovery[i] <- mean(l$eff_discovery)
	params$eff_replication[i] <- mean(l$eff_replication)
	params$simulation_power[i] <- sum(l$significant) / nrow(l)
	params$eff_total <- mean(l$eff_total)
	params$simulation_power_t <- sum(l$significant_t) / nrow(l)
	return(params[i,])
}

# do_sim_multi(params, 1, 20)
# do_sim_multi(params, 9, 20)


do_sim_multi_parallel <- function(params, nsim, ncores)
{
	a <- mclapply(1:nrow(params), function(x) {
		message(x)
		return(do_sim_multi(params, x, nsim))
	}, mc.cores=ncores)
	return(bind_rows(a))
}

params <- do_sim_multi_parallel(params, 100, 16)

save(params, file="../results/simulations.RData")



# N <- 105000
# fstat <- pwr.f2.test(u=1, v=N, sig=5e-8, power=0.75)$f2

# res <- array(0,100)
# res2 <- array(0,100)
# for(i in 1:100)
# {
# 	message(i)
# 	rsq <- fstat / (1 + fstat)
# 	x <- rnorm(N)
# 	ve <- (var(x) - rsq * var(x)) / rsq
# 	y <- x + rnorm(N, sd=sqrt(ve))
# 	res[i] <- fastAssoc(y,x)$pval < 5e-8
# 	res2[i] <- cor(y,x)^2
# }
# table(res)
# mean(res2)
# rsq

