library(pwr)
library(dplyr)
library(parallel)

n <- 500000
ndiscovery <- c(250000, 300000, 350000, 400000, 450000)
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
	params$rsq <- pwr.f2.test(u=1, v=n, sig.level=5e-8, power=params$power[i])$f2
	params$eff <- params$rsq / (2 * params$maf * (1-params$maf))
	params$discovery_power[i] <- pwr.f2.test(u=1, v=params$ndiscovery[i], sig.level=params$discovery_threshold[i], f2=params$rsq[i])$power
	params$replication_power[i] <- pwr.f2.test(u=1, v=params$nreplication[i], sig.level=params$replication_threshold[i], f2=params$rsq[i])$power
}

params$total_power <- params$discovery_power * params$replication_power


## Simulations
# For each row do 10 Simulations
# Keep the stats, including mean effect size


do_sim <- function(params, i)
{
	g1 <- scale(rbinom(params$ndiscovery[i], 2, params$maf[i]))
	g2 <- scale(rbinom(params$nreplication[i], 2, params$maf[i]))

	phen1 <- g1 * params$eff[i]
	ve1 <- var(phen1) / (params$rsq[i]) - var(phen1)
	phen1 <- scale(phen1 + rnorm(params$ndiscovery[i], 0, sqrt(ve1)))

	phen2 <- g2 * params$eff[i]
	ve2 <- var(phen2) / (params$rsq[i]) - var(phen2)
	phen2 <- scale(phen2 + rnorm(params$nreplication[i], 0, sqrt(ve2)))

	mod1 <- summary(lm(phen1 ~ g1))$coefficients
	mod2 <- summary(lm(phen2 ~ g2))$coefficients
	mod3 <- summary(lm(c(phen1,phen2) ~ c(g1,g2)))$coefficients

	return(as.data.frame(list(
		eff_discovery = mod1[2,1],
		se_discovery = mod1[2,2],
		p_discovery = mod1[2,4],
		eff_replication = mod2[2,1],
		se_replication = mod2[2,2],
		p_replication = mod2[2,4],
		significant = mod1[2,4] < params$discovery_threshold[i] & mod2[2,4] < params$replication_threshold[i],
		eff_total = mod3[2,1],
		se_total = mod3[2,2],
		p_total = mod3[2,4],
		significant_t = mod3[2,4] < 5e-8
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

do_sim_multi(params, 60, 5)

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

