---
title: Split biobank?
author: Gibran Hemani
date: "`r Sys.Date()`"
output: 
  pdf_document
---

```{r, echo=FALSE}
suppressPackageStartupMessages(suppressWarnings(library(knitr)))
opts_chunk$set(warning=FALSE, message=FALSE, echo=FALSE, cache=FALSE)
```

```{r }
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(lubridate))



library(ggplot2)
library(ggthemes)
library(dplyr)
library(reshape2)

load("../results/simulations.RData")


# Theoretical

plot_dat <- subset(params, maf==0.5)
plot_dat$power <- as.factor(round(plot_dat$power, 2))
plot_dat$ndiscovery <- plot_dat$ndiscovery / 1000

ggplot(plot_dat, aes(x=ndiscovery, y=total_power, group=power)) +
geom_line(aes(colour=power)) +
geom_point(aes(colour=power)) +
facet_grid(. ~ discovery_threshold) +
labs(x="Discovery sample size (x1000)", y="Power to detect", colour="Simulated\npower for\nn=500000") +
theme_bw() +
scale_colour_brewer(type="seq", palette = 2)


dev.new()
plot_dat <- subset(params, maf==0.5)
plot_dat$power_change <- plot_dat$total_power - plot_dat$power
plot_dat$power <- as.factor(round(plot_dat$power, 2))
plot_dat$ndiscovery <- plot_dat$ndiscovery / 1000

ggplot(plot_dat, aes(x=ndiscovery, y=power_change, group=power)) +
geom_line(aes(colour=power)) +
geom_point(aes(colour=power)) +
facet_grid(. ~ discovery_threshold) +
labs(x="Discovery sample size (x1000)", y="Change in power compared to using 500k discovery", colour="Simulated\npower for\nn=500000") +
theme_bw() +
scale_colour_brewer(type="seq", palette = 2)





# Simulations

plot_dat <- params %>% group_by(ndiscovery, power, discovery_threshold) %>%
	dplyr::summarise(simulation_power=mean(simulation_power), simulation_power_t=mean(simulation_power_t))
plot_dat$power <- as.factor(plot_dat$power)
plot_dat$ndiscovery <- plot_dat$ndiscovery / 1000

ggplot(plot_dat, aes(x=ndiscovery, y=simulation_power, group=power)) +
geom_line(aes(colour=power)) +
geom_point(aes(colour=power)) +
facet_grid(. ~ discovery_threshold) +
theme_minimal() +
labs(x="Discovery sample size (x1000)", y="Power to detect", colour="Simulated\npower for\nn=500000")




plot_dat <- params %>% group_by(ndiscovery, power, discovery_threshold) %>%
	dplyr::summarise(simulation_power=mean(simulation_power), simulation_power_t=mean(simulation_power_t))
plot_dat$power_change <- plot_dat$simulation_power - plot_dat$simulation_power_t
plot_dat$power <- as.factor(plot_dat$power)
plot_dat$ndiscovery <- plot_dat$ndiscovery / 1000

p <- plot_dat %>% group_by(power) %>% summarise(m=mean(simulation_power_t))

ggplot(plot_dat, aes(x=ndiscovery, y=power_change, group=power)) +
geom_line(aes(colour=power)) +
geom_point(aes(colour=power)) +
facet_grid(. ~ discovery_threshold) +
labs(x="Discovery sample size (x1000)", y="Change in power compared to using 500k discovery", colour="Simulated\npower for\nn=500000")


plot_dat <- params %>% group_by(power, discovery_threshold) %>%
	dplyr::summarise(
		True=mean(eff_replication),
		TwoStep=mean(eff_replication_sig),
		OneStep=mean(eff_total_sig)
	) %>% as.data.frame()


pd <- melt(plot_dat, measure.vars=c("True", "TwoStep", "OneStep"))
pd$power <- as.factor(round(pd$power,2))

ggplot(pd, aes(x=power,y=value)) +
geom_bar(stat="identity", aes(fill=variable), position="dodge") +
facet_grid(discovery_threshold ~ .) +
scale_fill_brewer(type="qual") +
labs(x="Simulated power for n=500000", y="Mean effect size", x="Method")
