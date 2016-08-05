library(ggplot2)
library(ggthemes)
library(dplyr)

load("../results/simulations.RData")


# Theoretical

plot_dat <- subset(params, maf==0.5)
plot_dat$power <- as.factor(round(plot_dat$power, 2))
plot_dat$ndiscovery <- plot_dat$ndiscovery / 1000

ggplot(plot_dat, aes(x=ndiscovery, y=total_power, group=power)) +
geom_line(aes(colour=power)) +
geom_point(aes(colour=power)) +
facet_grid(. ~ discovery_threshold) +
theme_minimal() +
labs(x="Discovery sample size (x1000)", y="Power to detect", colour="Simulated\npower for\nn=500000")


dev.new()
plot_dat <- subset(params, maf==0.5)
plot_dat$power_change <- plot_dat$total_power - plot_dat$power
plot_dat$power <- as.factor(plot_dat$power)
plot_dat$ndiscovery <- plot_dat$ndiscovery / 1000

ggplot(plot_dat, aes(x=ndiscovery, y=power_change, group=power)) +
geom_line(aes(colour=power)) +
geom_point(aes(colour=power)) +
facet_grid(. ~ discovery_threshold) +
theme_minimal() +
labs(x="Discovery sample size (x1000)", y="Power to detect", colour="Simulated\npower for\nn=500000")





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

