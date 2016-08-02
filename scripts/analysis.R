library(ggplot2)
library(ggthemes)

load("../results/simulations.RData")


# Theoretical

plot_dat <- subset(params, maf==0.5)
plot_dat$power <- as.factor(plot_dat$power)
plot_dat$ndiscovery <- plot_dat$ndiscovery / 1000

ggplot(plot_dat, aes(x=ndiscovery, y=total_power, group=power)) +
geom_line(aes(colour=power)) +
geom_point(aes(colour=power)) +
facet_grid(. ~ discovery_threshold) +
theme_minimal() +
labs(x="Discovery sample size (x1000)", y="Power to detect", colour="Simulated\npower for\nn=500000")



# Simulations

plot_dat <- subset(params, maf==0.5)
plot_dat$power <- as.factor(plot_dat$power)
plot_dat$ndiscovery <- plot_dat$ndiscovery / 1000

ggplot(plot_dat, aes(x=ndiscovery, y=simulation_power, group=power)) +
geom_line(aes(colour=power)) +
geom_point(aes(colour=power)) +
facet_grid(. ~ discovery_threshold) +
theme_minimal() +
labs(x="Discovery sample size (x1000)", y="Power to detect", colour="Simulated\npower for\nn=500000")


