# power-calcs.R
# 
# This script calculates power to estimating differences in microbe relative 
# abundances and number of simultaneous covariates that can be estimated for the
#  protocol "A Pilot Study of the Effect of the Microbiome on Immune Checkpoint 
#  Inhibitor Response in Melanoma (OSU-19125)"
# 

# Load packages
library(FDRsampsize)
library(ggplot2)
library(magrittr)
library(tidyr)
library(dplyr)
library(pwr)


# Define the false positive rate
sig.level <- 0.05
num.tests <- 70
typeIerror <- sig.level / num.tests

# Sample size
n.total <- 63

###### Power to detect differences in irAE vs no irAE ###########

# Define the fraction likely to develop an irAE
# This calculation is based on findings fromOwen DH, Wei L, Bertino EM, Edd T,
#  Villalona-Calero MA, He K, Shields PG, Carbone DP, Otterson GA. Incidence, 
#  Risk Factors, and Effect on Survival of Immune-related Adverse Events in 
#  Patients With Non–Small-cell Lung Cancer. Clinical lung cancer. 2018 Nov 1;
#  19(6):e893-900.
n.irae <- n.total * 0.16

# Compliance rate
compliance <- 0.8

# Define ranges of fold changes and coefficients of variation 
# to calculate the power
fold.change <- c(2, 2.5, 3, 3.8)
coef.variation <- c(0.3, 0.4, 0.5, 0.6)

power.l <- list()
for (c in 1:length(coef.variation)) {
  power.l[[as.character(coef.variation[c])]] <- 
    power.hart(n = n.irae * compliance, 
               alpha = typeIerror, 
               log.fc = log(fold.change),
               mu = rep(100, length(fold.change)), 
               sig = rep(coef.variation[c], length(fold.change)))
}

# Reformat the list for plotting
power.df <- bind_rows(power.l)
power.df$fold.change <- fold.change

# Plot as heatmap
power.df %>%
  gather(CV, Power, -fold.change) %>%
  ggplot(aes(x = factor(fold.change), y = CV)) +
  geom_tile(aes(fill = Power)) +
  scale_fill_gradient(low = "white", high = "blue") +
  geom_label(aes(factor(fold.change), CV, label = round(Power, 2))) +
  labs(x = "Fold Change",
       y = "Coefficient of Variation") +
  # Save as Figure 3 for the protocol
  ggsave("Fig3A_power-heatmap_cv-vs-fc.png", height = 4, width = 5)



###### Power to detect differences in response ###########

# Sample size
n.total <- 63
# This calculation assumes 16% of study participants will require steroids, which
# is defined by Owen DH, Wei L, Bertino EM, Edd T, Villalona-Calero MA, He K, 
# Shields PG, Carbone DP, Otterson GA. Incidence, Risk Factors, and Effect on 
# Survival of Immune-related Adverse Events in Patients With Non–Small-cell 
# Lung Cancer. Clinical lung cancer. 2018 Nov 1;19(6):e893-900.
n.respond <- n.total * 0.5
n.nr <- n.total - n.respond

# Compliance rate
compliance <- 0.8


# Define ranges of fold changes and coefficients of variation 
# to calculate the power
fold.change <- c(1.5, 2, 2.5, 3, 3.33)
coef.variation <- c(0.3, 0.4, 0.5, 0.6, 1)

power.l <- list()
for (c in 1:length(coef.variation)) {
  power.l[[as.character(coef.variation[c])]] <- 
    power.hart(n = n.respond * compliance, 
               alpha = typeIerror, 
               log.fc = log(fold.change),
               mu = rep(100, length(fold.change)), 
               sig = rep(coef.variation[c], length(fold.change)))
}

# Reformat the list for plotting
power.df <- bind_rows(power.l)
power.df$fold.change <- fold.change

# Plot as heatmap
power.df %>%
  gather(CV, Power, -fold.change) %>%
  ggplot(aes(x = factor(fold.change), y = CV)) +
  geom_tile(aes(fill = Power)) +
  scale_fill_gradient(low = "white", high = "blue") +
  geom_label(aes(factor(fold.change), CV, label = round(Power, 2))) +
  labs(x = "Fold Change",
       y = "Coefficient of Variation")


###### Power to simultaneously estimate covariates ###########

# Degrees of freedom for the tests  (X continuous variables - 1)
u <- seq(1, 20)

# Degrees of freedom of the error (v = n − u − 1)
v <- c()

effect.size <- 0.35

for (i in u) {
  v[i] <- pwr.f2.test(u = i, f2 = effect.size, sig.level = sig.level, power = 0.8)$v
}

n <- v + u + 1

df <- data.frame(u = u,
                 v = v,
                 n = n)
df %>%
  filter(n < 70 & n > 5) %>%
  ggplot(aes(n, u)) +
  geom_line(lwd = 1) +
  theme_bw() +
  labs(x = "Sample size",
       y = "Number of Covariates") +
  annotate("text", x = 25, y = 13, 
           label = "Power = 80%", size = 3, hjust = 0) +
  annotate("text", x = 25, y = 12,
           label = "alpha = 0.05", size = 3, hjust = 0) +
  annotate("text", x = 25, y = 11, 
           label = "Effect size = 0.35", size = 3, hjust = 0) +
  ggsave("Fig3B_f2-test.png", height = 3, width = 4)
