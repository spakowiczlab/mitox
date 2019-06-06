# Load packages
library(FDRsampsize)
library(ggplot2)
library(magrittr)
library(tidyr)
library(dplyr)


# Define the false positive rate
sig.level <- 0.05
num.tests <- 70
typeIerror <- sig.level / num.tests

# Sample size
n.total <- 52
# This calculation assumes 15% of study participants will require steroids
n.steroids <- n.total * 0.15
n.no.steroids <- n.total - n.steroids

# Compliance rate
steroid.compliance <- 0.9
other.compliance <- 0.9

# Define ranges of fold changes and coefficients of variation 
# to calculate the power
fold.change <- c(2, 2.5, 3, 4)
coef.variation <- c(0.3, 0.4, 0.5, 0.6)

power.l <- list()
for (c in 1:length(coef.variation)) {
  power.l[[as.character(coef.variation[c])]] <- 
    power.hart(n = n.steroids * steroid.compliance, 
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