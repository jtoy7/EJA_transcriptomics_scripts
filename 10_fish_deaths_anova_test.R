# fish deaths anova test
# Jason A. Toy
# 2021-07-29

library(tidyverse)


#load data
deaths <- read_tsv("C:/Users/jason/OneDrive/Documents/ucsc/projects/black_perch_experiment_2017/cumulative_fish_deaths.txt") %>% 
  mutate(treatment = as.factor(treatment))

str(deaths)

#calculate means and sds
deaths_sum <- 
  deaths %>% 
  group_by(treatment) %>% 
  summarise(mean_deaths_per_rep = mean(total_deaths), sd = sd(total_deaths))


#plot
ggplot(deaths_sum, aes(x = treatment, y = mean_deaths_per_rep)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_deaths_per_rep-sd, ymax = mean_deaths_per_rep+sd), width=.2) +
  theme_bw()


#run ANOVA (F = 2.088, p = 0.133)
deaths_anova <- aov(total_deaths ~ treatment, data = deaths)
summary(deaths_anova)
