# "Statistical modeling to adjust for time trends in adaptive platform trials utilising non-concurrent controls"
# Pavla Krotka, Martin Posch, Mohamed Gewily, Günter Höglinger, Marta Bofill Roig
# 2024

# This script contains all code to reproduce the simulations presented in the main paper and supplementary material

###########################################################################################################################################################################################
# devtools::install_github("pavlakrotka/NCC@v1.3", force = TRUE)
library(NCC)
library(tidyverse)

n_sim <- 10000

# Treatment effects to be used
# power.t.test(n=250, sd=1, sig.level = 0.025, power=0.8, type = "two.sample", alternative = "one.sided") -> treatment effect of 0.25
###########################################################################################################################################################################################

# SETTING 1

## - 10 treatment arms
## - Linear time trend
## - Varying d and lambda
## - Aim: evaluate the fixed effect model with period adjustment and spline regression

## Sample sizes:
# get_ss_matrix(num_arms = 10, n_arm = 250, d = 250*(0:9))

## Type I error rate:
set.seed(1)
setting_1_alpha <- data.frame(num_arms = 10, 
                              n_arm = 250, 
                              d1 = seq(0, 500, by = 25)*0,
                              d2 = seq(0, 500, by = 25)*1,
                              d3 = seq(0, 500, by = 25)*2,
                              d4 = seq(0, 500, by = 25)*3,
                              d5 = seq(0, 500, by = 25)*4,
                              d6 = seq(0, 500, by = 25)*5,
                              d7 = seq(0, 500, by = 25)*6,
                              d8 = seq(0, 500, by = 25)*7,
                              d9 = seq(0, 500, by = 25)*8,
                              d10 = seq(0, 500, by = 25)*9,
                              period_blocks = 2,
                              mu0 = 0,
                              sigma = 1,
                              theta1 = 0,
                              theta2 = 0,
                              theta3 = 0,
                              theta4 = 0,
                              theta5 = 0,
                              theta6 = 0,
                              theta7 = 0,
                              theta8 = 0,
                              theta9 = 0,
                              theta10 = 0,
                              lambda0 = rep(seq(-0.5, 0.5, length.out = 9), each = 21), 
                              lambda1 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                              lambda2 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                              lambda3 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                              lambda4 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                              lambda5 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                              lambda6 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                              lambda7 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                              lambda8 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                              lambda9 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                              lambda10 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                              trend = "linear",
                              unit_size = 100,
                              alpha = 0.025)

results_setting_1_alpha <- sim_study_par(nsim = n_sim, scenarios = setting_1_alpha, arms = c(2:10), models = c("fixmodel", "sepmodel", "poolmodel"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_setting_1_alpha, "results/results_setting_1_alpha.csv")



set.seed(2)
setting_1_splines_alpha <- data.frame(num_arms = 10, 
                                      n_arm = 250, 
                                      d1 = seq(0, 500, by = 25)*0,
                                      d2 = seq(0, 500, by = 25)*1,
                                      d3 = seq(0, 500, by = 25)*2,
                                      d4 = seq(0, 500, by = 25)*3,
                                      d5 = seq(0, 500, by = 25)*4,
                                      d6 = seq(0, 500, by = 25)*5,
                                      d7 = seq(0, 500, by = 25)*6,
                                      d8 = seq(0, 500, by = 25)*7,
                                      d9 = seq(0, 500, by = 25)*8,
                                      d10 = seq(0, 500, by = 25)*9,
                                      period_blocks = 2, 
                                      mu0 = 0,
                                      sigma = 1,
                                      theta1 = 0,
                                      theta2 = 0,
                                      theta3 = 0,
                                      theta4 = 0,
                                      theta5 = 0,
                                      theta6 = 0,
                                      theta7 = 0,
                                      theta8 = 0,
                                      theta9 = 0,
                                      theta10 = 0,
                                      lambda0 = rep(seq(-0.5, 0.5, length.out = 9), each = 21), 
                                      lambda1 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                                      lambda2 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                                      lambda3 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                                      lambda4 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                                      lambda5 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                                      lambda6 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                                      lambda7 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                                      lambda8 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                                      lambda9 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                                      lambda10 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                                      trend = c(rep("linear", 189), rep("stepwise_2", 189), rep("inv_u", 189), rep("seasonal", 189*2)),
                                      alpha = 0.025,
                                      ncc = TRUE,
                                      unit_size = 450,
                                      ci = FALSE,
                                      n_wave = c(rep(NA, 189*3), rep(1, 189), rep(2, 189)),
                                      N_peak = NA,
                                      bs_degree = rep(c(1, 2, 3), each = 189*5)) %>%
  rowwise() %>% mutate(N_peak = ifelse(trend=="inv_u", round(sum(get_ss_matrix(num_arms = 10, n_arm = 250, d = d2*(0:9)), na.rm = T)/2), NA))

results_setting_1_splines_alpha <- sim_study_par(nsim = n_sim, scenarios = setting_1_splines_alpha, arms = c(2:10), models = c("splines", "splines_cal", "fixmodel", "sepmodel"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_setting_1_splines_alpha, "results/results_setting_1_splines_alpha.csv")




## Power:
set.seed(3)
setting_1_pow <- data.frame(num_arms = 10, 
                            n_arm = 250, 
                            d1 = seq(0, 500, by = 25)*0,
                            d2 = seq(0, 500, by = 25)*1,
                            d3 = seq(0, 500, by = 25)*2,
                            d4 = seq(0, 500, by = 25)*3,
                            d5 = seq(0, 500, by = 25)*4,
                            d6 = seq(0, 500, by = 25)*5,
                            d7 = seq(0, 500, by = 25)*6,
                            d8 = seq(0, 500, by = 25)*7,
                            d9 = seq(0, 500, by = 25)*8,
                            d10 = seq(0, 500, by = 25)*9,
                            period_blocks = 2,
                            mu0 = 0,
                            sigma = 1,
                            theta1 = 0.25,
                            theta2 = 0.25,
                            theta3 = 0.25,
                            theta4 = 0.25,
                            theta5 = 0.25,
                            theta6 = 0.25,
                            theta7 = 0.25,
                            theta8 = 0.25,
                            theta9 = 0.25,
                            theta10 = 0.25,
                            lambda0 = rep(seq(-0.5, 0.5, length.out = 9), each = 21), 
                            lambda1 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                            lambda2 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                            lambda3 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                            lambda4 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                            lambda5 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                            lambda6 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                            lambda7 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                            lambda8 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                            lambda9 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                            lambda10 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                            trend = "linear",
                            unit_size = 100,
                            alpha = 0.025)

results_setting_1_pow <- sim_study_par(nsim = n_sim, scenarios = setting_1_pow, arms = c(2:10), models = c("fixmodel", "sepmodel", "poolmodel"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_setting_1_pow, "results/results_setting_1_pow.csv")



set.seed(4)
setting_1_splines_pow <- data.frame(num_arms = 10, 
                                    n_arm = 250, 
                                    d1 = seq(0, 500, by = 25)*0,
                                    d2 = seq(0, 500, by = 25)*1,
                                    d3 = seq(0, 500, by = 25)*2,
                                    d4 = seq(0, 500, by = 25)*3,
                                    d5 = seq(0, 500, by = 25)*4,
                                    d6 = seq(0, 500, by = 25)*5,
                                    d7 = seq(0, 500, by = 25)*6,
                                    d8 = seq(0, 500, by = 25)*7,
                                    d9 = seq(0, 500, by = 25)*8,
                                    d10 = seq(0, 500, by = 25)*9,
                                    period_blocks = 2, 
                                    mu0 = 0,
                                    sigma = 1,
                                    theta1 = 0.25,
                                    theta2 = 0.25,
                                    theta3 = 0.25,
                                    theta4 = 0.25,
                                    theta5 = 0.25,
                                    theta6 = 0.25,
                                    theta7 = 0.25,
                                    theta8 = 0.25,
                                    theta9 = 0.25,
                                    theta10 = 0.25,
                                    lambda0 = rep(seq(-0.5, 0.5, length.out = 9), each = 21), 
                                    lambda1 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                                    lambda2 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                                    lambda3 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                                    lambda4 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                                    lambda5 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                                    lambda6 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                                    lambda7 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                                    lambda8 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                                    lambda9 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                                    lambda10 = rep(seq(-0.5, 0.5, length.out = 9), each = 21),
                                    trend = c(rep("linear", 189), rep("stepwise_2", 189), rep("inv_u", 189), rep("seasonal", 189*2)),
                                    alpha = 0.025,
                                    ncc = TRUE,
                                    unit_size = 450,
                                    ci = FALSE,
                                    n_wave = c(rep(NA, 189*3), rep(1, 189), rep(2, 189)),
                                    N_peak = NA,
                                    bs_degree = rep(c(1, 2, 3), each = 189*5)) %>%
  rowwise() %>% mutate(N_peak = ifelse(trend=="inv_u", round(sum(get_ss_matrix(num_arms = 10, n_arm = 250, d = d2*(0:9)), na.rm = T)/2), NA))

results_setting_1_splines_pow <- sim_study_par(nsim = n_sim, scenarios = setting_1_splines_pow, arms = c(2:10), models = c("splines", "splines_cal", "fixmodel", "sepmodel"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_setting_1_splines_pow, "results/results_setting_1_splines_pow.csv")

###########################################################################################################################################################################################

# SETTING 2 

## - 4 treatment arms
## - Linear, stepwise, inverted-U and seasonal time trend
## - Varying lambda and calendar time unit size
## - Aim: evaluate the definition of time as calendar time intervals

## Sample sizes:
# get_ss_matrix(num_arms = 4, n_arm = 250, d = 250*(0:3))

## Type I error rate:
set.seed(5)
setting_2_alpha <- data.frame(num_arms = 4, 
                              n_arm = 250, 
                              d1 = 250*0,
                              d2 = 250*1,
                              d3 = 250*2,
                              d4 = 250*3,
                              period_blocks = 2, 
                              mu0 = 0,
                              sigma = 1,
                              theta1 = 0,
                              theta2 = 0,
                              theta3 = 0,
                              theta4 = 0,
                              lambda0 = rep(seq(-0.5, 0.5, length.out = 9), each = 30),
                              lambda1 = rep(seq(-0.5, 0.5, length.out = 9), each = 30),
                              lambda2 = rep(seq(-0.5, 0.5, length.out = 9), each = 30),
                              lambda3 = rep(seq(-0.5, 0.5, length.out = 9), each = 30),
                              lambda4 = rep(seq(-0.5, 0.5, length.out = 9), each = 30),
                              trend = c(rep("linear", 270), rep("stepwise_2", 270), rep("inv_u", 270), rep("seasonal", 270*2)),
                              alpha = 0.025,
                              ncc = TRUE,
                              unit_size = rep(seq(25, 750, by = 25), 9),
                              n_wave = c(rep(NA, 270*3), rep(1, 270), rep(2, 270)),
                              N_peak = c(rep(NA, 270*2), rep(750, 270), rep(NA, 270*2)))

results_setting_2_alpha <- sim_study_par(nsim = n_sim, scenarios = setting_2_alpha, arms = c(2:4), models = c("fixmodel", "fixmodel_cal", "sepmodel"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_setting_2_alpha, "results/results_setting_2_alpha.csv")




## Power:
set.seed(6)
setting_2_pow <- data.frame(num_arms = 4, 
                            n_arm = 250, 
                            d1 = 250*0,
                            d2 = 250*1,
                            d3 = 250*2,
                            d4 = 250*3,
                            period_blocks = 2, 
                            mu0 = 0,
                            sigma = 1,
                            theta1 = 0.25,
                            theta2 = 0.25,
                            theta3 = 0.25,
                            theta4 = 0.25,
                            lambda0 = rep(seq(-0.5, 0.5, length.out = 9), each = 30),
                            lambda1 = rep(seq(-0.5, 0.5, length.out = 9), each = 30),
                            lambda2 = rep(seq(-0.5, 0.5, length.out = 9), each = 30),
                            lambda3 = rep(seq(-0.5, 0.5, length.out = 9), each = 30),
                            lambda4 = rep(seq(-0.5, 0.5, length.out = 9), each = 30),
                            trend = c(rep("linear", 270), rep("stepwise_2", 270), rep("inv_u", 270), rep("seasonal", 270*2)),
                            alpha = 0.025,
                            ncc = TRUE,
                            unit_size = rep(seq(25, 750, by = 25), 9),
                            n_wave = c(rep(NA, 270*3), rep(1, 270), rep(2, 270)),
                            N_peak = c(rep(NA, 270*2), rep(750, 270), rep(NA, 270*2)))

results_setting_2_pow <- sim_study_par(nsim = n_sim, scenarios = setting_2_pow, arms = c(2:4), models = c("fixmodel", "fixmodel_cal", "sepmodel"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_setting_2_pow, "results/results_setting_2_pow.csv")

###########################################################################################################################################################################################

# SETTING 2

## - 4 treatment arms
## - Linear, stepwise, inverted-U and seasonal time trend
## - Varying lambda and calendar time unit size
## - Aim: evaluate the mixed models

## Sample sizes:
# get_ss_matrix(num_arms = 4, n_arm = 250, d = 250*(0:3))

## Type I error rate:
set.seed(7)
setting_2_mix_alpha <- data.frame(num_arms = 4, 
                                  n_arm = 250, 
                                  d1 = 250*0,
                                  d2 = 250*1,
                                  d3 = 250*2,
                                  d4 = 250*3,
                                  period_blocks = 2, 
                                  mu0 = 0,
                                  sigma = 1,
                                  theta1 = 0,
                                  theta2 = 0,
                                  theta3 = 0,
                                  theta4 = 0,
                                  lambda0 = rep(seq(-0.5, 0.5, length.out = 9), each = 3),
                                  lambda1 = rep(seq(-0.5, 0.5, length.out = 9), each = 3),
                                  lambda2 = rep(seq(-0.5, 0.5, length.out = 9), each = 3),
                                  lambda3 = rep(seq(-0.5, 0.5, length.out = 9), each = 3),
                                  lambda4 = rep(seq(-0.5, 0.5, length.out = 9), each = 3),
                                  trend = c(rep("linear", 27), rep("stepwise_2", 27), rep("inv_u", 27), rep("seasonal", 27*2)),
                                  alpha = 0.025,
                                  ncc = TRUE,
                                  unit_size = rep(c(25, 50, 100), 9),
                                  n_wave = c(rep(NA, 27*3), rep(1, 27), rep(2, 27)),
                                  N_peak = c(rep(NA, 27*2), rep(750, 27), rep(NA, 27*2)))

results_setting_2_mix_alpha <- sim_study_par(nsim = n_sim, scenarios = setting_2_mix_alpha, arms = c(2:4), models = c("fixmodel", "mixmodel", "mixmodel_AR1", "mixmodel_cal", "mixmodel_AR1_cal", "sepmodel"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_setting_2_mix_alpha, "results/results_setting_2_mix_alpha.csv")




## Power:
set.seed(8)
setting_2_mix_pow <- data.frame(num_arms = 4, 
                                n_arm = 250, 
                                d1 = 250*0,
                                d2 = 250*1,
                                d3 = 250*2,
                                d4 = 250*3,
                                period_blocks = 2, 
                                mu0 = 0,
                                sigma = 1,
                                theta1 = 0.25,
                                theta2 = 0.25,
                                theta3 = 0.25,
                                theta4 = 0.25,
                                lambda0 = rep(seq(-0.5, 0.5, length.out = 9), each = 3),
                                lambda1 = rep(seq(-0.5, 0.5, length.out = 9), each = 3),
                                lambda2 = rep(seq(-0.5, 0.5, length.out = 9), each = 3),
                                lambda3 = rep(seq(-0.5, 0.5, length.out = 9), each = 3),
                                lambda4 = rep(seq(-0.5, 0.5, length.out = 9), each = 3),
                                trend = c(rep("linear", 27), rep("stepwise_2", 27), rep("inv_u", 27), rep("seasonal", 27*2)),
                                alpha = 0.025,
                                ncc = TRUE,
                                unit_size = rep(c(25, 50, 100), 9),
                                n_wave = c(rep(NA, 27*3), rep(1, 27), rep(2, 27)),
                                N_peak = c(rep(NA, 27*2), rep(750, 27), rep(NA, 27*2)))

results_setting_2_mix_pow <- sim_study_par(nsim = n_sim, scenarios = setting_2_mix_pow, arms = c(2:4), models = c("fixmodel", "mixmodel", "mixmodel_AR1", "mixmodel_cal", "mixmodel_AR1_cal", "sepmodel"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_setting_2_mix_pow, "results/results_setting_2_mix_pow.csv")

###########################################################################################################################################################################################

# SETTING 3

## - 4 treatment arms
## - Different linear time trend
## - Varying lambda and trend case
## - Aim: evaluate mixed models with interaction

## Sample sizes:
# get_ss_matrix(num_arms = 4, n_arm = 250, d = 250*c(0,1,2,3))

## Type I error rate:
set.seed(9)
setting_3_alpha <- data.frame(num_arms = 4, 
                              n_arm = 250, 
                              d1 = 250*0,
                              d2 = 250*1,
                              d3 = 250*2,
                              d4 = 250*3,
                              period_blocks = 2, 
                              mu0 = 0,
                              sigma = 1,
                              theta1 = 0,
                              theta2 = 0,
                              theta3 = 0,
                              theta4 = 0,
                              lambda0 = c(rep(seq(-0.5, 0.5, length.out = 9), each = 3), rep(0, 27*3)),
                              lambda1 = c(rep(seq(-0.5, 0.5, length.out = 9), each = 3), rep(seq(-0.5, 0.5, length.out = 9), each = 3), rep(seq(-0.5, 0.5, length.out = 9), each = 3), rep(seq(-0.5, 0.5, length.out = 9), each = 3)),
                              lambda2 = c(rep(seq(-0.5, 0.5, length.out = 9), each = 3), rep(0, 27), rep(seq(-0.5, 0.5, length.out = 9), each = 3), rep(seq(-0.5, 0.5, length.out = 9), each = 3)),
                              lambda3 = c(rep(seq(-0.5, 0.5, length.out = 9), each = 3), rep(0, 27*3)),
                              lambda4 = c(rep(seq(-0.5, 0.5, length.out = 9), each = 3), rep(0, 27*2), rep(seq(-0.5, 0.5, length.out = 9), each = 3)),
                              trend = "linear",
                              alpha = 0.025,
                              ncc = TRUE,
                              unit_size = rep(c(25, 50, 100), 9),
                              trend_case = c(rep("eq", 27), rep("diff_1", 27), rep("diff_12", 27), rep("diff_124", 27)))

results_setting_3_alpha <- sim_study_par(nsim = n_sim, scenarios = setting_3_alpha, arms = c(3), models = c("mixmodel_int", "mixmodel_int_cal", "fixmodel", "sepmodel"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_setting_3_alpha, "results/results_setting_3_alpha.csv")




## Power:
set.seed(10)
setting_3_pow <- data.frame(num_arms = 4, 
                            n_arm = 250, 
                            d1 = 250*0,
                            d2 = 250*1,
                            d3 = 250*2,
                            d4 = 250*3,
                            period_blocks = 2, 
                            mu0 = 0,
                            sigma = 1,
                            theta1 = 0.25,
                            theta2 = 0.25,
                            theta3 = 0.25,
                            theta4 = 0.25,
                            lambda0 = c(rep(seq(-0.5, 0.5, length.out = 9), each = 3), rep(0, 27*3)),
                            lambda1 = c(rep(seq(-0.5, 0.5, length.out = 9), each = 3), rep(seq(-0.5, 0.5, length.out = 9), each = 3), rep(seq(-0.5, 0.5, length.out = 9), each = 3), rep(seq(-0.5, 0.5, length.out = 9), each = 3)),
                            lambda2 = c(rep(seq(-0.5, 0.5, length.out = 9), each = 3), rep(0, 27), rep(seq(-0.5, 0.5, length.out = 9), each = 3), rep(seq(-0.5, 0.5, length.out = 9), each = 3)),
                            lambda3 = c(rep(seq(-0.5, 0.5, length.out = 9), each = 3), rep(0, 27*3)),
                            lambda4 = c(rep(seq(-0.5, 0.5, length.out = 9), each = 3), rep(0, 27*2), rep(seq(-0.5, 0.5, length.out = 9), each = 3)),
                            trend = "linear",
                            alpha = 0.025,
                            ncc = TRUE,
                            unit_size = rep(c(25, 50, 100), 9),
                            trend_case = c(rep("eq", 27), rep("diff_1", 27), rep("diff_12", 27), rep("diff_124", 27)))

results_setting_3_pow <- sim_study_par(nsim = n_sim, scenarios = setting_3_pow, arms = c(3), models = c("mixmodel_int", "mixmodel_int_cal", "fixmodel", "sepmodel"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_setting_3_pow, "results/results_setting_3_pow.csv")

###########################################################################################################################################################################################

# SETTING 3 - Additional simulations

## - 4 treatment arms
## - Different linear time trend (differing also within arms 1, 2 and 4)
## - Varying lambda
## - Aim: evaluate mixed models with interaction

## Type I error rate:
set.seed(11)
setting_3_alpha_add <- data.frame(num_arms = 4, 
                                  n_arm = 250, 
                                  d1 = 250*0,
                                  d2 = 250*1,
                                  d3 = 250*2,
                                  d4 = 250*3,
                                  period_blocks = 2, 
                                  mu0 = 0,
                                  sigma = 1,
                                  theta1 = 0,
                                  theta2 = 0,
                                  theta3 = 0,
                                  theta4 = 0,
                                  lambda0 = 0,
                                  lambda1 = rep(seq(-0.5, 0.5, length.out = 9)*1, each = 3),
                                  lambda2 = rep(seq(-0.5, 0.5, length.out = 9)*2, each = 3),
                                  lambda3 = 0,
                                  lambda4 = rep(seq(-0.5, 0.5, length.out = 9)*3, each = 3),
                                  trend = "linear",
                                  alpha = 0.025,
                                  ncc = TRUE,
                                  unit_size = rep(c(25, 50, 100), 9),
                                  trend_case = "diff_124_add")

results_setting_3_alpha_add <- sim_study_par(nsim = n_sim, scenarios = setting_3_alpha_add, arms = c(3), models = c("mixmodel_int", "mixmodel_int_cal", "fixmodel", "sepmodel"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_setting_3_alpha_add, "results/results_setting_3_alpha_add.csv")




## Power:
set.seed(12)
setting_3_pow_add <- data.frame(num_arms = 4, 
                                n_arm = 250, 
                                d1 = 250*0,
                                d2 = 250*1,
                                d3 = 250*2,
                                d4 = 250*3,
                                period_blocks = 2, 
                                mu0 = 0,
                                sigma = 1,
                                theta1 = 0.25,
                                theta2 = 0.25,
                                theta3 = 0.25,
                                theta4 = 0.25,
                                lambda0 = 0,
                                lambda1 = rep(seq(-0.5, 0.5, length.out = 9)*1, each = 3),
                                lambda2 = rep(seq(-0.5, 0.5, length.out = 9)*2, each = 3),
                                lambda3 = 0,
                                lambda4 = rep(seq(-0.5, 0.5, length.out = 9)*3, each = 3),
                                trend = "linear",
                                alpha = 0.025,
                                ncc = TRUE,
                                unit_size = rep(c(25, 50, 100), 9),
                                trend_case = "diff_124_add")

results_setting_3_pow_add <- sim_study_par(nsim = n_sim, scenarios = setting_3_pow_add, arms = c(3), models = c("mixmodel_int", "mixmodel_int_cal", "fixmodel", "sepmodel"), endpoint = "cont", perc_cores = 0.99)
write_csv(results_setting_3_pow_add, "results/results_setting_3_pow_add.csv")

###########################################################################################################################################################################################