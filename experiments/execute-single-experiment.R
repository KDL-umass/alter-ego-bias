#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

source("experiments/test-suite.R")

configs <- read.csv("experiments/configs/all_adv_configurations.csv")
test.single.config(args[[1]], configs, trials=100, all=FALSE)