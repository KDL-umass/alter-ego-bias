#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

source("experiments/test-suite.R")

configs <- read.csv("experiments/configs/all_adv_configurations_rw.csv")
test.single.config(args[[1]], configs, trials=10, all=FALSE)
