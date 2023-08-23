source("/Users/kavery/workspace/non-cooperative-spillover/experiments/adversary-experiment.R")
source("/Users/kavery/workspace/non-cooperative-spillover/experiments/multiple-account-experiment.R")

test.mult.config <- function(idx, configs, trials, all=FALSE) { 
  cat("Running", idx, "\n")
  print(configs[idx,])
  
  results <- data.frame(index=numeric(), size.of.dom=logical(), method=character(), 
                        pt.uncovered=numeric(), adversary.influence=numeric(), ATE.true=numeric(), 
                        variable=numeric(), value=numeric(), pt.covered=numeric(), n=numeric(), 
                        graph.type=character(), power=numeric(), degree=numeric(), p=numeric(), 
                        mu=numeric(), ncoms=numeric(), maxc=numeric(), minc=numeric(), 
                        lambda_0=numeric(), lambda_1=numeric(), lambda_2=numeric(), stringsAsFactors=FALSE)
  
  graph.params <- build.graph.params(configs, idx)
  adversary.params <- list()
  adversary.params$model <- reduction.adv.model
  adversary.params$all <- all
  adversary.params$setting <- "dominating"
  adversary.params$weighting <- "inf"
  outcome.params <- build.outcome.params(configs[idx,"lambda_0"], configs[idx,"lambda_1"], configs[idx,"lambda_2"], configs[idx,"sd.noise"])
  clustering <- "infomap"
  
  for(i in 1:trials) {
    graph.params$ind <- i
    
    cat("trial", i, "\n")
    bias.behavior.ATE <- multiple.account.experiment(graph.params, clustering, adversary.params, outcome.params, adversary.params$setting)
    bias.behavior.ATE$adversary.influence <- as.numeric(bias.behavior.ATE$adversary.influence)
    bias.behavior.ATE$gui.beta <- as.numeric(bias.behavior.ATE$gui.beta)
    bias.behavior.ATE$gui.gamma <- as.numeric(bias.behavior.ATE$gui.gamma)
    
    bias.behavior.ATE <- add.graph.params(bias.behavior.ATE, graph.params)
    bias.behavior.ATE <- add.outcome.params(bias.behavior.ATE, outcome.params)
    bias.behavior.ATE$graph.id <- configs[idx,"graph.no"]
    bias.behavior.ATE$adv.bias <- bias.behavior.ATE$nonadv.ATE - bias.behavior.ATE$ATE.adv.gui
    
    bias.behavior.ATE$bias <- bias.behavior.ATE$ATE.true - bias.behavior.ATE$ATE.adv.gui
    bias.behavior.ATE$est.diff <- bias.behavior.ATE$nonadv.ATE - bias.behavior.ATE$ATE.adv.gui
    bias.behavior.ATE$diff.norm <- bias.behavior.ATE$est.diff / bias.behavior.ATE$nonadv.ATE
    bias.behavior.ATE$pt.adversaries <- (bias.behavior.ATE$index*2)/bias.behavior.ATE$n

    results <- rbind(results, bias.behavior.ATE)
    # write.table(results, paste0("/Users/kavery/workspace/non-cooperative-spillover/results/facebook-sybil-results-", graph.params$graph.type, "-", outcome.params["lambda_2"], "-", i, ".csv"), append = TRUE , col.names = FALSE,sep = ",")
    write.csv(results, paste0("/Users/kavery/workspace/non-cooperative-spillover/results/new-dominating-results-", graph.params$graph.type, "-", outcome.params["lambda_1"], "-", outcome.params["lambda_2"], "-", i, ".csv"))
  }
}


test <- function() { 
  test.all(100)  
}

test.all <- function(trials, all=FALSE) { 
  configs <- read.csv("/Users/kavery/workspace/non-cooperative-spillover/experiments/configs/all_adv_configurations.csv")
  
  for(idx in 1:length(configs[[1]])) { 
    test.single.config(idx, configs, trials, all)
  }
}

command.args.mult <- function(all=False) {
  args <- commandArgs(trailingOnly = TRUE)
  config <- args[1]
  trials <- args[2]
  test.mult.config(0, config, trials, all)
}

command.args.mult()