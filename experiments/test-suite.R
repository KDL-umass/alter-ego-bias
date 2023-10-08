source("/work/pi_jensen_umass_edu/kavery_umass_edu/non-cooperative-spillover/experiments/experiment-utils.R")
source("/work/pi_jensen_umass_edu/kavery_umass_edu/non-cooperative-spillover/experiments/multiple-account-experiment.R")

test.mult.config <- function(idx, setting, configs, trials, all=FALSE) { 
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
  adversary.params$setting <- setting
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
    
    results <- rbind(results, bias.behavior.ATE)
    write.csv(results, paste0("/work/pi_jensen_umass_edu/kavery_umass_edu/non-cooperative-spillover/results/new-",setting,"-results-", graph.params$graph.type, "-", outcome.params["lambda_1"], "-", outcome.params["lambda_2"], "-", i+20, ".csv"))
  }
}

test <- function() { 
  test.all(100)  
}

test.all <- function(trials, all=FALSE) { 
  configs <- read.csv("/work/pi_jensen_umass_edu/kavery_umass_edu/non-cooperative-spillover/experiments/configs/all_adv_configurations.csv")
  
  for(idx in 1:length(configs[[1]])) { 
    test.single.config(idx, configs, trials, all)
  }
}

test.all.mult <- function(trials, all=FALSE) { 
  configs <- read.csv("/work/pi_jensen_umass_edu/kavery_umass_edu/non-cooperative-spillover/experiments/configs/all_adv_configurations.csv")
  
  for(idx in 1:length(configs[[1]])) { 
    test.mult.config(idx, configs, trials, all)
  }
}

args <- commandArgs(trailingOnly = TRUE)
print(args)
idx <- as.integer(args[1])
configs <- read.csv("/work/pi_jensen_umass_edu/kavery_umass_edu/non-cooperative-spillover/experiments/configs/all_adv_configurations_rw.csv")
setting <- args[2]
trials <- as.integer(args[3])
test.mult.config(idx, setting, configs, trials, FALSE)

# test.all.mult(10)  
