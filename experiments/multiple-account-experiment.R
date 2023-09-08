library(plyr)
library(dplyr)
library(expm)
library(reshape2)
library(bit)

source("/work/pi_jensen_umass_edu/kavery_umass_edu/non-cooperative-spillover/graphs/graph-utils.R")
source("/work/pi_jensen_umass_edu/kavery_umass_edu/non-cooperative-spillover/experiments/experiment-utils.R")

select.adversaries <- function(adversaries, treatment.assignments, setting) {
  ads <- which(adversaries==1)
  treat <- which(treatment.assignments==1)
  ctrl <- which(treatment.assignments==0)
  if(setting == "random") {
    to.combine <- sample(ads, 2)
    return(to.combine)
  }
  else {
    treatment.ads <- intersect(ads, treat) 
    control.ads <- intersect(ads, ctrl) 

    if(length(treatment.ads) == 1){
      combine.treat <- as.integer(treatment.ads)
    }
    else{
      combine.treat <- as.integer(sample(treatment.ads, 1))
    }
    if(length(control.ads) == 1){
      combine.ctrl <- as.integer(control.ads)
    }
    else{
      combine.ctrl <- as.integer(sample(control.ads, 1))
    }
    return("selected"=c(combine.treat, combine.ctrl))
  }
}


multiple.account.experiment <- function(graph.params, clustering, ncp.params, outcome.params, setting="dominating") { 
  # generate graph structure
  g <- generate.graph(graph.params)
  graph.properties <- get.graph.properties(g)
  graph.params$n <- graph.properties$n
  V(g)$name <- 1:graph.properties$n
  
  avg.degree <- mean(graph.properties$degrees)
  
  # generate graph clustering
  clusters <- generate.clusters(graph.properties$g, clustering)
  if(sum(clusters==1)==graph.properties$n) stop("Only one cluster found")
  
  # assign treatment 
  treatment <- treatment.assignment(graph.properties$g, clusters)
  treatment.assignments <- treatment[clusters]
  
  # prepare outcome model parameters
  if(graph.params$graph.type=="facebook") { 
    noise <- FALSE
  } else { 
    noise <- TRUE
  }
  stochastic.vars <- get.stochastic.vars(graph.properties$n, 3, 0.1, noise)
  
  bias.behavior <- data.frame(index=numeric(), size.of.dom=logical(), method=character(), pt.uncovered=numeric(), adversary.influence=numeric(), ATE.true=numeric(), ATE.adv.gui=numeric(), gui.beta=numeric(), gui.gamma=numeric(), stringsAsFactors=FALSE)
  nonadv.ATE <- as.numeric(calculate.ATE.various(0, graph.properties, matrix(0,1,graph.properties$n), outcome.params, ncp.params, treatment.assignments, stochastic.vars, bias.behavior)$ATE.adv.gui[1])
  # print(nonadv.ATE)

  ncp.params$setting <- setting
  ncp.params$max <- TRUE
  ncp.params$weighting <- "inf"
  if(graph.params$graph.type=="facebook") { 
    adversaries <- matrix(0, 1, graph.properties$n)
    if(ncp.params$setting == "dominating"){
      ordered.adversaries <- c(108, 3438, 1, 1685, 1913, 349, 415, 3981, 687, 699)
    }
    else{
      ordered.adversaries <- sample(1:graph.properties$n,10,replace = FALSE)
    }
    ncp.params$num.adv <- length(ordered.adversaries)
    adversaries[,sample(ordered.adversaries, ncp.params$num.adv)] <- 1
    dominating.adversaries.deg <- adversaries
    print(dominating.adversaries.deg)
    print(ncp.params$num.adv)
  } 
  else { 
    adversary.list <- determine.adversaries(graph.properties, ncp.params)
    dominating.adversaries.deg <- adversary.list
    print(dominating.adversaries.deg)
  }
  
  ncp.params$max.dom.adv <- max(sum(dominating.adversaries.deg), ncp.params$max.dom.adv)

  ncp.params$max <- FALSE
  ads.left <- clone(dominating.adversaries.deg)
  adversaries <- matrix(0,1,graph.properties$n)

  treat <- which(treatment.assignments==1)
  ctrl <- which(treatment.assignments==0)
  print(graph.params$graph.type)
  all.selected <- list()

  # cycle through increasing numbers of adversaries
  while(sum(ads.left)>=2) { 
    ads <- which(ads.left==1)
    treat <- which(treatment.assignments==1)
    ctrl <- which(treatment.assignments==0)
    treatment.ads <- intersect(ads, treat) 
    control.ads <- intersect(ads, ctrl) 

    if(sum(treatment.ads)==0 | sum(control.ads)==0){ # check if there's no nodes in treatment or control
      break
    }

    selected <- select.adversaries(ads.left, treatment.assignments, ncp.params$setting)
    adversaries[selected] <- 1
    ads.left[selected] <- 0
    all.selected <- append(all.selected, list(selected))

    ncp.params$max.dom.adv <- ncp.params$max.dom.adv-1
    ncp.params$num.adv <- sum(adversaries)
    
    bias.behavior <- calculate.ATE.various(length(all.selected), graph.properties, adversaries, outcome.params, ncp.params, treatment.assignments, stochastic.vars, bias.behavior, selected=all.selected, benign=TRUE)
  }
  
  if(graph.params$graph.type != "facebook"){
    ads.left <- rep(0,graph.properties$n)
    ads.left[which(clone(dominating.adversaries.deg)==0)] <- 1
    while(sum(ads.left)>=2){
      ads <- which(ads.left==1)
      treat <- which(treatment.assignments==1)
      ctrl <- which(treatment.assignments==0)
      treatment.ads <- intersect(ads, treat) 
      control.ads <- intersect(ads, ctrl) 
    
      if(sum(treatment.ads)==0 | sum(control.ads)==0){ # check if there's no nodes in treatment or control
        break
      }
      
      selected <- select.adversaries(ads.left, treatment.assignments, ncp.params$setting)
      adversaries[selected] <- 1
      ads.left[selected] <- 0
      all.selected <- append(all.selected, list(selected))
      
      ncp.params$max.dom.adv <- ncp.params$max.dom.adv-1
      ncp.params$num.adv <- sum(adversaries)
      bias.behavior <- calculate.ATE.various(length(all.selected), graph.properties, adversaries, outcome.params, ncp.params, treatment.assignments, stochastic.vars, bias.behavior, selected=all.selected, benign=TRUE)
      print(bias.behavior)
    }
  }
  else{
    while(sum(ads.left)>=2){
      ads <- which(ads.left==1)
      treat <- which(treatment.assignments==1)
      ctrl <- which(treatment.assignments==0)
      treatment.ads <- intersect(ads, treat)
      control.ads <- intersect(ads, ctrl)
      
      selected <- select.adversaries(ads.left, treatment.assignments, "random")
      adversaries[selected] <- 1
      ads.left[selected] <- 0
      all.selected <- append(all.selected, list(selected))

      ncp.params$max.dom.adv <- ncp.params$max.dom.adv-1
      ncp.params$num.adv <- sum(adversaries)
      bias.behavior <- calculate.ATE.various(length(all.selected), graph.properties, adversaries, outcome.params, ncp.params, treatment.assignments, stochastic.vars, bias.behavior, selected=all.selected, benign=TRUE)
      print(bias.behavior)
    }
  }

  ncp.params$max.dom.adv <- max(sum(dominating.adversaries.deg), ncp.params$max.dom.adv)
  
  bias.behavior$index <- as.numeric(bias.behavior$index)
  bias.behavior$pt.uncovered <- as.numeric(bias.behavior$pt.uncovered)
  bias.behavior$ATE.true <- as.numeric(bias.behavior$ATE.true)
  bias.behavior$ATE.adv.gui <- as.numeric(bias.behavior$ATE.adv.gui)
  
  bias.behavior$pt.covered <- 1 - bias.behavior$pt.uncovered
  bias.behavior$nonadv.ATE <- nonadv.ATE
  bias.behavior$avg.degree <- avg.degree
  return(bias.behavior)
}
