library(plyr)
library(dplyr)
library(expm)
library(reshape2)
library(bit)

source("/work/pi_jensen_umass_edu/kavery_umass_edu/non-cooperative-spillover/graphs/graph-utils.R")
source("/work/pi_jensen_umass_edu/kavery_umass_edu/non-cooperative-spillover/experiments/experiment-utils.R")

build.outcome.params <- function(lambda_0, lambda_1, lambda_2, sd.noise) { 
  outcome.params <- list()
  outcome.params$lambda_0 <- lambda_0
  outcome.params$lambda_1 <- lambda_1
  outcome.params$lambda_2 <- lambda_2
  outcome.params$sd.noise <- sd.noise
  
  return(outcome.params)
}

get.stochastic.vars <- function(num, steps, sdnoise, noise) { 
  stochastic.vars <- lapply(1:steps, function(x) { 
    if(noise) return(rnorm(num, 0, sdnoise))
    else return(matrix(0,num,1))
  })
  names(stochastic.vars) <- paste0("t", 1:steps)
  
  return(stochastic.vars)
}

specify.sybils <- function(g,i){
  sybil.g <- make_full_graph(i, directed = FALSE, loops = FALSE)
  
  V(g)$name <- 1:length(V(g))
  V(sybil.g)$name <- (length(V(g))+1):(length(V(g))+i)
  attrs <- rbind(as_data_frame(g, "vertices"), as_data_frame(sybil.g, "vertices")) 
  el <- rbind(as_data_frame(g), as_data_frame(sybil.g))
  combined.g <- graph_from_data_frame(el, directed = FALSE, vertices = attrs)
  return(combined.g)
}

add.sybil.edges <- function(g, i, sybil.conn, adversaries){
  g.conn <- clone(g)
  ads <- which(adversaries==1)
  
  for(j in 1:sybil.conn){
    if(i == 1){
      rand.sybil <- length(V(g.conn))
    }
    else{
      rand.sybil <- sample((length(V(g.conn))-(i-1)):(length(V(g.conn))), 1)
    }
    g.conn <- g.conn + edge(c(rand.sybil, ads[j]))
  }
  return(g.conn)
}


sybil.experiment <- function(graph.params, clustering, ncp.params, outcome.params, setting="dominating") { 
  # generate graph structure
  g <- generate.graph(graph.params)
  graph.properties <- get.graph.properties(g)
  graph.params$n <- graph.properties$n
  if(graph.params$graph.type=="facebook") { 
    sybilnum <- 808
  }
  else{
    sybilnum <- 200
  }
  
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
  } 
  else { 
    noise <- TRUE
  }
  stochastic.vars <- get.stochastic.vars(graph.properties$n+sybilnum, 3, 0.1, noise)
  stochastic.vars.0 <- clone(stochastic.vars)
  stochastic.vars.0$t1 <- stochastic.vars.0$t1[1:graph.properties$n]
  stochastic.vars.0$t2 <- stochastic.vars.0$t2[1:graph.properties$n]
  stochastic.vars.0$t3 <- stochastic.vars.0$t3[1:graph.properties$n]
  
  bias.behavior <- data.frame(index=numeric(), size.of.dom=logical(), method=character(), pt.uncovered=numeric(), adversary.influence=numeric(), ATE.true=numeric(), ATE.adv.gui=numeric(), gui.beta=numeric(), gui.gamma=numeric(), stringsAsFactors=FALSE)
  nonadv.ATE <- as.numeric(calculate.ATE.various(0, graph.properties, matrix(0,1,graph.properties$n), outcome.params, ncp.params, treatment.assignments, stochastic.vars.0, bias.behavior)$ATE.adv.gui[1])

  ncp.params$setting <- "dominating"
  ncp.params$max <- TRUE
  ncp.params$weighting <- "inf"
  if(graph.params$graph.type=="facebook") { 
    adversaries <- matrix(0, 1, graph.properties$n)
    dominating.adversaries.load <- c(108, 3438, 1, 1685, 1913, 349, 415, 3981, 687, 699)
    ncp.params$num.adv <- length(dominating.adversaries.load)
    nondominating <- 1:graph.properties$n
    nondominating <- nondominating[!nondominating %in% dominating.adversaries.load]
    adversaries[,sample(dominating.adversaries.load, ncp.params$num.adv)] <- 1
    adversaries[,sample(nondominating, 80-ncp.params$num.adv)] <- 1
    ncp.params$num.adv <- 80
    dominating.adversaries.deg <- adversaries
  } else { 
    adversary.list <- determine.adversaries(graph.properties, ncp.params)
    dominating.set <- adversary.list$dominating.set[1:20]
    dominating.adversaries.deg <- matrix(0,1,graph.properties$n)
    dominating.adversaries.deg[dominating.set] <- 1
  }
  
  sybil.assignment <- sample(c(0,1), 1)

  if(setting == "dominating"){
    for(i in 1:sybilnum){
      # print(i)

      ncp.params$max.dom.adv <- i
      
      g.sybil <- specify.sybils(g,i)
      adversaries <- matrix(0,1,graph.properties$n+i)
      if(i==1){
        adversaries[length(adversaries)] <- 1
      }
      else{
        adversaries[(length(adversaries)-(i-1)):length(adversaries)] <- 1
      }

      g.sybil.conn <- add.sybil.edges(g.sybil, i, 10, dominating.adversaries.deg)

      sybil.assignments <- rep(sybil.assignment, i)
      combined.assignments <- c(treatment.assignments, sybil.assignments)

      stochastic.vars.i <- clone(stochastic.vars)
      stochastic.vars.i$t1 <- stochastic.vars.i$t1[1:(graph.properties$n+i)]
      stochastic.vars.i$t2 <- stochastic.vars.i$t2[1:(graph.properties$n+i)]
      stochastic.vars.i$t3 <- stochastic.vars.i$t3[1:(graph.properties$n+i)]

      ncp.params$max <- FALSE
      ncp.params$num.adv <- i
      sybil.graph.properties <- get.graph.properties(g.sybil.conn)
      
      # compare to Gui estimator
      bias.behavior <- calculate.ATE.various(i, sybil.graph.properties, adversaries, outcome.params, ncp.params, combined.assignments, stochastic.vars.i, bias.behavior)
    }
  }
  
  bias.behavior$index <- as.numeric(bias.behavior$index)
  bias.behavior$pt.uncovered <- as.numeric(bias.behavior$pt.uncovered)
  bias.behavior$ATE.true <- as.numeric(bias.behavior$ATE.true)
  bias.behavior$ATE.adv.gui <- as.numeric(bias.behavior$ATE.adv.gui)
  
  bias.behavior$pt.covered <- 1 - bias.behavior$pt.uncovered
  bias.behavior$nonadv.ATE <- nonadv.ATE
  bias.behavior$avg.degree <- avg.degree
  return(bias.behavior)
}
