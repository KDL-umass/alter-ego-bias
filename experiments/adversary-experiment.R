library(plyr)
library(dplyr)
library(expm)
library(reshape2)

#source("../graphs/graph-utils.R")


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

ncp.experiment <- function(graph.params, clustering, ncp.params, outcome.params) { 
  # generate graph structure
  g <- generate.graph(graph.params)
  graph.properties <- get.graph.properties(g)
  graph.params$n <- graph.properties$n
  
  avg.degree <- mean(graph.properties$degrees)
  
  # generate graph clustering
  clusters <- generate.clusters(graph.properties$g, clustering)
  if(sum(clusters==1)==graph.properties$n) stop("Only one cluster found")
  
  # assign treatment 
  treatment <- treatment.assignment(graph.properties$g, clusters)
  treatment.assignments <- treatment[clusters]
  
  # prepre outcome model parameters
  if(graph.params$graph.type=="facebook") { 
    noise <- FALSE
  } else { 
    noise <- TRUE
  }
  stochastic.vars <- get.stochastic.vars(graph.properties$n, 3, 0.1, noise)
  
  bias.behavior <- data.frame(index=numeric(), size.of.dom=logical(), method=character(), pt.uncovered=numeric(), ncp.influence=numeric(), ATE.true=numeric(), ATE.adv.est=numeric(), ATE.adv.gui=numeric(), gui.beta=numeric(), gui.gamma=numeric(), stringsAsFactors=FALSE)
  nonadv.ATE <- as.numeric(calculate.ATE.various(0, graph.properties, matrix(0,1,graph.properties$n), outcome.params, ncp.params, treatment.assignments, stochastic.vars, bias.behavior)$ATE.adv.gui[1])
  
  ncp.params$setting <- "dominating"
  ncp.params$max <- TRUE
  ncp.params$weighting <- "degree"
  if(graph.params$graph.type=="facebook") { 
    adversaries <- matrix(0, 1, graph.properties$n)
    dominating.adversaries.load <- c(108, 3438, 1, 1685, 1913, 349, 415, 3981, 687, 699)
    ncp.params$num.adv <- length(dominating.adversaries.load)
    adversaries[,sample(dominating.adversaries.load, ncp.params$num.adv)] <- 1
    dominating.adversaries.deg <- adversaries
  } else { 
    dominating.adversaries.deg <- unlist(list(determine.adversaries(graph.properties, ncp.params)))
  }
  ncp.params$max.dom.adv <- max(sum(dominating.adversaries.deg), ncp.params$max.dom.adv)
  
  ncp.params$max <- FALSE
  # cycle through increasing numbers of adversaries
  for(x in 1:sum(dominating.adversaries.deg == 1)) { 
    # decide adversaries
    ncp.params$num.adv <- x
    adversaries <- matrix(0,1,graph.properties$n)
    adversaries[which(dominating.adversaries.deg==1)[1:x]] <- 1
    
    # compare to Gui estimator
    bias.behavior <- calculate.ATE.various(x, graph.properties, adversaries, outcome.params, ncp.params, treatment.assignments, stochastic.vars, bias.behavior)
  }
  
  ncp.params$max <- TRUE
  ncp.params$setting <- "random"
  
  #random.adversaries <- unlist(list(determine.adversaries(graph.properties, ncp.params)))
  rand.size <- sum(dominating.adversaries.deg==1)
  random.selection <- sample(1:graph.properties$n, rand.size)
  random.adversaries <- matrix(0,1,graph.properties$n)
  random.adversaries[random.selection] <- 1
  
  for(x in 1:rand.size) {
    # decide adversaries
    ncp.params$num.adv <- x
    adversaries <- matrix(0,1,graph.properties$n)
    adversaries[which(random.adversaries==1)[1:x]] <- 1
    
    # compare to Gui estimator
    bias.behavior <- calculate.ATE.various(x, graph.properties, adversaries, outcome.params, ncp.params, treatment.assignments, stochastic.vars, bias.behavior)
  }
  
  bias.behavior$index <- as.numeric(bias.behavior$index)
  bias.behavior$pt.uncovered <- as.numeric(bias.behavior$pt.uncovered)
  bias.behavior$ATE.adv.est <- as.numeric(bias.behavior$ATE.adv.est)
  bias.behavior$ATE.true <- as.numeric(bias.behavior$ATE.true)
  bias.behavior$ATE.adv.gui <- as.numeric(bias.behavior$ATE.adv.gui)
  
  #bias.behavior.ATE <- melt(bias.behavior, id.vars=c("index", "size.of.dom", "method", "pt.uncovered", "ncp.influence", "ATE.true"))
  bias.behavior$pt.covered <- 1 - bias.behavior$pt.uncovered
  bias.behavior$nonadv.ATE <- nonadv.ATE
  bias.behavior$avg.degree <- avg.degree
  #ggplot(bias.behavior, aes(ncp.influence, ATE.adv.gui, color=method)) + geom_point()
  return(bias.behavior)
}

treatment.assignment <- function(g, clusters, prob=0.5) { 
  return(rbinom(length(unique(clusters)), 1, prob))
}

determine.adversaries <- function(graph.properties, ncp.params) {
  adversaries <- matrix(0, 1, graph.properties$n)
  if(ncp.params$setting == "random") { 
    if(ncp.params$weighting == "inf") { 
      infs <- colSums(graph.properties$transition)
      if(ncp.params$max) { 
          ii <- 1
          while(!check.dominating.set(graph.properties, adversaries)) { 
            adversaries[infs[order(infs, decreasing = TRUE)]]  
          }
      }
    } else { 
      rand.order <- sample(1:graph.properties$n, graph.properties$n, replace = FALSE)
      if(ncp.params$max) { 
        idx <- 1
        while(!check.dominating.set(graph.properties, adversaries)) { 
          adversaries[rand.order[idx]] <- 1
          idx <- idx + 1
        }  
    }
      else adversaries[sample(1:n, ncp.params$num.adv, replace=FALSE)] <- 1
    }
  }
  if(ncp.params$setting == "dominating") {
    if(ncp.params$weighting == "degree"){
      dominating.set <- dominate.greedy(graph.properties)
      if(ncp.params$max) ncp.params$num.adv <- length(dominating.set)
      adversaries[,sample(dominating.set, ncp.params$num.adv)] <- 1
    } else {
      dominating.set <- dominate.greedy.inf(graph.properties)
      if(ncp.params$max) ncp.params$num.adv <- length(dominating.set)
      adversaries[,sample(dominating.set, ncp.params$num.adv)] <- 1
    }
  }
  return(adversaries)
}

calculate.ATE.various <- function(idx, graph.properties, adversaries, outcome.params, ncp.params, treatment.assignments, stochastic.vars, bias.behavior) { 
  uncovered.vertices <- 1 - adversaries %*% graph.properties$adj - adversaries
  pt.uncovered <- sum(uncovered.vertices == 1)/graph.properties$n
  #prepare.for.plots(g, adversaries, ncp.params$ncp.exposure, treatment.assignments)
  
  # calculate true outcome without adversaries
  ATE.true <- outcome.params$lambda_1 + outcome.params$lambda_2
  
  # calculate outcome with adversaries
  ncp.params <- exposure.probs(ncp.params, graph.properties, treatment.assignments, adversaries)
  outcome.adv <- outcome.model(outcome.params, treatment.assignments, graph.properties, adversaries, ncp.params, stochastic.vars)
  
  # estimate ATE allowing adversaries
  lm.estimator.adv <- lam.I.adv(treatment.assignments, ncp.params, outcome.adv)
  ATE.adv.est <- lm.estimator.adv$coefficients[2] + lm.estimator.adv$coefficients[3]
  ATE.bias.adv <- lm.estimator.adv$coefficients[4] - lm.estimator.adv$coefficients[5]
  if(is.na(ATE.bias.adv)) ATE.bias.adv <- 0
  
  # estimate ATE using the Gui framework
  lm.estimator.gui <- lam.I(graph.properties, treatment.assignments, outcome.adv)
  gui.beta <- lm.estimator.gui$coefficients[2]
  gui.gamma <- lm.estimator.gui$coefficients[3]
  ATE.adv.gui <- gui.beta + gui.gamma
  
  over.dom.max <- ifelse(ncp.params$setting == "dominating", FALSE, ncp.params$max.dom.adv < sum(adversaries))
  if(idx == 0) over.dom.max <- FALSE
  ad.inf <- sum(ncp.params$influence.as.ncp[which(adversaries==1)])/graph.properties$n
  
  method <- ifelse(ncp.params$setting=="dominating", ncp.params$weighting, ncp.params$setting)
  if(is.null(ncp.params$setting)) method <- "none"
  
  # determine bias in estimate due to adversaries
  bias.behavior[nrow(bias.behavior)+1,] <- c(idx, over.dom.max, method, pt.uncovered, ad.inf, ATE.true, ATE.bias.adv, ATE.adv.gui, gui.beta, gui.gamma)  
  return(bias.behavior)
}

lam.I <- function(graph.properties, treatment.assignments, outcome) {
  frac.treated <- as.numeric((graph.properties$adj %*% treatment.assignments) / graph.properties$degrees)
  
  outcome.model <- lm(outcome ~ treatment.assignments + frac.treated)
  return(outcome.model)
}

lam.I.adv <- function(treatment.assignments, ncp.params, outcome) {
  #outcome.model <- lm(outcome ~ treatment.assignments + ncp.params$nonadv.treat.exposure + ncp.params$ncp.treat.exposure + ncp.params$ncp.control.exposure)
  #outcome.model <- lm(outcome ~ treatment.assignments + ncp.params$treatment.exposure.total.rw + ncp.params$ncp.treat.exposure + ncp.params$ncp.control.exposure)
  #outcome.model <- lm(outcome ~ treatment.assignments + ncp.params$treatment.exposure.neighbors + ncp.params$ncp.treat.exposure + ncp.params$ncp.control.exposure)
  outcome.model <- lm(outcome ~ treatment.assignments + ncp.params$nonadv.treat.exposure.neighbors + ncp.params$ncp.treat.exposure.neighbors + ncp.params$ncp.control.exposure.neighbors)
  
  return(outcome.model)
}

exposure.probs <- function(ncp.params, graph.properties, treatment.assignments, adversaries, lambda=0.1, p=2) { 
  nonadv <- 1 - adversaries
  treated.nonadv <- treatment.assignments * nonadv
  control.nonadv <- (1 - treatment.assignments) * nonadv
  treated.adv <- treatment.assignments * adversaries
  control.adv <- adversaries - treated.adv
  
  ncp.params$empty <- as.vector(matrix(0, 1, graph.properties$n))
  ncp.params$ncp.exposure.neighbors <- as.vector(t(graph.properties$adj %*% t(adversaries) / graph.properties$degrees))
  ncp.params$ncp.treat.exposure.neighbors <- as.vector(t(graph.properties$adj %*% t(treated.adv) / graph.properties$degrees))
  ncp.params$ncp.control.exposure.neighbors <- as.vector(t(graph.properties$adj %*% t(control.adv) / graph.properties$degrees))
  ncp.params$nonadv.treat.exposure.neighbors <- as.vector(t(graph.properties$adj %*% t(treated.nonadv) / graph.properties$degrees))
  ncp.params$nonadv.control.exposure.neighbors <- as.vector(t(graph.properties$adj %*% t(control.nonadv) / graph.properties$degrees))
  
  ncp.params$treatment.exposure.neighbors <-as.vector( t(graph.properties$adj %*% treatment.assignments / graph.properties$degrees))
  ncp.params$influence.as.ncp <- colSums(graph.properties$transition)
  
  return(ncp.params)
}

outcome.model <- function(outcome.params, treat, graph.properties, adversaries, ncp.params, stochastic.vars) { 
  treated.adv <- treat * adversaries
  control.adv <- adversaries - treated.adv
  
  out.t0 <- matrix(0, 1, graph.properties$n)
  out.t1 <- outcome.params$lambda_0 + outcome.params$lambda_1 * treat + outcome.params$lambda_2 * rowSums(graph.properties$adj %*% diag(as.numeric(out.t0)) / graph.properties$degrees) + stochastic.vars$t1
  #out.t1 <- (out.t1 > 0) + 0
  out.t1[which(adversaries==1)] <- ncp.params$model(treated.adv, control.adv, outcome.params)[which(adversaries==1)]
  
  out.t2 <- outcome.params$lambda_0 + outcome.params$lambda_1 * treat + outcome.params$lambda_2 * rowSums(graph.properties$adj %*% diag(as.numeric(out.t1)) / graph.properties$degrees) + stochastic.vars$t2
  #out.t2 <- (out.t2 > 0) + 0
  out.t2[which(adversaries == 1)] <- ncp.params$model(treated.adv, control.adv, outcome.params)[which(adversaries == 1)]
  
  out.t3 <- outcome.params$lambda_0 + outcome.params$lambda_1 * treat + outcome.params$lambda_2 * rowSums(graph.properties$adj %*% diag(as.numeric(out.t2)) / graph.properties$degrees) + stochastic.vars$t3
  #out.t3 <- (out.t3 > 0) + 0
  out.t3[which(adversaries == 1)] <- ncp.params$model(treated.adv, control.adv, outcome.params)[which(adversaries == 1)]
  
  return(out.t3) 
}

add.graph.params <- function(bias.behavior, graph.params) { 
  params <- c("n", "graph.type", "power", "degree", "p", "mu", "ncoms", "maxc", "minc")
  
  for(pa in params) { 
    bias.behavior[[pa]] <- ifelse(!is.null(graph.params[[pa]]), graph.params[[pa]], "")  
  }
  
  return(bias.behavior)
}

add.outcome.params <- function(bias.behavior, outcome.params) { 
  params <- c("lambda_0", "lambda_1", "lambda_2")
  
  for(pa in params) { 
    bias.behavior[[pa]] <- outcome.params[[pa]]
  }
  
  return(bias.behavior)
}

reduction.adv.model <- function(treated.adv, control.adv, outcome.params) { 
  n <- length(treated.adv)
  out <- matrix(0, 1, n)  
  out[1,which(treated.adv==1)] <- outcome.params$lambda_0
  out[1,which(control.adv==1)] <- outcome.params$lambda_0 + outcome.params$lambda_1
  
  #out[1,which(treated.adv==1)] <- 0
  #out[1,which(control.adv==1)] <- 1
  return(out) 
}