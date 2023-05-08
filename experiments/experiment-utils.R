library(plyr)
library(dplyr)
library(expm)
library(reshape2)
library(bit)

source("/Users/kavery/workspace/non-cooperative-spillover/graphs/graph-utils.R")

treatment.assignment <- function(g, clusters, prob=0.5) { 
  return(rbinom(length(unique(clusters)), 1, prob))
}

determine.adversaries <- function(graph.properties, ncp.params) {
  adversaries <- matrix(0, 1, graph.properties$n)
  if(ncp.params$setting == "random") { 
    rand.order <- sample(1:graph.properties$n, graph.properties$n, replace = FALSE)
    if(ncp.params$max) { 
      idx <- 1
      while(!check.dominating.set(graph.properties, adversaries)) { 
        adversaries[rand.order[idx]] <- 1
        idx <- idx + 1
      }  
    }
    else adversaries[sample(1:graph.properties$n, ncp.params$num.adv, replace=FALSE)] <- 1
  }
  if(ncp.params$setting == "dominating") {
    dominating.set <- dominate.greedy.inf(graph.properties)
    if(ncp.params$max) ncp.params$num.adv <- length(dominating.set)
    adversaries[,sample(dominating.set, ncp.params$num.adv)] <- 1
  }
  return(adversaries)
}


calculate.ATE.various <- function(idx, graph.properties, adversaries, outcome.params, ncp.params, treatment.assignments, stochastic.vars, bias.behavior, selected=NULL, benign=FALSE) { 
  tryCatch(
    expr = {
      uncovered.vertices <- 1 - adversaries %*% graph.properties$adj - adversaries
      pt.uncovered <- sum(uncovered.vertices == 1)/graph.properties$n
    },
    error = function(e){
      print(e)
      print("adversaries")
      print(adversaries)
    }
  )
  
  # calculate true outcome without adversaries
  ATE.true <- outcome.params$lambda_1 + outcome.params$lambda_2
  
  # calculate outcome with adversaries
  ncp.params <- exposure.probs(ncp.params, graph.properties, treatment.assignments, adversaries)
  outcome.adv <- outcome.model(outcome.params, treatment.assignments, graph.properties, adversaries, ncp.params, stochastic.vars, selected, benign)

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
  bias.behavior[nrow(bias.behavior)+1,] <- c(idx, over.dom.max, method, pt.uncovered, ad.inf, ATE.true, ATE.adv.gui, gui.beta, gui.gamma)  
  return(bias.behavior)
}


lam.I <- function(graph.properties, treatment.assignments, outcome) {
  frac.treated <- as.numeric((graph.properties$adj %*% treatment.assignments) / graph.properties$degrees)
  outcome.model <- lm(outcome ~ treatment.assignments + frac.treated)
  
  return(outcome.model)
}


exposure.probs <- function(ncp.params, graph.properties, treatment.assignments, adversaries, lambda=0.1, p=2) { 
  nonadv <- 1 - adversaries
  treated.nonadv <- treatment.assignments * nonadv
  control.nonadv <- (1 - treatment.assignments) * nonadv
  treated.adv <- treatment.assignments * adversaries
  control.adv <- adversaries - treated.adv
  
  ncp.params$empty <- as.vector(matrix(0, 1, graph.properties$n))
  # print(as.vector(adversaries))
  ncp.params$ncp.exposure.neighbors <- as.vector(t(graph.properties$adj %*% t(adversaries) / graph.properties$degrees))
  ncp.params$ncp.treat.exposure.neighbors <- as.vector(t(graph.properties$adj %*% t(treated.adv) / graph.properties$degrees))
  ncp.params$ncp.control.exposure.neighbors <- as.vector(t(graph.properties$adj %*% t(control.adv) / graph.properties$degrees))
  ncp.params$nonadv.treat.exposure.neighbors <- as.vector(t(graph.properties$adj %*% t(treated.nonadv) / graph.properties$degrees))
  ncp.params$nonadv.control.exposure.neighbors <- as.vector(t(graph.properties$adj %*% t(control.nonadv) / graph.properties$degrees))
  
  ncp.params$treatment.exposure.neighbors <-as.vector( t(graph.properties$adj %*% treatment.assignments / graph.properties$degrees))
  ncp.params$influence.as.ncp <- colSums(graph.properties$transition)
  
  return(ncp.params)
}

outcome.model <- function(outcome.params, treat, graph.properties, adversaries, ncp.params, stochastic.vars, selected=NULL, benign=FALSE) { 
  treated.adv <- treat * adversaries
  control.adv <- adversaries - treated.adv
  
  out.t0 <- matrix(0, 1, graph.properties$n)
  
  out.t1 <- outcome.params$lambda_0 + outcome.params$lambda_1 * treat + outcome.params$lambda_2 * rowSums(graph.properties$adj %*% diag(as.numeric(out.t0)) / graph.properties$degrees) + stochastic.vars$t1
  if(!benign) out.t1[which(adversaries==1)] <- ncp.params$model(treated.adv, control.adv, outcome.params)[which(adversaries==1)]
  if(!is.null(selected)){
      for(sel in selected){
          out.t1[sel[2]] = out.t1[sel[1]]
      }
  }
  
  out.t2 <- outcome.params$lambda_0 + outcome.params$lambda_1 * treat + outcome.params$lambda_2 * rowSums(graph.properties$adj %*% diag(as.numeric(out.t1)) / graph.properties$degrees) + stochastic.vars$t2
  if(!benign) out.t2[which(adversaries == 1)] <- ncp.params$model(treated.adv, control.adv, outcome.params)[which(adversaries == 1)]
  if(!is.null(selected)){
      for(sel in selected){
          out.t2[sel[2]] = out.t2[sel[1]]
      }
  }

  out.t3 <- outcome.params$lambda_0 + outcome.params$lambda_1 * treat + outcome.params$lambda_2 * rowSums(graph.properties$adj %*% diag(as.numeric(out.t2)) / graph.properties$degrees) + stochastic.vars$t3
  if(!benign) out.t3[which(adversaries == 1)] <- ncp.params$model(treated.adv, control.adv, outcome.params)[which(adversaries == 1)]
  if(!is.null(selected)){
      for(sel in selected){
          out.t3[sel[2]] = out.t3[sel[1]]
      }
  }

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

reduction.adv.model <- function(treated.adv, control.adv, outcome.params, clusters=NULL) { 
  n <- length(treated.adv)
  out <- matrix(0, 1, n)  
  out[1,which(treated.adv==1)] <- outcome.params$lambda_0 + outcome.params$lambda_1
  out[1,which(control.adv==1)] <- outcome.params$lambda_0 
  
  return(out) 
}



