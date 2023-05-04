library(plyr)
library(dplyr)
library(expm)
library(reshape2)
library(bit)

source("/Users/kavery/workspace/non-cooperative-spillover/graphs/graph-utils.R")


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


combine.adversaries <- function(g, adversaries, treatment.assignments, setting) {
  ads <- which(adversaries==1)
  treat <- which(treatment.assignments==1)
  ctrl <- which(treatment.assignments==0)
  if(setting == "random") {
    mapping <- array( V(g) )
    to.combine <- sample(ads, 2)
    mapping[to.combine[1]] <- to.combine[2]
    # print("to.combine")
    # print(to.combine)
    g <- contract(g, mapping, vertex.attr.comb="first")
    # print(g)
    tryCatch(
      expr = {
        g <- delete_vertices(g, c(to.combine[1]))
      },
      error = function(e){
        print(e)
      }
    )
    return(list("graph"=g, "selected"=to.combine))
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
    # print(combine.treat)
    # print(combine.ctrl)
    mapping <- array( V(g) )
    mapping[as.integer(combine.ctrl)] <- as.integer(combine.treat)
    g <- contract(g, mapping, vertex.attr.comb="first")
    # print(g)
    tryCatch(
      expr = {
        g <- delete_vertices(g, combine.ctrl)
      },
      error = function(e){
        print(e)
      }
    )
    return(list("graph"=g, "selected"=c(combine.treat, combine.ctrl)))
  }
}

specify.sybils <- function(g,i){
  sybil.g <- make_full_graph(i, directed = FALSE, loops = FALSE)
  
  #combine the two graphs
  # V(g)$name <- V(g)$label
  V(g)$name <- 1:length(V(g))
  V(sybil.g)$name <- (length(V(g))+1):(length(V(g))+i)
  # V(sybil.g)$name <- V(sybil.g)$label
  attrs <- rbind(as_data_frame(g, "vertices"), as_data_frame(sybil.g, "vertices")) # %>% unique()
  # print(attrs)
  el <- rbind(as_data_frame(g), as_data_frame(sybil.g))
  # print(el)
  combined.g <- graph_from_data_frame(el, directed = FALSE, vertices = attrs)
  # print(V(combined.g))
  # print(combined.g)
  return(combined.g)
}

add.sybil.edges <- function(g, i, sybil.conn, adversaries){
  g.conn <- clone(g)
  ads <- which(adversaries==1)
  # ads <- adversary.order
  # print((length(V(g.conn))-(i-1)))
  # print((length(V(g.conn))))
  # print(adversaries)
  
  for(j in 1:sybil.conn){
    if(i == 1){
      rand.sybil <- length(V(g.conn))
    }
    else{
      rand.sybil <- sample((length(V(g.conn))-(i-1)):(length(V(g.conn))), 1)
    }
    # print(rand.sybil)
    # print(ads[j])
    g.conn <- g.conn + edge(c(rand.sybil, ads[j]))
  }
  return(g.conn)
}


sybil.experiment <- function(graph.params, clustering, ncp.params, outcome.params, setting="dominating") { 
  # generate graph structure
  g <- generate.graph(graph.params)
  # g <- read.graph("/Users/kavery/workspace/non-cooperative-spillover/example.txt", directed=FALSE)
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
  # clusters <- c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2)
  # print(clusters)
  if(sum(clusters==1)==graph.properties$n) stop("Only one cluster found")
  
  # assign treatment 
  treatment <- treatment.assignment(graph.properties$g, clusters)
  # treatment <- c(1,0)
  # print(treatment)
  treatment.assignments <- treatment[clusters]
  # print(treatment.assignments)
  
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
  print(nonadv.ATE)

  ncp.params$setting <- "dominating"
  ncp.params$max <- TRUE
  ncp.params$weighting <- "inf"
  if(graph.params$graph.type=="facebook") { 
    adversaries <- matrix(0, 1, graph.properties$n)
    dominating.adversaries.load <- c(108, 3438, 1, 1685, 1913, 349, 415, 3981, 687, 699)
    ncp.params$num.adv <- length(dominating.adversaries.load)
    nondominating <- 1:graph.properties$n
    nondominating <- nondominating[!nondominating %in% dominating.adversaries.load]
    # print(nondominating)
    adversaries[,sample(dominating.adversaries.load, ncp.params$num.adv)] <- 1
    adversaries[,sample(nondominating, 80-ncp.params$num.adv)] <- 1
    # print(adversaries)
    # print(which(adversaries==1))
    ncp.params$num.adv <- 80
    dominating.adversaries.deg <- adversaries
  } else { 
    adversary.list <- determine.adversaries(graph.properties, ncp.params)
    # dominating.adversaries.deg <- adversary.list$adversaries
    dominating.set <- adversary.list$dominating.set[1:20]
    dominating.adversaries.deg <- matrix(0,1,graph.properties$n)
    dominating.adversaries.deg[dominating.set] <- 1
    # print(dominating.set)
    # print(dominating.adversaries.deg)
  }
  
  # num.conn <- c(2, 3, 5, 10, 20, 30, 40, 50, 100)
  sybil.assignment <- sample(c(0,1), 1)
  # x <- 497:808
  # x <- x[x %% 10 == 0]
  
  if(setting == "dominating"){
    for(i in 1:sybilnum){
      # print(i)

      ncp.params$max.dom.adv <- i
      
      g.sybil <- specify.sybils(g,i)
      adversaries <- matrix(0,1,graph.properties$n+i)
      # print("length(adversaries)")
      # print(length(adversaries))
      # print((length(adversaries)-(i-1)):length(adversaries))
      if(i==1){
        adversaries[length(adversaries)] <- 1
      }
      else{
        adversaries[(length(adversaries)-(i-1)):length(adversaries)] <- 1
      }

      g.sybil.conn <- add.sybil.edges(g.sybil, i, 10, dominating.adversaries.deg)
      # print(adversaries)
      # print(length(adversaries))

      sybil.assignments <- rep(sybil.assignment, i)
      combined.assignments <- c(treatment.assignments, sybil.assignments)
      # print(combined.assignments)
      # print(length(combined.assignments))

      stochastic.vars.i <- clone(stochastic.vars)
      stochastic.vars.i$t1 <- stochastic.vars.i$t1[1:(graph.properties$n+i)]
      stochastic.vars.i$t2 <- stochastic.vars.i$t2[1:(graph.properties$n+i)]
      stochastic.vars.i$t3 <- stochastic.vars.i$t3[1:(graph.properties$n+i)]
      # print("stochastic.vars.i")
      # print(length(stochastic.vars.i$t1))

      ncp.params$max <- FALSE
      ncp.params$num.adv <- i
      sybil.graph.properties <- get.graph.properties(g.sybil.conn)
      # print(graph.properties$degrees)
      
      # compare to Gui estimator
      bias.behavior <- calculate.ATE.various(i, sybil.graph.properties, adversaries, outcome.params, ncp.params, combined.assignments, stochastic.vars.i, bias.behavior)
      # write.table(bias.behavior, paste0("/Users/kavery/workspace/non-cooperative-spillover/results/facebook-sybil-test-results-", graph.params$graph.type, "-", outcome.params["lambda_2"], ".csv"), append = TRUE , col.names = FALSE,sep = ",")
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


multiple.account.experiment <- function(graph.params, clustering, ncp.params, outcome.params, setting="dominating") { 
  # generate graph structure
  g <- generate.graph(graph.params)
  # print(E(g))
  # g <- read.graph("/Users/kavery/workspace/non-cooperative-spillover/example.txt", directed=FALSE)
  graph.properties <- get.graph.properties(g)
  graph.params$n <- graph.properties$n
  V(g)$name <- 1:graph.properties$n
  
  avg.degree <- mean(graph.properties$degrees)
  
  # generate graph clustering
  clusters <- generate.clusters(graph.properties$g, clustering)
  # clusters <- c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2)
  if(sum(clusters==1)==graph.properties$n) stop("Only one cluster found")
  
  # assign treatment 
  treatment <- treatment.assignment(graph.properties$g, clusters)
  # treatment <- c(1,0)
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

  ncp.params$setting <- "dominating"
  ncp.params$max <- TRUE
  ncp.params$weighting <- "inf"
  if(graph.params$graph.type=="facebook") { 
    adversaries <- matrix(0, 1, graph.properties$n)
    ordered.adversaries <- c(108, 3438, 1, 1685, 1913, 349, 415, 3981, 687, 699)
    ncp.params$num.adv <- length(ordered.adversaries)
    adversaries[,sample(ordered.adversaries, ncp.params$num.adv)] <- 1
    dominating.adversaries.deg <- adversaries
  } 
  else { 
    adversary.list <- determine.adversaries(graph.properties, ncp.params)
    dominating.adversaries.deg <- adversary.list$adversaries
  }
  
  ncp.params$max.dom.adv <- max(sum(dominating.adversaries.deg), ncp.params$max.dom.adv)

  ncp.params$max <- FALSE
  ads.left <- clone(dominating.adversaries.deg)
  g.dom <- clone(g)
  dom.assignments <- clone(treatment.assignments)
  dom.stoc.vars <- clone(stochastic.vars)
  adversaries <- matrix(0,1,graph.properties$n)

  treat <- which(dom.assignments==1)
  ctrl <- which(dom.assignments==0)
  print(graph.params$graph.type)
  # V(g.dom)$color <- ifelse(dom.assignments==1, "lightblue", "orange")
  # lay <- layout_with_lgl(g)
  # LO2 <- lay
  # LO2[V(g.dom)$color == "lightblue"] = LO2[V(g.dom)$color == "lightblue"] + 5000
  # LO2[V(g.dom)$color == "orange"] = LO2[V(g.dom)$color == "orange"] - 5000
  # plot(g.dom,layout=LO2,vertex.size=2,vertex.label=NA)
  # print_all(g.dom)
  # print(dom.assignments==1)
  # print(dom.assignments==0)

  # cycle through increasing numbers of adversaries
  while(sum(ads.left)>=2) { 
    ads <- which(ads.left==1)
    treat <- which(dom.assignments==1)
    ctrl <- which(dom.assignments==0)
    treatment.ads <- intersect(ads, treat) 
    control.ads <- intersect(ads, ctrl) 

    if(sum(treatment.ads)==0 | sum(control.ads)==0){ # check if there's no nodes in treatment or control
      write_graph(g.dom, paste0("/Users/kavery/workspace/non-cooperative-spillover/results/checkpoint/",graph.params$graph.type,"/",length(V(g.dom)),".txt"), "edgelist")
      break
    }

    combined <- combine.adversaries(g.dom, ads.left, dom.assignments, ncp.params$setting)
    g.dom <- combined$graph
    selected <- combined$selected
    adversaries[selected] <- 1
    ads.left[selected] <- 0

    graph.properties.dom <- get.graph.properties(g.dom)
    
    dom.assignments <- dom.assignments[-selected[2]]
    adversaries <- t(as.matrix(adversaries[-selected[2]]))
    ads.left <- ads.left[-selected[2]]
    dom.stoc.vars$t1 <- dom.stoc.vars$t1[-selected[2]]
    dom.stoc.vars$t2 <- dom.stoc.vars$t2[-selected[2]]
    dom.stoc.vars$t3 <- dom.stoc.vars$t3[-selected[2]]

    ncp.params$max.dom.adv <- ncp.params$max.dom.adv-1
    ncp.params$num.adv <- sum(adversaries)
    
    bias.behavior <- calculate.ATE.various(length(V(g.dom)), graph.properties.dom, adversaries, outcome.params, ncp.params, dom.assignments, dom.stoc.vars, bias.behavior, benign=TRUE)
    # write.table(bias.behavior, paste0("/Users/kavery/workspace/non-cooperative-spillover/results/facebook-test-mult-results-", graph.params$graph.type, "-", outcome.params["lambda_2"], ".csv"), append = TRUE , col.names = FALSE,sep = ",")
    # if((length(V(g.dom)) %% 10) == 0){
    #   write_graph(g.dom, paste0("/Users/kavery/workspace/non-cooperative-spillover/results/checkpoint/",length(V(g.dom)),".txt"), "edgelist")
    # }
  }
  ads.left <- matrix(1,1,length(V(g.dom)))
  while(sum(ads.left)>=2){
    ads <- which(ads.left==1)
    treat <- which(dom.assignments==1)
    ctrl <- which(dom.assignments==0)
    treatment.ads <- intersect(ads, treat) 
    control.ads <- intersect(ads, ctrl) 

    if(sum(treatment.ads)==0 | sum(control.ads)==0){ # check if there's no nodes in treatment or control
      print(treatment.ads)
      print(control.ads)
      break
    }
    
    combined <- combine.adversaries(g.dom, ads.left, dom.assignments, ncp.params$setting)
    g.dom <- combined$graph
    selected <- combined$selected
    adversaries[selected] <- 1
    ads.left[selected] <- 0
    
    graph.properties.dom <- get.graph.properties(g.dom)
    
    dom.assignments <- dom.assignments[-selected[2]]
    adversaries <- t(as.matrix(adversaries[-selected[2]]))
    ads.left <- ads.left[-selected[2]]
    dom.stoc.vars$t1 <- dom.stoc.vars$t1[-selected[2]]
    dom.stoc.vars$t2 <- dom.stoc.vars$t2[-selected[2]]
    dom.stoc.vars$t3 <- dom.stoc.vars$t3[-selected[2]]

    ncp.params$max.dom.adv <- ncp.params$max.dom.adv-1
    ncp.params$num.adv <- sum(adversaries)
    
    bias.behavior <- calculate.ATE.various(length(V(g.dom)), graph.properties.dom, adversaries, outcome.params, ncp.params, dom.assignments, dom.stoc.vars, bias.behavior, benign=TRUE)
    # print(bias.behavior)
  }

  # print_all(g.dom)
  # print(dom.assignments==1)
  # print(dom.assignments==0)
  # V(g.dom)$color <- ifelse(dom.assignments==1, "lightblue", "orange")
  # plot(g.dom,layout=LO2,vertex.size=2,vertex.label=NA)
  # write_graph(g.dom, paste0("/Users/kavery/workspace/non-cooperative-spillover/results/checkpoint/test-",length(V(g.dom)),".txt"), "edgelist")
  
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
    else adversaries[sample(1:n, ncp.params$num.adv, replace=FALSE)] <- 1
  }
  if(ncp.params$setting == "dominating") {
    dominating.set <- dominate.greedy.inf(graph.properties)
    if(ncp.params$max) ncp.params$num.adv <- length(dominating.set)
    adversaries[,sample(dominating.set, ncp.params$num.adv)] <- 1
  }
  return(list("adversaries" = adversaries, "dominating.set" = dominating.set))
}


calculate.ATE.various <- function(idx, graph.properties, adversaries, outcome.params, ncp.params, treatment.assignments, stochastic.vars, bias.behavior, benign=FALSE) { 
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
  #prepare.for.plots(g, adversaries, ncp.params$ncp.exposure, treatment.assignments)
  
  # calculate true outcome without adversaries
  ATE.true <- outcome.params$lambda_1 + outcome.params$lambda_2
  
  # calculate outcome with adversaries
  ncp.params <- exposure.probs(ncp.params, graph.properties, treatment.assignments, adversaries)
  outcome.adv <- outcome.model(outcome.params, treatment.assignments, graph.properties, adversaries, ncp.params, stochastic.vars, benign)
  # outcome.adv <- c(outcome.adv)
  # print(outcome.adv)
  # t <- !treatment.assignments
  # print(outcome.adv[t])
  # print(outcome.adv[!t])
  # print(mean(outcome.adv[t]))
  # print(mean(outcome.adv[!t]))
  # print(length(outcome.adv))
  # print(outcome.adv[2])
  
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
  # print(bias.behavior)
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

outcome.model <- function(outcome.params, treat, graph.properties, adversaries, ncp.params, stochastic.vars, benign=FALSE) { 
  treated.adv <- treat * adversaries
  control.adv <- adversaries - treated.adv
  
  out.t0 <- matrix(0, 1, graph.properties$n)
  
  out.t1 <- outcome.params$lambda_0 + outcome.params$lambda_1 * treat + outcome.params$lambda_2 * rowSums(graph.properties$adj %*% diag(as.numeric(out.t0)) / graph.properties$degrees) + stochastic.vars$t1
  if(!benign) out.t1[which(adversaries==1)] <- ncp.params$model(treated.adv, control.adv, outcome.params)[which(adversaries==1)]
  # print(length(out.t1))
  
  out.t2 <- outcome.params$lambda_0 + outcome.params$lambda_1 * treat + outcome.params$lambda_2 * rowSums(graph.properties$adj %*% diag(as.numeric(out.t1)) / graph.properties$degrees) + stochastic.vars$t2
  if(!benign) out.t2[which(adversaries == 1)] <- ncp.params$model(treated.adv, control.adv, outcome.params)[which(adversaries == 1)]
  
  out.t3 <- outcome.params$lambda_0 + outcome.params$lambda_1 * treat + outcome.params$lambda_2 * rowSums(graph.properties$adj %*% diag(as.numeric(out.t2)) / graph.properties$degrees) + stochastic.vars$t3
  if(!benign) out.t3[which(adversaries == 1)] <- ncp.params$model(treated.adv, control.adv, outcome.params)[which(adversaries == 1)]

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



