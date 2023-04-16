mdl <- function(q, p, M){
  # q is of length M; p is of length V
  mdl <- sum(q)*log(sum(q)) - 2*sum(q*lapply(q,log)) - sum(p*lapply(p,log)) 
  for(m in M){
    sum_p_alpha <- 0
    for(alpha in m){
      sum_p_alpha <- sum_p_alpha + p[alpha]
    }
    mdl <- mdl + (q[m] + sum_p_alpha)*log(q[m] + sum_p_alpha)
  }
  return(mdl)
}

find.mdl.edges <- function(g, ncp.params){
  ncp.params$model() # get the responses

}

heterogeneous.adv.model <- function(treated.adv, control.adv, outcome.params, clusters) { 
  n <- length(treated.adv)
  treated.diff.response <- sample(which(treated.adv==1), length(which(treated.adv==1))/2) 
  control.diff.response <- sample(which(control.adv==1), length(which(control.adv==1))/2) 
  out <- matrix(0, 1, n)  

  out[1,which(treated.adv==1)] <- outcome.params$lambda_0 + outcome.params$lambda_1
  out[1,which(control.adv==1)] <- outcome.params$lambda_0
  out[1,treated.diff.response==1] <- outcome.params$lambda_0 + outcome.params$lambda_1 - 0.25
  out[1,control.diff.response==1] <- outcome.params$lambda_0 + 0.25
  
  return(out) 
}

clustering.experiment <- function(graph.params, clustering, ncp.params, outcome.params, setting="dominating") { 
  # generate graph structure
  g <- generate.graph(graph.params)
  # g <- read.graph("/Users/kavery/workspace/non-cooperative-spillover/example.txt", directed=FALSE)
  graph.properties <- get.graph.properties(g)
  graph.params$n <- graph.properties$n
  
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

  #TODO cycle through and find ordering of edges that minimizes mdl. returns list of 25 vectors in order of lowest-highest mdl. vectors contain 2 nodes that form an edge
  mdl.edges <- list(c(1,2), c(3,4))
  
  ncp.params$max.dom.adv <- 25
  
  if(setting == "dominating"){
    ncp.params$max <- FALSE
    # cycle through increasing numbers of adversaries
    for(x in 1:length(mdl.edges)) { 
      # decide adversaries
      ncp.params$num.adv <- x
      # TODO add 1 node
      adversaries <- matrix(0,1,graph.properties$n)
      adversaries[which(dominating.adversaries.deg==1)[1:x]] <- 1
      
      # compare to Gui estimator
      bias.behavior <- calculate.ATE.various(x, graph.properties, adversaries, outcome.params, ncp.params, treatment.assignments, stochastic.vars, bias.behavior)
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
