library(igraph)

source('../../experiments/adversary-experiment.R')
source('../../graphs/graph-utils.R')
source('../../experiments/plot-utils.R')

g <- read_graph(file = "../../graphs/realworld/lesmis/lesmis.gml", format = "gml")
graph.properties <- get.graph.properties(g)

graph.params <- list()
graph.params$graph.type <- "lesmis"
graph.params$n <- graph.properties$n

avg.degree <- mean(graph.properties$degrees)

# generate graph clustering
clusters <- generate.clusters(graph.properties$g, clustering)
g$cluster_id <- clusters

ncp.params <- list()
ncp.params$setting <- "dominating"
ncp.params$max <- TRUE
ncp.params$weighting <- "degree"