library(igraph)

source('experiments/adversary-experiment.R')
source('graphs/graph-utils.R')
source('experiments/plot-utils.R')

g <- read_graph(file = "graphs/realworld/lesmis/lesmis.gml", format = "gml")
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

adversaries.deg <- unlist(list(determine.adversaries(graph.properties, ncp.params)))
adversaries <- matrix(0,1,graph.properties$n)
adversaries[which(adversaries.deg==1)] <- 1

# random NCP
# n_adversaries <- sum(adversaries.deg)
# random.adversaries <- sample(1:graph.properties$n, n_adversaries, replace=FALSE)
# adversaries <- matrix(0,1,graph.properties$n)
# adversaries[random.adversaries] <- 1

print(adversaries)

treatment <- treatment.assignment(graph.properties$g, clusters)
treatment.assignments <- treatment[clusters]

uncovered.vertices <- 1 - adversaries %*% graph.properties$adj - adversaries
pt.uncovered <- sum(uncovered.vertices == 1)/graph.properties$n
prepare.for.plots(g, adversaries, ncp.params$ncp.exposure, treatment.assignments, labels=TRUE)


exposure.params <- exposure.probs(ncp.params, graph.properties, treatment.assignments, adversaries)
print(exposure.params$influence.as.ncp)
#V(g)$color <- adversaries.deg

g$palette <- grey.colors(100)
V(g)$color <- exposure.params$ncp.exposure.neighbors * 100
#V(g)$color[which(exposure.params$ncp.exposure.neighbors) == 0] <- "lightblue"
V(g)$color[which(adversaries > 0)] <- "red"


lo <- layout_with_kk(g) # create a layout
lo <- norm_coords(lo, ymin=-1, ymax=1, xmin=-1, xmax=1)
plot(g, layout=lo*2, vertex.color=V(g)$color)
print('done')

