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

adversary.params <- list()
adversary.params$setting <- "dominating"
adversary.params$max <- TRUE
adversary.params$weighting <- "degree"

adversaries.deg <- unlist(list(determine.adversaries(graph.properties, adversary.params)))
adversaries <- matrix(0,1,graph.properties$n)
adversaries[which(adversaries.deg==1)[1:x]] <- 1

random.adversaries <- sample(1:graph.properties$n, 11, replace=FALSE)
adversaries <- matrix(0,1,graph.properties$n)
adversaries[random.adversaries] <- 1

print(adversaries)

treatment <- treatment.assignment(graph.properties$g, clusters)
treatment.assignments <- treatment[clusters]

uncovered.vertices <- 1 - adversaries %*% graph.properties$adj - adversaries
pt.uncovered <- sum(uncovered.vertices == 1)/graph.properties$n
prepare.for.plots(g, adversaries, adversary.params$adversary.exposure, treatment.assignments, labels=TRUE)


exposure.params <- exposure.probs(adversary.params, graph.properties, treatment.assignments, adversaries)
print(exposure.params$influence.as.adversary)
#V(g)$color <- adversaries.deg

g$palette <- grey.colors(100)
V(g)$color <- exposure.params$adversary.exposure.neighbors * 100
#V(g)$color[which(exposure.params$adversary.exposure.neighbors) == 0] <- "lightblue"
V(g)$color[which(adversaries > 0)] <- "red"


lo <- layout_with_kk(g) # create a layout
lo <- norm_coords(lo, ymin=-1, ymax=1, xmin=-1, xmax=1)
plot(g, layout=lo*2, vertex.color=V(g)$color)
print('done')

