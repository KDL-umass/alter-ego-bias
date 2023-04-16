library(ggplot2)

prepare.for.plots <- function(g, adversaries, adversary.exposure, treatment.assignments="", labels=FALSE) { 
  V(g)$color <- "lightblue"
  V(g)$color[which(adversary.exposure[1,] > 0)] <- "orange"
  V(g)$color[which(adversaries > 0)] <- "red"
  
  bdr <- rep("black", length(V(g)))
  if(length(treatment.assignments) > 2) bdr <- ifelse(treatment.assignments, "green", "black")
  
  if(labels == FALSE) plot(g, vertex.frame.color=bdr, vertex.label=NA, layout=layout_with_fr) 
  else plot(g, vertex.frame.color=bdr) 
}

plot.realworld.ATE.bias <- function() { 
  res <- read.csv("/Users/kavery/workspace/non-cooperative-spillover/results/results-forest-fire-0.25-1.csv") 
  res$n <- 4039
  
  res$bias <- res$ATE.true - res$ATE.adv.gui
  res$est.diff <- res$nonadv.ATE - res$ATE.adv.gui
  res$bias.norm <- res$bias / res$ATE.true
  res$diff.norm <- res$est.diff / res$ATE.adv.gui #res$nonadv.ATE
  
  # res$method <- ifelse(res$method == "random", "random", "dominating")
  
  res$pt.adversaries <- res$index / res$n
  
  plot1 <- ggplot(res, aes(pt.adversaries, abs(diff.norm))) + 
    geom_smooth() + 
    xlab("Adversarial fraction of network") + ylab("Bias in Estimated ATE / Estimated ATE") +  
    theme_bw()+ theme(text = element_text(size = 15)) + theme(legend.position="bottom") +
    guides(color=guide_legend(override.aes=list(fill=NA))) +  ylim(c(0,1)) + xlim(c(0,0.2)) +
    theme(axis.text.x = element_text(angle = 70, hjust = 1))
  plot(plot1)
}



plot.increase.ATE.bias <- function() {
  #res <- read.csv("adversary-results-revised.csv")
  #res2 <- read.csv("adversary-results-revised-sbm.csv")
  #res <- rbind(res, res2)
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
  res <- read.csv("/Users/kavery/workspace/non-cooperative-spillover/results/mult-results-forest-fire-0.25-2.csv")
  
  res$bias <- res$ATE.true - res$ATE.adv.gui
  res$est.diff <- res$nonadv.ATE - res$ATE.adv.gui
  # res$bias.norm <- res$bias / res$ATE.true
  res$diff.norm <- res$est.diff / res$nonadv.ATE
  
  res$graph.type <- ifelse(res$graph.type == "sbm", "SBM", as.character(res$graph.type))
  
  res$pt.adversaries <- 1 - (res$index / res$n)
  # res$pt.adversaries <- res$index / res$n
  
  # index <- append(res$index, c(res$index[length(res$index)]:1))
  # print(index)
  # pt.adversaries <- 1 - (index / 1000)
 
  # print(pt.adversaries)
  # diff.norm <- append(res$diff.norm, rep(pt.adversaries[length(pt.adversaries)], 1000-nrow(res)))
  # print(diff.norm)

  # print(1000-nrow(res))
  # res[nrow(res)+(1000-nrow(res)),] <- NA
  # res$index <- index
  # res$pt.adversaries <- pt.adversaries
  # res$diff.norm <- diff.norm
  
  # #remove late indices
  # res <- subset(res, !(graph.type == "forest-fire" & index > 350 & res$n==1000))
  # res <- subset(res, !(graph.type == "small-world" & index > 233 & res$n==1000))
  # res <- subset(res, !(graph.type == "SBM" & index > 245 & res$n==1000))
  
  # #res <- subset(res, !(graph.type == "forest-fire" & index > 350))
  # res <- subset(res, !(graph.type == "small-world" & index > 160 & res$n==5000))
  # res <- subset(res, !(graph.type == "SBM" & index > 51 & res$n==5000))
  
  # res <- subset(res, !(graph.type == "forest-fire" & index > 180 & res$n==500))
  # res <- subset(res, !(graph.type == "small-world" & index > 143 & res$n==500))
  # res <- subset(res, !(graph.type == "SBM" & index > 143 & res$n==500))
  
  # df <- subset(res, size.of.dom==FALSE & graph.type == "small-world")
  # df <- subset(res, size.of.dom==FALSE & graph.type == "facebook")
  
  df <- subset(res, size.of.dom==FALSE & graph.type == "forest-fire")

  plot3 <- ggplot(df, aes(pt.adversaries, abs(diff.norm))) + geom_smooth(color="#56B4E9") + geom_point() +
    xlab("Sybil fraction of network") + ylab("Bias in Estimated ATE / Estimated non-Sybil ATE") + 
    geom_abline(slope=0) + theme_bw() + theme(text = element_text(size = 20)) + ylim(c(0,1)) + #xlim(c(0,0.2)) + 
    theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(fill=NA))) + 
    theme(axis.text.x = element_text(angle = 70, hjust = 1)) + scale_colour_manual(values=cbPalette)
  plot(plot3) 
  
  # df <- subset(res, size.of.dom==FALSE & graph.type == "SBM")
  
  # plot5 <- ggplot(df, aes(pt.adversaries, abs(diff.norm))) + geom_smooth(color="#E69F00") + 
  #   xlab("Sybil fraction of network") + ylab("Bias in Estimated ATE / Estimated non-Sybil ATE") + 
  #   geom_abline(slope=0) + theme_bw() + theme(text = element_text(size = 20)) + ylim(c(0,1)) + xlim(c(0,0.5)) +
  #   theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(fill=NA))) + 
  #   theme(axis.text.x = element_text(angle = 70, hjust = 1)) + scale_colour_manual(values=cbPalette)
  # plot(plot5) 
  
  # df <- subset(res, size.of.dom==FALSE & graph.type == "forest-fire")

  # plot8 <- ggplot(df, aes(pt.adversaries, abs(diff.norm))) + geom_smooth(color="#009E73") + 
  #   xlab("Sybil fraction of network") + ylab("Bias in Estimated ATE / Estimated non-Sybil ATE") + 
  #   geom_abline(slope=0) + theme_bw() + theme(text = element_text(size = 20)) + ylim(c(0,1)) + xlim(c(0,0.5)) +
  #   theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(fill=NA))) + 
  #   theme(axis.text.x = element_text(angle = 70, hjust = 1)) + scale_colour_manual(values=cbPalette)
  # plot(plot8) 
  
  # res$graph.type <- factor(res$graph.type, levels=c("SBM", "small-world", "scale-free", "forest-fire"))
  
  # plot5 <- ggplot(df, aes(pt.adversaries, abs(diff.norm))) + geom_smooth(color="#009E73") + 
  #   xlab("Adversarial fraction of network") + ylab("Bias in Estimated ATE / Estimated non-Sybil ATE") + 
  #   geom_abline(slope=0) + theme_bw() + theme(text = element_text(size = 20)) + ylim(c(0,1)) + 
  #   theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(fill=NA))) + 
  #   theme(axis.text.x = element_text(angle = 70, hjust = 1)) + scale_colour_manual(values=cbPalette)
  # plot(plot5)
  
  # plot6 <- ggplot(subset(res, size.of.dom==FALSE), aes(pt.adversaries, adversary.influence, color=graph.type)) + 
  #   geom_smooth() + xlab("Adversarial fraction of network") + ylab("Adversary influence") + 
  #    theme_bw() + theme(text = element_text(size = 20)) + 
  #   theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(fill=NA))) + 
  #   theme(axis.text.x = element_text(angle = 70, hjust = 1))
  # plot(plot6)
  
}

plot.increase.ATE.bias()
# plot.realworld.ATE.bias()