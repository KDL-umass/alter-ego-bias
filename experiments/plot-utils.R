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
  res <- read.csv("/Users/kavery/workspace/non-cooperative-spillover/results/all-results-facebook-zero.txt") 
  res$n <- 3732
  
  res$bias <- res$ATE.true - res$ATE.adv.gui
  res$est.diff <- res$nonadv.ATE - res$ATE.adv.gui
  res$bias.norm <- res$bias / res$ATE.true
  res$diff.norm <- res$est.diff / res$ATE.adv.gui #res$nonadv.ATE
  
  # res$method <- ifelse(res$method == "random", "random", "dominating")
  
  res$pt.adversaries <- res$index / res$n
  
  plot1 <- ggplot(res, aes(pt.adversaries, diff.norm)) + 
    geom_smooth() + 
    xlab("Adversarial fraction of network") + ylab("Bias in Estimated ATE / Estimated ATE") +  
    theme_bw()+ theme(text = element_text(size = 15)) + theme(legend.position="bottom") +
    guides(color=guide_legend(override.aes=list(fill=NA))) + 
    theme(axis.text.x = element_text(angle = 70, hjust = 1))
  plot(plot1)
}



plot.increase.ATE.bias <- function() {
  #res <- read.csv("adversary-results-revised.csv")
  #res2 <- read.csv("adversary-results-revised-sbm.csv")
  #res <- rbind(res, res2)
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
  res <- read.csv("/Users/kavery/workspace/non-cooperative-spillover/results/sybil/adversary-results-small-world-0.25-4.csv")
  
  res$bias <- res$ATE.true - res$ATE.adv.gui
  res$est.diff <- res$nonadv.ATE - res$ATE.adv.gui
  res$bias.norm <- res$bias / res$ATE.true
  res$diff.norm <- res$est.diff / res$nonadv.ATE
  
  #res <- subset(res, method != "degree")
  res$graph.type <- ifelse(res$graph.type == "sbm", "SBM", as.character(res$graph.type))
  
  res$pt.adversaries <- res$index / res$n
  
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
  
  df <- subset(res, size.of.dom==FALSE & graph.type == "small-world")
  # plot2 <- ggplot(df, aes(pt.adversaries, abs(diff.norm))) + 
  #   geom_smooth() + xlab("Adversarial fraction of network (small world)") + ylim(c(0,1)) + xlim(c(0,0.175)) +
  #   ylab("Bias in Estimated ATE / Estimated ATE") + geom_abline(slope=0) + 
  #   theme_bw() + theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(fill=NA))) + 
  #   theme(axis.text.x = element_text(angle = 70, hjust = 1))
  # plot(plot2)
  
  plot3 <- ggplot(df, aes(pt.adversaries, abs(diff.norm))) + geom_smooth(color="#56B4E9") + 
    xlab("Adversarial fraction of network (small-world)") + ylab("Bias in Estimated ATE / Estimated ATE") + ylim(c(0,1)) + xlim(c(0,1)) +
    geom_abline(slope=0) + theme_bw() + theme(text = element_text(size = 15)) + 
    theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(fill=NA))) + 
    theme(axis.text.x = element_text(angle = 70, hjust = 1))
  plot(plot3) 
  
  df <- subset(res, size.of.dom==FALSE & graph.type == "SBM")

  # plot7 <- ggplot(df, aes(pt.adversaries, abs(diff.norm))) + 
  #   geom_smooth() + geom_ribbon(aes(y = abs(diff.norm), ymin = abs(diff.norm) - sd, ymax = abs(diff.norm) + sd, fill = method), alpha = .2) + 
  #   xlab("Adversarial fraction of network (SBM)") + ylim(c(0,1)) + xlim(c(0,0.175)) +
  #   ylab("Bias in Estimated ATE / Estimated ATE") + geom_abline(slope=0) + 
  #   theme_bw() + theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(fill=NA))) + 
  #   theme(axis.text.x = element_text(angle = 70, hjust = 1))
  # plot(plot7)
  
  # plot8 <- ggplot(df, aes(pt.adversaries, abs(diff.norm))) + geom_smooth(color="#E69F00") + 
  #   xlab("Adversarial fraction of network (SBM)") + ylab("Bias in Estimated ATE / Estimated ATE") + 
  #   geom_abline(slope=0) + theme_bw() + theme(text = element_text(size = 15)) + 
  #   theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(fill=NA))) + 
  #   theme(axis.text.x = element_text(angle = 70, hjust = 1)) + scale_colour_manual(values=cbPalette)
  # plot(plot8) 
  
  df <- subset(res, size.of.dom==FALSE & graph.type == "forest-fire")
  # plot7 <- ggplot(df, aes(pt.adversaries, diff.norm)) + 
  #   geom_smooth() + xlab("Adversarial fraction of network (forest-fire)") + ylim(c(0,1)) + xlim(c(0,0.175)) +
  #   ylab("Bias in Estimated ATE / Estimated ATE") + geom_abline(slope=0) + 
  #   theme_bw() + theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(fill=NA))) + 
  #   theme(axis.text.x = element_text(angle = 70, hjust = 1))
  # plot(plot7)

  # plot8 <- ggplot(df, aes(pt.adversaries, diff.norm)) + geom_smooth(color="#009E73") + 
  #   xlab("Adversarial fraction of network (forest fire)") + ylab("Bias in Estimated ATE / Estimated ATE") + 
  #   geom_abline(slope=0) + theme_bw() + theme(text = element_text(size = 15)) + 
  #   theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(fill=NA)), linetype=guide_legend(override.aes=list(fill=NA))) + 
  #   theme(axis.text.x = element_text(angle = 70, hjust = 1)) + scale_colour_manual(values=cbPalette)
  # plot(plot8) 
  
  # res$graph.type <- factor(res$graph.type, levels=c("SBM", "small-world", "scale-free", "forest-fire"))
  
  # plot5 <- ggplot(subset(res, size.of.dom==FALSE), aes(pt.adversaries, adversary.influence)) + 
  #   geom_smooth() + xlab("Adversarial fraction of network") + ylab("Adversary influence") + 
  #   facet_wrap(~graph.type, scales="free_x") + theme_bw() + theme(text = element_text(size = 15)) + 
  #   theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(fill=NA))) + 
  #   theme(axis.text.x = element_text(angle = 70, hjust = 1))
  # plot(plot5)
  
  # plot6 <- ggplot(subset(res, size.of.dom==FALSE), aes(pt.adversaries, adversary.influence, color=graph.type)) + 
  #   geom_smooth() + xlab("Adversarial fraction of network") + ylab("Adversary influence") + 
  #    theme_bw() + theme(text = element_text(size = 15)) + 
  #   theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(fill=NA))) + 
  #   theme(axis.text.x = element_text(angle = 70, hjust = 1))
  # plot(plot6)
  
}

plot.increase.ATE.bias()