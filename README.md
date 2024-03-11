# The Effect of Alter Ego Accounts on A/B Tests in Social Networks
This repo contains the code for the paper "The Effect of Alter Ego Accounts on A/B Tests in Social Networks." 

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.10806591.svg)](http://dx.doi.org/10.5281/zenodo.10806591)

## Abstract

>Social network users often maintain multiple active accounts, sometimes referred to as *alter egos*. Examples of alter egos include personal and professional accounts or named and anonymous accounts. If alter egos are common on a platform, they can affect the results of A/B testing because a user’s alter egos can influence each other. For a single user, one account may be assigned treatment, while another is assigned control. Alter-ego bias is relevant when the treatment affects the individual user rather than the account. 

>Through experimentation and theoretical analysis, we examine the worst and expected case bias for different numbers of alter egos and for a variety of network structures and peer effect strengths. We show that alter egos moderately bias the results of simulated A/B tests on several network structures, including a real-world Facebook subgraph and several types of synthetic networks: small world networks, forest fire networks, stochastic block models, and a worst-case structure. We also show that bias increases with the number of alter egos and that different network structures have different upper bounds on bias.  

## Run the Experiments

[REPLICATION.md](REPLICATION.md) contains installation, configuration, and run instructions for experiments in the paper.

## Overview
Given a network graph and simulation parameters identified in a configuration file, 
* Cluster-randomize treatment assignment
* Generate a set of alter egos
* Simulate the outcomes per N alter egos 
* Estimate the outcomes with a linear estimator
* Calculate the bias

#### Experiment Configs
Options: 
* Network structure for A/B test
* Alter ego network placement (random, greedy)
* Outcome simulation parameters

More information is provided in  [this README](experiments/configs/README.md).

#### Network Graph Resources
* **Synthetic Corpus**:
We provide a corpus of graph types used for the simulation studies in the paper, generated under the same set of parameters. 
* **Real-World**:
We use a real-world network released in the SNAP Library, which can be found here: https://snap.stanford.edu/data/ego-Facebook.html.  

Alternatively, run the scripts used to generate graphs in the corpus to generate new graphs. 
More details at [this README](REPLICATION.md).
