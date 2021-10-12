# Network A/B Testing with Non-Cooperative Behavior 
This repo contains the code for our work on A/B testing in networks with non-cooperative behavior. 

Cite our paper at USENIX 2022: 
```
@inproceedings{clary2021noncooperative, 
    author={Clary, Kaleigh and Tosch, Emma and Onaolapo, Jeremiah and Jensen, David D.}, 
    title={Stick it to {T}he {M}an: Correcting for Non-Cooperative Behavior in {A}/{B} Network Testing},
    journal={{USENIX} Security},
    year={2022}
}
```

## Run the Experiments

[REPLICATION.md](REPLICATION.md) contains installation, configuration, and run instructions for experiments and figures in the paper.

## Overview
Given a network graph and simulation parameters identified in a configuration file, 
* Cluster-randomized treatment assignment
* Non-cooperative participant network set construction
* Outcome simulation per n adversaries up to dominating set
* Model fit and outcome estimation with linear estimator
* Bias calculation in sample

#### Experiment Configs
Options: 
* Network structure for A/B test
* Non-cooperative participant network placement (random RNCP, greedy GNCP)
* Outcome simulation parameters
* Non-cooperative outcome models 

More information is provided in  [this README](experiments/configs/README.md).

#### Network Graph Resources
We provide some existing graph structures for the experiments:
* **Synthetic Corpus**:
We provide a corpus of graph types used for the simulation studies in the paper, generated under the same set of parameters. 
* **Real-World**:
We use networks released in the SNAP Library. 
* **Custom**:
Instructions for adding custom graph source.

Alternatively, run the scripts used to generate graphs in the corpus to generate new graphs. 
More details at [this README](graphs/README.md).

