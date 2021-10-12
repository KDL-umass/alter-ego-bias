This markdown provides instructions for running the experiments in Stick it to The Man: Correcting for Non-Cooperative Behavior in Network A/B Testing. 


## Setup 
1. Download R for your system: [https://cran.r-project.org/](https://cran.r-project.org/)

2. Clone the repo:
 `git clone git@github.com:KDL-umass/Non-cooperative-spillover.git`

3. Set up R environment. Run install script and verify setup: `Rscript install/install.R`

## Configure Experiment Settings

1. Set up experiment configs. More details at [this README](experiments/configs/README.md).

4. Set up graph dir. 
* Synthetic Corpus: unzip graph corpus
* Real-World: run load scripts for SNAP library graphs: 
* Or alteratively, generate synthetic graphs with graph generation algorithms. More info and additional options at [this README](graphs/README.md).


## Run the Experiments

## Plots from Experiment Data
### Topology and non-cooperative influence
### Non-adversarial bias in ATE estimation