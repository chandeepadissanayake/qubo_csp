# Quantum Annealing Solutions for the Closest String Problem with D-Wave Systems - Codebase

This repository contains the codebase for the paper titled "Quantum Annealing Solutions for the Closest String Problem with D-Wave Systems" published on arXiv (ID: [2310.12852](https://arxiv.org/abs/2310.12852)).

## Abstract

The Closest String Problem is an NP-complete problem which appears more commonly in bioinformatics and coding theory. Less surprisingly, classical approaches have been pursued with two prominent algorithms being the genetic algorithm and simulated annealing. Latest improvements to quantum computing devices with a specialization in optimization tasks such as DWave systems, suggest that an attempt to embed the problem in a model accepted by such systems is worthwhile. In this work, two QUBO formulations have been proposed, with one being a slight modification over the other. Subsequently, an evaluation based on a few simple test cases had been carried out on both formulations. In this regard, the D-Wave annealers have been used, while providing guidelines for optimality on certain platform-specific concerns. For evaluation purposes, a metric termed Occurrence Ratio (OR) has been defined. With minimal hyperparameter tuning, the expected solutions were obtained for every test case and the optimality was guaranteed. To address practical and implementation issues, an inherent decomposition strategy based on the possibility of having substrings has been elucidated to accommodate the restricted qubit count. Conclusively, the need for further investigation on tuning the hyperparameters is emphasized.

## Getting Started

You simply have to install the requirements and run main.py. Just note that you need to have D-Wave SDK setup locally.

```bash
$ git clone https://github.com/chandeepadissanayake/qubo_csp.git
$ cd qubo_csp
$ pip install -r requirements.txt
```

### Setting up the D-Wave SDK

Installing the requirements will automatically build the SDK in the local environment. In order to setup the environment with D-Wave credentials follow the instructions [here](https://docs.ocean.dwavesys.com/en/latest/overview/install.html#set-up-your-environment)

### Configuring and Submitting the Job to D-Wave Leap Cloud

Modify the parameters as necessary in main.py. Afterward, you can simply run the file as follows to submit the job and get the results.

```bash
$ python3 main.py
```

# Citation

If you use this codebase or are inspired by the research presented in the paper, please consider citing the original paper:

```
@misc{dissanayake2023quantum,
      title={Quantum Annealing Solutions for the Closest String Problem with D-Wave Systems}, 
      author={Chandeepa Dissanayake},
      year={2023},
      eprint={2310.12852},
      archivePrefix={arXiv},
      primaryClass={quant-ph}
}
```
