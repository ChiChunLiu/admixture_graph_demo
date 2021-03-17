
Note that this module is still under development, please use it at your own risk.

This module attempts to perform Markov chain Monte Carlo for the topology of admixture graphs, along with their parameters such as branch lengths, event times etc. The underlying likelihood uses the framework in momi. Please refer to demo for visualization and examples.

![An example MCMC step](illustration_figure.png)


To start:
conda env create -f environment.yml  
conda activate agraph

agraph is a module containing the admixture_graph class built upon the networkx package that stores parameters of past demographic events. mcmc contains utility that perform metropolis-Hasting iterations to obtain posterior distribution of topologies.

