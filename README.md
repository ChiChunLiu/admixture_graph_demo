
Note that this module is still under development, please use it at your own risk.

Repo structure:
.
├── agraph
│   └── __pycache__
├── data
├── demo
│   └── __pycache__
└── logs

To start:
conda env create -f environment.yml  
conda activate agraph

agraph is a module containing the admixture_graph class built upon the networkx package that stores parameters of past demographic events. mcmc contains utility that perform metropolis-Hasting iterations to obtain posterior distribution of topologies.

Please refer to demo for visualization and examples.

