# geneticAlgorithmNetworkHamiltonian
Genetic algorithm for automated network Hamiltonian model discovery

This code implements the genetic algorithm described in the article:

"A Genetic Algorithm for Automated Parameterization of Network Hamiltonian 
                    Models of Amyloid Fibril Formation"

by: Gianmarc Grazioli, Andy Tao, Inika Bhatia, and Patrick Regan 

currently under review at Digital Discovery - Royal Society of Chemistry
preprint: https://chemrxiv.org/engage/chemrxiv/article-details/64a3d048ba3e99daef7c671a

This main function carries out automated discovery of network Hamiltonian models
that can self-assemble into the 2-ribbon amyloid fibril topological structure.
The example used in the main function is very similar to the code used to generate 
the figure in the aforementioned article showing the evolution of different 
generations of models as the genetic algorithm converges on a region of parameter 
space that produces maximal fibril yield for 2-ribbon type amyloid fibril structures.
To give an idea of run time, this main function runs in under 10 minutes on a 
2019 MacBook Pro with a 2.4 GHz 8-Core Intel Core i9 processor and 64 GB of RAM.

Code written by Gianmarc Grazioli 
