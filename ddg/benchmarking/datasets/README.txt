This directory contains datasets from various sources:

(1) Kellogg et al 2011 dataset:
- Obtained from the Kortemme Lab ddg benchmark repo at https://github.com/Kortemme-Lab/ddg/, downloaded on 20200709
- Paper doi: 10.1002/prot.22921
- file stored as kellogg.csv

(2) Protherm dataset
- Also taken from the Kortemme lab ddg benchmark repo
- compiled by Kortemme Lab (see ddg/input/README.rst and the paper by Ó Conchúir)
- removed the comment character before the column names
- file stored as curatedprotherm.csv

(3) Hahnboem's dataset
- This is the dataset used for the ddG benchmarking in the score function paper (Park et al 2016 JCTC). Hahnboem sent it to me and noted: "The benchmark set is brought from Liz Kellog's paper, which is originally from ProTherm database. I am pretty sure the identical database should be in Tanja's group. Please confirm the entries against mine -- I dropped some of X-to-P mutations in my benchmark." 
- Note that he suggests running with the current master (ref2015) without the -beta_cart option, as it is now standard.
- file stored as ref2015.Lizddg.txt

(4) Brandon's dataset
- This is a further curated protherm dataset from Frenz et al 2020, which tries to balance the number of entries for each type of mutation (i.e. X to Ala, small to large, etc).
- file stored as balanced_benchmark.csv

Kortemme Lab ddg benchmark repo:
- repo at https://github.com/Kortemme-Lab/ddg/
- downloaded on 20200709 `git clone git@github.com:Kortemme-Lab/ddg.git`

References
(1) Kellogg, E. H., Leaver-Fay, A., & Baker, D. (2011). Role of conformational sampling in computing mutation-induced changes in protein structure and stability. Proteins, 79(3), 830–838. 
(2) Ó Conchúir, S., Barlow, K. A., Pache, R. A., Ollikainen, N., Kundert, K., O’Meara, M. J., Smith, C. A., & Kortemme, T. (2015). A Web Resource for Standardized Benchmark Datasets, Metrics, and Rosetta Protocols for Macromolecular Modeling and Design. PloS One, 10(9), e0130433.
(3) Park, H., Bradley, P., Greisen, P., Jr, Liu, Y., Mulligan, V. K., Kim, D. E., Baker, D., & DiMaio, F. (2016). Simultaneous Optimization of Biomolecular Energy Functions on Features from Small Molecules and Macromolecules. Journal of Chemical Theory and Computation, 12(12), 6201–6212.
(4) Frenz, B., Lewis, S., King, I., Park, H., DiMaio, F., & Song, Y. (2020). Prediction of protein mutational free energy: benchmark and sampling improvements increase classification accuracy. In bioRxiv (p. 2020.03.18.989657). https://doi.org/10.1101/2020.03.18.989657
