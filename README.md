# Inverse_Design_of_Polymers
Inverse Design of Polymers via Evolutionary Algorithm and Bayesian optimization
## Description
We propose a robust machine learning workflow for the inverse design of high thermal conductivity (TC) polymers. Our work starts from a computational database containing 1144 polymers with known TCs. Using those data, we construct a surrogate deep neural network model for TC calculation and extract a polymer-unit library with 32 sequences. Two state-of-the-art algorithms of unified non-dominated sorting genetic algorithm III (U-NSGA-III) and q-noisy expected hypervolume improvement (qNEHVI) are employed for sequence-controlled polymer design. They are the multi-objective evolutionary algorithm and multi-objective Bayesian optimization algorithm, respectively, since the synthesizability of the emerging polymers is also evaluated using the synthetic accessibility score. The approach proposed is flexible and universal, and can be extended to the design of polymers with other property targets. Please refer to our work "Machine learning-assisted inverse design of sequence-controlled high intrinsic thermal conductivity polymers" for additional details.
![Workflow](https://github.com/SJTU-MI/Inverse_Design_of_Polymers/blob/main/Workflow.png)
## Installation

### Files loading and environment setup:

To download, clone this repository:<br>
````
git clone https://github.com/SJTU-MI/Inverse_Design_of_Polymers.git
````

To run most code in this repository, the relevant anaconda environment can be installed from environment.yml. To build this environment, run：<br>
````
cd ./Inverse_Design_of_Polymers
conda env create -f environment.yml
conda activate IDPoly
````
### Try the desired parts of the project:
#### Code in the sources folder
**Cal_TC.py**: A script for calculating thermal conductivity of polymers, here predicted by a DNN surrogate model <br>
**Cal_SA.py**: A script for evaluating the synthetic accessibility scores of polymers <br>
**utility.py**: Some utility functions <br>
#### Code in the results folder
**MOEA_candidates.csv**: Polymers designed by MOEA <br>
**MOEA_HV.csv**: Convergence of the MOEA based on hypervolumes <br>
**MOBO_candidates.csv**: Polymers designed by MOBO <br>
**MOBO_HV.csv**: Convergence of the MOBO based on hypervolumes <br>
#### Case Studies
**MOEA_Case.ipynb**: A case for MOEA <br>
**MOBO_Case.ipynb**: A case for MOBO <br>
**visualization.ipynb**: Visualization of results for MOEA and MOBO cases <br>
## Authors

| **AUTHORS** |Xiang Huang, Shenghong Ju            |
|-------------|--------------------------------------------------|
| **VERSION** | V1.0 / October,2023                               |
| **EMAILS**  | shenghong.ju@sjtu.edu.cn                         |

## Related projects
- pymoo(Multi-objective Optimization in Python)[[Link]](https://github.com/anyoptimization/pymoo)
- botorch(Bayesian optimization in PyTorch)[[Link]](https://github.com/pytorch/botorch)
- RadonPy(automated physical property calculation using all-atom classical molecular dynamics simulations for polymer informatic)[[Link]](https://github.com/RadonPy/RadonPy)

## Attribution

This work is under BSD-2-Clause License. Please, acknowledge use of this work with the appropiate citation to the repository and research article.
