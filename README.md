
[![License: GPL-3.0](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fpbagos%2Fmage%2Fedit%2Fmaster%2FREADME.md&count_bg=%2379C83D&title_bg=%23555555&icon=latex.svg&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)


# MAGE
MAGE :: Meta-Analysis of Gene Expression

## Web tool
MAGE (v2.0) web tool : http://rs.dib.uth.gr:3839/MAGE/ and http://195.251.108.185:3838/MAGE/

## Citation
Tamposis, I.A.; Manios, G.A.; Charitou, T.; Vennou, K.E.; Kontou, P.I.; Bagos, P.G. MAGE: An Open-Source Tool for Meta-Analysis of Gene Expression Studies. Biology 2022, 11, 895. https://doi.org/10.3390/biology11060895

## Introduction
MAGE, is an acronym for Meta-Analysis of Gene Expression.

The overall aim of this work has been to develop a software tool that would offer a large collection of meta-analysis options, as well as several extensions to evaluate the software applied to various biological problems. The MAGE framework is characterized by:
Speed: It takes a small amount of time to perform the functions which are included 
Effectiveness: It gives reliable estimations and results, thanks to the mathematical models which are implemented.
Compatibility: It can be executed either from a Windows or a UNIX operating system.             

Also, there is no need for the user to be expert of any programming or computer science knowledge to run MAGE.

MAGE is a Python package that can be run from the command line. MAGE is written in Python (ver. 3.7.9) and requires the following Python libraries and packages to run:

Requirements
- Pandas~=1.2.4
- Numpy~=1.19.5
- Matplotlib~=3.3.4
- scipy~=1.6.2
- statsmodels~=0.12.2
- PythonMeta~=1.23
- requests~=2.25.1
- statistics~=1.0.3.5
- Seaborn~=0.11.1


## Installation guide

1)	Download MAGE from: https://github.com/pbagos/mage

    Otherwise, you can run MAGE from its online infrastructure at: http://www.compgen.org/tools/mage (Mozila Firefox browser is suggested)
  	MAGE (v2.0) web tool : http://rs.dib.uth.gr:3839/MAGE/ and http://195.251.108.185:3838/MAGE/

3)	After downloading the .zip folder of MAGE from GitHub, extract it to a working directory. 

4)	Το install the requirements, pip needs to be installed. Download the script for pip, from: https://bootstrap.pypa.io/get-pip.py.

5)	Open a terminal/command prompt, cd to the folder containing the get-pip.py file and run:
    
    python get-pip.py

6)	To install the mentioned requirements with pip, open a terminal/command prompt and run:
    
    pip install -r /path/to/requirements.txt

7)	To execute MAGE, execute with:
   

    python mage.py  -c conf.txt  -o results/

## Arguments and options
MAGE provides the following command-line arguments:
  
  -c: The configuration (.txt) file which contains the settings selected from the user.
  
  -o: The output file where the user wants to store the results extracted from MAGE


## Methods
MAGE is consisted of three basic functions. 

#### GISU (Gene ID / Symbol update)
MAGE uses an optional component called GISU to transform the platform's probe identifiers to gene symbols identifiers. These can be helpful when one is comparing datasets aris-ing from different platforms, then the probe identifiers must be con-verted to gene identifiers. Considering that multiple probes may corre-spond to the same gene in a microarray experiment [1], the multiple entries of the same gene can be combined into one using the minimum, maximum or arithmetic mean (average) [1,2,3]. If the experiment's platform is not included in the list, the user can upload the platform file in order to proceed to the transformation.
#### Meta-analysis: 
For meta-analysis, the package supports the standard meta-analysis, bootstrap meta-analysis and multivariate meta-analysis functions.

#### Enrichment analysis:  
Furthermore, the software uses g: Profiler tool [4] to perform functional enrichment analysis with a given gene list produced from the meta-analysis by using the implemented python module (http://biit.cs.ut.ee/gprofiler/). The software returns multiple files con-taining an output file of gene definitions, a file with statistically signifi-cant enriched GO terms, biological pathways, regulatory motifs in DNA, or phenotype ontologies that these genes are highly enriched and provides to the user the option to visualize results with a Manhattan or a heatmap plot.

## References
1. Ramasamy A, Mondry A, Holmes CC, Altman DG. Key issues in conducting a meta-analysis of gene expression microarray datasets. PLoS Med. 2008;5(9):e184. doi:10.1371/journal.pmed.005018
2. Li, Q., N. J. Birkbak, B. Gyorffy, Z. Szallasi and A. C. Eklund (2011). "Jetset: selecting the optimal microarray probe set to represent a gene." BMC Bioinformatics 12: 474.
3. Warnat, P., R. Eils and B. Brors (2005). "Cross-platform analysis of cancer microarray data improves gene expression based classification of phenotypes." BMC Bioinformatics 6: 265.
4. Uku Raudvere, Liis Kolberg, Ivan Kuzmin, Tambet Arak, Priit Adler, Hedi Peterson, Jaak Vilo: g:Profiler: a web server for functional enrichment anal-ysis and conversions of gene lists (2019 update) Nucleic Acids Research 2019; doi:10.1093/nar/gkz369
