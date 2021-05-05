# StrainXpress
## Description
StrainXpress is a de novo assembly method which base on overlap-layout-consensus (OLC) paradigm and can fast and accurately assemble high complexity metagenome sequencing data at strain resolution. 
## Installation and dependencies
Please note that StrainXpress is built for linux-based systems and python3 only. StrainXpress relies on the following dependencies:
StrainXpress relies on the following dependencies:
- [minimap2](https://github.com/lh3/minimap2)
- *overlap graph construction module from [HaploConduct](https://github.com/HaploConduct/HaploConduct)*
- g++ >=5.5.0 and with boost libraries

To install phasebook, firstly, it is recommended to intall the dependencies through [Conda](https://docs.conda.io/en/latest/):
```
conda create -n strainxpress
conda activate strainxpress
conda install -c bioconda python=3.6 minimap2
```
Subsequently, pull down the code to the directory where you want to install, and compile the code:
```
git clone https://github.com/kangxiongbin/StrainXpress.git
cd StrainXpress
sh install.sh
```
