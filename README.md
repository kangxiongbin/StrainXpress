# StrainXpress
## Description
StrainXpress is a de novo assembly method which base on overlap-layout-consensus (OLC) paradigm and can fast and accurately assemble high complexity metagenome sequencing data at strain resolution. 
## Installation and dependencies
Please note that StrainXpress is built for linux-based systems and python3 only. StrainXpress relies on the following dependencies:
StrainXpress relies on the following dependencies:
- [minimap2](https://github.com/lh3/minimap2)
- [HaploConduct](https://github.com/HaploConduct/HaploConduct)*
- g++ >=5.5.0 and with boost libraries

To install strainxpress, firstly, it is recommended to intall the dependencies through [Conda](https://docs.conda.io/en/latest/):
```
conda create -n strainxpress
conda activate strainxpress
conda install -c bioconda python=3.6 scipy pandas minimap2
```
Subsequently, pull down the code to the directory where you want to install, and compile the code:
```
git clone https://github.com/kangxiongbin/StrainXpress.git
cd StrainXpress
sh install.sh
```
## Examples
```
- Illumina miseq (reads length 2X250bp)
python ../scripts/strainxpress.py -fq all_reads.fq

- Data set is big

python ../scripts/strainxpress.py -fq all_reads.fq -fast

```
## Possible issues during installation (optional)

- If `g++` version of the system is not satisfied, one could try this to install:
```
conda install -c conda-forge gxx_linux-64=7.3.0
# replace the /path/to/ with your own path
ln -s /path/to/miniconda3/envs/strainxpress/bin/x86_64-conda-cos6-linux-gnu-g++ /path/to/miniconda3/envs/strainxpress/bin/g++
ln -s /path/to/miniconda3/envs/strainxpress/bin/x86_64-conda-cos6-linux-gnu-gcc /path/to/miniconda3/envs/strainxpress/bin/gcc
```
- If `boost` library is not installed, you could try this to install:
```
conda install -c conda-forge boost
# set envionment variables
export LD_LIBRARY_PATH=/path/to/miniconda3/envs/strainxpress/lib/:$LD_LIBRARY_PATH
export CPATH=/path/to/miniconda3/envs/strainxpress/include/:$CPATH
```

- If compile error occurs something like `/path/to/miniconda3/envs/strainxpress/x86_64-conda_cos6-linux-gnu/bin/ld: cannot find -lboost_timer `
or `cannot find -lgomp`, 
 which means it fails to link `boost` or `libgomp` library, one could try this to solve:
```
ln -s /path/to/miniconda3/envs/strainxpress/lib/libboost_* /path/to/miniconda3/envs/strainxpress/x86_64-conda_cos6-linux-gnu/lib/.
ln -s /path/to/miniconda3/envs/strainxpress/lib/libgomp* /path/to/miniconda3/envs/strainxpress/x86_64-conda_cos6-linux-gnu/lib/.
# then re-complile and install
sh install.sh
```
