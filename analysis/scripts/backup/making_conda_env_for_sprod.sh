# make a new conda environment in the current directory - in my external ssd
conda create --prefix sprod

# activate the conda environment in the current directory - or specify full path
conda activate ./sprod

# I added the conda-forge as channel source for installing tools
conda config --append channels conda-forge

# Install the required tools
conda install r-base=4.0 r-optparse r-distances r-dplyr r-gtools r-mvtnorm r-ggplot2 r-umap
conda install python=3.7 pandas scikit-image=0.17 numpy scipy scikit-learn umap-learn

