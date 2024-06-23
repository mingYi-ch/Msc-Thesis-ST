
mamba create -n c2l_2 ipykernel python=3.9

source activate c2l_2
ipython kernel install --user --name c2l_2

# in notebook
!pip install cell2location[tutorials]

data/.local/share/jupyter/kernels/c2l_2/kernel.json