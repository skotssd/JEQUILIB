# JEQUILIB
Jupyter file using matlab kernel and PHREEQC (optional) to solve equilibrium problems formulated in tableau notation.

Install anaconda or install miniconda

to add matlab kernel

FOR MATLAB (based on http://www.jmlilly.net/jupyter-matlab)

conda create -vv -n jmatlab python=3.8.8 jupyter
conda activate jmatlab
conda install -c conda-forge jupyterlab
pip install matlab_kernel
python -m matlab_kernel install
jupyter kernelspec list (to verify)

still in jmatlab …

cd /Applications/MATLAB_R2021a.app/extern/engines/python
python setup.py install

but find the right directory for your installation

and might get an error about privileges.  If you do then type this instead after changing the directory

python setup.py build -b C:\Temp install


