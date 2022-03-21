This directory contains the source code (in Python) of FORGE model and its output files (in NetCDF format) for this study, including the three sets of global simulations (S0, S1, S2). The corresponding input files for the FORGE model can be downloaded at https://doi.org/10.6084/m9.figshare.14995320.v2 (total file size about 10 GB). See detailed descriptions in the Methods in Zhu et al., 2021: https://www.nature.com/articles/s41559-021-01548-3 

The cpu-time for a global (2-degree spatial resolution) simulation with 300 years (loop over the 10-years' inputs) is about 6 hours with 8 cores in parallel. To facilitate a quick single-point test, there is also an input file for one site (41.62N, 122.7W).

To run the model, change the path to your own directory in the file run.py, which calls main.py.

Contact: zhudan@pku.edu.cn
