higgsino
==============

Plotting and tabulating utilities for the RA2b Higgsino analysis based on 
[`manuelfs/analysis_code`](https://github.com/manuelfs/analysis_code) and 
[`ald77/ra4_macros`](https://github.com/ald77/ra4_macros).


#### Plotting
The plot engine is in `src/utilities_macros.cpp`, and the interface in files like `sr/plot_basic.cxx`.
To produce the plots, make sure the file is pointing to the ntuples that you have, and run

    ./compile.sh && ./run/plot_basic.exe

The .pdf files will appear in subfolders of `plots`


#### Tables
The two basic scripts to make tables with background and signal yields are `src/table_cutflow.cxx`,
which produces a cutflow and row with N-1 cuts, and `src/table_regions.cxx`, which produces a table
with the yields in the different sideband and signal regions.

These scripts are run with

    ./compile.sh && ./run/table_cutflow.exe
    ./compile.sh && ./run/table_regions.exe
