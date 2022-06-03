# PlottingTools

This is a repository for plotting ROOT histograms from TTrees and compare them.
The purpose is to 

Common utility functions and dictionaries are defined in `plottingUtils.py`.

`ProduceRootFileFromTop.py` takes a comma separated list of input text files (`-f`) and creates an output ROOT file (`-o`) with a `TList()` object that contains all corresponding histograms.

`PlotComparison.py` is a script that plots parton level comparisons of fixed-order calculations and Monte-Carlo simulations.
The scripts produces .pdf files as output
