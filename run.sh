#!/bin/zsh
echo 'Compiling markers data and computing features ...'
Rscript compile.R
echo 'Visualizing markers and marker features ...'
Rscript visualize_marx.R