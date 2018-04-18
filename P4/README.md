# Project4 Computational Physics
The program used for the project is main.cpp while data_analysis.cpp is a class used to generate data file. All these codes are in the folder program.

main.cpp is a parallelized code used to analyze a phase transition in a two dimensional lattice in a range of temperature



it is possible to change in the program the temperature step size and the temperature range, monte carlo trials and lattice dimension

main.cpp can be compiled as:
mpic++ -O3 -o p4.x *.cpp
and run with 2 processors as:
mpirun -n 2 p4.x

In addition there are also commented parts usful for:
a histogram of the energy configurations
a plot of the accepted configurations
compare the 2x2 solutions with the theoretical ones
