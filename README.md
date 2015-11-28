# Markov-Ising-System
This repository contains a basic ising system. The simulation is done using the metropolis method.

COmpile using : g++ -Wall -std=c++11 -o3 -o ... Markov.cpp

Some code in the project is not compatible with windows (I use the system() ) command to call GNUplot.
Also, in order for the plots to be made this project assumes that GNUplot is installed.
In order to make the movie this program assumes that ffmpeg is installed.

If either of those two programs is not installed comment out that section and make the plots, movie with your own favorite plotting programs.
