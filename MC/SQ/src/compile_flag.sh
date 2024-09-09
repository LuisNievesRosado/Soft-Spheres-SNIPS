#!/bin/sh

gfortran -Wall -fcheck=all -fbacktrace -c energy.f95
gfortran -Wall -fcheck=all -fbacktrace -c celllists.f95
gfortran -Wall -fcheck=all -fbacktrace -c output.f95
gfortran -Wall -fcheck=all -fbacktrace -c verletlists.f95
gfortran -Wall -fcheck=all -fbacktrace -c moves.f95 
gfortran -Wall -fcheck=all -fbacktrace -c main.f95
gfortran -Wall -fcheck=all -fbacktrace -g -o MC celllists.o verletlists.o output.o energy.o moves.o main.f95