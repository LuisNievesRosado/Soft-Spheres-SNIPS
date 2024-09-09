#!/bin/sh

gfortran -c energy.f95
gfortran -c celllists.f95
gfortran -c output.f95
gfortran -c verletlists.f95
gfortran -c moves.f95 
gfortran -c main.f95
gfortran -g -o MC celllists.o verletlists.o output.o energy.o moves.o main.f95