CC=mpicxx
CXXFLAGS=-std=c++11 -Wall -O2

default: ldcs

LidDrivenCavitySolver.o: LidDrivenCavitySolver.cpp
	mpicxx LidDrivenCavitySolver.o LidDrivenCavitySolver.cpp    
LidDrivenCavity.o: LidDrivenCavity.cpp
	mpicxx LidDrivenCavity.o LidDrivenCavity.cpp

ldcs: LidDrivenCavitySolver.o LidDrivenCavity.o
	mpicxx -o ldcs LidDrivenCavitySolver.o LidDrivenCavity.o
