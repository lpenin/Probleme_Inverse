#include<iostream>
#include<mpi.h>
#include<cmath>
#include<vector>


double dot(std::vector<double> u, std::vector<double> v);

std::vector<double> CG(std::vector<std::vector<double> > Aloc, std::vector<double> bloc, std::vector<double> x0loc , double err, int kmax, int nx, int ny);
