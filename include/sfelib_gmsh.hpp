#ifndef SFELIB_GMSH
#define SFELIB_GMSH

#include <gmsh.h>
#include <cassert>
#include <iostream>

using namespace gmsh::model;
using namespace std;

void sfelib_gmsh_create_hexes(const vector<double> &bounds, int naxis_points, vector<size_t> &hexes, vector<double> &coords);

void sfelib_gmsh_create_quads(const vector<double> &bounds, int naxis_points, vector<size_t> &quads, vector<double> &coords);

#endif
