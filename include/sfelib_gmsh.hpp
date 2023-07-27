#ifndef SFELIB_GMSH
#define SFELIB_GMSH

#include <gmsh.h>
#include <cassert>

using namespace gmsh::model;
using namespace std;

void sfelib_gmsh_create_hexes(const vector<double> &bounds, const int &naxis_points, vector<size_t> &hexes, vector<double> &coords);


#endif
