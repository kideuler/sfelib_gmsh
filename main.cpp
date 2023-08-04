#include <set>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include "sfelib_gmsh.hpp"

int main(int argc, char **argv){
    gmsh::initialize();

    gmsh::model::add("circle");

    gmsh::merge("../models/circle");
    vector<size_t> quads;
    vector<double> coords;

    sfelib_gmsh_create_quads({0,0,1,1}, 20, quads, coords);
    cout << "finished creating geometry" << endl;

    gmsh::model::geo::synchronize();
    //gmsh::model::mesh::generate(2);
    gmsh::write("circtest.geo_unrolled");
    gmsh::write("circtest.msh");
    gmsh::finalize();
    return 0;
}