#include <set>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include "sfelib_gmsh.hpp"

int main(int argc, char **argv){
    gmsh::initialize();

    gmsh::model::add("sphere");

    gmsh::merge("../sphere.msh");

    vector<size_t> hexes;
    vector<double> coords;
    vector<double> bounds = {0.0,0.0,0.0,1.0,1.0,1.0};
    sfelib_gmsh_create_hexes(bounds, 20, hexes, coords);
    int nv = coords.size() / 3;
    vector<size_t> nodetags(nv);
    for (int i = 0; i<nv; ++i){
        nodetags[i] = i+1;
    }
    gmsh::model::mesh::addNodes(3,1, nodetags, coords);
    gmsh::model::mesh::addElementsByType(1, 5, {}, hexes);
    gmsh::model::mesh::createGeometry();
    
    gmsh::option::setNumber("Mesh.MeshSizeMax", 0.1);
    gmsh::model::geo::synchronize();
    double xmin, ymin, zmin, xmax, ymax, zmax;

    gmsh::model::getBoundingBox(-1,-1, xmin, ymin, zmin,xmax, ymax, zmax);

    cout << xmin << ", " << ymin << ", " << zmin << " | " << xmax << ", " << ymax << ", " << zmax << endl;
    //gmsh::model::mesh::generate(2);

    //gmsh::write("../sphere.msh");
    gmsh::finalize();
    return 0;
}