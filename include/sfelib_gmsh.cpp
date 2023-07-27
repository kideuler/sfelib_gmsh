#include "sfelib_gmsh.hpp"

void sfelib_gmsh_create_hexes(const vector<double> &bounds, const int &naxis_points, vector<size_t> &hexes, vector<double> &coords){
    assert(bounds.size() == 6);
    
    size_t nv = naxis_points*naxis_points*naxis_points;
    size_t nelems = (naxis_points-1)*(naxis_points-1)*(naxis_points-1);

    coords.resize(3*nv);
    hexes.resize(8*nelems);

    int n = 0;
    double nh = ((double)naxis_points-1.0);
    double dx = (bounds[3]-bounds[0])/nh;
    double dy = (bounds[4]-bounds[1])/nh;
    double dz = (bounds[5]-bounds[2])/nh;

    
    int np = 0;
    for (int k = 0; k<naxis_points; ++k){
        for(int j = 0; j<naxis_points; ++j){
            for(int i = 0; i<naxis_points; ++i){
                coords[3*np] = bounds[0] + ((double)i)*dx; 
                coords[3*np+1] = bounds[1] + ((double)j)*dy; 
                coords[3*np+2] = bounds[2] + ((double)k)*dz; 
                ++np;
            }
        }
    }

    int eid = 0;
    int ii,jj,kk;
    for (int k = 0; k<naxis_points-1; ++k){
        for(int j = 0; j<naxis_points-1; ++j){
            ii = k*n*n + j*n + 1;
            jj = ii + n;
            kk = ii + n*n;
            for(int i = 0; i<naxis_points-1; ++i){
                hexes[8*eid] = ii;
                hexes[8*eid+1] = ii+1;
                hexes[8*eid+2] = jj+1;
                hexes[8*eid+3] = jj;
                hexes[8*eid+4] = kk;
                hexes[8*eid+5] = kk+1;
                hexes[8*eid+6] = kk+n+1;
                hexes[8*eid+7] = kk+n;
                ++eid;
                ++ii; ++jj; ++kk;
            }
        }
    }
}