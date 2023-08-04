#include "sfelib_gmsh.hpp"

void sfelib_gmsh_create_hexes(const vector<double> &bounds, int naxis_points, vector<size_t> &hexes, vector<double> &coords){
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
            ii = k*naxis_points*naxis_points + j*naxis_points + 1;
            jj = ii + naxis_points;
            kk = ii + naxis_points*naxis_points;
            for(int i = 0; i<naxis_points-1; ++i){
                hexes[8*eid] = ii;
                hexes[8*eid+1] = ii+1;
                hexes[8*eid+2] = jj+1;
                hexes[8*eid+3] = jj;
                hexes[8*eid+4] = kk;
                hexes[8*eid+5] = kk+1;
                hexes[8*eid+6] = kk+naxis_points+1;
                hexes[8*eid+7] = kk+naxis_points;
                ++eid;
                ++ii; ++jj; ++kk;
            }
        }
    }
    return;
}

void sfelib_gmsh_create_quads(const vector<double> &bounds, int naxis_points, vector<size_t> &quads, vector<double> &coords){
    assert(bounds.size() == 4);
    
    size_t nv = naxis_points*naxis_points;
    size_t nelems = (naxis_points-1)*(naxis_points-1);

    coords.resize(3*nv);
    quads.resize(4*nelems);
    vector<bool> inside(nv);

    double nh = ((double)naxis_points-1.0);
    double dx = (bounds[2]-bounds[0])/nh;
    double dy = (bounds[3]-bounds[1])/nh;

    // create nodes
    int np = 0;
    int tag;
    for (int j = 0; j<naxis_points; ++j){
        for(int i = 0; i<naxis_points; ++i){
            coords[3*np] = bounds[0] + ((double)i)*dx; 
            coords[3*np+1] = bounds[1] + ((double)j)*dy; 
            coords[3*np+2] = 0.0; 
            tag = gmsh::model::geo::addPoint(coords[3*np],coords[3*np+1],coords[3*np+2]);
            inside[np] = gmsh::model::isInside(2,1,{coords[3*np],coords[3*np+1],coords[3*np+2]});
            ++np;
        }
    }
    int tag_offset = 1 + tag - np;
    gmsh::model::geo::synchronize();
    cout << "done creating nodes" << endl;


    // create lines
    int v1, v2, v3, v4, nlines = 0;
    for (int j = 0; j<naxis_points; ++j){
        for(int i = 0; i<naxis_points-1; ++i){
            v1 = i + naxis_points*(j) + tag_offset;
            v2 = v1+1;
            tag = gmsh::model::geo::addLine(v1,v2);
            ++nlines;
        }
    }
    int vert = nlines;
    for (int j = 0; j<naxis_points; ++j){
        for(int i = 0; i<naxis_points-1; ++i){
            v1 = j + naxis_points*(i) + tag_offset;
            v2 = v1 + naxis_points;
            tag = gmsh::model::geo::addLine(v1,v2);
            ++nlines;
        }
    }
    int line_offset = 1 + tag - nlines;
    gmsh::model::geo::synchronize();

    for (int i = 0; i< nlines; ++i){
        gmsh::model::geo::mesh::setTransfiniteCurve(i+line_offset,2);
    }
    cout << "done creating lines" << endl;

    // create curve loops
    int cl = 0;
    int L1,L2,L3,L4;
    for (int j = 0; j<naxis_points-1; ++j){
        for (int i = 0; i<naxis_points-1; ++i){
            L1 = i + line_offset + j*(naxis_points-1);
            L2 = vert + line_offset + j + (i+1)*(naxis_points-1);
            L3 = i + line_offset + (j+1)*(naxis_points-1);
            L4 = vert + line_offset + j + i*(naxis_points-1);
            tag = gmsh::model::geo::addCurveLoop({L1,L2,L3,L4},-1,true);
            ++cl;
        }
    }

    int CL_offset = 1 + tag - cl;
    gmsh::model::geo::synchronize();
    cout << "done creating curve loops" << endl;

    // create Plane Surfaces
    int ns=0;
    for (int i = 0; i<cl; ++i){
        tag = gmsh::model::geo::addPlaneSurface({i+CL_offset});
        ++ns;
    }
    int Surf_offset = 1 + tag - ns;
    gmsh::model::geo::synchronize();
    for (int i = 0; i<ns; ++i){
        gmsh::model::geo::mesh::setTransfiniteSurface(i+Surf_offset);
        gmsh::model::geo::mesh::setRecombine(2,i+Surf_offset);
    }
    gmsh::model::geo::synchronize();
    cout << "done creating plane surfaces" << endl;


    // finding initial boundary to compare to later
    gmsh::vectorpair pairs, bdypairs;
    gmsh::vectorpair bdy;
    bdypairs.resize(nelems);
    for (int i = 0; i<nelems; ++i){
        bdypairs[i] = make_pair(2,Surf_offset+i);
    }
    gmsh::model::getBoundary(bdypairs, bdy, true, true, false);
    vector<bool> LineMask(4*nelems);
    for (int i = 0; i<bdy.size(); ++i){
        LineMask[bdy[i].second] = true;
    }
    

    // remove quads that are outside computational domain
    pairs.resize(nelems);
    bdypairs.resize(nelems);
    tag = 1;
    int nout = 0;
    int nin = 0;
    for (int j = 0; j<naxis_points-1; ++j){
        for (int i = 0; i<naxis_points-1; ++i){
            v1 = i + j*(naxis_points);
            v2 = i+1 + (j)*(naxis_points);
            v3 = i+1 + (j+1)*(naxis_points);
            v4 = i + (j+1)*(naxis_points);
            if (inside[v1] || inside[v2] || inside[v3] || inside[v4]){
                pairs[nout] = (make_pair(2,tag+1));
                ++nout;
            } else {
                bdypairs[nin] = (make_pair(2,tag+1));
                ++nin;
            }
            tag++;
        }
    }
    pairs.resize(nout);
    bdypairs.resize(nin);
    cout << "done tagging quads outside" << endl;
    gmsh::model::geo::remove(pairs, true);
    gmsh::model::geo::synchronize();
    cout << "done removing quads outside domain" << endl;
    
    
    gmsh::model::getBoundary(bdypairs, bdy, true, true, false);
    int nbdy = bdy.size();
    vector<int> loop(nbdy+line_offset-1);
    int nb = 0;
    for (int i = 0; i<nbdy; i++){
        cout << "(" << bdy[i].first << "," << bdy[i].second << ")" << endl;
        if (LineMask[bdy[i].second]){
            loop[nb] = bdy[i].second;
            ++nb;
        }
    }
    for (int i = 0; i<(line_offset-1); i++){
        loop[nb] = -(i+1);
        ++nb;
    }
    cout << tag_offset << endl;
    cout << "done finding boundary" << endl;

    tag = gmsh::model::geo::addCurveLoop(loop,-1,true);
    gmsh::model::geo::addPlaneSurface({tag});
}