#include <string>
#include <sstream>
#include <armadillo>
#include "3Dsavers.h"


std::string cubeToDf3(const arma::cube &m) {
    std::stringstream ss(std::stringstream::out | std::stringstream::binary);
    uint nx = m.n_rows;
    uint ny = m.n_cols;
    uint nz = m.n_slices;
    ss.put(nx >> 8u);
    ss.put(nx & 0xffu);
    ss.put(ny >> 8u);
    ss.put(ny & 0xffu);
    ss.put(nz >> 8u);
    ss.put(nz & 0xffu);
    double theMin = 0.0;
    double theMax = m.max();
    for (uint k = 0; k < m.n_slices; k++) {
        for (uint j = 0; j < m.n_cols; j++) {
            for (uint i = 0; i < m.n_rows; i++) {
                uint v = 255 * (fabs(m(i, j, k)) - theMin) / (theMax - theMin);
                ss.put(v);
            }
        }
    }
    return ss.str();
}


std::string cubeToRaw(const arma::cube &m) {
    std::stringstream ss(std::stringstream::out | std::stringstream::binary);
    double theMin = 0.0;
    double theMax = m.max();
    for (uint k = 0; k < m.n_slices; k++) {
        for (uint j = 0; j < m.n_cols; j++) {
            for (uint i = 0; i < m.n_rows; i++) {
                uint v = 255 * (fabs(m(i, j, k)) - theMin) / (theMax - theMin);
                ss.put(v);
            }
        }
    }
    return ss.str();
}