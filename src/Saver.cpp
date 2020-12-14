#include "Saver.h"

#include <utility>

/**
 * @param d the z values where functions were evaluated
 * @param filename the matrix containing a function on each row
 */
void Saver::saveToCSV(const arma::mat& d, std::string filename)
{
  d.save(std::move(filename), arma::csv_ascii);
}

/**
 *
 */
void Saver::cubeToDf3(const arma::cube &m, const std::string& filename) {
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
    std::ofstream file;
    file.open(filename);
    file << ss.str();
    file.close();
}