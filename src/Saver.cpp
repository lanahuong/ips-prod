#include "Saver.h"

/**
 * @param z the z values where functions were evaluated
 * @param f the matrix containing a function on each row
 */
void Saver::saveToCSV(arma::rowvec z, arma::mat f)
{
  // Namefile
  std::string filename = "tmp/schrodinger_solutions.csv";

  arma::mat z_and_functions = arma::mat(f.n_rows + 1, f.n_cols);
  z_and_functions.row(0) = z;
  z_and_functions.rows(1, z_and_functions.n_rows - 1) = f;

  z_and_functions.save(filename, arma::csv_ascii);
}
