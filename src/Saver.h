/**
 * @file Saver.h
 *
 * This file contains the Saver class that saves a function to CSV.
 */
#ifndef SAVER_H
#define SAVER_H

#include <armadillo>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

/**
 * @class Saver
 *
 * This class saves a function to a CSV file
 */
class Saver {
 public:
  /**
   * @brief Save matrix \a d to CSV
   */
  static void saveToCSV(const arma::mat&, std::string);

  /**
   *
   * @param m the cube to save to df3
   * @param filename the name of the df3 file (must end with .df3)
   */
  static void cubeToDf3(const arma::cube &m, const std::string& filename);
};

#endif // !SAVER_H
