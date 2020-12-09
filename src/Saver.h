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

  static void cubeToDf3(const arma::cube &m, std::string filename);
};

#endif // !SAVER_H
