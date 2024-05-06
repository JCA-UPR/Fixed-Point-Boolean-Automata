#ifndef AUTOMATA_HPP
#define AUTOMATA_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <string>

#include "Primes.hpp"

using Eigen::MatrixXd;
using namespace Eigen;

// #define DEBUG 1;

class Automata {
 private:
  MatrixXd power_matrix;

  double modulo_field;

 public:
  // Getters

  // Returns an Eigen::MatrixXd style matrix
  MatrixXd get_power_matrix() const;

  // Returns the field (% operation or Z_p) that the automata is in
  double get_modulo_field() const;

  // Gets the rows from the powers matrix
  int get_rows() const;

  // Gets the columns from the powers matrix
  int get_columns() const;

  // Setters

  ////////////////////////////////////////////////////////
  // Replace the current matrix with a given one        //
  // Input: An Eigen::MatrixXd style matrix             //
  ////////////////////////////////////////////////////////
  void set_power_matrix(const MatrixXd new_power_matrix);

  ////////////////////////////////////////////////////////
  // Replaces the current modulo field with a given one //
  // Input: an int > 2                                  //
  ////////////////////////////////////////////////////////
  void set_modulo_field(const double new_modulo_field);

  // Constructors

  // Builds an empty automata with a 0x0 matrix and a 0 modulo field. Matrix
  // elements are 0.
  Automata() : power_matrix(), modulo_field(2){};

  // Builds an automata with a rows x columns matrix and a 0 field. Matrix
  // elements are 0.
  Automata(int rows, int columns)
      : power_matrix(rows, columns), modulo_field(2){/* Something goes here */};

  // Builds an automata with a rows x columns matrix and a given field. Matrix
  // elements are 0.
  Automata(int rows, int columns, int field)
      : power_matrix(rows, columns),
        modulo_field(field){/* Something goes here */};

  Automata(const Automata& source_automata, double new_modulo);

  // Functions
  bool populate_from_user();
  bool populate_from_csv();
  double& operator()(int row, int column);
  const double& operator()(int row, int column) const;

  friend Automata& operator*(const Automata& automata_1,
                             const Automata& automata_2);
};

double& mod(const double& dividend, const double& divisor);
std::ostream& operator<<(std::ostream& out, const Automata& automata);
bool subset_sum(int find, std::vector<int>& list, std::vector<int>& gen_list);

std::ostream& operator<<(std::ostream& out, const std::vector<double>& vector);
std::ostream& operator<<(std::ostream& out, const std::vector<int>& vector);
Automata& operator-(const Automata& auto_1, const Automata& auto_2);

//////////////////////////////////////////////////////////////////////////
// This function determines wether a given number is a prime or not.    //
//                                                                      //
// Input: a const double& number                                        //
// Output: 1 if true, 0 if false
//////////////////////////////////////////////////////////////////////////
bool is_prime(const double& num);

///////////////////////////////////////////////////////////////////////////
// This function returns whether any given number is a boricua prime.   //
//                                                                      //
// Input: a const double& number                                        //
// Output: a boolean 1 or 0                                             //
///////////////////////////////////////////////////////////////////////////
bool is_boricua_prime(const double& num);

///////////////////////////////////////////////////////////////////////////
// This function takes any number and returns all of its factors.       //
// If the number is prime, returns a vector containing only itself.     //
//                                                                      //
// Input: a double number                                               //
// Output: An std::vector of doubles containing all of the factors      //
// of the given number.                                                 //
///////////////////////////////////////////////////////////////////////////
const std::vector<double>& decompose_boricua(double field);

// std::unordered_map<int,Automata*>& is_fixed_point_matrix(const Automata&
// my_automata);
bool has_lone_factor(const Automata& my_automata, int& factor);

bool is_all_zeros(const Automata& my_automata);

int binary_search_for_real_zero(
    int left_power, int right_power, int midpoint,
    std::vector<int>& calculated_powers_list,
    std::unordered_map<int, Automata*>& calculated_powers_map);

int brute_force(const Automata& my_automata);

std::vector<int>& brute_force_composite(const Automata& automata_1);
Automata& readMatrixFromCSV(const std::string& filePath);

#endif