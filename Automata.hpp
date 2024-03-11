#ifndef AUTOMATA_HPP
#define AUTOMATA_HPP

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
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

    //////////////////////////////////////////////////////////
    // Replace the current matrix with a given one          //
    // Input: An Eigen::MatrixXd style matrix               //
    //////////////////////////////////////////////////////////
    void set_power_matrix(const MatrixXd new_power_matrix);

    ////////////////////////////////////////////////////////
    // Replaces the current modulo field with a given one //
    // Input: an int > 2                                  //
    ////////////////////////////////////////////////////////
    void set_modulo_field(const double new_modulo_field);

    // Constructors

    // Builds an empty automata with a 0x0 matrix and a 0 modulo field. Matrix elements are 0.
    Automata() : power_matrix(), modulo_field() {};

    // Builds an automata with a rows x columns matrix and a 0 field. Matrix elements are 0.
    Automata(int rows, int columns) : power_matrix(rows, columns), modulo_field() {/* Something goes here */};

    // Builds an automata with a rows x columns matrix and a given field. Matrix elements are 0.
    Automata(int rows, int columns, int field) : power_matrix(rows, columns), modulo_field(field) {/* Something goes here */};

    // Functions
    bool populate_from_user();
    bool populate_from_csv();
    double& operator () (int row, int column);
};

#endif