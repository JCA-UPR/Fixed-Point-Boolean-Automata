// #ifndef AUTOMATA_HPP
// #define AUTOMATA_HPP

#include <iostream>
#include <Eigen/Dense>
#include "Automata.hpp"

using Eigen::MatrixXd;
using namespace Eigen;

#define DEBUG 0


MatrixXd Automata::get_power_matrix() const {
    return power_matrix;
}
double Automata::get_modulo_field() const {
    return modulo_field;
}
int Automata::get_rows() const {
    return this->power_matrix.rows();
}
int Automata::get_columns() const {
    return this->power_matrix.cols();
}

// Setters
void Automata::set_power_matrix(const MatrixXd new_power_matrix) {
    #ifdef DEBUG
        std::cout << "Power Matrix:" << std::endl 
        << power_matrix << 
        std::endl << std::endl 

        << "New Power Matrix:" << std::endl 
        << new_power_matrix << 
        std::endl << std::endl;
    #endif
    power_matrix = new_power_matrix;
}
void Automata::set_modulo_field(const double new_modulo_field) {
    #ifdef DEBUG
        std::cout 
        << "Modulo Field: " << modulo_field << std::endl
        << "New Modulo field: " << new_modulo_field << std::endl;
    #endif
    modulo_field = new_modulo_field;
}

bool Automata::populate_from_user() {
    int rows {}, columns {}, tmp {};

    std::cout << "How many rows or dimentions does the system have?" << std::endl;
    std::cin >> rows;

    if (rows < 0) {
        std::cout << "The number of rows or dimentions should be 1 or larger." << std::endl;
        return false;
    }

    std::cout << "How many columns or variables does the system have?" << std::endl;
    std::cin >> columns;

    if (columns < 0) {
        std::cout << "The number of columnsor variables should be 1 or larger." << std::endl;
        return false;
    } 

    std::cout << "Field(%) of the system (if none type 0): " << std::endl;
    std::cin >> tmp;
    this->set_modulo_field(tmp);

    // Resize the matrix to the new dimentions
    this->power_matrix.resize(rows, columns);

    // Populate the matrix with given inputs
    for (auto row {0}; row < rows; row++) {
        for (auto col {0}; col < columns; col++) {
            std::cout << "Input value for element (" << row << "," << col << "): ";
            std::cin >> this->power_matrix(row,col);
        }
    }

    #if DEBUG == 1 
        std::cout << this->power_matrix << std::endl; 
    #endif

    return true;
}
// bool Automata::populate_from_csv() {

// }

double& Automata::operator () (int row, int column) {
    return this->power_matrix(row,column);
}


