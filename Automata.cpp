// #ifndef AUTOMATA_HPP
// #define AUTOMATA_HPP

#include "Automata.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>

using Eigen::MatrixXd;
using namespace Eigen;

#define DEBUG 0

std::ostream& operator<<(std::ostream& out, const std::vector<double>& vector) {
  for (const auto& num : vector) {
    out << num << ", ";
  }
  return out;
}

std::ostream& operator<<(std::ostream& out, const std::vector<int>& vector) {
  for (const auto& num : vector) {
    out << num << ", ";
  }
  return out;
}

//////////////////////////////////////////////////////////////////////////
// This function determines wether a given number is a prime or not.    //
//                                                                      //
// Input: a const double& number                                        //
// Output: 1 if true, 0 if false                                        //
//////////////////////////////////////////////////////////////////////////
bool is_prime(const double& num) {
  // Calculate √num
  double sq{sqrt(num)};

#if DEBUG == 2
  std::cout << "Square root of " << num << " is : " << sq << std::endl;
#endif

  // Check if num is divisible by primes up to √num or we pass first 10,000
  // primes
  for (int i{0}; (i < 10000) && (primes[i] <= sq); i++) {
    // std::cout << "Here at " << i << std::endl;

#if DEBUG == 2
    std::cout << "Mod of " << num << " and " << primes[i] << " is "
              << mod(num, primes[i]) << std::endl;
#endif

    if (mod(num, primes[i]) == 0) return false;
  }

  return true;
}

///////////////////////////////////////////////////////////////////////////
// This function returns whether any given number is a boricua prime.   //
//                                                                      //
// Input: a const double& number                                        //
// Output: a boolean 1 or 0                                             //
///////////////////////////////////////////////////////////////////////////
bool is_boricua_prime(const double& num) {
  // if is prime
  if (!is_prime(num)) return false;

  // and if num-1 is p1*p2*p3...
  auto minus_one{num - 1};

  // None of the factors of a number P should be larger than √P
  auto sq{sqrt(minus_one)};

  for (auto i{0}; i < 10000 && primes[i] <= sq; i++) {
    // Make sure all factors are at most to the power of 1
    // if (num/x) % x == 0, is not boricua prime
    if (mod((minus_one / primes[i]), primes[i]) == 0) return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////
// This function takes any number and returns all of its factors.       //
// If the number is prime, returns a vector containing only itself.     //
//                                                                      //
// Input: a double number                                               //
// Output: An std::vector of doubles containing all of the factors      //
// of the given number.                                                 //
///////////////////////////////////////////////////////////////////////////
const std::vector<double>& decompose_boricua(double field) {
  std::vector<double>* factors{new std::vector<double>()};

  // If the number is prime, return a list of just itself
  if (is_prime(field - 1)) {
    factors->emplace_back(field);
    return *factors;
  } else {  // If the number is not prime
    // double sq {sqrt(field)+1};
    // Go through each prime ≤ √field in the list of Primes
    double sq{sqrt(field - 1)};
    for (auto i{0}; (i < 10000) && primes[i] <= (field - 1) / 2;
         i++) {  // DANGER HERE!!!!!!!!!! DANGER!!!!!
      if (mod(field - 1, primes[i]) == 0) {  // If the current prime is a factor
        factors->emplace_back(primes[i]);    // Add it to the vector
      }
    }

    // Return the address of the vector containing all factors
    return *factors;
  }
}

MatrixXd Automata::get_power_matrix() const { return power_matrix; }
double Automata::get_modulo_field() const { return modulo_field; }
int Automata::get_rows() const { return this->power_matrix.rows(); }
int Automata::get_columns() const { return this->power_matrix.cols(); }

// Setters
void Automata::set_power_matrix(const MatrixXd new_power_matrix) {
#if DEBUG == 1
  std::cout << "Power Matrix:" << std::endl
            << power_matrix << std::endl
            << std::endl

            << "New Power Matrix:" << std::endl
            << new_power_matrix << std::endl
            << std::endl;
#endif

  power_matrix = new_power_matrix;
}
void Automata::set_modulo_field(const double new_modulo_field) {
#if DEBUG == 1
  std::cout << "Modulo Field: " << modulo_field << std::endl
            << "New Modulo field: " << new_modulo_field << std::endl;
#endif

  // If the number is not a boricua prime, reject it.
  if (!is_boricua_prime(new_modulo_field))
    std::cerr
        << "The modulo field of an Automata has to be a boricua prime number."
        << std::endl;

  // Assign new modulo field
  modulo_field = new_modulo_field;

  // Constrain each element to the new modulo field
  for (int i{0}; i < this->get_rows(); i++) {
    for (int j{0}; j < this->get_columns(); j++) {
      (*this)(i, j) = mod((*this)(i, j), new_modulo_field);
    }
  }
}

bool Automata::populate_from_user() {
  int rows{}, columns{}, tmp{};

  std::cout << "How many rows or dimentions does the system have?" << std::endl;
  std::cin >> rows;

  if (rows < 0) {
    std::cout << "The number of rows or dimentions should be 1 or larger."
              << std::endl;
    return false;
  }

  std::cout << "How many columns or variables does the system have?"
            << std::endl;
  std::cin >> columns;

  if (columns < 0) {
    std::cout << "The number of columnsor variables should be 1 or larger."
              << std::endl;
    return false;
  }

  std::cout << "Field(%) of the system (if none type 0): " << std::endl;
  std::cin >> tmp;
  this->set_modulo_field(tmp);

  // Resize the matrix to the new dimentions
  this->power_matrix.resize(rows, columns);

  // Populate the matrix with given inputs
  for (auto row{0}; row < rows; row++) {
    for (auto col{0}; col < columns; col++) {
      std::cout << "Input value for element (" << row << "," << col << "): ";
      std::cin >> this->power_matrix(row, col);
    }
  }

#if DEBUG == 1
  std::cout << this->power_matrix << std::endl;
#endif

  return true;
}
// bool Automata::populate_from_csv() {

// }

double& Automata::operator()(int row, int column) {
  return this->power_matrix(row, column);
}

const double& Automata::operator()(int row, int column) const {
  return this->power_matrix(row, column);
}

double& mod(const double& dividend, const double& divisor) {
  // if the divisor is 0, throw an error
  if (!divisor) std::cerr << "Cannot divide by 0." << std::endl;

  // Generate the address to be returned
  double* remainder{new double};

  // If divisor is larger than dividend, return dividend
  if (divisor > dividend) {
    *remainder = dividend;
    return *remainder;
  }

  // Remainder = Dividend – (Divisor x Quotient)
  *remainder = dividend - (divisor * floor(dividend / divisor));
  // Using floor() function because integer division is a requirment of Galois
  // Fields

  // Return the remainder
  return *remainder;
}

std::ostream& operator<<(std::ostream& out, const Automata& automata) {
  out << automata.get_power_matrix();
  return out;
}

Automata& operator*(const Automata& automata_1, const Automata& automata_2) {
  // Create memory location for automata
  Automata* new_automata{new Automata};

  // Set the matrix to be the multiplication of both given matrices
  new_automata->set_power_matrix(automata_1.get_power_matrix() *
                                 automata_2.get_power_matrix());

  // Constrain all values to automata_1's modulo field, and assign said modulo
  // field to the automata
  new_automata->set_modulo_field(automata_1.get_modulo_field());

  // return the memory address where the resulting automata is
  return *new_automata;
}

Automata& operator-(const Automata& auto_1, const Automata& auto_2) {
  Automata* tmp{new Automata};
  tmp->set_power_matrix(auto_1.get_power_matrix() - auto_2.get_power_matrix());
  tmp->set_modulo_field(auto_1.get_modulo_field());
  return *tmp;
}

Automata::Automata(const Automata& source_automata, double new_modulo)
    : power_matrix(source_automata.get_power_matrix()),
      modulo_field(new_modulo) {
  // After values have been copied from source_automata

  for (auto i{0}; i < this->get_rows(); i++) {       // for each column
    for (auto j{0}; j < this->get_columns(); j++) {  // and for each row
      (*this)(i, j) =
          mod((*this)(i, j),
              new_modulo);  // Aply the modulo operation under the new_modulo
    }
  }
};

bool is_all_zeros(const Automata& my_automata) {
  for (auto i{0}; i < my_automata.get_rows(); i++) {
    for (auto j{0}; j < my_automata.get_columns(); j++) {
      if (my_automata(i, j) != 0) return false;
    }
  }
  return true;
}

bool subset_sum(int find, std::vector<int>& list, std::vector<int>& gen_list) {
  std::sort(list.begin(), list.end(), std::greater<int>());
// Base case:
#if DEBUG == 5
  std::cout << "Finding: " << find << std::endl;
  std::cout << "list: " << list << std::endl;
  std::cout << "gen_list: " << gen_list << std::endl;
#endif
  if (find == 0) return true;
  if (find < 0) return false;

  // Recursive case
  if (find > 0) {
    for (const auto& num : list) {
      // std::cout << "find-num: "<< (find-num) << std::endl << std::endl;
      if (find - num > 0) {
        gen_list.emplace_back(num);
        return subset_sum(find - num, list, gen_list);
      } else if (find - num == 0) {
        gen_list.emplace_back(num);
        return true;
      }
    }
    return false;
  }
}

// bool subset_sum(int number_to_calculate, std::vector<int>& list_of_numbers,
// std::vector<int>& adds_to_number_list) {
//     std::cout << "Fuck you1"<< std::endl;
//     // Make sure return array is empty
//     // adds_to_number_list.clear();

//     // std::cout << "Fuck you3"<< std::endl;

//     // Make a copy of the given array
//     // std::vector<int> list_of_numbers {list_of_numbers};
//     // Sort it in descending order
//     std::cout << "Fuck you2"<< std::endl;
//     std::sort(list_of_numbers.begin(), list_of_numbers.end(),
//     std::greater<int>());

//     auto tmp = number_to_calculate;
//     auto working = tmp;

//     // While the desired number has not been reached
//     std::cout << "Fuck you"<< std::endl;
//     while (tmp > 0) {
//         std::cout << "tmp: " << tmp << std::endl;
//         // Itereat through the list of numbers
//         for (int i {0}; i < list_of_numbers.size(); i++) {
//             std::cout << "Here" << std::endl;
//             // If you can get closer to the desired number, do so
//             if (tmp - list_of_numbers[i] >= 0) {
//                 // Add the new component to the list
//                 adds_to_number_list.emplace_back(list_of_numbers[i]);
//                 // Update the new desired number
//                 tmp-=list_of_numbers[i];
//                 // Break and go find the new desired number
//                 break;
//             }
//         }
//         // If a full cycle is completed without any progress, there is no way
//         to generate the desired number with the provided list if (tmp ==
//         working) {
//             // Clear the list that was generated so far
//             adds_to_number_list.clear();
//             // Return false
//             return false;
//         }
//     }

//     return true;
// }

int binary_search_for_real_zero(
    int left_power, int right_power, int midpoint,
    std::vector<int>& calculated_powers_list,
    std::unordered_map<int, Automata*>& calculated_powers_map) {
  if (midpoint == left_power || midpoint == right_power) return false;

  std::vector<int>* midpoint_generators{new std::vector<int>};
  // std::cout << "Got here 0" << std::endl;
  // If the number can be generated
  if (subset_sum(midpoint, calculated_powers_list, *midpoint_generators)) {
    // Generate the base matrix
    // std::cout << "Got here 0.3" << std::endl;
    auto mata{*(calculated_powers_map[(*midpoint_generators)[0]])};
    // std::cout << "Got here 0.4" << std::endl;
    auto mods{(calculated_powers_map[1])->get_modulo_field()};
    // std::cout << "Got here 0.4.1" << std::endl;
    auto* tmp_automata{new Automata(mata, mods)};
    // std::cout << "Got here 0.5" << std::endl;
    // For every element in the list
    for (int i{1}; i < midpoint_generators->size(); i++) {
      // Multiply it by the next factor to generate the desired power
      *tmp_automata =
          *tmp_automata * *(calculated_powers_map[(*midpoint_generators)[i]]);
    }

    // Add generated midpoint to map and list
    calculated_powers_map[midpoint] = &*tmp_automata;
    calculated_powers_list.emplace_back(midpoint);
  }

  std::cout << "Got here 1" << std::endl;
  // Empty list
  midpoint_generators->clear();
  // Clear pointer memory
  delete midpoint_generators;

  // If the value is all zero
  if (is_all_zeros(*(calculated_powers_map[midpoint]))) {
    // If midpoint-1 exists
    std::cout << "Surprise motherfucke_3r" << std::endl;
    if (calculated_powers_map.find(midpoint - 1) !=
        calculated_powers_map.cend()) {
      // if it is NOT all Zeros
      if (!is_all_zeros(*calculated_powers_map[midpoint - 1])) return midpoint;

    } else {
      std::cout << "Surprise motherfucker" << std::endl;
      return binary_search_for_real_zero(
          left_power, midpoint, std::midpoint(left_power, midpoint),
          calculated_powers_list, calculated_powers_map);
    }
  } else {
    std::cout << "Surprise motherfucker_2" << std::endl;
    return binary_search_for_real_zero(
        midpoint, right_power, std::midpoint(left_power, midpoint),
        calculated_powers_list, calculated_powers_map);
  }
}

bool has_lone_factor(const Automata& my_automata, int& factor) {
  // vector to hold the powers I have
  std::vector<int>* power_matrices_list{new std::vector<int>};

  // Hashmap to hold the actual matrices I have calculated
  std::unordered_map<int, Automata*>* power_matrices_map{
      new std::unordered_map<int, Automata*>};
  // auto top {my_automata.get_modulo_field()};
  auto top{534};
  auto field{my_automata.get_modulo_field()};
  Automata* tmp{nullptr};

  // Add the original matrix to the map
  (*power_matrices_map)[1] = new Automata(my_automata);
  power_matrices_list->emplace_back(1);
#if DEBUG == 3
  std::cout << "Got to point 1" << std::endl;
#endif
  for (int power{2}, flag{0}; power <= top || flag == 0; power += power) {
    if (power == 2) {
      tmp = &(my_automata * my_automata);
      // std::cout << "Matrix mult: " << std::endl << tmp << std::endl;
      (*power_matrices_map)[power] = new Automata(*tmp, field);
      power_matrices_list->emplace_back(power);

      if (is_all_zeros(*tmp)) {
#if DEBUG == 3
        std::cout << "Got to point 1.5. Power: " << power << std::endl;
#endif
        factor = power;
        return true;
      }
#if DEBUG == 3
      std::cout << "Got to point 2. Modulo: " << my_automata.get_modulo_field()
                << power << ", Matrix: " << std::endl
                << *tmp << std::endl;
#endif

    } else if (power <= top) {
      tmp = &((*(*power_matrices_map)[power / 2]) *
              (*(*power_matrices_map)[power / 2]));

      (*power_matrices_map)[power] = new Automata(*tmp, field);
      power_matrices_list->emplace_back(power);

      if (is_all_zeros(*tmp)) {
        factor = power;
        // binary_search_for_real_zero(power/2, power,
        // std::midpoint(power/2,power), *power_matrices_list,
        // *power_matrices_map);
        return true;
      }
#if DEBUG == 3
      std::cout << "Got to point 3. Modulo: " << my_automata.get_modulo_field()
                << power << ", Matrix: " << std::endl
                << *tmp << std::endl;
#endif
    } else if (power > top) {
      flag = 1;
      // Create a list to hold the numbers that add up to the endpoint
      // (dimension of th automata)
      std::vector<int> end_point;
      // Fill the list with numbers that add up to the endpoint
      subset_sum(top, *power_matrices_list, end_point);

      // Make the base automata out of the highest power automata (first element
      // in the list)
      Automata* endpoint{
          new Automata(*((*power_matrices_map)[end_point[0]]), field)};

      // Multiply that automata for every other element in the list
      for (auto i{1}; i < end_point.size(); i++) {
        *endpoint = *endpoint * (*((*power_matrices_map)[end_point[i]]));
      }

      // Save the endpoint automata to the hashmap
      (*power_matrices_map)[top] = &*endpoint;
      // Free the memory form the ponter that was holding the endpoint automata
      delete endpoint;

      // If the endpoint has non-zero entries, return false and factor = -1
      if (!is_all_zeros((*(*power_matrices_map)[top]))) {
        factor = -1;
        return false;
      } else {
#if DEBUG == 3
        std::cout << "Got to point 4. Modulo: "
                  << my_automata.get_modulo_field() << power
                  << ", Matrix: " << std::endl
                  << *tmp << std::endl;
#endif
      }
    }
  }
// is power <, =, > than top? Calculate this Jimmy!!!
#if DEBUG == 3
  std::cout << "Finished lone factor." << std::endl;
#endif
  return false;
}

int brute_force(const Automata& my_automata) {
  Automata tmp{my_automata};
  int power{1};

  for (; power <= tmp.get_columns();) {
    tmp = tmp * my_automata;
    power++;
    if (is_all_zeros(tmp)) {
      return power;
    }
  }
  return power;
}

Automata& generate_identity_matrix(int n) {
  Automata* tmp{new Automata(n, n, 2)};

  for (int i{0}; i < n; i++) {
    (*tmp)(i, i) = 1;
  }
  return *tmp;
}

std::vector<int>& brute_force_composite(const Automata& automata_1) {
  auto X{automata_1};
  auto B{automata_1 - generate_identity_matrix(automata_1.get_columns())};
  auto power{1};
  std::vector<int>* tmp_vec{new std::vector<int>};

  for (; power <= X.get_columns();) {
    X = X * automata_1;
    power++;
    if (is_all_zeros(X)) {
      tmp_vec->emplace_back(power);
      // std::cout << *tmp_vec << std::endl;
      return *tmp_vec;
    } else if (is_all_zeros(X * B)) {
      tmp_vec->emplace_back(power);
      tmp_vec->emplace_back(0);

      // std::cout << *tmp_vec << std::endl;

      return *tmp_vec;
    }
  }
  return *tmp_vec;
}

// Function to open the file, read its contents, and extract the comma-separated
// string
std::string readCSVFile(const std::string& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return "";
  }

  std::string line;
  std::getline(file, line);
  file.close();

  return line;
}

// Function to parse the comma-separated string and return numbers as a vector
std::vector<int>& parseNumbers(const std::string& line) {
  std::vector<int>* numbers{new std::vector<int>};
  std::istringstream iss(line);
  std::string token;
  while (std::getline(iss, token, ',')) {
    numbers->push_back(std::stoi(token));
  }
  return *numbers;
}

MatrixXd& matrix_from_linear_vector(MatrixXd& my_automata,
                                    const std::vector<int>& line, int N) {
  // Resize the matrix
  my_automata.resize(N, N);

  // Go through the vector
  for (int i{0}, row{0}, column{0}; i < line.size(); i++) {
    // Every N elements, add a row
    if (column < N) {
      my_automata(row, column++) = line.at(i);
    } else if (column == N) {
      column = 0;
      my_automata(++row, column++) = line.at(i);
    }
  }
  return my_automata;
}

bool Automata::populate_from_csv() {
  // Get P and N
  const std::string Path{"data.csv"};
  const std::string Path_2{"data2.csv"};

  std::ifstream file(Path);
  // std::vector<Data> data;

  if (!file.is_open()) {
    std::cerr << "Error opening file: " << Path << std::endl;
    return false;
  }

  std::string line;
  // Skip the header line if present
  std::getline(file, line);

  std::getline(file, line);
  // std::cout << "got here 2" << std::endl;
  std::istringstream ss(line);
  std::string token;

  // Read P number
  std::getline(ss, token, ',');
  this->set_modulo_field(std::stoi(token));

  // Read N number
  std::getline(ss, token, ',');

  this->power_matrix.resize(std::stoi(token), std::stoi(token));

  file.close();

  // Get the matrix

  auto vec{parseNumbers(readCSVFile(Path_2))};

  auto mat{MatrixXd()};

  this->power_matrix =
      matrix_from_linear_vector(mat, vec, this->power_matrix.rows());

  // std::cout << power_matrix << std::endl;

  return true;
}

// MatrixXd readMatrixFromCSV(const std::string& filePath) {
Automata& readMatrixFromCSV(const std::string& filePath) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filePath << std::endl;
        // Return an empty matrix
        // return Automata();
    }

    // Read P and N from the first two lines
    int P, N;
    file >> P >> N;

    // Initialize matrix with size N x N
    auto matrix {new Automata(N, N, P)};
    // MatrixXd matrix(N, N);

    // Read matrix values from the file
    std::string line;
    getline(file, line); // Consume the end-of-line character after reading N
    for (int i = 0; i < N; ++i) {
        getline(file, line); // Read a line from the file
        std::stringstream ss(line);
        std::string cell;
        for (int j = 0; j < N; ++j) {
            if (!getline(ss, cell, ',')) {
                std::cerr << "Error reading value from file." << std::endl;
                // Return partially filled matrix
                return *matrix;
            }
            int value = std::stoi(cell);
            (*matrix)(i, j) = value;
        }
    }

    file.close();

    return *matrix;
}
