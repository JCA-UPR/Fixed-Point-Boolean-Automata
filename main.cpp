#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>

#include "Automata.hpp"

using Eigen::MatrixXd;
using namespace Eigen;
using namespace std;

#define DEBUG 0
// #define SHOW

std::vector<Automata*>& create_field_prime_factor_automatas(
    const Automata& my_automata) {
  std::vector<double> prime_factors{
      decompose_boricua(my_automata.get_modulo_field())};
  auto size{prime_factors.size()};

  // Create a vector that will hold all of the matrices
  std::vector<Automata*>* components{new std::vector<Automata*>};

  // For each prime factor of the given automatas boricua number (so n-1, not
  // n):
  for (auto i{0}; i < size; i++) {
    // Copy the given automata, but assign it the prime_factor[i] as the new
    // modulo field
    components->emplace_back(new Automata(my_automata, prime_factors[i]));
  }

  return *components;
}

bool is_fixed_point(const Automata& MyAutomata) {
  if (is_boricua_prime(MyAutomata.get_modulo_field())) {
    auto factors{decompose_boricua(MyAutomata.get_modulo_field())};

    for (const auto& factor : factors) {
      auto tmp{brute_force_composite(Automata(MyAutomata, factor))};

      // If size is 0, no solution was found
      if (tmp.size() == 0) {
        return false;
      }
#ifdef SHOW
      else if (tmp.size() == 1) {
        std::cout << "The min-poly of the automata in GF(" << factor
                  << ") is X^" << tmp[0] << std::endl;
      } else if (tmp.size() == 2) {
        std::cout << "The min-poly of the automata in GF(" << factor
                  << ") is X^" << tmp[0] << "(X-1)" << std::endl;
      }
#endif
    }
    return true;
  } else {
    return false;
  }
}

// Function to read a boolean value from a file
bool read_boolean_from_file(const std::string& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return false;  // Or handle the error as appropriate
  }

  char value;
  file >> value;
  file.close();

  // Interpret the character as a boolean value
  return (value == '1');
}

int main() {
  Automata MyAutomata;
  // MyAutomata.populate_from_user();
  MyAutomata.populate_from_csv();
  std::ofstream ofs("mat.csv", std::ofstream::trunc);
  ofs << MyAutomata;
  ofs.close();

  system("python3 main.py mat.csv");

  auto discrete{read_boolean_from_file("demofile2.txt")};

  std::ofstream result("result.txt", std::ofstream::trunc);
  if (discrete) {
    std::cout << "The discretization of the system is fixed point."
              << std::endl;
    result << "The discretization of the system is fixed point." << std::endl;
  } else {
    std::cout << "The discretization of the system is not fixed point."
              << std::endl;
    result << "The discretization of the system is not fixed point."
           << std::endl;
  }

  if (is_fixed_point(MyAutomata) && discrete) {
    std::cout << "The decomposition of the system is of fixed point."
              << std::endl;
    result << "The decomposition of the system is of fixed point." << std::endl;
  } else {
    std::cout << "The decomposition of the system is of not fixed point."
              << std::endl;
    result << "The decomposition of the system is of not fixed point."
           << std::endl;
  }

  result.close();

  return 0;
}