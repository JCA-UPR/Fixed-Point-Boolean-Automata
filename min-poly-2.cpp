#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "Automata.hpp"

using Eigen::MatrixXd;
using namespace Eigen;
using namespace std;

#define DEBUG 0

std::ostream& operator << (std::ostream& out, const std::vector<double>& vector) {
    for (const auto& num : vector) {
        out << num << ", ";
    }
    return out;
}

std::ostream& operator << (std::ostream& out, const Automata& automata) {
    out << automata.get_power_matrix();
    return out;
}

Automata& operator * (const Automata& automata_1, const Automata& automata_2) {

    // Check that the dimentions of the matrices match
    // if (automata_1.get_columns() != automata_2.get_rows()) return Automata&(nullptr);

    Automata* new_automata {new Automata(   automata_1.get_rows(), 
                                            automata_1.get_columns(), 
                                            automata_1.get_modulo_field())};

    new_automata->set_power_matrix(automata_1.get_power_matrix() * automata_2.get_power_matrix());

    // Apply modulo operation to elements if the automata has a modulo > 1
    if (new_automata->get_modulo_field() > 1) {
        for (auto i {0}; i < new_automata->get_rows(); i++) {
            for (auto j {0}; j < new_automata->get_columns(); j++) {
                (*new_automata)(i,j) = int((*new_automata)(i,j)) % int(new_automata->get_modulo_field());
            }
        }
    } else { // If it does not, return the multiplication of the matrices
        return *new_automata;
    }

    return *new_automata;
}

double& mod(const double& dividend, const double& divisor) {
    // if the divisor is 0
    if (!divisor) std::cerr << "Cannot divide by 0." << std::endl;
    // If divisor is larger than dividend
    if (divisor > dividend) std::cerr << "Divisor should be larger than the dividend." << std::endl;

    double* remainder {new double};

    // Remainder = Dividend – (Divisor x Quotient)
    *remainder = dividend - (divisor * floor(dividend/divisor));

    // while (*remainder - divisor >= divisor)  {
    //     #if DEBUG == 2
    //         std::cout << *remainder << " - " << divisor << " >= " << divisor << "?" << std::endl; 
    //     #endif
    //     *remainder -= divisor;
    // }
    // *remainder -=divisor;

    return *remainder;
}

bool is_prime(const double& num) {
    // Calculate √num
    double sq {sqrt(num)};

    #if DEBUG == 2
        std::cout << "Square root of " << num << " is : " << sq << std::endl;
    #endif

    // Check if num is divisible by primes up to √num or we pass first 10,000 primes
    for (int i {0}; (i < 10000) && (primes[i] <= sq); i++) {
        // std::cout << "Here at " << i << std::endl;

        #if DEBUG == 2
            std::cout << "Mod of " << num << " and " << primes[i] << " is " << mod(num,primes[i]) << std::endl;
        #endif

        if (mod(num, primes[i]) == 0) return false;
    }


    return true;
}

const std::vector<double>& decompose_field(double field) {

    std::vector<double>* factors {new std::vector<double>()};

    if (is_prime(field)) {
        factors->emplace_back(field);
        return *factors;
    } else {
        // double sq {sqrt(field)+1};
        for (auto i {0}; (i < 10000) && primes[i] < field; i++) {
            if (mod(field, primes[i]) == 0)  {
                factors->emplace_back(primes[i]);
            }
        }

        return *factors;
    }

}


bool create_field_prime_factor_automatas(const Automata& automata) {
    std::vector<double> prime_factors {decompose_field(automata.get_modulo_field())};
    int n {automata.get_columns()};

    // Create a vector that will hold all of the matrices
    // std::vector<std::vector<automata&>> components(n, std::vector<automata&>(n, automata));

    // If below log-algorthm optimal point
    // if (n <= 6) {

    // } 


    // If above log-algorithm optimal point

    return 1;
}


int main() {
    // Automata MyAutomata, Identity;
    // MyAutomata.populate_from_user();
    // Identity.populate_from_user();

    // std::cout << std::endl << (MyAutomata * Identity) << std::endl;

    int var {};

    std::cout << "Prime?" << std::endl;
    std::cin >> var;
    std::cout << decompose_field(var) << std::endl;

    return 0;
}