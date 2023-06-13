#ifndef SOLVE4_SOLVE_H
#define SOLVE4_SOLVE_H

#include <vector>

const double EPSILON = 1e-13;

struct EquationRoot4
{
    int numRoots;
    std::vector<double> roots;
};

EquationRoot4 Solve4(double a, double b, double c, double d, double e);

#endif //SOLVE4_SOLVE_H
