#include "Solve.h"
#include <cmath>
#include <limits>
#include <iostream>

using namespace std;

struct EquationRoot2
{
    int numRoots;
    vector<double> roots;
};

struct EquationRoot3
{
    int numRoots;
    vector<double> roots;
};

EquationRoot2 Solve2(double b, double c)
{
    double discriminant = b * b - 4 * c;
    if (discriminant < 0)
    {
        return EquationRoot2{0, {}};
    }
    if (discriminant == 0)
    {
        return EquationRoot2{2, {-b / (2), -b / (2)}};
    }
    return EquationRoot2{
        2,
        {
            (-b + sqrt(discriminant)) / (2),
            (-b - sqrt(discriminant)) / (2)
        }
    };
}

int Sign(double x)
{
    if (x < 0)
    {
        return -1;
    }
    return 1;
}

struct CoefficientsVieta
{
    double q;
    double r;
    double s;
};

CoefficientsVieta GetVeitaCoefs(double a, double b, double c)
{
    double q = (a * a - 3 * b) / 9;
    double r = ((2 * a * a * a) - (9 * a * b) + (27 * c)) / 54;
    double s = q * q * q - r * r;

    return CoefficientsVieta{q, r, s};
}

EquationRoot3 Solve3(double a, double b, double c)
{
    CoefficientsVieta vieta = GetVeitaCoefs(a, b, c);
    vector<double> v;
    if (vieta.s > 0)
    {
        double phi = acos(vieta.r / sqrt(vieta.q * vieta.q * vieta.q)) / 3;
        v.emplace_back(-2 * sqrt(vieta.q) * cos(phi) - a / 3);
        v.emplace_back(-2 * sqrt(vieta.q) * cos(phi + M_PI * 2 / 3) - a / 3);
        v.emplace_back(-2 * sqrt(vieta.q) * cos(phi - M_PI * 2 / 3) - a / 3);
    }
    if (vieta.s == 0)
    {
        v.emplace_back(-2 * pow(vieta.r, 1.0 / 3.0) - a / 3);
        v.emplace_back(pow(vieta.r, 1.0 / 3.0) - a / 3);
    }
    if (vieta.s < 0)
    {
        if (vieta.q > 0)
        {
            double sq3 = sqrt(abs(vieta.q * vieta.q * vieta.q));
            double phi = (acosh(abs(vieta.r) / sq3)) / 3;
            v.emplace_back(-2 * Sign(vieta.r) * sqrt(vieta.q) * cosh(phi) - a / 3);
        }
        if (vieta.q < 0)
        {
            double sq3 = sqrt(abs(vieta.q * vieta.q * vieta.q));
            double phi = asinh(abs(vieta.r) / sq3) / 3;
            v.emplace_back(-2 * Sign(vieta.r) * sqrt(abs(vieta.q)) * sinh(phi) - a / 3.0);
        }
        if (vieta.q == 0)
        {
            v.emplace_back(-pow(c - a * a * a / 27, 1.0/3.0) - a / 3);
        }
    }
    return EquationRoot3{static_cast<int>(v.size()), v};
}

struct CubeResolventeSolutions
{
    double p1;
    double p2;
    double q1;
    double q2;
};

EquationRoot4 GetRootsOfSquareEq(double p1, double p2, double q1, double q2, double a1)
{
    double tmp = p1 * q2 + p2 * q1;
    double delta = abs(tmp - a1);
    if (isnan(tmp))
    {
        throw domain_error("No real roots");
    }
    if (delta > EPSILON)
    {
        swap(q1, q2);
    }

    vector<double> v;
    EquationRoot2 eq = Solve2(p1, q1);
    if (eq.numRoots > 0)
    {
        v.push_back(eq.roots[0]);
        v.push_back(eq.roots[1]);
    }
    eq = Solve2(p2, q2);
    if (eq.numRoots > 0)
    {
        v.push_back(eq.roots[0]);
        v.push_back(eq.roots[1]);
    }
    return EquationRoot4{static_cast<int>(v.size()), v};
}

EquationRoot4 Solve4(double a, double b, double c, double d, double e)
{
    if (a == 0)
    {
        throw invalid_argument("First argument must not be zero");
    }
    double a3 = b / a, a2 = c / a, a1 = d / a, a0 = e / a;

    double cubeQ2 = -a2;
    double cubeQ1 = a1 * a3 - 4 * a0;
    double cubeQ0 = -(a1 * a1 + a0 * a3 * a3 - 4 * a0 * a2);

    EquationRoot3 roots = Solve3(cubeQ2, cubeQ1, cubeQ0);
    double solution = roots.roots[0];

    double p1 = (a3 / 2) + sqrt(a3 * a3 / 4 + solution - a2);
    double p2 = (a3 / 2) - sqrt(a3 * a3 / 4 + solution - a2);
    double q1 = (solution / 2) + sqrt(solution * solution / 4 - a0);
    double q2 = (solution / 2) - sqrt(solution * solution / 4 - a0);

    // x^4+a_3x^3+a_2x^2+a_1x+a_0=(x^2+p_1x+q_1)(x^2+p_2x+q_2)
    return GetRootsOfSquareEq(p1, p2, q1, q2, a1);
}