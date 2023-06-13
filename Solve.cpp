#include "Solve.h"
#include <cmath>
#include <limits>
#include <iostream>

using namespace std;
vector<double> GetRootsWithBettaZero(double a, double b, double alpha, double gamma)
{
    vector<double> roots;
    roots.reserve(4);
    double rootInside = sqrt(alpha * alpha - 4 * gamma);
    double outsideRoot1 = sqrt((-alpha + rootInside) / 2);
    double outsideRoot2 = sqrt((-alpha - rootInside) / 2);
    double minusBDiv4A = - b / (4 * a);
    roots.emplace_back(minusBDiv4A + outsideRoot1);
    roots.emplace_back(minusBDiv4A - outsideRoot1);
    roots.emplace_back(minusBDiv4A + outsideRoot2);
    roots.emplace_back(minusBDiv4A - outsideRoot2);
    return roots;
}

EquationRoot4 Solve4_1(double a, double b, double c, double d, double e)
{
    // решение взято из статьи https://ru.wikipedia.org/wiki/Метод_Феррари
    if (a == 0)
    {
        throw invalid_argument("First argument must not be negative");
    }

    double alpha = -((3 * b * b) / (8 * a * a)) + (c / a);
    double betta = ((b * b * b) / (8 * a * a * a)) - ((b * c) / (2 * a * a)) + (d / a);
    double gamma = -((3 * b * b * b * b) / (256 * a * a * a * a)) + ((b * b * c) / (16 * a * a * a))
                    - ((b * d) / (4 * a * a)) + (e / a);

    if (abs(betta) < numeric_limits<double>::epsilon()) // betta == 0
    {
        vector<double> roots = GetRootsWithBettaZero(a, b, alpha, betta);
        return EquationRoot4{4, roots};
    }

    double p = - (alpha * alpha / 12) - gamma;
    double q = - (alpha * alpha * alpha / 108) + (alpha * gamma / 3) - (betta * betta / 8);
    double sqr = sqrt((q * q / 4.0) + ((p * p * p) / 27.0));
    /*double p3 = p * p * p;
    double sqr = sqrt(p3);*/
    double r = - (q / 2) + sqr;
    double u = pow(r , 1.0 / 3.0);
    double y = - (5 * alpha / 6) + u;
    if (abs(u) < numeric_limits<double>::epsilon())
    {
        double negCubeRootQ = - pow(q, 1.0 / 3.0);
        y -= negCubeRootQ;
    }
    else
    {
        y -= (- p / (3 * u));
    }
    double w = sqrt(alpha + 2 * y);
    double minusBDiv4A = - b / (4 * a);

    vector<double> res;
    res.reserve(4);
    res.emplace_back(minusBDiv4A + (w + sqrt(- (3 * alpha + 2 * y + (2 * betta / w)))) / 2);
    res.emplace_back(minusBDiv4A + (w - sqrt(- (3 * alpha + 2 * y + (2 * betta / w)))) / 2);
    res.emplace_back(minusBDiv4A + (-w + sqrt(- (3 * alpha + 2 * y - (2 * betta / w)))) / 2);
    res.emplace_back(minusBDiv4A + (-w - sqrt(- (3 * alpha + 2 * y - (2 * betta / w)))) / 2);
    return EquationRoot4{4, res};
}

struct EquationRoot2
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

struct EquationRoot3
{
    int numRoots;
    vector<double> roots;
};

int sign(double x)
{
    if (x < 0)
    {
        return -1;
    }
    return 1;
}

EquationRoot3 Solve3(double a, double b, double c)
{
    double q = (a * a - 3 * b) / 9;
    double r = ((2 * a * a * a) - (9 * a * b) + (27 * c)) / 54;
    double s = q * q * q - r * r;

    if (s > 0)
    {
        double phi = acos(r / sqrt(q * q * q)) / 3;
        vector<double> v;
        v.emplace_back(-2 * sqrt(q) * cos(phi) - a / 3);
        v.emplace_back(-2 * sqrt(q) * cos(phi + M_PI * 2 / 3) - a / 3);
        v.emplace_back(-2 * sqrt(q) * cos(phi - M_PI * 2 / 3) - a / 3);
        return EquationRoot3 {3, v};
    }
    if (s == 0)
    {
        vector<double> v;
        v.emplace_back(-2 * pow(r, 1.0/3.0) - a / 3);
        v.emplace_back(pow(r, 1.0/3.0) - a / 3);
        return EquationRoot3{2, v};
    }
    if (s < 0)
    {
        vector<double> v;
        if (q > 0)
        {
            double phi = acos(abs(r) / sqrt(q * q * q)) / 3;
            v.emplace_back(-2 * pow(r, 1.0/3.0) * cosh(phi) - a / 3);
            return EquationRoot3{1, v};
        }
        if (q < 0)
        {
            double sq3 = sqrt(abs(q * q * q));
            double asVal = abs(r) / sq3;
            double phi = asinh(asVal) / 3;
            double val = -2 * sign(r) * sqrt(abs(q)) * sinh(phi) - a / 3.0;
            v.push_back(val);
            return EquationRoot3{1, v};
        }
        v.emplace_back(-pow(c - a * a * a / 27, 1.0/3.0) - a / 3);
        return EquationRoot3{1, v};
    }
    return EquationRoot3{0, {}};
}

EquationRoot4 Solve4(double a, double b, double c, double d, double e)
{
    if (a == 0)
    {
        throw invalid_argument("First argument must not be negative");
    }
    double a3 = b / a;
    double a2 = c / a;
    double a1 = d / a;
    double a0 = e / a;

    double cubeQ2 = -a2;
    double cubeQ1 = a1 * a3 - 4 * a0;
    double cubeQ0 = -(a1 * a1 + a0 * a3 * a3 - 4 * a0 * a2);

    EquationRoot3 roots = Solve3(cubeQ2, cubeQ1, cubeQ0);
    double solution = roots.roots[0];

    double p1 = (a3 / 2) + sqrt(a3 * a3 / 4 + solution - a2);
    double p2 = (a3 / 2) - sqrt(a3 * a3 / 4 + solution - a2);
    double q1 = (solution / 2) + sqrt(solution * solution / 4 - a0);
    double q2 = (solution / 2) - sqrt(solution * solution / 4 - a0);

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