#include <iostream>
#include <cmath>
#include <complex>

#include "Solve.h"

using namespace std;
// 1 -5 6 5 -6 два действительных корня

void HandleInput(double& a, double& b, double& c, double& d, double& e, string const& inputLine)
{
    stringstream ss(inputLine);
    ss >> a >> b >> c >> d >> e;
    if (ss.fail())
    {
        throw runtime_error("Too few args");
    }
    if (!ss.eof())
    {
        throw runtime_error("Too many args");
    }
}

int main()
{
    double a, b, c, d, e;
    string inputLine, startMessage = "Enter coefficients of the equation ax^4 + bx^3 + cx^2 + dx + e = 0:";
    cout << startMessage << endl;
    while(getline(cin, inputLine))
    {
        try
        {
            HandleInput(a, b, c, d, e, inputLine);
        }
        catch (runtime_error const& ex)
        {
            cout << ex.what() << endl << startMessage << endl;
            continue;
        }
        try
        {
            EquationRoot4 res = Solve4(a, b, c, d, e);
            cout << "Solutions count: " << res.numRoots << endl;
            copy(res.roots.begin(), res.roots.end(), ostream_iterator<double>(cout, " "));
            cout << endl;
        }
        catch (invalid_argument const& ex)
        {
            cout << ex.what() << endl;
        }
        catch (domain_error const& ex)
        {
            cout << ex.what() << endl;
        }
        cout << startMessage << endl;
    }
    return EXIT_SUCCESS;
}