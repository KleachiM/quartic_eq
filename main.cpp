#include <iostream>
#include <cmath>
#include <complex>

#include "Solve.h"

using namespace std;

//std::invalid_argument a == 0
//std::domain_error - нет действительных корней

int main()
{
    double a, b, c, d, e;
    string inputLine, msg = "Enter coefficients of the equation ax^4 + bx^3 + cx^2 + dx + e = 0:";
    cout << msg << endl;
    while(getline(cin, inputLine))
    {
        stringstream ss(inputLine);
        ss >> a >> b >> c >> d >> e;
        if (ss.fail())
        {
            cout << "Too few args" << endl << msg << endl;
            continue;
        }
        if (!ss.eof())
        {
            cout << "Too many args" << endl << msg << endl;
            continue;
        }
        try
        {
            EquationRoot4 res = Solve4(a, b, c, d, e);
            cout << "res: " << res.numRoots << endl;
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
        cout << msg << endl;
    }
    return EXIT_SUCCESS;
}