#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/catch_approx.hpp>
#include "../Solve.h"

TEST_CASE("First arg 0")
{
    double a = 0, b = 1, c = 1, d = 1, e = 1;
    REQUIRE_THROWS_AS(Solve4(a, b, c, d, e), std::invalid_argument);
    REQUIRE_THROWS_WITH(Solve4(a, b, c, d, e), "First argument must not be zero");
}

TEST_CASE("Real roots")
{
    double a = 1, b = 4, c = -4, d = -20, e = -5;
    EquationRoot4 eq = Solve4(a, b, c, d, e);
    REQUIRE(eq.numRoots == 4);
    REQUIRE(eq.roots[0] == Catch::Approx(-0.27).margin(0.01));
    REQUIRE(eq.roots[1] == Catch::Approx(-2.24).margin(0.01));
    REQUIRE(eq.roots[2] == Catch::Approx(2.24).margin(0.01));
    REQUIRE(eq.roots[3] == Catch::Approx(-3.73).margin(0.01));
}

TEST_CASE("2 real roots")
{
    double a = 1, b = 4, c = -4, d = 20, e = -5;
    EquationRoot4 eq = Solve4(a, b, c, d, e);
    REQUIRE(eq.numRoots == 2);
    REQUIRE(eq.roots[0] == Catch::Approx(0.26).margin(0.01));
    REQUIRE(eq.roots[1] == Catch::Approx(-5.44).margin(0.01));
}

TEST_CASE("0 real roots")
{
    double a = 1, b = 0, c = 0, d = 0, e = 5;
    REQUIRE_THROWS_AS(Solve4(a, b, c, d, e), std::domain_error);
    REQUIRE_THROWS_WITH(Solve4(a, b, c, d, e), "No real roots");
}

TEST_CASE("2 real roots 1 -5 6 5 -6")
{
    double a = 1, b = -5, c = 6, d = 5, e = -6;
    EquationRoot4 eq = Solve4(a, b, c, d, e);
    REQUIRE(eq.numRoots == 2);
    REQUIRE(eq.roots[0] == Catch::Approx(0.84).margin(0.01));
    REQUIRE(eq.roots[1] == Catch::Approx(-0.96).margin(0.01));
}