#include "integration.h"

#include <math.h>
#include <stdlib.h>

double integrateSimpson(Function f, double a, double b, unsigned resolution) {
    resolution += resolution % 2;  // round up to even
    double h = (b - a) / resolution;

#define X(I) (a + (I) * h)

    double sum = 0;
    for (unsigned i = 1; i <= resolution; i += 2) {
        sum += f(X(i - 1)) + 4.0 * f(X(i)) + f(X(i + 1));
    }

#undef X

    return sum * h / 3.0;
}

static constexpr double PI = 3.1415926535897931;
static constexpr unsigned NEWTON_STEPS = 6;

static double exchange(double* x, double y) {
    double tmp = *x;
    *x = y;
    return tmp;
}

static double square(double x) {
    return x * x;
}

static double legendrePolynomial(unsigned n, double x) {
    if (n == 0)
        return 1;

    // https://en.wikipedia.org/wiki/Legendre_polynomials#Recurrence_relations
    double p1 = 1, p2 = x;
    for (unsigned i = 1; i < n; ++i) {
        p1 = exchange(&p2, ((2 * i + 1) * x * p2 - i * p1) / (i + 1));
    }
    return p2;
}

static double legendrePolynomialDerivativeAtRoot(unsigned n, double x) {
    if (n == 0)
        return 0;

    // https://en.wikipedia.org/wiki/Legendre_polynomials#Recurrence_relations
    // simplified considering P_n(x) = 0 (x is a root)
    return n * legendrePolynomial(n - 1, x) / (1 - square(x));
}

static double* legendrePolynomialRoots(unsigned n) {
    double* res = malloc(sizeof(double) * n);

    for (unsigned i = 0; i < n / 2 + n % 2; ++i) {
        // initial guess derived from Bruns' theorem:
        // https://math.stackexchange.com/revisions/1280800/2
        double root = cos(PI * (4 * (i + 1) - 1) / (4 * n + 2));

        // a few steps of Newton's algorithm
        for (unsigned step = 0; step < NEWTON_STEPS; ++step) {
            root = root - legendrePolynomial(n, root) /
                              legendrePolynomialDerivativeAtRoot(n, root);
        }

        // roots are symmetrical
        res[i] = -root;
        res[n - i - 1] = root;
    }

    return res;
}

double integrateGaussLegendre(Function f,
                              double a,
                              double b,
                              unsigned resolution) {
    double res = 0.0;

    double* roots = legendrePolynomialRoots(resolution);

    for (unsigned i = 0; i < resolution; ++i) {
        double t = roots[i];

        double dpn = legendrePolynomialDerivativeAtRoot(resolution, t);
        double weight = 2.0 / ((1 - square(t)) * square(dpn));
        double x = (b - a) / 2 * t + (a + b) / 2;

        res += weight * f(x);
    }

    free(roots);

    return (b - a) / 2 * res;
}
