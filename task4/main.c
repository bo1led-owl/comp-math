#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "integration.h"

static constexpr double PI = 3.1415926535897931;

static void integrateSomeFunctions();
static void testConvergence();

static double square(double x) {
    return x * x;
}

static double f1(double x) {
    double x5 = x * x * x * x * x;
    return (sin(PI * x5)) / (x5 * (1 - x));
}

static double f2(double t) {
    double x = t / (1 - t);
    return exp(-sqrt(x) + sin(x / 10)) / square(1 - t);
}

int main() {
    // integrateSomeFunctions();
    // putchar('\n');

    testConvergence();
    putchar('\n');

    puts("computed integrals:");

    // wolfram gave 8.03491060...
    printf("I1 = %.8lf\n", integrateGaussLegendre(f1, 0, 1, 16));
    printf("I2 = %.8lf\n", integrateGaussLegendre(f2, 0, 1, 32));

    for (unsigned i = 8; i <= 256; i *= 2) {
    printf("%.8lf\n", fabs(integrateGaussLegendre(f2, 0, 1, i) - 2.981));
    }
}

static void testConvergence() {
    static constexpr double A = 0;
    static constexpr double B = 1;
    double exact = exp(1) - 1;

    static constexpr unsigned N = 6;

    double* simpsonErrors = calloc(N, sizeof(double));
    double* gaussLegendreErrors = calloc(N, sizeof(double));

    for (unsigned i = 0; i < N; ++i) {
        unsigned resolution = 2 << i;

        simpsonErrors[i] =
            fabs(exact - integrateSimpson(exp, A, B, resolution));
        gaussLegendreErrors[i] =
            fabs(exact - integrateGaussLegendre(exp, A, B, resolution));
    }

    printf("Simpson order of convergence:\n");
    for (unsigned i = 1; i < N; ++i) {
        double order = log2(simpsonErrors[i - 1] / simpsonErrors[i]);
        printf("%.16lf %.16lf\n", simpsonErrors[i - 1], simpsonErrors[i]);
        printf("%lf\n", order);
    }

    printf("\nGauss-Legendre order of convergence:\n");
    for (unsigned i = 1; i < N; ++i) {
        double order =
            log2(gaussLegendreErrors[i - 1] / gaussLegendreErrors[i]);
        printf("%.16lf %.16lf\n", gaussLegendreErrors[i - 1], gaussLegendreErrors[i]);
        printf("%lf\n", order);
    }

    free(simpsonErrors);
    free(gaussLegendreErrors);
}

static double f(double x) {
    return 2 * x;
}

[[maybe_unused]] static void integrateSomeFunctions() {
    typedef struct {
        char* name;
        Function f;
        double a, b;
    } Integral;

    static constexpr Integral integrals[] = {
        (Integral){.name = "2x", .f = f, .a = 0, .b = 1},
        (Integral){.name = "sin x", .f = sin, .a = 0, .b = PI},
        (Integral){.name = "cos x", .f = cos, .a = 0, .b = PI},
        (Integral){.name = "sin x", .f = sin, .a = 0, .b = PI / 2},
    };

    static constexpr unsigned GAUSS_RESOLUTION = 8;
    static constexpr unsigned SIMPSON_RESOLUTION = 64;

    static constexpr size_t n = sizeof(integrals) / sizeof(Integral);
    for (size_t i = 0; i < n; ++i) {
        Integral integral = integrals[i];

        printf("int_%.2lf^%.2lf %s dx\n",
               integral.a,
               integral.b,
               integral.name);
        printf("gauss-legendre:\t%lf\nsimpson:\t%lf\n",
               integrateGaussLegendre(integral.f,
                                      integral.a,
                                      integral.b,
                                      GAUSS_RESOLUTION),
               integrateSimpson(integral.f,
                                integral.a,
                                integral.b,
                                SIMPSON_RESOLUTION));

        if (i + 1 < n) {
            putchar('\n');
        }
    }
}
