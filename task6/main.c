#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static constexpr double PI = 3.1415926535897931;

static double square(double x) {
    return x * x;
}

static void thomas(size_t n,
                   const double* restrict a,
                   const double* restrict b,
                   const double* restrict c,
                   const double* restrict d,
                   double* restrict x) {
    double* alpha = malloc(2 * n * sizeof(double));
    double* beta = alpha + n;

    alpha[0] = c[0] / b[0];
    beta[0] = d[0] / b[0];

    for (size_t i = 1; i < n; ++i) {
        alpha[i] = c[i] / (b[i] - a[i] * alpha[i - 1]);
        beta[i] = (d[i] - a[i] * beta[i - 1]) / (b[i] - a[i] * alpha[i - 1]);
    }

    x[n - 1] = beta[n - 1];
    for (size_t i = n - 1; i-- > 0;) {
        x[i] = beta[i] - alpha[i] * x[i + 1];
    }

    free(alpha);
}

// function left-right
static void solveLfRf(size_t n, double* u) {
    double h = PI / n;

    n -= 1;
    double* a = malloc(4 * n * sizeof(double));
    double* b = a + n;
    double* c = b + n;
    double* d = c + n;

    for (size_t i = 0; i < n; ++i) {
        double x = (i + 1) * h;
        a[i] = 1;
        b[i] = -2;
        c[i] = 1;
        d[i] = square(h) * sin(x);
    }

    thomas(n, a, b, c, d, u + 1);
    u[0] = 0;
    u[n + 1] = 0;

    free(a);
}

// function left, derivative right
static void solveLfRd(size_t n, double* u) {
    double h = PI / n;
    double* a = malloc(4 * n * sizeof(double));
    double* b = a + n;
    double* c = b + n;
    double* d = c + n;

    for (size_t i = 0; i < n - 1; ++i) {
        double x = (i + 1) * h;
        a[i] = 1;
        b[i] = -2;
        c[i] = 1;
        d[i] = square(h) * sin(x);
    }

    double xN = PI;
    a[n - 1] = 2;
    b[n - 1] = -2;
    c[n - 1] = 0;
    d[n - 1] = square(h) * sin(xN) + 2 * h * cos(xN);

    thomas(n, a, b, c, d, u + 1);
    u[0] = 0;

    free(a);
}

// derivative left, function right
static void solveLdRf(size_t n, double* u) {
    double h = PI / n;
    double* a = malloc(4 * n * sizeof(double));
    double* b = a + n;
    double* c = b + n;
    double* d = c + n;

    a[0] = 0;
    b[0] = -2;
    c[0] = 2;
    d[0] = square(h) * sin(0) - 2 * h * cos(0);

    for (size_t i = 1; i < n - 1; ++i) {
        double x = i * h;
        a[i] = 1;
        b[i] = -2;
        c[i] = 1;
        d[i] = h * h * sin(x);
    }

    a[n - 1] = 1;
    b[n - 1] = -2;
    c[n - 1] = 0;
    d[n - 1] = square(h) * sin((n - 1) * h);

    thomas(n, a, b, c, d, u);
    u[n] = 0;

    free(a);
}

static void dft(const double complex* restrict in,
                double complex* restrict out,
                size_t n) {
    for (size_t i = 0; i < n; ++i) {
        out[i] = 0;
        for (size_t j = 0; j < n; ++j) {
            double angle = -2 * PI * j * i / n;
            out[i] += in[j] * cexp(I * angle);
        }
    }
}

static void idft(const double complex* restrict in,
                 double complex* restrict out,
                 size_t n) {
    for (size_t i = 0; i < n; ++i) {
        out[i] = 0;
        for (size_t j = 0; j < n; ++j) {
            double angle = 2 * PI * i * j / n;
            out[i] += in[j] * cexp(I * angle);
        }
        out[i] /= n;
    }
}

static void solveDft(size_t n, double* u) {
    double complex* f = malloc(4 * n * sizeof(double complex));
    double complex* F = f + n;
    double complex* Y = F + n;
    double complex* yComplex = Y + n;

    double h = 2 * PI / n;

    for (size_t i = 0; i < n; ++i) {
        double x = i * h;
        f[i] = sin(x);
    }

    dft(f, F, n);

    for (size_t i = 0; i < n; ++i) {
        size_t m = (i <= n / 2) ? i : (i - n);

        if (m == 0) {
            Y[i] = 0;
        } else {
            Y[i] = -F[i] / (m * m);
        }
    }

    idft(Y, yComplex, n);

    for (size_t i = 0; i < n; ++i) {
        u[i] = creal(yComplex[i]);
    }

    free(f);
}

static double exact(double x) {
    return -sin(x);
}

static double maxError(size_t n, const double* y, double a, double b) {
    double h = (b - a) / n;
    double err = 0;
    for (size_t i = 0; i <= n; ++i) {
        double x = a + i * h;
        double diff = fabs(y[i] - exact(x));
        err = fmax(diff, err);
    }
    return err;
}

static double maxErrorPeriodic(size_t n, const double* y, double a, double b) {
    double err = 0;
    double h = (b - a) / n;
    for (size_t i = 0; i < n; ++i) {
        double x = a + i * h;
        double diff = fabs(y[i] - exact(x));
        err = fmax(err, diff);
    }
    return err;
}

int main() {
    {
        printf("thomas:\n");

        constexpr size_t N_BOUNDARY_CONDITIONS = 3;

        const char* names[N_BOUNDARY_CONDITIONS] = {
            "function left-right",
            "function left, derivative right",
            "derivative left, function right"};

        void (*solvers[N_BOUNDARY_CONDITIONS])(size_t, double*) = {solveLfRf,
                                                                   solveLfRd,
                                                                   solveLdRf};

        constexpr size_t ORDERS = 5;
        constexpr size_t INITIAL_SIZE = 16;
        double* orders = malloc((ORDERS - 1) * sizeof(double));

        for (size_t i = 0; i < N_BOUNDARY_CONDITIONS; ++i) {
            double prevErr;

            for (size_t j = 0; j < ORDERS; ++j) {
                size_t n = INITIAL_SIZE << j;
                double* y = malloc((n + 1) * sizeof(double));

                solvers[i](n, y);
                double err = maxError(n, y, 0, PI);

                if (j > 0) {
                    orders[j - 1] = log2(prevErr / err);
                }

                prevErr = err;
                free(y);
            }

            printf("%s\n", names[i]);
            printf("N\tOrder\n");
            for (size_t j = 0; j < ORDERS - 1; ++j) {
                size_t n = INITIAL_SIZE << (j + 1);
                printf(" %zu\t %lf\n", n, orders[j]);
            }
            putchar('\n');
        }
        free(orders);
    }

    {
        printf("periodic + DFT:\n");

        constexpr size_t ORDERS = 5;
        constexpr size_t INITIAL_SIZE = 2;

        printf("N\tMax error\n");
        for (size_t j = 0; j < ORDERS; ++j) {
            size_t n = INITIAL_SIZE << j;
            double* y = malloc((n + 1) * sizeof(double));

            solveDft(n, y);
            double err = maxErrorPeriodic(n, y, 0, 2 * PI);
            printf(" %zu\t %e\n", n, err);

            free(y);
        }
    }

    return 0;
}
