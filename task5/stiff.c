#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static constexpr double PERT = 1e-8;

static void vadd(size_t n, double* dest, const double* a, const double* b) {
    for (size_t i = 0; i < n; ++i) {
        dest[i] = a[i] + b[i];
    }
}

static void vsub(size_t n, double* dest, const double* a, const double* b) {
    for (size_t i = 0; i < n; ++i) {
        dest[i] = a[i] - b[i];
    }
}

static void vscale(size_t n, double* dest, double x, const double* y) {
    for (size_t i = 0; i < n; ++i) {
        dest[i] = x * y[i];
    }
}

static void vneg(size_t n, double* dest, const double* x) {
    for (size_t i = 0; i < n; ++i) {
        dest[i] = -x[i];
    }
}

static void eulerExplicit(size_t n,
                          double* dest,
                          void (*f)(double*, double, const double*),
                          const double* y0,
                          size_t gridSize,
                          const double* grid) {
    double* yn = dest;
    memcpy(yn, y0, sizeof(double) * n);

    double tn = grid[0];

    for (size_t i = 1; i < gridSize; ++i) {
        double h = grid[i] - tn;

        f(yn + n, tn, yn);
        vscale(n, yn + n, h, yn + n);
        vadd(n, yn + n, yn, yn + n);

        tn = grid[i];
        yn += n;
    }
}

static double magnitude(size_t n, const double* v) {
    double res = 0;
    for (size_t i = 0; i < n; ++i) {
        double a = fabs(v[i]);
        res = fmax(a, res);
    }
    return res;
}

static void swap(double* a, double* b) {
    double tmp = *a;
    *a = *b;
    *b = tmp;
}

static void gauss(size_t n, double* A, double* b) {
    for (size_t i = 0; i < n; ++i) {
        size_t maxRow = i;
        double maxVal = fabs(A[i * n + i]);
        for (size_t j = i + 1; j < n; ++j) {
            double val = fabs(A[j * n + i]);
            if (val > maxVal) {
                maxVal = val;
                maxRow = j;
            }
        }

        if (maxRow != i) {
            for (size_t j = i; j < n; ++j) {
                swap(A + i * n + j, A + maxRow * n + j);
            }

            swap(b + i, b + maxRow);
        }

        for (size_t j = i + 1; j < n; ++j) {
            double c = A[j * n + i] / A[i * n + i];
            A[j * n + i] = 0;
            for (size_t k = j + 1; k < n; ++k) {
                A[j * n + k] -= c * A[i * n + k];
            }
            b[j] -= c * b[i];
        }
    }

    for (size_t i = n; i-- > 0;) {
        double sum = b[i];
        for (size_t j = i + 1; j < n; ++j) {
            sum -= A[i * n + j] * b[j];
        }
        b[i] = sum / A[i * n + i];
    }
}

void eulerImplicit(size_t n,
                   size_t steps,
                   double delta,
                   double* dest,
                   void (*f)(double*, double, const double*),
                   const double* y0,
                   size_t gridSize,
                   const double* grid) {
    double* yn = dest;
    memcpy(yn, y0, sizeof(double) * n);

    double tn = grid[0];

    double* mem = malloc((6 * n + (n * n)) * sizeof(double));

    double* yGuess = mem;
    double* f0 = yGuess + n;
    double* fPert = f0 + n;
    double* fVal = fPert + n;
    double* yPert = fVal + n;
    double* rhs = yPert + n;
    double* J = rhs + n;

    for (size_t step = 1; step < gridSize; ++step) {
        double t_next = grid[step];
        double h = t_next - tn;

        memcpy(yGuess, yn, n * sizeof(double));

        for (size_t step = 0; step < steps; ++step) {
            f(f0, t_next, yGuess);

            vscale(n, fVal, h, f0);
            vsub(n, fVal, yGuess, fVal);
            vsub(n, fVal, fVal, yn);

            if (magnitude(n, fVal) < delta) {
                break;
            }

            // J = E - h * (partial f)/(partial y)

            memset(J, 0, n * n * sizeof(double));
            for (size_t i = 0; i < n; ++i) {
                J[i * n + i] = 1;
            }

            for (size_t i = 0; i < n; ++i) {
                memcpy(yPert, yGuess, n * sizeof(double));

                yPert[i] += PERT;
                f(fPert, t_next, yPert);

                for (size_t j = 0; j < n; ++j) {
                    J[j * n + i] -= h * (fPert[j] - f0[j]) / PERT;
                }
            }

            vneg(n, rhs, fVal);

            gauss(n, J, rhs);

            for (size_t i = 0; i < n; ++i) {
                yGuess[i] += rhs[i];
            }
        }

        memcpy(yn + n, yGuess, n * sizeof(double));

        tn = t_next;
        yn += n;
    }

    free(mem);
}

static void f(double* res, double, const double* y) {
    double u = y[0], v = y[1];

    res[0] = 998.0 * u + 1998.0 * v;
    res[1] = -999.0 * u - 1999.0 * v;
}

static void exact(double* res, double t) {
    res[0] = 2 * exp(-t) - exp(-1000 * t);
    res[1] = exp(-1000 * t) - exp(-t);
}

static double maxError(size_t n,
                       size_t gridSize,
                       const double* res,
                       const double* grid) {
    double* tmp = malloc(n * sizeof(double));

    double err = 0;
    for (size_t i = 0; i < gridSize; ++i) {
        double t = grid[i];
        const double* y = res + i * n;

        exact(tmp, t);

        vsub(n, tmp, tmp, y);
        err = fmax(err, magnitude(n, tmp));
    }

    free(tmp);
    return err;
}

static constexpr size_t GRID_SIZE = 10;
static constexpr double A = 0;
static constexpr double B = 0.1;

int main() {
    double* grid = malloc(GRID_SIZE * sizeof(double));
    for (size_t i = 0; i < GRID_SIZE; ++i) {
        grid[i] = A + i * ((B - A) / GRID_SIZE);
    }

    double* dest = malloc(GRID_SIZE * 2 * sizeof(double));
    double y0[2] = {1, 0};

    eulerExplicit(2, dest, f, y0, GRID_SIZE, grid);

    printf("euler's explicit max error magnitude:\t%e\n",
           maxError(2, GRID_SIZE, dest, grid));

    eulerImplicit(2, 10, 1e-8, dest, f, y0, GRID_SIZE, grid);
    printf("euler's implicit max error magnitude:\t%e\n",
           maxError(2, GRID_SIZE, dest, grid));

    free(grid);
    free(dest);
}
