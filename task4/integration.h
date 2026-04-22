#pragma once

typedef double (*Function)(double);

double integrateSimpson(Function f, double a, double b, unsigned resolution);
double integrateGaussLegendre(Function f, double a, double b, unsigned resolution);
