#ifndef CONSTITUTIVE_H
#define CONSTITUTIVE_H

#include "utils.h"

double f(double lambda);

double energy_f(double lambda);

double cohen_func(double x, double lambda_lim);

double relative_energy_cohen_func(double x, double lambda_lim);

#endif