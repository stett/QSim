#include "gsl/gsl_constants.h"

double gsl_fdiv (const double x, const double y) {
  return x / y;
}
double gsl_nan (void) {
  return gsl_fdiv (0.0, 0.0);
}
double gsl_posinf (void) {
  return gsl_fdiv (+1.0, 0.0);
}
double gsl_neginf (void) {
  return gsl_fdiv (-1.0, 0.0);
}