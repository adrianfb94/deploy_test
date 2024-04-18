#ifndef PERCENTILES_H
#define PERCENTILES_H

double percentile(double *array, size_t len, double pn);
void percentile_set_method(const char *methodstr);
void percentile_check_number(double pn);

#endif /* PERCENTILES_H */
