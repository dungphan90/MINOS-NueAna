#ifndef TYWEIGHT_C
#define TYWEIGHT_C

#include "NueAna/NueRecord.h"

double TYHadronWeight(double w2, double totpt, int npi0, int res, double *par);
double TYHadronWeight(NueRecord *nr, double *par);
double TYHadronWeight(NueRecord *nr);
double TYHadronWeight(double w2, double totpt, int npi0, int res);

#endif
