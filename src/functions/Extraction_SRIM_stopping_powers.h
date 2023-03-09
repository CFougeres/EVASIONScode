#ifndef EXTRACTION_SRIM_STOPPING_POWERS_H
#define EXTRACTION_SRIM_STOPPING_POWERS_H

#include "main.h"

#define MAX_LINES 5000
#define NUCLEI 5

int Extraction_SRIM_stopping_powers(Double_t IonEnergy[NUCLEI][MAX_LINES], Double_t dEdx_e[NUCLEI][MAX_LINES], Double_t dEdx_n[NUCLEI][MAX_LINES], int Nfile_energy[NUCLEI]);

#endif
