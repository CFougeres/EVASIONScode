#ifndef TOOLS_H
#define TOOLS_H

#include "main.h"


#define MAX_LINES 5000

int integral_cumulative(TH1* h, int noiseRange[2]);

double integral_counts(TH2F* h2, int Range[2]);

double ponderation(TH1F* values);

double chi2_ndf_1D(TH1D* Sobs, TH1D* Sfit, int number_fit_param);

double log_likelihood(TH1F* signalObservation, TH1F* signalFit , double fit_range[2]);

int extraction_noise(TH2F* h2, Double_t angle[], Double_t ponderation_angle[]);

double conversion_angle_random_distribution(double r, Double_t angle[], Double_t ponderation_angle[]);

double interpol(double x[MAX_LINES],double y[MAX_LINES],int const n,double p);

double loss_E(double Ei, double distance, Double_t EnergyPart[MAX_LINES], Double_t SPel_Part[MAX_LINES], Double_t SPnu_Part[MAX_LINES], int N_energy);

Double_t pol2_func(Double_t* x, Double_t* par) ;

double Brho_Ekin_conv(double input, double Mnuc, double Anuc, double Znuc, int conv);
//double errorIntegral(TH1F* h, int bin1, int bin2);
//Int_t binConversion1D(TH1* signal, double value);

#endif
