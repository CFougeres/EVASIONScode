#include "Tools.h"
/*
 #include "Riostream.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <TRandom1.h>
 */

#define MAX_LINES 5000
//using namespace std;

int integral_cumulative(TH1* h, int noiseRange[2]) {
    TAxis* axis = h->GetXaxis();
    int bin_down = axis->FindBin(noiseRange[0]);
    int bin_up = axis->FindBin(noiseRange[1]);
    int S = 0;    int s;
    for (int k = bin_down; k < bin_up + 1; k++) {
        s = h->GetBinContent(k);
        S = S + s;
    }
    return S;
}

double integral_counts(TH2F* h2, int Range[2]) {
    TAxis* Xaxis = h2->GetXaxis();
    TAxis* Yaxis = h2->GetYaxis();
    int bminX = Xaxis->FindBin(Range[0]);
    int bmaxX = Xaxis->FindBin(Range[1]);
    double integ = h2->Integral(bminX, bmaxX, 0,-1);
    return integ;
}

double ponderation(TH1F* values) {
    int nX = values->GetNbinsX();
    double pond = 0;
    double I = values->Integral();
    double x, b;
    for (int i = 0; i < nX; i++) {
        b = values->GetBinContent(i);
        x = values->GetBinCenter(i);
        pond = pond + x * (b / I);
    }
    return pond;
}

double chi2_ndf_1D(TH1D* Sobs, TH1D* Sfit, int number_fit_param) {
    double chi = 0.0;
    double o, s, ero, ers;
    int nbins = 0;
    int Nbins = Sobs->GetNbinsX();
    for (int i = 1; i < Nbins; i++) {
        o = Sobs->GetBinContent(i);        s = Sfit->GetBinContent(i);
        ero = Sobs->GetBinError(i);      ers = Sfit->GetBinError(i);
        if (o > 0) {
            chi = chi + (o - s) * (o - s) / (ero*ero + ers*ers);
            nbins = nbins + 1;
        }
    }
    chi = chi / (nbins -1. - number_fit_param);
    return chi;
}

double log_likelihood(TH1F* signalObservation, TH1F* signalFit , double fit_range[2]) {
    double ll = 0.0;
    double o, s;
    TAxis* Xaxis = signalObservation->GetXaxis();
    int fit_bin[2] = { Xaxis->FindBin(fit_range[0]),Xaxis->FindBin(fit_range[1]) };
    int nbins = signalObservation->GetNbinsX();
    for (int i = fit_bin[0]; i < fit_bin[1]; i++) {
        o = signalObservation->GetBinContent(i);
        s = signalFit->GetBinContent(i);
        if (s > 0 && o > 0) {
            ll = ll + (s - o + o * log(o / s));
        }
        /*    if (o == 0) {
                ll = ll + (s - o);
            }*/
    }
    ll = 2 * ll;
    return ll;
}

int extraction_noise(TH2F* h2, Double_t angle[], Double_t ponderation_angle[]){
    TAxis* Xaxis = h2->GetXaxis();
    TAxis* Yaxis = h2->GetYaxis();
    int NbinsY = h2->GetNbinsY();
    TH1D* projX;
    int sizeArray = sizeof(angle)/sizeof(angle[0]);
    double noiseTotal, noiseCount, noiseFactor;
    double noiseSumUp = 0;
    noiseTotal = h2->Integral();
    projX = h2->ProjectionX("",0,2);
    noiseCount = projX->Integral();
    noiseFactor = noiseCount / noiseTotal;
    noiseSumUp = noiseFactor + noiseSumUp;
    angle[0] = Yaxis->GetBinCenter(0);
    ponderation_angle[0] = noiseSumUp;
    for (int j = 1; j < NbinsY; j++) {
        projX = h2->ProjectionX("", j - 1, j + 1);
        noiseCount  = projX->Integral();
        noiseFactor = noiseCount / noiseTotal;
        noiseSumUp = noiseFactor + noiseSumUp;
        angle[j] = Yaxis->GetBinCenter(j);
        ponderation_angle[j] = noiseSumUp;
    }
    for (int i = 0; i < sizeArray; i++) {
        ponderation_angle[i] = ponderation_angle[i] / ponderation_angle[sizeArray-1];
    }
    return 1;
}

double conversion_angle_random_distribution(double r, Double_t angle[], Double_t ponderation_angle[]) {
    if (r < ponderation_angle[0]) { return angle[0]; }
    int sizeArray = sizeof(angle)/sizeof(angle[0]);
    for (int i = 1; i <sizeArray; i++) {
        if (r > ponderation_angle[i - 1] && r < ponderation_angle[i]) { return angle[i]; break; }
    }
}

double interpol(double x[MAX_LINES],double y[MAX_LINES],int const n,double p){
  int i=1;
  double val,b1,c1,d1,c2,d2;
// b1 false if something is wrong
  bool bb1=true;
// loop breaking
  bool bb2=true;

  if (p<=x[1]){      bb1=false;}
  else{
      do{
          i++;
          if (p<=x[i])bb2=false;
      }while (bb2&&(i<n-2));
    }
  if (p>x[n-2])bb1=false;
  
  if (bb1==true){
      b1=(y[i-1]-y[i-2])/(x[i-1]-x[i-2]);
      b1=b1*p-b1*x[i-2]+y[i-2];
      
      c1=(y[i]-y[i-2])/(x[i]-x[i-2]);
      c1=c1*p-c1*x[i-2]+y[i-2];
      
      d1=(y[i+1]-y[i-2])/(x[i+1]-x[i-2]);
      d1=d1*p-d1*x[i-2]+y[i-2];
      
      c2=(c1-b1)/(x[i]-x[i-1]);
      c2=c2*p-c2*x[i-1]+b1;

      d2=(d1-b1)/(x[i+1]-x[i-1]);
      d2=d2*p-d2*x[i-1]+b1;

      val=(d2-c2)/(x[i]-x[i-1]);
      val=val*p-val*x[i]+c2;
    }
  else{ val=0.;}
  return val;
}

double loss_E(double Ei, double distance, Double_t EnergyPart[MAX_LINES], Double_t SPel_Part[MAX_LINES], Double_t SPnu_Part[MAX_LINES], int N_energy){

    double dx = 1 * 0.0001;
    double Etemp = Ei;
    double dpart = 0;
    while (dpart < distance) {
        Etemp = Etemp - 0.001 * interpol(EnergyPart, SPel_Part, N_energy, Etemp * 1000.0) * dx - 0.001 * interpol(EnergyPart, SPnu_Part, N_energy, Etemp * 1000.0) * dx;
        dpart = dpart + dx;
    }
    return Ei - Etemp;
}
Double_t pol2_func(Double_t* x, Double_t* par) {
    return  par[0] + par[1] * x[0] + par[2] * x[0] * x[0];
}
double Brho_Ekin_conv(double input, double Mnuc, double Anuc, double Znuc, int conv){
    double resultat;
    if(conv==1){resultat = Mnuc * (1 / sqrt(1 - (input* input * 0.3204 * 0.3204 * Znuc* Znuc ) / ( Anuc * Anuc )) - 1);}
    if(conv==0){resultat = Anuc * (sqrt(1-pow(1+input/Mnuc,-2))/(0.3204*Znuc));}
    return resultat;
}

/*
 double errorIntegral(TH1F* h, int bin1, int bin2) {
     double e = 0;
     for (int i = bin1; i < bin2 + 1; i++) {
         e = e + sqrt(h->GetBinContent(i));
     }
     return e;
 }
 Int_t binConversion1D(TH1* signal, double value) {
     TAxis* xaxis = signal->GetXaxis();
     Int_t x = xaxis->FindBin(value);
     return x;
 }*/
