#include "Extraction_SRIM_stopping_powers.h"
#define MAX_LINES 5000
#define NUCLEI 5
using namespace std;

//********************************************
// Extraction function of stopping powers along energy
// from a table generated with SRIM code
//********************************************
int Extraction_SRIM_stopping_powers(Double_t IonEnergy[NUCLEI][MAX_LINES], Double_t dEdx_e[NUCLEI][MAX_LINES], Double_t dEdx_n[NUCLEI][MAX_LINES], int Nfile_energy[NUCLEI]){
    
    string path = gSystem->pwd();
    
    string fileSPbeam = path+"/src/functions/SP_SRIM/beam_in_matter.dat";
    string fileSPrecoil = path+"/src/functions/SP_SRIM/recoil_in_matter.dat";
    string fileSPejectil = path+"/src/functions//SP_SRIM/ejectil_in_matter.dat";
    string fileSPparticle = path+"/src/functions//SP_SRIM/particle_in_matter.dat";
    string fileSPparticle_si = path+"/src/functions//SP_SRIM/particle_in_silicon.dat";
    
    int LoadSRIMFile(string FileName, Double_t IonEnergy[MAX_LINES], Double_t dEdx_e[MAX_LINES], Double_t dEdx_n[MAX_LINES], int N_energy);

    
    int extractSPSRIM= LoadSRIMFile(fileSPbeam, IonEnergy[0], dEdx_e[0], dEdx_n[0],  Nfile_energy[0]);
    extractSPSRIM= LoadSRIMFile(fileSPparticle, IonEnergy[1], dEdx_e[1],dEdx_n[1],Nfile_energy[1]);
    extractSPSRIM= LoadSRIMFile(fileSPejectil, IonEnergy[2], dEdx_e[2], dEdx_n[2],Nfile_energy[2]);
    extractSPSRIM= LoadSRIMFile(fileSPparticle, IonEnergy[3], dEdx_e[3], dEdx_n[3], Nfile_energy[3]);
    extractSPSRIM= LoadSRIMFile(fileSPparticle_si, IonEnergy[4], dEdx_e[4], dEdx_n[4],Nfile_energy[4]);

    return 1;
}
int LoadSRIMFile(string FileName, Double_t IonEnergy[], Double_t dEdx_e[], Double_t dEdx_n[], int N_energy)
{
    int GoodELossFile; int   Energy_in_range;
 // double IonEnergy, dEdx_e, dEdx_n;
  string aux, unit;
  int str_counter=0;
  ifstream Read(FileName.c_str());
 // this->FileName = FileName;
  int last_point = 0;
  if(!Read.is_open()) {
    cout << "*** EnergyLoss Error: SRIM file " << FileName << " was not found." << endl;
    GoodELossFile = 0;
  }
  else {
    GoodELossFile = 1;
    Energy_in_range = 1;
    // Read all the string until you find "Straggling", then read the next 7 strings.
    // (this method is not elegant at all but there is no time to make it better)
    do
      Read >> aux;
    while (aux!="Straggling");
    for (int i=0; i<7; i++)
      Read >> aux;

    do {
      Read >> aux;
      str_counter++;
    } while (aux!="-----------------------------------------------------------");
    Read.close();
   
    str_counter--;
    int points = str_counter/10;
   
   /* // Create the arrays depending on the number rows in the file.
    this->IonEnergy = new double[points];
    this->dEdx_e = new double[points];
    this->dEdx_n = new double[points];
 */
    // Go to the begining of the file and read it again to now save the info in the
    // newly created arrays.
    Read.open(FileName.c_str());
    do
      Read >> aux;
    while (aux!="Straggling");

    for (int i=0; i<7; i++)
      Read >> aux;
    
    for (int p=0; p<points; p++) {
      Read >> IonEnergy[p] >> unit >> dEdx_e[p]>> dEdx_n[p]>> aux >> aux >> aux >> aux >> aux >> aux;
      if (unit=="eV")
    IonEnergy[p] *= 1E-6;
      else if (unit=="keV")
    IonEnergy[p] *= 0.001;
      else if (unit=="GeV")
    IonEnergy[p] *= 1000;
      //cout << p << " " << IonEnergy << " " << unit << " " << dEdx_e << " " << dEdx_n << endl;
     /* this->IonEnergy[p] = IonEnergy;
      this->dEdx_e[p] = dEdx_e;
      this->dEdx_n[p] = dEdx_n;*/
    }
      N_energy=points;
  }
  return GoodELossFile;
}
