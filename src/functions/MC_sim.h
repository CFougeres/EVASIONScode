#ifndef MC_SIM_H
#define MC_SIM_H

#include "main.h"
#include "Extraction_SRIM_stopping_powers.h"
#include "Tools.h"
#include "kinematics.h"

#define MAX_LINES 5000
#define NUCLEI 5
#define NUCLEI_REACTION 4
#define TARGET_PARAMETERS 4
#define SPECTRO_PARAMETERS 5

double MCsim_kinematics( double Ebeam[2], double Ex,Double_t MReaction[NUCLEI_REACTION], Double_t AZReaction[2][NUCLEI_REACTION], Double_t target_dim[TARGET_PARAMETERS], double BrhoMeas, double Brholim[2], double Eloss, double theta_VAMOS_max,Double_t IonEnergy[NUCLEI][MAX_LINES], Double_t dEdx_e[NUCLEI][MAX_LINES], Double_t dEdx_n[NUCLEI][MAX_LINES], int Nfile_energy[NUCLEI]);

double MCsim_emission(TTree* build_back_data_sim,  double Ebeam[2],  Double_t MReaction[NUCLEI_REACTION],  Double_t AZReaction[2][NUCLEI_REACTION], Double_t target_dim[TARGET_PARAMETERS], Double_t level_spectro[SPECTRO_PARAMETERS], double theta_VAMOS_max, Double_t AGATA_resolution[2], Double_t AGATA_angle[2],  Double_t IonEnergy[NUCLEI][MAX_LINES], Double_t dEdx_e[NUCLEI][MAX_LINES], Double_t dEdx_n[NUCLEI][MAX_LINES], int Nfile_energy[NUCLEI]);

double MCsim_Brho_to_ElossEx(double Ebeam[2],  Double_t MReaction[NUCLEI_REACTION], Double_t AZReaction[2][NUCLEI_REACTION], Double_t level_spectro[SPECTRO_PARAMETERS], Double_t target_dim[TARGET_PARAMETERS],  double pol2Eloss[3], double pol2Ex[3], double theta_VAMOS_max, Double_t IonEnergy[NUCLEI][MAX_LINES], Double_t dEdx_e[NUCLEI][MAX_LINES], Double_t dEdx_n[NUCLEI][MAX_LINES], int Nfile_energy[NUCLEI]);
    
int MCsim_noise(TChain* data, TTree* build_back_data_sim, double Ebeam, double theta_VAMOS_max, double pol2Eloss[3], Double_t MReaction[NUCLEI_REACTION], Double_t AZreaction[2][NUCLEI_REACTION], Double_t AGATA_E_range_noise[2], Double_t AGATA_angle[2], Double_t AGATA_position[3],  Double_t AGATA_resolution[2], Double_t level_spectro[SPECTRO_PARAMETERS]);
    
#endif
