#ifndef CODEA_GAMMA_H
#define CODEA_GAMMA_H

#include "main.h"
#include "Extraction_SRIM_stopping_powers.h"
#include "Tools.h"
#include "kinematics.h"
#include "MC_sim.h"

#define NUCLEI_REACTION 4
#define TARGET_PARAMETERS 4
#define SPECTRO_PARAMETERS 5

int CodeA_gamma(TChain* data,  TTree* build_back_data_meas,  TTree* build_back_data_sim[2],  double Ebeam[2], Double_t MReaction[NUCLEI_REACTION], Double_t AZreaction[2][NUCLEI_REACTION], Double_t target_dim[TARGET_PARAMETERS], Double_t level_spectro[SPECTRO_PARAMETERS], Double_t AGATA_position[3], Double_t AGATA_angle[2], Double_t AGATA_resolution[2],Double_t AGATA_E_range_noise[2], double theta_VAMOS_max, int extract_width,  double Loop_width[2][3]);

#endif
