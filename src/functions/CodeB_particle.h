#ifndef CODEB_PARTICLE_H
#define CODEB_PARTICLE_H

#include "main.h"
#include "Extraction_SRIM_stopping_powers.h"
#include "Tools.h"
#include "kinematics.h"
#include "MC_sim.h"

#define NUCLEI_REACTION 4
#define TARGET_PARAMETERS 4
#define DECAY_SPECTRO_PARAMETERS 6
#define MAX_LINES 5000
#define NUCLEI 5

int CodeB_particle(TChain* data,  TTree* build_back_data_meas,   TTree* build_back_data_sim,  double Ebeam[2], Double_t MReaction[NUCLEI_REACTION], Double_t AZReaction[2][NUCLEI_REACTION], Double_t MDecay[NUCLEI_REACTION], Double_t target_dim[TARGET_PARAMETERS], Double_t level_spectro[SPECTRO_PARAMETERS], Double_t decay_spectro[DECAY_SPECTRO_PARAMETERS], double theta_VAMOS_max, Double_t SPIDER_geometry[3], int extract_width);

#endif


