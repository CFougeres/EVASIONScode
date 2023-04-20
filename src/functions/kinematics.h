#ifndef KINEMATICS_H
#define KINEMATICS_H

#include "main.h"
#include "Tools.h"
#include "CodeA_gamma.h"
#include "CodeB_particle.h"

#define NUCLEI_REACTION 4
#define NUCLEI_DECAY 2
#define DECAY_SPECTRO_PARAMETERS 6
#define SPECTRO_PARAMETERS 5
#define TARGET_PARAMETERS 4

double VAMOS_kinematics(double angleEjectil[2], double angleRecoil[2], double betaRecoil, double Brho, double  pol2Eloss[3], double Ebeam, Double_t MReaction[NUCLEI_REACTION], Double_t AZReaction[2][NUCLEI_REACTION]);

double lab_to_Doppler_angle(double angleRecoil[2], double angleGamma[2]);

double Brho_to_Eloss(double Brho, double pol2Eloss[3]) ;

double VAMOS_Ex(double angleEjectil[2], double Brho, double  pol2Eloss[3], double Ebeam, Double_t MReaction[NUCLEI_REACTION], Double_t AZReaction[2][NUCLEI_REACTION]);

double Ep_to_Ex(double Ep, double thetap, double betaRecoil,  Double_t MDecay[NUCLEI_REACTION], Double_t level_spectro[SPECTRO_PARAMETERS], Double_t decay_spectro[DECAY_SPECTRO_PARAMETERS]) ;

double Si_Ep_to_Ep(double Si_Ep, double angle_particle, Double_t target_dim[TARGET_PARAMETERS], Double_t SPIDER_geometry[3]);
#endif
