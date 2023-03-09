#ifndef KINEMATICS_H
#define KINEMATICS_H

#include "main.h"
#include "Tools.h"
#include "CodeA_gamma.h"

#define NUCLEI_REACTION 4

double VAMOS_kinematics(double angleEjectil[2], double angleRecoil[2], double betaRecoil, double Brho, double  pol2Eloss[3], double Ebeam, Double_t MReaction[NUCLEI_REACTION], Double_t AZReaction[2][NUCLEI_REACTION]);

double lab_to_Doppler_angle(double angleRecoil[2], double angleGamma[2]);

double Brho_to_Eloss(double Brho, double pol2Eloss[3]) ;

#endif
