#include "kinematics.h"

#define NUCLEI_REACTION 4

using namespace std;

double Brho_to_Eloss(double Brho, double pol2Eloss[3]) {
    double Eloss = pol2Eloss[0] + pol2Eloss[1] * Brho + pol2Eloss[2] * Brho * Brho;
    return Eloss;
}

double VAMOS_kinematics(double angleEjectil[2], double angleRecoil[2], double betaRecoil, double Brho, double  pol2Eloss[3], double Ebeam, Double_t MReaction[NUCLEI_REACTION], Double_t AZReaction[2][NUCLEI_REACTION]){
    
    double c = 3.0 * pow(10, 8); // speed velocity m.s-1

    double gammabeam = Ebeam / MReaction[0]+ 1.;
    double impulsionbeam = sqrt(Ebeam * Ebeam + 2 * Ebeam * MReaction[0]) / c;

    int conv=1;
    double Eejectil = Brho_Ekin_conv(Brho, MReaction[3], AZReaction[0][3], AZReaction[1][3], conv);

    double thetaRecoil, ERecoil, gammaRecoil, impulsionRecoil;
    double gammaEjectil, impulsionEjectil;
    
    double ElossEjectil = Brho_to_Eloss(Brho, pol2Eloss);
    double Eejectil_Reaction = Eejectil+ElossEjectil;
    
    gammaEjectil = Eejectil_Reaction  / MReaction[3] + 1;
    impulsionEjectil = sqrt(Eejectil_Reaction * Eejectil_Reaction  + 2 * MReaction[3] * Eejectil_Reaction ) / c;
    
    gammaRecoil = sqrt(1. + (MReaction[3] * MReaction[3] * (gammaEjectil * gammaEjectil - 1) + MReaction[0] * MReaction[0] * (gammabeam * gammabeam - 1) - 2. * cos(angleEjectil[0] * TMath::Pi() / 180.) * MReaction[0]* MReaction[3] * sqrt(gammabeam * gammabeam - 1) * sqrt(gammaEjectil * gammaEjectil - 1)) / (MReaction[2] * MReaction[2]));
    ERecoil = MReaction[2]* (gammaRecoil - 1.);
    impulsionRecoil = sqrt(ERecoil * ERecoil + 2 * MReaction[2] * ERecoil) / c;
    betaRecoil= sqrt(1-pow(gammaRecoil,-2));
    
    angleRecoil[0] = acos((impulsionbeam - impulsionEjectil * cos(angleEjectil[0] * TMath::Pi() / 180)) / impulsionRecoil) * 180. / TMath::Pi();
    angleRecoil[1] = angleEjectil[1] * 180. / TMath::Pi();
    if (cos(angleEjectil[0] * TMath::Pi() / 180) > 0) { angleRecoil[1] = angleRecoil[1] + 180.; }
    if (angleRecoil[1] < 0) angleRecoil[1] = angleRecoil[1] + 360.0;
    
    return 1.0;
}

double lab_to_Doppler_angle(double angleRecoil[2], double angleGamma[2]) {
    double stheta23Mg = sin(angleRecoil[0] * TMath::Pi() / 180.);
    double ctheta23Mg = cos(angleRecoil[0] * TMath::Pi() / 180.);
    double sphi23Mg = sin(angleRecoil[1] * TMath::Pi() / 180.);
    double cphi23Mg = cos(angleRecoil[1] * TMath::Pi() / 180.);
    double cthetagamma = cos(angleGamma[0] * TMath::Pi() / 180.);
    double sthetagamma = sin(angleGamma[0] * TMath::Pi() / 180.);
    double cphigamma = cos(angleGamma[1] * TMath::Pi() / 180);
    double sphigamma = sin(angleGamma[1] * TMath::Pi() / 180);
    double D_angle = stheta23Mg * cphi23Mg * sthetagamma * cphigamma + sphi23Mg * stheta23Mg * sthetagamma * sphigamma + ctheta23Mg * cthetagamma;
    D_angle = 180.0 / TMath::Pi() * acos(D_angle);
    return  D_angle;
}


