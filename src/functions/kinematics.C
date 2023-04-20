#include "kinematics.h"

#define NUCLEI_REACTION 4
#define NUCLEI_DECAY 2
#define DECAY_SPECTRO_PARAMETERS 6
#define SPECTRO_PARAMETERS 5
#define TARGET_PARAMETERS 4

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

double VAMOS_Ex(double angleEjectil[2], double Brho, double  pol2Eloss[3], double Ebeam, Double_t MReaction[NUCLEI_REACTION], Double_t AZReaction[2][NUCLEI_REACTION]){

    double c = 3.0 * pow(10, 8); // speed velocity m.s-1
    
    double thetaRecoil, ERecoil, gammaRecoil, impulsionRecoil, Ex;
    double gammaEjectil, impulsionEjectil, Eejectil, ElossEjectil,Eejectil_Reaction  ;
    
    double Qreaction=MReaction[3]+MReaction[2]-MReaction[1]-MReaction[0];
    
    double gammabeam = Ebeam / MReaction[0]+ 1.;
    double impulsionbeam = sqrt(Ebeam * Ebeam + 2 * Ebeam * MReaction[0]) / c;

    int conv=1;
    Eejectil = Brho_Ekin_conv(Brho, MReaction[3], AZReaction[0][3], AZReaction[1][3], conv);
    ElossEjectil = Brho_to_Eloss(Brho, pol2Eloss);
    Eejectil_Reaction = Eejectil+ElossEjectil;
    gammaEjectil = Eejectil_Reaction  / MReaction[3] + 1;
    impulsionEjectil = sqrt(Eejectil_Reaction * Eejectil_Reaction  + 2 * MReaction[3] * Eejectil_Reaction ) / c;

    gammaRecoil = sqrt(1. + (MReaction[3] * MReaction[3] * (gammaEjectil * gammaEjectil - 1) + MReaction[0] * MReaction[0] * (gammabeam * gammabeam - 1) - 2. * cos(angleEjectil[0] * TMath::Pi() / 180.) * MReaction[0]* MReaction[3] * sqrt(gammabeam * gammabeam - 1) * sqrt(gammaEjectil * gammaEjectil - 1)) / (MReaction[2] * MReaction[2]));
    ERecoil = MReaction[2]* (gammaRecoil - 1.);
    
    Ex = Ebeam -   ERecoil  - (Eejectil_Reaction) + Qreaction;

    return Ex;
}
double Ep_to_Ex(double Ep, double thetap, double betaRecoil,  Double_t MDecay[NUCLEI_REACTION], Double_t level_spectro[SPECTRO_PARAMETERS], Double_t decay_spectro[DECAY_SPECTRO_PARAMETERS] ) {
    double Ex;
    double rM = MDecay[1] / (MDecay[0] * (MDecay[1]  + MDecay[0] ));
    double b = 1 - rM * ( decay_spectro[2]+level_spectro[2]);
    double vp = sqrt(1 - 1 / pow((Ep / MDecay[0] + 1), 2));
    double cthetap = cos(thetap * TMath::Pi()/ 180.);
    double sthetap = sin(thetap * TMath::Pi()/ 180.);
    double Af = (betaRecoil - vp * cthetap) / (betaRecoil* vp * cthetap - 1.);
    double v1_2 = Af * Af + vp * vp * sthetap * sthetap * (1. - betaRecoil* Af) * (1. - betaRecoil * Af) / (1 - betaRecoil* betaRecoil);
    Ex = (1. / sqrt(1. - v1_2) - b) / rM;
    return Ex;
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
double Si_Ep_to_Ep(double Si_Ep, double angle_particle, Double_t target_dim[TARGET_PARAMETERS], Double_t SPIDER_geometry[3]) {
    int angle_index= int((tan(angle_particle*TMath::Pi()/190.)*SPIDER_geometry[0] - SPIDER_geometry[1] - 1.75)*16./(SPIDER_geometry[1]+1));
    double Ep;
    double p[6];
    int target_thickness=0;
    if(target_dim[1]>2){target_thickness ==0;}//thick
    if(target_dim[1]<1){target_thickness ==1;}//thin
    if (target_thickness ==0 ) {//thick
        if (angle_index == 0) { p[0] = 2.77473; p[1] = 0.357285; p[2] = 0.119577; p[3] = -0.0140169; p[4] = 0.000883156; p[5] = -2.25265e-05; }
        if (angle_index == 1) { p[0] = 2.77875; p[1] = 0.358112; p[2] = 0.119032; p[3] = -0.0139127; p[4] = 0.000874611; p[5] = -2.22705e-05; }
        if (angle_index == 2) { p[0] = 2.78272; p[1] = 0.359275; p[2] = 0.118337; p[3] = -0.0137811; p[4] = 0.000863835; p[5] = -2.19475e-05; }
        if (angle_index == 3) { p[0] = 2.78698; p[1] = 0.36043; p[2] = 0.117633; p[3] = -0.0136474; p[4] = 0.00085288; p[5] = -2.16189e-05; }
        if (angle_index == 4) { p[0] = 2.7914; p[1] = 0.361695; p[2] = 0.116873; p[3] = -0.0135034; p[4] = 0.000841078; p[5] = -2.12649e-05; }
        if (angle_index == 5) { p[0] = 2.79573; p[1] = 0.363448; p[2] = 0.115894; p[3] = -0.0133191; p[4] = 0.00082597; p[5] = -2.08111e-05; }
        if (angle_index == 6) { p[0] = 2.80034; p[1] = 0.365105; p[2] = 0.114948; p[3] = -0.0131406; p[4] = 0.000811341; p[5] = -2.03717e-05; }
        if (angle_index == 7) { p[0] = 2.80488; p[1] = 0.367167; p[2] = 0.11382; p[3] = -0.0129286; p[4] = 0.000793962; p[5] = -1.98492e-05; }
        if (angle_index == 8) { p[0] = 2.80966; p[1] = 0.369236; p[2] = 0.11268; p[3] = -0.012714; p[4] = 0.000776368; p[5] = -1.93201e-05; }
        if (angle_index == 9) { p[0] = 2.81455; p[1] = 0.371472; p[2] = 0.111459 ; p[3] = -0.0124847; p[4] = 0.000757562; p[5] = -1.87544e-05; }
        if (angle_index == 10) { p[0] = 2.81955; p[1] = 0.373896; p[2] = 0.110149; p[3] = -0.0122386; p[4] = 0.000737375; p[5] = -1.81468e-05; }
        if (angle_index == 11) { p[0] = 2.82442; p[1] = 0.376758; p[2] = 0.108643; p[3] = -0.0119566; p[4] = 0.000714222; p[5] = -1.74492e-05; }
        if (angle_index == 12) { p[0] = 2.82957; p[1] = 0.379614; p[2] = 0.107129; p[3] = -0.0116727; p[4] = 0.000690923; p[5] = -1.67472e-05; }
        if (angle_index == 13) { p[0] = 2.83475; p[1] = 0.382625; p[2] = 0.105544; p[3] = -0.0113758; p[4] = 0.000666544; p[5] = -1.60123e-05; }
        if (angle_index == 14) { p[0] = 2.84004; p[1] = 0.385831; p[2] = 0.103867; p[3] = -0.0110616; p[4] = 0.000640739; p[5] = -1.5234e-05; }
        if (angle_index == 15) { p[0] = 2.84524; p[1] = 0.389459; p[2] = 0.101998; p[3] = -0.010712; p[4] = 0.000612005; p[5] = -1.43666e-05; }
    }
    if (target_thickness == 1) {//thin
        if (angle_index == 0) { p[0] = 2.45824; p[1] = 0.38364; p[2]= 0.122151; p[3] = -0.0149843; p[4] = 0.000976974; p[5] = -2.56245e-05; }
        if (angle_index == 1) { p[0] = 2.4617; p[1] = 0.384549; p[2] = 0.121563; p[3] = -0.0148693; p[4] = 0.000967333; p[5] = -2.53294e-05; }
        if (angle_index == 2) { p[0] = 2.46541; p[1] = 0.385386; p[2] = 0.120998; p[3] = -0.0147583; p[4] = 0.000958018; p[5] = -2.50444e-05;}
        if (angle_index == 3) { p[0] = 2.46907; p[1] = 0.38656; p[2] = 0.120276; p[3] = -0.014618; p[4] = 0.00094626; p[5] = -2.46842e-05;}
        if (angle_index == 4) { p[0] = 2.47304; p[1] = 0.387674; p[2] = 0.11957; p[3] = -0.0144802; p[4] = 0.000934695; p[5] = -2.433e-05; }
        if (angle_index == 5) { p[0] = 2.47696; p[1] = 0.389141; p[2] = 0.1187; p[3] = -0.0143116; p[4] = 0.000920567; p[5] = -2.38969e-05; }
        if (angle_index == 6) { p[0] = 2.4811; p[1] = 0.390557; p[2] = 0.117844; p[3] = -0.0141455; p[4] = 0.000906637; p[5] = -2.34699e-05; }
        if (angle_index == 7) { p[0] = 2.48524; p[1] = 0.392315; p[2] = 0.116829; p[3] = -0.0139493; p[4] = 0.000890195; p[5] = -2.29656e-05; }
        if (angle_index == 8) { p[0] = 2.48957; p[1] = 0.39407; p[2] = 0.115807; p[3] = -0.0137516; p[4] = 0.000873617; p[5] = -2.2457e-05; }
        if (angle_index == 9) { p[0] = 2.49402 ; p[1] = 0.39595; p[2] = 0.114721; p[3] = -0.0135419; p[4] = 0.000856035; p[5] = -2.19175e-05; }
        if (angle_index == 10) { p[0] = 2.49841; p[1] = 0.398221; p[2] = 0.113455; p[3] = -0.0132979; p[4] = 0.00083558; p[5] = -2.12894e-05; }
        if (angle_index == 11) { p[0] = 2.50302; p[1] = 0.400475; p[2] = 0.112188; p[3] = -0.0130538; p[4] = 0.000815108; p[5] = -2.06606e-05; }
        if (angle_index == 12) { p[0] = 2.5077 ; p[1] = 0.402846; p[2] = 0.110863; p[3] = -0.0127986; p[4] = 0.000793707; p[5] = -2.00032e-05; }
        if (angle_index == 13) { p[0] = 2.51232; p[1] = 0.405614; p[2] = 0.109354; p[3] = -0.0125088; p[4] = 0.0007694; p[5] = -1.92559e-05; }
        if (angle_index == 14) { p[0] = 2.51714; p[1] = 0.408373; p[2] = 0.107841; p[3] = -0.012218; p[4] = 0.000745002; p[5] = -1.85057e-05; }
        if (angle_index == 15) { p[0] = 2.52201; p[1] = 0.411275; p[2] = 0.10626; p[3] = -0.0119142; p[4] = 0.00071952; p[5] = -1.7722e-05; }
    }
    Ep = p[0] + p[1] * Si_Ep + p[2] * pow(Si_Ep, 2) + p[3] * pow(Si_Ep,3) + p[4] * pow(Si_Ep,4) + p[5] * pow(Si_Ep, 5);
    return Ep;
}
