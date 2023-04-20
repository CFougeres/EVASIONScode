#include "CodeB_particle.h"

#define NUCLEI_REACTION 4
#define NUCLEI_DECAY 2
#define TARGET_PARAMETERS 4
#define DECAY_SPECTRO_PARAMETERS 6
#define SPECTRO_PARAMETERS 5
#define MAX_LINES 5000
#define NUCLEI 5

//************************************************
//   CodeB_particle
//   Analysis and simultion of particle emissions
//   Reconstruction of state excitation energies
// ***********************************************  
int CodeB_particle(TChain* data,  TTree* build_back_data_meas,   TTree* build_back_data_sim,  double Ebeam[2], Double_t MReaction[NUCLEI_REACTION], Double_t AZReaction[2][NUCLEI_REACTION], Double_t MDecay[NUCLEI_REACTION], Double_t target_dim[TARGET_PARAMETERS], Double_t level_spectro[SPECTRO_PARAMETERS], Double_t decay_spectro[DECAY_SPECTRO_PARAMETERS], double theta_VAMOS_max, Double_t SPIDER_geometry[3], int extract_width){
    //SRIM
    Double_t IonEnergy[NUCLEI][MAX_LINES]; Double_t dEdx_e[NUCLEI][MAX_LINES] ; Double_t dEdx_n[NUCLEI][MAX_LINES];int Nfile_energy[NUCLEI];//sp keV/micron
    int  SRIMextra = Extraction_SRIM_stopping_powers(IonEnergy, dEdx_e,  dEdx_n,Nfile_energy);
  
    int extraction_dataExp(TChain* data, TTree* build_back_data, double Ebeam, double theta_VAMOS_max, double pol2Eloss[3], Double_t MReaction[NUCLEI_REACTION], Double_t AZReaction[2][NUCLEI_REACTION],  Double_t MDecay[NUCLEI_REACTION], Double_t level_spectro[SPECTRO_PARAMETERS], Double_t decay_spectro[DECAY_SPECTRO_PARAMETERS], Double_t target_dim[TARGET_PARAMETERS], Double_t SPIDER_geometry[3]) ;
   
    double pol2Eloss[3]={0};          double pol2Ex[3]={0};
    double applied_MCsim_Brho_to_ElossEx= MCsim_Brho_to_ElossEx(Ebeam, MReaction, AZReaction, level_spectro, target_dim , pol2Eloss, pol2Ex, theta_VAMOS_max, IonEnergy, dEdx_e,dEdx_n,Nfile_energy);
    if(extract_width==0){
        //Conversion from raw measurement
        int applied_extraction_dataExp=extraction_dataExp(data, build_back_data_meas, Ebeam[0], theta_VAMOS_max,  pol2Eloss,  MReaction, AZReaction, MDecay, level_spectro, decay_spectro, target_dim, SPIDER_geometry);
        //MC sim
    }
    
    if(extract_width==1){
        //quantification Ip & Ig
        //to be added
    }
    
    return 1;
}

int extraction_dataExp(TChain* data, TTree* build_back_data, double Ebeam, double theta_VAMOS_max, double pol2Eloss[3], Double_t MReaction[NUCLEI_REACTION], Double_t AZReaction[2][NUCLEI_REACTION],  Double_t MDecay[NUCLEI_REACTION], Double_t level_spectro[SPECTRO_PARAMETERS], Double_t decay_spectro[DECAY_SPECTRO_PARAMETERS], Double_t target_dim[TARGET_PARAMETERS], Double_t SPIDER_geometry[3]) {
    
    double Ex_ejectile, Ex_particle, E_particle;
    Float_t Brho, ThetaLdeg, Theta, Phi, Si_dE, Si_Eres,Si_Ep, Si_theta;
    int nbTrack;
    data->SetBranchAddress("Brho", &Brho);
    data->SetBranchAddress("Phi", &Phi); data->SetBranchAddress("Theta", &Theta); data->SetBranchAddress("ThetaLdeg", &ThetaLdeg);
    data->SetBranchAddress("Si_dE", &Si_dE); data->SetBranchAddress("Theta", &Theta); data->SetBranchAddress("Si_Eres", &Si_Eres);
    int nevent = data->GetEntries();
    
    build_back_data->SetBranchAddress("Ex_ejectile", &Ex_ejectile);
    build_back_data->SetBranchAddress("Ex_particle", &Ex_particle);
    //Angles
    double  Doppler_angle;
    double  angleRecoil[2]; double  angleEjectil[2]; double betaRecoil;
    int derivation_kinematics;
    for (int i = 0; i < nevent; i++) {
        data->GetEntry(i);
        if (Brho > 0 && Si_dE>0) {
            // ejectile angles
            angleEjectil[0] =  ThetaLdeg ;
            angleEjectil[1] = atan(sin(Phi * 0.001) / (-cos(Phi * 0.001) * sin(Theta * 0.001))) * 180. / TMath::Pi();
            if (-sin(Theta * 0.001) < 0) { angleEjectil[1]= angleEjectil[1]+ 180.; }
            if (angleEjectil[1] < 0) { angleEjectil[1]= angleEjectil[1] + 360.; }
            if (angleEjectil[0] < theta_VAMOS_max) {
                derivation_kinematics = VAMOS_kinematics(angleEjectil, angleRecoil, betaRecoil, Brho, pol2Eloss, Ebeam, MReaction, AZReaction);
                Ex_ejectile = VAMOS_Ex(angleEjectil,Brho, pol2Eloss,Ebeam, MReaction, AZReaction);
                if (Ex_ejectile > level_spectro[2]*0.9  && Ex_ejectile <(level_spectro[2]+decay_spectro[5])*1.1 ) {
                    Si_Ep = Si_dE+Si_Eres;
                    E_particle = Si_Ep_to_Ep(Si_Ep, Si_theta, target_dim, SPIDER_geometry);
                    Ex_particle = Ep_to_Ex( E_particle, Si_theta, betaRecoil, MDecay,level_spectro, decay_spectro);
                    build_back_data->Fill();
                }
            }
        }
    }
    return 1;
}
