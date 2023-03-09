#include "MC_sim.h"

#define MAX_LINES 5000
#define NUCLEI 5
#define NUCLEI_REACTION 4
#define TARGET_PARAMETERS 4
#define SPECTRO_PARAMETERS 5

using namespace std;

double MCsim_kinematics(double Ebeam[2], double Ex, Double_t MReaction[NUCLEI_REACTION], Double_t AZReaction[2][NUCLEI_REACTION], Double_t target_dim[TARGET_PARAMETERS], double BrhoMeas, double Brholim[2], double Eloss, double theta_VAMOS_max, Double_t IonEnergy[NUCLEI][MAX_LINES], Double_t dEdx_e[NUCLEI][MAX_LINES], Double_t dEdx_n[NUCLEI][MAX_LINES], int Nfile_energy[NUCLEI]){
    //RANDOM FCTs
    TRandom* rEnergy_Beam = new TRandom();//gaussian reponse of beam
    srand(time(NULL));
    gRandom->SetSeed(0);
    //MC statistics
    int Nbeam = 100;
    int Nreactions = 250;
    //REACTION
    double Efais = Ebeam[0];
    double Pfais = sqrt(Ebeam[0]* Ebeam[0]+ 2. * Ebeam[0]*MReaction[0]);
    Double_t masses[2] = { 0.001 * (MReaction[2]+ Ex), 0.001*MReaction[3] };
    TLorentzVector target(0.0, 0.0, 0.0, MReaction[1] * 0.001);
    TLorentzVector beam(0.0, 0.0, Pfais * 0.001, (Ebeam[0] + MReaction[0]) * 0.001);
    TLorentzVector Tr = beam + target;
    TLorentzVector * pRecoil;
    TLorentzVector * pEjectil;
    TGenPhaseSpace event;
    event.SetDecay(Tr, 2, masses);
    // Variables
    double dt = 0.1;
    double dbeam = 0.0; //Angstrom
    double dx = 5.0; //Angstrom
    double distance, dd; //micron
    int stop;
    double vbeamlab, vEjectil_lab, eEjectil_lab;
    double angleEjectil_lab[2]={0};
    int compteur =0;
    Eloss=0;double tempEloss=0; double tempBrho=0;
    Brholim[0]=100;       Brholim[1]=0;
    int conv=0;
    for (int i = 0; i < Nbeam; i++) {
        Efais = rEnergy_Beam->Gaus(Ebeam[0], Ebeam[1]);
        while (dbeam < target_dim[0]*1000 && Efais>0) {
            Efais = Efais - 0.001 * interpol(IonEnergy[0], dEdx_e[0], Nfile_energy[0], Efais * 1000) * dx * 0.0001 - 0.001 * interpol(IonEnergy[0], dEdx_n[0], Nfile_energy[0], Efais * 1000) * dx * 0.0001;
            dbeam = dbeam + dx;
            if (Efais < 0) { Efais = 0.; }
            vbeamlab = sqrt(1 - 1 / pow((Efais / MReaction[0] + 1), 2));
            Pfais = sqrt(Efais * Efais + 2. * Efais * MReaction[0]);
            beam.SetPxPyPzE(0, 0, Pfais * 0.001, (Efais + MReaction[0]) * 0.001);
            Tr = beam + target;
            event.SetDecay(Tr, 2, masses);
            for (int ir = 0; ir < Nreactions; ir++) {
                event.Generate();
                pEjectil = event.GetDecay(1);
                //Ejectil
                angleEjectil_lab[0] = pEjectil->Theta() * 180. / TMath::Pi();
                angleEjectil_lab[1]= pEjectil->Phi() * 180. / TMath::Pi();
                if (angleEjectil_lab[1]< 0) angleEjectil_lab[1] = angleEjectil_lab[1] + 360.0;
                vEjectil_lab = pEjectil->Beta();
                eEjectil_lab = (1 / sqrt(1 - pow(vEjectil_lab, 2)) - 1) * MReaction[3];
                if ( angleEjectil_lab[0] < theta_VAMOS_max && eEjectil_lab>15) {
                    tempEloss = eEjectil_lab - loss_E(eEjectil_lab , target_dim[1] + target_dim[2] - dbeam * 0.0001, IonEnergy[2], dEdx_n[2], dEdx_n[2], Nfile_energy[2]);
                    if (tempEloss < 0) {tempEloss=0;}
                    tempBrho = Brho_Ekin_conv(tempEloss, MReaction[3], AZReaction[0][3], AZReaction[1][3], conv);
                    Eloss = Eloss + tempEloss;
                    BrhoMeas = BrhoMeas + tempEloss;
                    compteur+=1;
                    if(tempBrho>Brholim[1]){Brholim[1]=tempBrho;}
                    if(tempBrho<Brholim[0]){Brholim[0]=tempBrho;}
                }
            }
        }
    }
    Eloss=Eloss/double(compteur);
    BrhoMeas=BrhoMeas/double(compteur);
    return 1.0;
}



double MCsim_emission(TTree* build_back_data_sim,  double Ebeam[2],  Double_t MReaction[NUCLEI_REACTION], Double_t AZReaction[2][NUCLEI_REACTION],  Double_t target_dim[TARGET_PARAMETERS], Double_t level_spectro[SPECTRO_PARAMETERS], double theta_VAMOS_max, Double_t AGATA_resolution[2], Double_t AGATA_angle[2],  Double_t IonEnergy[NUCLEI][MAX_LINES], Double_t dEdx_e[NUCLEI][MAX_LINES], Double_t dEdx_n[NUCLEI][MAX_LINES], int Nfile_energy[NUCLEI]){
    
    
    Double_t Beta_emission, Beta_reaction, DeltaBeta, Egamma_DS,theta_doppler;

    build_back_data_sim->SetBranchAddress("Beta_emission", &Beta_emission);
    build_back_data_sim->SetBranchAddress("Beta_reaction", &Beta_reaction);
    build_back_data_sim->SetBranchAddress("DeltaBeta", &DeltaBeta);
    build_back_data_sim->SetBranchAddress("Egamma_DS", &Egamma_DS);
    build_back_data_sim->SetBranchAddress("theta_doppler", &theta_doppler);
    
    //RANDOM FCTs
    TRandom* r1 = new TRandom1();//deexcitation, Poisson
    TRandom* rEnergy_Gamma = new TRandom();//gaussian reponse of AGATA
    TRandom* rAngle_Gamma = new TRandom();//gaussian reponse of AGATA
    TRandom* rEnergy_Beam = new TRandom();//gaussian reponse of beam
    srand(time(NULL));
    gRandom->SetSeed(0);
    //MC statistics
    int Nbeam = 100;
    int Nreactions = 250;
    double Ndexci = 5000; //such as Ndexci x dt > 500
    //REACTION
    double Efais = Ebeam[0];
    double Pfais = sqrt(Ebeam[0]* Ebeam[0]+ 2. * Ebeam[0]*MReaction[0]);
    Double_t masses[2] = { 0.001 * (MReaction[2]+ level_spectro[0]), 0.001*MReaction[3] };
    TLorentzVector target(0.0, 0.0, 0.0, MReaction[1] * 0.001);
    TLorentzVector beam(0.0, 0.0, Pfais * 0.001, (Ebeam[0] + MReaction[0]) * 0.001);
    TLorentzVector Tr = beam + target;
    TLorentzVector * pRecoil;
    TLorentzVector * pEjectil;
    TGenPhaseSpace event;
    event.SetDecay(Tr, 2, masses);
    // Variables
    double dt = 0.1;
    double dbeam = 0.0; //Angstrom
    double dx = 5.0; //Angstrom
    double distance, dd; //micron
    double time, Nt, Ntdt, bibi_deexcitation;
    int stop;
    double vbeamlab, vEjectil_lab, vRecoil_lab, eRecoil_lab, eEjectil_lab;
    double angleRecoil_lab[2]={0};
    double angleEjectil_lab[2]={0};
    double angleGamma[2]={0};
    float lembda;
    double R, cangle, thetagammaDoppler;
    double Energy_Gamma_Reponse, Angle_Gamma_Reponse;
    double beta_recoil_mean=0; int compteur =0;
    
    for (int i = 0; i < Nbeam; i++) {
        Efais = rEnergy_Beam->Gaus(Ebeam[0], Ebeam[1]);
        while (dbeam < target_dim[0]*1000 && Efais>0) {
            Efais = Efais - 0.001 * interpol(IonEnergy[0], dEdx_e[0], Nfile_energy[0], Efais * 1000) * dx * 0.0001 - 0.001 * interpol(IonEnergy[0], dEdx_n[0], Nfile_energy[0], Efais * 1000) * dx * 0.0001;
            dbeam = dbeam + dx;
            if (Efais < 0) { Efais = 0.; }
            vbeamlab = sqrt(1 - 1 / pow((Efais / MReaction[0] + 1), 2));
            Pfais = sqrt(Efais * Efais + 2. * Efais * MReaction[0]);
            beam.SetPxPyPzE(0, 0, Pfais * 0.001, (Efais + MReaction[0]) * 0.001);
            Tr = beam + target;
            event.SetDecay(Tr, 2, masses);
            for (int ir = 0; ir < Nreactions; ir++) {
                event.Generate();
                pRecoil = event.GetDecay(0);
                pEjectil = event.GetDecay(1);
                //Ejectil
                angleEjectil_lab[0] = pEjectil->Theta() * 180. / TMath::Pi();
                angleEjectil_lab[1]= pEjectil->Phi() * 180. / TMath::Pi();
                if (angleEjectil_lab[1]< 0) angleEjectil_lab[1] = angleEjectil_lab[1] + 360.0;
                vEjectil_lab = pEjectil->Beta();
                eEjectil_lab = (1 / sqrt(1 - pow(vEjectil_lab, 2)) - 1) * MReaction[3];
                if ( angleEjectil_lab[0] < theta_VAMOS_max && eEjectil_lab>20) {
                    //RECOIL
                    stop = 0;
                    angleRecoil_lab[0] = pRecoil->Theta() * 180. / TMath::Pi();
                    angleRecoil_lab[1]  = pRecoil->Phi() * 180. / TMath::Pi();
                    if (angleRecoil_lab[1] < 0) angleRecoil_lab[1]  = angleRecoil_lab[1]  + 360.0;
                    Beta_reaction = pRecoil->Beta();
                    vRecoil_lab=Beta_reaction;
                    eRecoil_lab = (1 / sqrt(1 - pow(Beta_reaction , 2)) - 1) * MReaction[2];
                    Nt = Ndexci;
                    time = 0;
                    lembda = 1.0 / level_spectro[3];
                    dd = 0;
                    distance = 0.0001 * dbeam;
                    while (Nt > 0 && distance < target_dim[1] + target_dim[2]+target_dim[3]){
                        dd = vRecoil_lab * 0.3 * dt;
                        if (distance < target_dim[1]|| distance >= target_dim[1] + target_dim[3])
                        {
                            eRecoil_lab = eRecoil_lab - 0.001 * interpol(IonEnergy[1], dEdx_e[1],  Nfile_energy[1], eRecoil_lab  * 1000) * dd - 0.001 * interpol(IonEnergy[1], dEdx_n[1],  Nfile_energy[1], eRecoil_lab  * 1000) * dd; //dd en micron
                            if (eRecoil_lab  > 0) {   vRecoil_lab = sqrt(1 - 1 / pow((eRecoil_lab  / MReaction[2] + 1), 2)); }
                            if (eRecoil_lab  < 0) {   vRecoil_lab = 0.; eRecoil_lab  = 0.; }
                        }
                        if (distance >= target_dim[1] && distance < target_dim[1] + target_dim[3]) {
                            eRecoil_lab  = eRecoil_lab ;
                            if (eRecoil_lab  > 0) {   vRecoil_lab = sqrt(1 - 1 / pow((eRecoil_lab  / MReaction[2] + 1), 2)); }
                            if (eRecoil_lab  < 0) {   vRecoil_lab = 0.; eRecoil_lab  = 0.; }
                        }
                        distance = distance + dd;
                        Ntdt = Ndexci * exp(-(time + dt) * lembda);
                        if (eRecoil_lab > 0) {
                            bibi_deexcitation = Nt - Ntdt;
                            if (bibi_deexcitation > 1) {
                                for (int kk = 0; kk < bibi_deexcitation; kk++) {
                                    angleGamma[0] = -1.0 + 2.0 * 0.0001 * (rand() % 10001);
                                    angleGamma[0]= 180.0 / TMath::Pi() * acos( angleGamma[0] );
                                    angleGamma[1] = 360.0 * 0.0001 * (rand() % 10001);
                                    thetagammaDoppler = lab_to_Doppler_angle(angleRecoil_lab, angleGamma);
                                    Egamma_DS = level_spectro[1] * sqrt(1 -   vRecoil_lab *   vRecoil_lab) / (1 -   vRecoil_lab * cos(TMath::Pi() * thetagammaDoppler / 180.0));
                                    Angle_Gamma_Reponse = rAngle_Gamma->Gaus(thetagammaDoppler, AGATA_resolution[1]);
                                    Energy_Gamma_Reponse = rEnergy_Gamma->Gaus(Egamma_DS * 1000.0, AGATA_resolution[0]);
                                    if( Angle_Gamma_Reponse>=AGATA_angle[0]&& Angle_Gamma_Reponse<=AGATA_angle[1]){
                                        R = Energy_Gamma_Reponse / (1000. * level_spectro[1]);
                                        cangle = cos(TMath::Pi() * Angle_Gamma_Reponse/ 180.0);
                                        Beta_emission = (cangle * R * R + sqrt(1 + R * R * cangle * cangle - R * R)) / (R * R * cangle * cangle + 1);
                                        DeltaBeta=Beta_reaction-Beta_emission;
                                        theta_doppler = Angle_Gamma_Reponse;
                                        Egamma_DS = Energy_Gamma_Reponse;
                                        build_back_data_sim->Fill();
                                        beta_recoil_mean+= vRecoil_lab; compteur +=1;
                                    }
                                }
                            }
                        }
                        if (eRecoil_lab <= 0) {
                            stop = 1;
                            if (Ntdt > 1) {
                                for (int kk = 0; kk < Ntdt; kk++) {
                                    DeltaBeta=Beta_reaction;
                                    Beta_emission=0;
                                    angleGamma[0] = -1.0 + 2.0 * 0.0001 * (rand() % 10001);
                                    angleGamma[0]= 180.0 / TMath::Pi() * acos( angleGamma[0] );
                                    angleGamma[1] = 360.0 * 0.0001 * (rand() % 10001);
                                    Angle_Gamma_Reponse  = rAngle_Gamma->Gaus(angleGamma[0], AGATA_resolution[1]);
                                    if( Angle_Gamma_Reponse>=AGATA_angle[0]&& Angle_Gamma_Reponse<=AGATA_angle[1]){
                                        Energy_Gamma_Reponse = rEnergy_Gamma->Gaus(level_spectro[1] * 1000.0, AGATA_resolution[0]);
                                        theta_doppler = Angle_Gamma_Reponse;
                                        Egamma_DS = Energy_Gamma_Reponse;
                                        build_back_data_sim->Fill();
                                        beta_recoil_mean+= 0; compteur +=1;
                                    }
                                }
                            }
                        }
                        time = time + dt;
                        Nt = Ntdt;
                        if (stop == 1) { Nt = 0; }
                    }
                }
            }
        }
    }
    return beta_recoil_mean/double(compteur);
}



double MCsim_Brho_to_ElossEx(double Ebeam[2],  Double_t MReaction[NUCLEI_REACTION], Double_t AZReaction[2][NUCLEI_REACTION], Double_t level_spectro[SPECTRO_PARAMETERS], Double_t target_dim[TARGET_PARAMETERS],  double pol2Eloss[3], double pol2Ex[3], double theta_VAMOS_max,  Double_t IonEnergy[NUCLEI][MAX_LINES], Double_t dEdx_e[NUCLEI][MAX_LINES], Double_t dEdx_n[NUCLEI][MAX_LINES], int Nfile_energy[NUCLEI]){
    int nStates = 20; double Ex[nStates];
    double Eloss[nStates];   double BrhoMeas[nStates];
    double appliedMCsim_kinematics;
    double Brholim[2]={0};
    for (int i = 0; i < nStates; i++) {
        Ex[i]=double(i)*1.5*level_spectro[0]/double(nStates);
        appliedMCsim_kinematics= MCsim_kinematics(Ebeam, Ex[i], MReaction, AZReaction,  target_dim, BrhoMeas[i],  Brholim,  Eloss[i], theta_VAMOS_max, IonEnergy, dEdx_e, dEdx_n, Nfile_energy);

    }
    TGraph* grloss= new TGraph(nStates, BrhoMeas, Eloss);
    TGraph* grEx= new TGraph(nStates, BrhoMeas, Ex);
    TF1* funcLoss = new TF1("funcLoss", pol2_func, 0.5, 1.5, 3);
    TF1* funcEx = new TF1("funcEx", pol2_func, 0.5, 1.5, 3);
    funcLoss->SetParameter(0, 6);    funcLoss->SetParameter(1, -0.0898);  funcLoss->SetParameter(2, 5.78*pow(10,-4));
    funcEx->SetParameter(0, 20);  funcEx->SetParameter(1, -0.0898); funcEx->SetParameter(2, 5.78*pow(10,-4));
    grloss->Fit("funcLoss", "N");
    grEx->Fit("funcEx", "N");
    for(int p=0;p<3;p++){
        pol2Eloss[p]=funcLoss->GetParameter(p);
        pol2Ex[p]=funcEx->GetParameter(p);
        
    }
    return 1.0;
}



int MCsim_noise(TChain* data, TTree* build_back_data_sim, double Ebeam, double theta_VAMOS_max, double pol2Eloss[3], Double_t MReaction[NUCLEI_REACTION], Double_t AZReaction[2][NUCLEI_REACTION], Double_t AGATA_E_range_noise[2], Double_t AGATA_angle[2], Double_t AGATA_position[3],  Double_t AGATA_resolution[2], Double_t level_spectro[SPECTRO_PARAMETERS]) {
    
    double Noise_simu = 50000;
    TRandom* r3 = new TRandom();//uniform in energy, angle distribution derived from exp. noise experimental
    TRandom* rEnergy_Gamma = new TRandom();//gaussian reponse of AGATA
    TRandom* rAngle_Gamma = new TRandom();//gaussian reponse of AGATA
    
    //Experimental noise distribution
    TH1F* ProfileBrho = new TH1F("ProfileBrho", ";B#rho^{VAMOS} (MeV)", 2000, 0, 2);
    TH1F* ProfileBeta = new TH1F("ProfileBeta", ";#beta^{VAMOS} (MeV)", 1000, 0, 1);
    TH2F* GammaNOISE = new TH2F("GammaNOISE", ";E_{#gamma} (keV);#theta (deg)",int((AGATA_E_range_noise[1]-AGATA_E_range_noise[0])/1.5)+1, AGATA_E_range_noise[0]-1, AGATA_E_range_noise[1]+1, int(AGATA_angle[1]-AGATA_angle[0])+5, AGATA_angle[0]-2, AGATA_angle[1]+2);
    int N_gamma_angle=AGATA_angle[1]-AGATA_angle[0] +1;
    double gamma_angle_list[N_gamma_angle], ponderation_gamma_angle[N_gamma_angle];
    
    Double_t Beta_emission, Beta_reaction, DeltaBeta, Egamma_DS, theta_doppler;
    Float_t trackE[31]; Float_t trackX1[31]; Float_t trackY1[31]; Float_t trackZ1[31];
    Float_t Brho, ThetaLdeg, Theta, Phi;
    int nbTrack;
    
    build_back_data_sim->SetBranchAddress("Beta_emission", &Beta_emission);
    build_back_data_sim->SetBranchAddress("Beta_reaction", &Beta_reaction);
    build_back_data_sim->SetBranchAddress("DeltaBeta", &DeltaBeta);
    build_back_data_sim->SetBranchAddress("Egamma_DS", &Egamma_DS);
    build_back_data_sim->SetBranchAddress("theta_doppler", &theta_doppler);
    data->SetBranchAddress("Brho", &Brho);
    data->SetBranchAddress("Phi", &Phi); data->SetBranchAddress("Theta", &Theta); data->SetBranchAddress("ThetaLdeg", &ThetaLdeg);
    data->SetBranchAddress("nbTrack", &nbTrack);
    data->SetBranchAddress("trackE", &trackE);  data->SetBranchAddress("trackX1", &trackX1); data->SetBranchAddress("trackY1", &trackY1); data->SetBranchAddress("trackZ1", &trackZ1);
    int nevent = data->GetEntries();

    double Z, Y, X, R, cA, betaRecoil;
    double  Doppler_angle;          double  angleRecoil[2]; double  angleEjectil[2]; double angleGamma[2];
    int derivation_kinematics;
    
    for (int i = 0; i < nevent; i++) {
        data->GetEntry(i);
        //  ProfileBrho->Fill(Brho);
        for (int nt = 0; nt < nbTrack; nt++) {
            if(trackE[nt]>=AGATA_E_range_noise[0] && trackE[nt]<=AGATA_E_range_noise[1]){
                Z = trackZ1[nt] + AGATA_position[2];
                Y = trackY1[nt] + AGATA_position[1];
                X = trackX1[nt] + AGATA_position[0];
                //g-ray angles
                angleGamma[0] = (180. / TMath::Pi()) * acos(Z / sqrt(X * X + Y * Y + Z * Z));
                angleGamma[1] = atan(-X / Y) * 180. / TMath::Pi();
                if (Y < 0) angleGamma[1]= angleGamma[1] + 180.0;
                if (angleGamma[1] < 0) angleGamma[1] = angleGamma[1] + 360.0;
                // ejectile angles
                angleEjectil[0] =  ThetaLdeg ;
                angleEjectil[1] = atan(sin(Phi * 0.001) / (-cos(Phi * 0.001) * sin(Theta * 0.001))) * 180. / TMath::Pi();
                if (-sin(Theta * 0.001) < 0) { angleEjectil[1]= angleEjectil[1]+ 180.; }
                if (angleEjectil[1] < 0) { angleEjectil[1]= angleEjectil[1] + 360.; }
                if (angleEjectil[0] < theta_VAMOS_max) {
                    derivation_kinematics = VAMOS_kinematics(angleEjectil, angleRecoil, betaRecoil, Brho, pol2Eloss, Ebeam, MReaction, AZReaction);
                    Doppler_angle = lab_to_Doppler_angle(angleRecoil, angleGamma);
                    GammaNOISE->Fill(trackE[nt],Doppler_angle);
                    ProfileBeta->Fill(betaRecoil);
                }
            }
        }
    }
    int looked_for_AGATA_angle_noise_distrobution= extraction_noise(GammaNOISE, gamma_angle_list, ponderation_gamma_angle);
    for (int k = 0; k < Noise_simu; k++) {
        theta_doppler = conversion_angle_random_distribution(r3->Rndm(), gamma_angle_list, ponderation_gamma_angle);
        theta_doppler = rAngle_Gamma->Gaus(theta_doppler, AGATA_resolution[1]);
        Egamma_DS = r3->Uniform(AGATA_E_range_noise[0], AGATA_E_range_noise[1]);
        Egamma_DS = rEnergy_Gamma->Gaus(Egamma_DS, AGATA_resolution[0]);
        R =  Egamma_DS / (level_spectro[1]*1000.);
        cA = cos(theta_doppler * TMath::Pi() / 180.);
        Beta_emission = (cA * R * R + sqrt(1 + R * R * cA * cA - R * R)) / (R * R * cA * cA + 1);
        Beta_reaction=ProfileBeta->GetRandom();
        DeltaBeta= Beta_reaction - Beta_emission;
        build_back_data_sim->Fill();
    }
    return 1;
}
