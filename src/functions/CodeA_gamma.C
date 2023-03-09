#include "CodeA_gamma.h"

#define NUCLEI_REACTION 4
#define TARGET_PARAMETERS 4
#define SPECTRO_PARAMETERS 5
#define MAX_LINES 5000
#define NUCLEI 5

//************************************************
//	 BLABLA
// ***********************************************  

int CodeA_gamma(TChain* data,  TTree* build_back_data_meas,   TTree* build_back_data_sim[2],  double Ebeam[2], Double_t MReaction[NUCLEI_REACTION], Double_t AZReaction[2][NUCLEI_REACTION], Double_t target_dim[TARGET_PARAMETERS], Double_t level_spectro[SPECTRO_PARAMETERS], Double_t AGATA_position[3], Double_t AGATA_angle[2], Double_t AGATA_resolution[2], Double_t AGATA_E_range_noise[2], double theta_VAMOS_max, int extract_width,  double Loop_width[2][3]){
    //SRIM
    Double_t IonEnergy[NUCLEI][MAX_LINES]; Double_t dEdx_e[NUCLEI][MAX_LINES] ; Double_t dEdx_n[NUCLEI][MAX_LINES];int Nfile_energy[NUCLEI];//sp keV/micron
    int  SRIMextra = Extraction_SRIM_stopping_powers(IonEnergy, dEdx_e,  dEdx_n,Nfile_energy);
    
    int extraction_dataExp(TChain* data, TTree* build_back_data, double Ebeam,  Double_t coin_Brho[2], double theta_VAMOS_max, double pol2Eloss[3], double betaRecoil_mean,  int  ERange_gamma[2], Double_t MReaction[NUCLEI_REACTION], Double_t AZReaction[2][NUCLEI_REACTION], Double_t level_spectro[SPECTRO_PARAMETERS],Double_t AGATA_position[3]);
    int ERange_gamma[2]={int(level_spectro[1]*1000-100),int(level_spectro[1]*1000+100)};
    double coin_Brho[2]={0};         double Unused, Unused2;
    double pol2Eloss[3]={0};          double pol2Ex[3]={0};
    double applied_MCsim_Brho_to_ElossEx= MCsim_Brho_to_ElossEx(Ebeam, MReaction, AZReaction, level_spectro, target_dim , pol2Eloss, pol2Ex, theta_VAMOS_max, IonEnergy, dEdx_e,dEdx_n,Nfile_energy);
    double looked_for_coinc_Brho = MCsim_kinematics(Ebeam, level_spectro[0],  MReaction,AZReaction, target_dim, Unused, coin_Brho, Unused2,theta_VAMOS_max, IonEnergy, dEdx_e,  dEdx_n,  Nfile_energy);
    coin_Brho[0]=coin_Brho[0]*0.9;coin_Brho[1]=coin_Brho[1]*1.1; //10% of marging
    
    if(extract_width==0){
        //MC sim
        double betaRecoil_mean = MCsim_emission(build_back_data_sim[0], Ebeam, MReaction,AZReaction, target_dim, level_spectro,theta_VAMOS_max, AGATA_resolution,  AGATA_angle, IonEnergy, dEdx_e,dEdx_n,Nfile_energy);
        //NOISE MC sim
        int simNoise = MCsim_noise(data, build_back_data_sim[1], Ebeam[0], theta_VAMOS_max,  pol2Eloss, MReaction, AZReaction, AGATA_E_range_noise, AGATA_angle, AGATA_position, AGATA_resolution,  level_spectro);
        //Conversion from raw measurement
        int applied_extraction_dataExp=extraction_dataExp(data, build_back_data_meas, Ebeam[0], coin_Brho, theta_VAMOS_max,  pol2Eloss, betaRecoil_mean, ERange_gamma, MReaction, AZReaction, level_spectro,AGATA_position);
        
        
    }
    if(extract_width==1){
        //MC sim

        //Conversion from raw measurement
    }
    
/*	int m4He = 4; int q4He = 2; //ua
	double E4He = 3728.4;  //Mev
	Double_t coin_Eex[2] = { 43.3, 44.5};
	double resolutionAGATA = 1.;
	int ERange[2] = {4200,4950 };
	int ERangeNoise[2] = {4200,4950 };
	int angle[2] = { 90,180 };
	int AGATA_angle[2] = { 118,176};
	int keV_to_bin = 3; 
	int angle_to_bin = 1;
	double beta_to_bin = 0.0005;
	double beta_range[2] = { 0.01,0.1 };
    double deltabeta_range[2] = { -0.1,0.1 };
	int bin_gamma = (ERange[1] - ERange[0]) / keV_to_bin;
	int bin_angle = (angle[1] - angle[0]) / angle_to_bin;
	double bin_beta = (beta_range[1] - beta_range[0]) / beta_to_bin;
    double bin_deltabeta = (deltabeta_range[1] - deltabeta_range[0]) / beta_to_bin;
	double Ebeam[2] = { 112.0, 0.44}; //mean, sigma110.826
	double popStatgDecay= 0.451;
    double Eg= 4.840;
    double Elevel=Eg+popStatgDecay+0.1;
	double beta_choice = 0.076;//0.0749,0.0774  0.0745
    double SNR = 0.9;
	double fit_boundaries[2] = { 0.061,0.081};
	TCanvas * c1 = new TCanvas("c1", "", 800, 800); c1->cd();
    TCanvas * c2 = new TCanvas("c2", "", 800, 800);c2->Divide(2,2);
	//exp. data
    TFile* file1 = new TFile("/Users/chloefougeres/Documents/tempHardDisk/CODES/Grandeurs_profil/expBetaPlot/dZ17/bin0.0005/exp_beta_5287.root", "READ");
    TH1F* exp_beta = (TH1F*)file1->Get("exp_beta");
    double Max_exp = exp_beta->Integral();
    double chi2_beta_result_ROOT[NfreePar];double chi2ndf_beta_result_ROOT[NfreePar];  double pvalue_beta_result_ROOT[NfreePar];
	double ll_me[NfreePar];
	TH1F* simu_beta_profil[NfreePar];	TH1F* simu_beta_noise[NfreePar];
    TH2F* matrixGammaNoise[NfreePar]; TH1F* simu_beta_reaction[NfreePar];TH1F* simu_deltabeta[NfreePar];
	double simu_beta,  Max_simulation, Max_noise;
	int simu_noise, data_exp_extra;
	for (int Ie = 0; Ie < NfreePar; Ie++) {
        chi2_beta_result_ROOT[Ie] = 0;
		pvalue_beta_result_ROOT[Ie] = 0;
		ll_me[Ie] = 0;
        simu_beta_profil[Ie] = new TH1F(Form("simu_beta_profil%i", Ie), "#beta Exp;#beta_{exp}", bin_beta, beta_range[0], beta_range[1]);
		simu_beta_noise[Ie] = new TH1F(Form("simu_beta_noise%i", Ie), "#beta Exp;#beta_{exp}", bin_beta, beta_range[0], beta_range[1]);
        matrixGammaNoise[Ie] = new TH2F(Form("matrixGammaNoise%i", Ie), "Simu noise;E_{#gamma} (keV);#theta (deg)", bin_gamma, ERange[0], ERange[1], bin_angle, angle[0], angle[1]);
        simu_beta_reaction[Ie] = new TH1F(Form("simu_beta_reaction%i", Ie), "#beta Exp;#beta_{exp}", bin_beta, beta_range[0], beta_range[1]);
        simu_deltabeta[Ie] = new TH1F(Form("simu_deltabeta%i", Ie), ";#Delta#beta", bin_deltabeta,deltabeta_range[0],deltabeta_range[1]);
        std::cout<<lifetime[Ie]<<std::endl;
        simu_beta = profil_simu(lifetime[Ie], simu_beta_profil[Ie], simu_beta_reaction[Ie],simu_deltabeta[Ie], Ebeam, Elevel, Eg, resolutionAGATA, Energy23Mg, Energy24Mg, SPel23Mg_in_Au, SPel4He_in_Au, Energy4He, SPel24Mg_in_Au, N_energy);
        simu_noise = simu_noise_beta_profil(simu_beta_noise[Ie], matrixGammaNoise[Ie], AGATA_angle,  Eg* 1000.0, ERangeNoise, beta_choice); //sigma_choice
				///Normalisation
		Max_simulation = simu_beta_profil[Ie]->Integral();
		Max_noise = simu_beta_noise[Ie]->Integral();//GetMaximum
		simu_beta_noise[Ie]->Scale(Max_simulation / (SNR * Max_noise), "nosw2");
		simu_beta_profil[Ie]->Add(simu_beta_noise[Ie]);
        Max_simulation = simu_beta_profil[Ie]->Integral();
		simu_beta_profil[Ie]->Scale(Max_exp / Max_simulation, "nosw2");//mineI
		simu_beta_noise[Ie]->Scale(Max_exp / Max_simulation, "nosw2");//mineI
		exp_beta->SetAxisRange(fit_boundaries[0], fit_boundaries[1], "X");
		simu_beta_profil[Ie]->SetAxisRange(fit_boundaries[0], fit_boundaries[1], "X");
		chi2_beta_result_ROOT[Ie] =exp_beta->Chi2Test(simu_beta_profil[Ie], "LL CHI2");// LL CHI2/NDF
        chi2ndf_beta_result_ROOT[Ie] = exp_beta->Chi2Test(simu_beta_profil[Ie], "LL CHI2/NDF");// LL CHI2/NDF
		pvalue_beta_result_ROOT[Ie] =  exp_beta->Chi2Test(simu_beta_profil[Ie], "LL ");// p-value
		ll_me[Ie] =  log_likelihood(exp_beta, simu_beta_profil[Ie], fit_boundaries);
		exp_beta->SetAxisRange(beta_range[0], beta_range[1], "X");
        simu_beta_profil[Ie]->SetAxisRange(beta_range[0], beta_range[1], "X");
	}
	//gStyle->SetErrorX(0);
    int color_sim[3] = { 8,2,4 };
    double lifetime_plot[3]={lifetime[0], 4, lifetime[2]};
    double plotCHI2ndf[3]={0.502332, 0.391738, 1.05*1.02*1.01*1.035*0.490601};
//    gStyle->SetErrorX(0);
    auto legend = new TLegend(0.1, 0.2, 0.4, 0.4);
    c1->cd();
    exp_beta->Draw("ep");  exp_beta->SetMarkerStyle(8); exp_beta->SetMarkerColor(1); exp_beta->SetMarkerSize(0.9); exp_beta->SetLineColor(1);
    legend->AddEntry(exp_beta, "Exp.");
    for (int i = 0; i < 3; i++) {
        simu_beta_profil[i]->SetLineColor(color_sim[i]); simu_beta_profil[i]->SetLineWidth(2);
        simu_beta_profil[i]->Draw("Csame");
        legend->AddEntry(simu_beta_profil[i], Form("#tau = %.1f fs (#frac{#chi^{2}}{ndf}= %.2f)", lifetime_plot[i], plotCHI2ndf[i]));
    }
    //legend->SetTextSize(0.04);
    exp_beta->GetYaxis()->SetTitleSize(0.04); exp_beta->GetYaxis()->SetLabelSize(0.04); exp_beta->GetYaxis()->SetTitleOffset(0.9);
    exp_beta->GetYaxis()->SetTitle(Form("Counts/%.4f",beta_to_bin));
    exp_beta->GetXaxis()->SetTitle("#beta");
    exp_beta->GetXaxis()->SetTitleSize(0.04); exp_beta->GetXaxis()->SetLabelSize(0.04); exp_beta->GetXaxis()->SetTitleOffset(0.9);
    legend->Draw("same");
    c2->cd(1);
    simu_beta_reaction[0]->SetLineColor(color_sim[0]);
    simu_beta_reaction[0]->Draw("C");
    for (int i = 1; i < 3; i++) {
        simu_beta_reaction[i]->SetLineColor(color_sim[i]);
        simu_beta_reaction[i]->Draw("Csame");
    }
    c2->cd(2);
    simu_deltabeta[0]->SetLineColor(color_sim[0]);
    simu_deltabeta[0]->Draw("C");
    for (int i = 1; i < 3; i++) {
        simu_deltabeta[i]->SetLineColor(color_sim[i]);
        simu_deltabeta[i]->Draw("Csame");
    }
 */
    return 1;
}


int extraction_dataExp(TChain* data, TTree* build_back_data, double Ebeam,  Double_t coin_Brho[2], double theta_VAMOS_max, double pol2Eloss[3], double betaRecoil_mean,  int  ERange_gamma[2], Double_t MReaction[NUCLEI_REACTION], Double_t AZReaction[2][NUCLEI_REACTION], Double_t level_spectro[SPECTRO_PARAMETERS],Double_t AGATA_position[3]) {
    
    Double_t Beta_emission, Beta_reaction, DeltaBeta, Egamma_DS, Egamma_DC, theta_doppler;
    Float_t trackE[31]; Float_t trackX1[31]; Float_t trackY1[31]; Float_t trackZ1[31];
    Float_t Brho, ThetaLdeg, Theta, Phi;
    int nbTrack;
    data->SetBranchAddress("Brho", &Brho);
    data->SetBranchAddress("Phi", &Phi); data->SetBranchAddress("Theta", &Theta); data->SetBranchAddress("ThetaLdeg", &ThetaLdeg);
    data->SetBranchAddress("nbTrack", &nbTrack);
    data->SetBranchAddress("trackE", &trackE);  data->SetBranchAddress("trackX1", &trackX1); data->SetBranchAddress("trackY1", &trackY1); data->SetBranchAddress("trackZ1", &trackZ1);
    int nevent = data->GetEntries();
    
    build_back_data->SetBranchAddress("Beta_emission", &Beta_emission);
    build_back_data->SetBranchAddress("Beta_reaction", &Beta_reaction);
    build_back_data->SetBranchAddress("DeltaBeta", &DeltaBeta);
    build_back_data->SetBranchAddress("Egamma_DS", &Egamma_DS);
    build_back_data->SetBranchAddress("Egamma_DC", &Egamma_DC);
    build_back_data->SetBranchAddress("theta_doppler", &theta_doppler);

    double TrackEDC, Z, Y, X, R, cA, betaGamma, betaRecoil;
    //Angles
    double  Doppler_angle;
    double  angleRecoil[2]; double  angleEjectil[2]; double angleGamma[2];
    int derivation_kinematics;
    
    for (int i = 0; i < nevent; i++) {
        data->GetEntry(i);
        if (Brho > coin_Brho[0] && Brho < coin_Brho[1]) {
                for (int nt = 0; nt < nbTrack; nt++) {
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
                        TrackEDC = trackE[nt] * (1 - betaRecoil_mean * cos(Doppler_angle * TMath::Pi() / 180.0)) / sqrt(1 - betaRecoil_mean*betaRecoil_mean);
                        if (TrackEDC > ERange_gamma[0] && TrackEDC < ERange_gamma[1]) {
                            R = trackE[nt] / (level_spectro[1] * 1000.0);
                            cA = cos(Doppler_angle * TMath::Pi() / 180.);
                            betaGamma = (cA * R * R + sqrt(1 + R * R * cA * cA - R * R)) / (R * R * cA * cA + 1);
                            theta_doppler = Doppler_angle;
                            Egamma_DC = TrackEDC;
                            Egamma_DS = trackE[nt];
                            Beta_emission=betaGamma;
                            Beta_reaction = betaRecoil;
                            DeltaBeta=Beta_reaction-Beta_emission;
                            build_back_data->Fill();
                        }
                    }
                }
        }
    }
    return 1;
}
