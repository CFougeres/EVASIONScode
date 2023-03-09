//********************************************
//
// ROOT macro to run the EVASIONS code
//
//Author: C. Fougères (2022)
//any remark(s) to be sent at cfougeres@anl.gov
//********************************************

#include "functions/main.h"
#include "functions/CodeA_gamma.h"
#include "functions/CodeB_particle.h"
#define NUCLEI_REACTION 4
#define NUCLEI_DECAY 2
#define TARGET_PARAMETERS 4
#define SPECTRO_PARAMETERS 5

using namespace std;

double c = 3.0 * pow(10, 8); // speed velocity m.s-1
TH1F * difference_beta_profil;
TH2F * unbound_particle_recoil_coincEx;
TCanvas* C_beta;
TCanvas* C_unbound_particle;

int main(int argc, char* argv[]){
    printf("=====================================\n");
    printf("===          EVASION              ===\n");
    printf("=====================================\n");
    if (argc < 5)    {
          printf("Incorrect number of arguments (%s):\n",argv[0]);
          printf("Applied CodeA for gamma emission? 1=yes (else =0) \n");
          printf("Applied CodeB for particle emission? 1=yes (else =0) \n");
          printf("Width determination? 1=yes (else =0)\n");
          printf("Saving Trees? 1=yes (else =0)\n");
          printf("Plotting? 1=yes (else =0)\n");
      return 1;
    }
    if (argc > 4)    {
        int codeA = atoi(argv[1]); int codeB =  atoi(argv[2]); int extract_width = atoi(argv[3]);
        int saving = atoi(argv[4]); int plotting =  atoi(argv[5]);
        
        TApplication theApp("App",&argc, argv);
        
        printf("EVASION processed \n");
        
        string path = gSystem->pwd();    cout<<path<<endl;
        
        //USER INPUTS
        int extraction_inputs(string file, double Ebeam[2], Double_t MReaction[NUCLEI_REACTION], Double_t AZReaction[2][NUCLEI_REACTION], Double_t target_dim[TARGET_PARAMETERS], Double_t level_spectro[SPECTRO_PARAMETERS], Double_t  MDecay[NUCLEI_DECAY], Double_t decay_spectro[SPECTRO_PARAMETERS], Double_t AGATA_position[3], Double_t AGATA_angle[2], Double_t AGATA_resolution[2],Double_t AGATA_E_range_noise[2], double theta_VAMOS_max, Double_t SPIDER_geometry[3], Double_t SPIDER_resolution[3], double Loop_width[2][3]);
        string path_input = path+"/inputs.dat";     cout<<path_input<<endl;
        string path_save = path+"/tree_results/";   cout<<path_save<<endl;
        
        //DATA TO BE ANALYZED
        TChain* chain = new TChain("TreeMaster");
        string data_path =path+"/data/data.root"; cout<<data_path<<endl;
        chain->Add(data_path.c_str());
        printf("Amount of experimental entries %i \n", chain->GetEntries());
        
        double Ebeam[2]={0};
        Double_t MReaction[4]={0};        Double_t AZReaction[2][4];
        Double_t target_dim[4]={0};
        Double_t level_spectro[5]={0};
        Double_t  MDecay[2]={0};
        Double_t decay_spectro[5]={0};
        Double_t AGATA_position[3]={0};         Double_t AGATA_resolution[2]={0};           Double_t AGATA_angle[2]={0};    Double_t AGATA_E_range_noise[2]={0};
        Double_t SPIDER_geometry[3]={0};        Double_t SPIDER_resolution[3]={0};
        double theta_VAMOS_max=0;
        double Loop_width[2][3];
        int applied_extraction_inputs=extraction_inputs(path_input, Ebeam, MReaction, AZReaction,target_dim,level_spectro, MDecay, decay_spectro, AGATA_position, AGATA_angle, AGATA_resolution,AGATA_E_range_noise, theta_VAMOS_max, SPIDER_geometry, SPIDER_resolution, Loop_width);
        //  if(SPIDER_geometry[0]==0){theta_ejectil_lim = theta_VAMOS_max;}
        if(SPIDER_geometry[0]>0){theta_VAMOS_max = atan(SPIDER_geometry[1]/SPIDER_geometry[0])*(180./TMath::Pi());}
        
        
        //DATA CONVERTED && SIMULATED
        double Beta_emission, Beta_reaction, DeltaBeta, Egamma_DS, Egamma_DC, theta_doppler;
        TTree* build_back_data_meas = new TTree("tree", "tree");
        TTree* build_back_data_sim[2];  build_back_data_sim[0]  = new TTree("tree", "tree");        build_back_data_sim[1]  = new TTree("tree", "tree");
        
        //PLOTS DEFINITIONS
        double bining_Ex=0.05;  double Ex_range[2]={0, level_spectro[0]*1.5}; int NbinsEx = (Ex_range[1]-Ex_range[0])/bining_Ex;
        double bining_beta=2.*pow(10,-4);  double beta_range[2]={-1, 1}; int NbinsBeta = (beta_range[1]- beta_range[0])/bining_beta;
        if(codeA==1){
            C_beta= new TCanvas("C_beta","C_beta",400,400);
            difference_beta_profil = new TH1F(" difference_beta_profil",Form(";#Delta#beta=#beta_{reaction}-#beta_{emission};Counts/%.5f",bining_beta), NbinsBeta,  beta_range[0],beta_range[1]);
        }
        if(codeB==1){
            C_unbound_particle= new TCanvas("C_unbound_particle","C_unbound_particle",400,400);
            unbound_particle_recoil_coincEx = new TH2F("unbound_particle_recoil_coincEx",";Ex^{VAMOS} (MeV); Ex^{SPIDER} (MeV)", NbinsEx,  Ex_range[0],Ex_range[1],NbinsEx,  Ex_range[0],Ex_range[1]);
        }
            
        int applied_sim_analysis;
        if(codeA==1){
            //TREE EXPERIMENT
            build_back_data_meas->Branch("Beta_emission", &Beta_emission,"Beta_emission/D");        build_back_data_meas->SetBranchAddress("Beta_emission", &Beta_emission);
            build_back_data_meas->Branch("Beta_reaction", &Beta_reaction,"Beta_reaction/D");        build_back_data_meas->SetBranchAddress("Beta_reaction", &Beta_reaction);
            build_back_data_meas->Branch("DeltaBeta", &DeltaBeta,"DeltaBeta/D");                    build_back_data_meas->SetBranchAddress("DeltaBeta", &DeltaBeta);
            build_back_data_meas->Branch("Egamma_DS", &Egamma_DS,"Egamma_DS/D");                    build_back_data_meas->SetBranchAddress("Egamma_DS", &Egamma_DS);
            build_back_data_meas->Branch("Egamma_DC", &Egamma_DC,"Egamma_DC/D");                    build_back_data_meas->SetBranchAddress("Egamma_DC", &Egamma_DC);
            build_back_data_meas->Branch("theta_doppler", &theta_doppler,"theta_doppler/D");        build_back_data_meas->SetBranchAddress("theta_doppler", &theta_doppler);
            //TREE MC SIMULATION
            for(int i=0;i<2;i++){
                build_back_data_sim[i]->Branch("Beta_emission", &Beta_emission,"Beta_emission/D");        build_back_data_sim[i]->SetBranchAddress("Beta_emission", &Beta_emission);
                build_back_data_sim[i]->Branch("Beta_reaction", &Beta_reaction,"Beta_reaction/D");        build_back_data_sim[i]->SetBranchAddress("Beta_reaction", &Beta_reaction);
                build_back_data_sim[i]->Branch("DeltaBeta", &DeltaBeta,"DeltaBeta/D");                    build_back_data_sim[i]->SetBranchAddress("DeltaBeta", &DeltaBeta);
                build_back_data_sim[i]->Branch("Egamma_DS", &Egamma_DS,"Egamma_DS/D");                    build_back_data_sim[i]->SetBranchAddress("Egamma_DS", &Egamma_DS);
                build_back_data_sim[i]->Branch("theta_doppler", &theta_doppler,"theta_doppler/D");        build_back_data_sim[i]->SetBranchAddress("theta_doppler", &theta_doppler);
            }
            if(extract_width==0){
                applied_sim_analysis= CodeA_gamma(chain, build_back_data_meas, build_back_data_sim,Ebeam,MReaction, AZReaction,target_dim, level_spectro, AGATA_position,AGATA_angle, AGATA_resolution,AGATA_E_range_noise, theta_VAMOS_max, extract_width,Loop_width);
                if(saving==1){
                    gSystem->cd(path_save.c_str());
                    TFile* outputFileExp = new TFile("dataExp.root", "RECREATE");
                    build_back_data_meas->Write();      outputFileExp->Close();
                    TFile* outputFileMC = new TFile("dataMC_events.root", "RECREATE");
                    build_back_data_sim[0]->Write();      outputFileMC->Close();
                    TFile* outputFileMCnoise = new TFile("dataMC_noise.root", "RECREATE");
                    build_back_data_sim[1]->Write();      outputFileMCnoise->Close();
                }
                if(plotting==1){
                    printf("Enter Ctrl+C to stop ...");
                    C_beta->cd();
                    difference_beta_profil->Draw("C");difference_beta_profil->GetXaxis()->CenterTitle();difference_beta_profil->GetYaxis()->CenterTitle();
                    theApp.Run();
                }
              
            }
            if(extract_width==1){
                //LOOP
                if(plotting==1){
                    printf("Enter Ctrl+C to stop ...");
                }
            }
        }
        
        if(codeB==1){
            applied_sim_analysis=CodeB_particle();
            if(saving==1){
                gSystem->cd(path_save.c_str());
                TFile* outputFileExp = new TFile("dataExp.root", "RECREATE");
                build_back_data_meas->Write();      outputFileExp->Close();
                TFile* outputFileMC = new TFile("dataMC_events.root", "RECREATE");
                build_back_data_sim[0]->Write();      outputFileMC->Close();
                TFile* outputFileMCnoise = new TFile("dataMC_noise.root", "RECREATE");
                build_back_data_sim[1]->Write();      outputFileMCnoise->Close();
            }
            if(plotting==1){
                printf("Enter Ctrl+C to stop ...");
                C_unbound_particle->cd();
                gStyle->SetPalette(kThermometer);
                unbound_particle_recoil_coincEx->Draw("colz");unbound_particle_recoil_coincEx->GetXaxis()->CenterTitle();unbound_particle_recoil_coincEx->GetYaxis()->CenterTitle();
                theApp.Run();
            }
        }
        printf(" EVASIONS finished \n");
        return 1;
    }
}

//********************************************
// Extraction function of user inputs in
// inputs#.dat file
//********************************************
int extraction_inputs(string file, double Ebeam[2], Double_t MReaction[NUCLEI_REACTION], Double_t AZReaction[2][NUCLEI_REACTION], Double_t target_dim[TARGET_PARAMETERS], Double_t level_spectro[SPECTRO_PARAMETERS], Double_t  MDecay[NUCLEI_DECAY], Double_t decay_spectro[SPECTRO_PARAMETERS], Double_t AGATA_position[3], Double_t AGATA_angle[2], Double_t AGATA_resolution[2], Double_t AGATA_E_range_noise[2], double theta_VAMOS_max, Double_t SPIDER_geometry[3], Double_t SPIDER_resolution[3],double Loop_width[2][3]) {
    double c = 3.0 * pow(10, 8); // speed velocity m.s-1
    double mu = 931.5; //MeV
    
    char variable[500][100]={0};
    Double_t param[500]={0}; Double_t index[500]={0};
    Double_t param2[500]={0};     Double_t param3[500]={0};
    ifstream in;
    in.open(file.c_str());
    Int_t nlines = 0;
    while (1) {
        if(nlines!=4 && nlines!=2 && nlines!=29 && nlines!=0 && nlines!=33 && nlines!=36){
            in >> index[nlines] >> variable[nlines] >>param[nlines];
            cout<<index[nlines]<<" "<<variable[nlines] <<" "<<param[nlines] <<endl;
        }
        if(nlines==4 || nlines ==2 || nlines ==29){
            in >> index[nlines] >> variable[nlines] >> param[nlines] >> param2[nlines] >>param3[nlines];
            cout<<index[nlines]<<" "<<variable[nlines] <<" "<<param[nlines] <<" "<<param2[nlines] <<" "<<param3[nlines] <<endl;
        }
        if(nlines==0 || nlines ==33  || nlines ==36){
            in >> index[nlines] >> variable[nlines] >> param[nlines] >> param2[nlines];
            cout<<index[nlines]<<" "<<variable[nlines] <<" "<<param[nlines]<<" "<<param2[nlines] <<endl;
        }
        if (!in.good()) break;
        nlines++;
    }
    in.close();
    
    Ebeam[0]=param[0];   Ebeam[1]=param2[0];
    
    MReaction[0] = param[10] *mu+param[12]; //beam
    MReaction[1] = param[13] *mu+param[15]; //reactif
    MReaction[2] = param[16] *mu+param[18]; //recoil
    MReaction[3] = param[19] *mu+param[21]; //ejectil
    
    AZReaction[0][0] = param[10] ; //beam
    AZReaction[0][1] = param[13] ; //reactif
    AZReaction[0][2] = param[16] ; //recoil
    AZReaction[0][3] = param[19] ; //ejectil
    AZReaction[1][0] = param[11] ; //beam
    AZReaction[1][1] = param[14] ; //reactif
    AZReaction[1][2] = param[17] ; //recoil
    AZReaction[1][3] = param[20] ; //ejectil
    
    //recoil state populated
    level_spectro[0] = param[1]; //Ex
    level_spectro[1] = param[2]; //Eg
    level_spectro[2] = param[3]; //particle threshold
    level_spectro[3] = param[4]; //lifetime
    level_spectro[4] = param[6]; //spin
    Loop_width[0][0]= param[4];   Loop_width[0][1] = param2[4];  Loop_width[0][2] = param3[4];
    Loop_width[1][0]= param[2];    Loop_width[1][1] = param2[2];  Loop_width[1][2] = param3[2];

    //particle decay
    decay_spectro[0] = level_spectro[0]- level_spectro[2]; //particle c.o.m E
    decay_spectro[1] = param[7]; //particle spin
    decay_spectro[2] = param[5]; //daughter Ex
    decay_spectro[3] = param[8]; //daughter spin
    decay_spectro[4] = param[9]; //Particle Momentum Transfered
    MDecay[0] = param[22] *mu+param[24]; //particle
    MDecay[1] = param[25] *mu+param[27]; //daughter

    //Target structure
    target_dim[0]= param[28];//active target thickness (micron)
    target_dim[1]= param[29];//total target thickness (micron)
    target_dim[2]= param2[29];//foil thickness (micron)
    target_dim[3]= param3[29];//distance target-foil (micron)
    
    //AGATA Position
    AGATA_position[0]=param[30];//shift along x-axis ww.r.t nominal position
    AGATA_position[1]=param[31];//shift along x-axis ww.r.t nominal position
    AGATA_position[2]=param[32];//shift along Z-axis ww.r.t nominal position
    AGATA_angle[0]=param[33]; //theta lab mininal
    AGATA_angle[1]=param2[33];//theta lab maximal
    //AGATA Resolution
    AGATA_resolution[0]= param[34];//sigma energy resolution (kev)
    AGATA_resolution[1]= param[35];//sigma angle resolution (deg)
    //NOISE energy range AGATA
    AGATA_E_range_noise[0] =param[36]; // mininal
    AGATA_E_range_noise[1] =param2[36]; //maximal

    //VAMOS geometry
    theta_VAMOS_max = param[37];

    //SPIDER geometry
    SPIDER_geometry[0]=param[38];//distance target->SPIDER along beam axis
    SPIDER_geometry[1]=param[39];//SPIDER_inner_Radius hole
    SPIDER_geometry[2]=param[40];//SPIDER_outer_Radius
    //SPIDER Resolution
    SPIDER_resolution[0]= param[41];//sigma Si_Eres energy resolution (kev)
    SPIDER_resolution[1]= param[42];//sigma Si_dE energy resolution (kev)
    SPIDER_resolution[2]= param[43];//sigma theta resolution (kev)
    
   return 1;
}
