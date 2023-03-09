#include <iostream>
#include <stdlib.h>
string path = "/Volumes/TOSHIBA EXT/Thèse/Run/";
//string path = "...";

double Pi = TMath::Pi();

void extractionData() {
    TChain* chain = new TChain("TreeMaster");
    string f=path + "target5/*.root";                      chain->Add(f.c_str());
    f=path + "target2/*.root";                      chain->Add(f.c_str());
    f=path + "target3/*.root";                      chain->Add(f.c_str());
    f=path + "target4/*.root";                      chain->Add(f.c_str());
    f=path + "target1/*.root";                      chain->Add(f.c_str());
    std::cout<<chain->GetEntries()<<std::endl;
    
    UShort_t GATCONF ;
    Float_t Brho, ThetaLdeg, Theta, Phi, SIE_C;
    int nbTrack;
    Float_t DL_StrV[16], UL_StrV[16], DR_StrV[16], UR_StrV[16];
    Float_t DL_SecV[16], UL_SecV[16], DR_SecV[16], UR_SecV[16];
    short DL_StrVN[16], UL_StrVN[16], DR_StrVN[16], UR_StrVN[16];
    short DL_SecVN[16], UL_SecVN[16], DR_SecVN[16], UR_SecVN[16];
    Int_t  DL_StrVM, DR_StrVM, UL_StrVM, UR_StrVM, DL_SecVM, DR_SecVM, UL_SecVM, UR_SecVM;
    Float_t trackE[31];    Float_t trackX1[31]; Float_t trackY1[31]; Float_t trackZ1[31];
    Float_t Si_dE, Si_Eres, Si_theta; //Float_t Si_phi;

    chain->SetBranchAddress("GATCONF", &GATCONF); chain->SetBranchAddress("Brho", &Brho);
    chain->SetBranchAddress("Phi", &Phi); chain->SetBranchAddress("Theta", &Theta); chain->SetBranchAddress("ThetaLdeg", &ThetaLdeg);
    chain->SetBranchAddress("trackE", &trackE);   chain->SetBranchAddress("nbTrack", &nbTrack);
    chain->SetBranchAddress("trackX1", &trackX1); chain->SetBranchAddress("trackY1", &trackY1); chain->SetBranchAddress("trackZ1", &trackZ1);
    chain->SetBranchAddress("SIE_C", &SIE_C);
    chain->SetBranchAddress("UL_StrVN", &UL_StrVN);        chain->SetBranchAddress("UR_StrVN", &UR_StrVN);
    chain->SetBranchAddress("DL_StrVN", &DL_StrVN);        chain->SetBranchAddress("DR_StrVN", &DR_StrVN);
    chain->SetBranchAddress("UL_StrVM", &UL_StrVM);        chain->SetBranchAddress("UR_StrVM", &UR_StrVM);
    chain->SetBranchAddress("DL_StrVM", &DL_StrVM);        chain->SetBranchAddress("DR_StrVM", &DR_StrVM);
    chain->SetBranchAddress("UL_SecVN", &UL_SecVN);        chain->SetBranchAddress("UR_SecVN", &UR_SecVN);
    chain->SetBranchAddress("DL_SecVN", &DL_SecVN);        chain->SetBranchAddress("DR_SecVN", &DR_SecVN);
    chain->SetBranchAddress("UL_SecVM", &UL_SecVM);        chain->SetBranchAddress("UR_SecVM", &UR_SecVM);
    chain->SetBranchAddress("DL_SecVM", &DL_SecVM);        chain->SetBranchAddress("DR_SecVM", &DR_SecVM);
    chain->SetBranchAddress("UL_StrV", &UL_StrV);          chain->SetBranchAddress("UR_StrV", &UR_StrV);
    chain->SetBranchAddress("DL_StrV", &DL_StrV);          chain->SetBranchAddress("DR_StrV", &DR_StrV);
    chain->SetBranchAddress("UL_SecV", &UL_SecV);          chain->SetBranchAddress("UR_SecV", &UR_SecV);
    chain->SetBranchAddress("DL_SecV", &DL_SecV);          chain->SetBranchAddress("DR_SecV", &DR_SecV);
    

    TFile* outputFile = new TFile("data.root", "RECREATE");
    TTree* GANILData = new TTree("TreeMaster","TreeMaster");

    GANILData->Branch("nbTrack",&nbTrack,"nbTrack/I");            GANILData->SetBranchAddress("nbTrack", &nbTrack);
    GANILData->Branch("trackE", &trackE,"trackE[nbTrack]/F");          GANILData->SetBranchAddress("trackE", &trackE);
    GANILData->Branch("trackX1", &trackX1,"trackX1[nbTrack]/F");       GANILData->SetBranchAddress("trackX1", &trackX1);
    GANILData->Branch("trackY1", &trackY1,"trackY1[nbTrack]/F");       GANILData->SetBranchAddress("trackY1", &trackY1);
    GANILData->Branch("trackZ1", &trackZ1,"trackZ1[nbTrack]/F");       GANILData->SetBranchAddress("trackZ1", &trackZ1);
    GANILData->Branch("Brho", &Brho,"Brho/F");                                  GANILData->SetBranchAddress("Brho", &Brho);
    GANILData->Branch("Theta", &Theta,"Theta/F");                               GANILData->SetBranchAddress("Theta", &Theta);
    GANILData->Branch("Phi", &Phi,"Phi/F");                                     GANILData->SetBranchAddress("Phi", &Phi);
    GANILData->Branch("ThetaLdeg", &ThetaLdeg,"ThetaLdeg/F");                   GANILData->SetBranchAddress("ThetaLdeg", &ThetaLdeg);
    GANILData->Branch("Si_dE", &Si_dE,"Si_dE/F");                               GANILData->SetBranchAddress("Si_dE", &Si_dE);
    GANILData->Branch("Si_theta", &Si_theta,"Si_theta/F");                      GANILData->SetBranchAddress("Si_theta", &Si_theta);
  //  GANILData->Branch("Si_phi", &Si_phi,"Si_phi/F");                        GANILData->SetBranchAddress("Si_phi", &Si_phi);
    GANILData->Branch("Si_Eres", &Si_Eres,"Si_Eres/F");                         GANILData->SetBranchAddress("Si_Eres", &Si_Eres);


    double piedestal[2] = {0.3,.728};
    double Target_SPIDER =  119.8;//mm

    for(int i = 0;i<chain->GetEntries();i++){
        chain->GetEntry(i);
        if(GATCONF==1  && Brho>0){
            //POSSIBLE ADDITIONAL SELECTION: VTS - TStrack, QPLD_C  ,  T_HF_PLG_caOf
            nbTrack=nbTrack;
            Brho=Brho;
            Theta=Theta;    Phi=Phi;
            ThetaLdeg = ThetaLdeg;
            Si_Eres=SIE_C+piedestal[1];
            for (int nt = 0; nt < nbTrack; nt++) {
                if (TMath::Abs(trackE[nt] - 6875) > 8) {
                    trackX1[nt]=trackX1[nt];
                    trackY1[nt]=trackY1[nt];
                    trackZ1[nt]=trackZ1[nt];
                    trackE[nt]=trackE[nt];
                }
            }
            Si_dE=0;  Si_theta=0;
            if (UR_StrVM == 1 && UR_SecVM == 1 && UL_StrVM + DR_StrVM + DL_StrVM + UL_SecVM + DR_SecVM + DL_SecVM < 1) {
                if (TMath::Abs(UR_SecV[0] - UR_StrV[0]) < 0.1) {
                    Si_dE = UR_StrV[0] - piedestal[0];
                    Si_theta = atan((24. + 0.75 + 24. / 16. * UR_StrVN[0]) / Target_SPIDER) * 180. / Pi;
                }
            }
            if (UL_StrVM == 1 && UL_SecVM == 1 && UR_StrVM + DR_StrVM + DL_StrVM + UR_SecVM + DR_SecVM + DL_SecVM < 1) {
                if (TMath::Abs(UL_SecV[0] - UL_StrV[0]) < 0.1) {
                    Si_dE= UL_StrV[0] - piedestal[0];
                    Si_theta = atan((24. + 0.75 + 24. / 16. * UL_StrVN[0]) / Target_SPIDER) * 180. / Pi;
                }
            }
            if (DR_StrVM == 1 && DR_SecVM == 1 && UL_StrVM + UR_StrVM + DL_StrVM + UL_SecVM + UR_SecVM + DL_SecVM < 1) {
                if (TMath::Abs(DR_SecV[0] - DR_StrV[0]) < 0.1) {
                    Si_dE= DR_StrV[0] - piedestal[0];
                    Si_theta  = atan((24. + 0.75 + 24. / 16. * DR_StrVN[0]) / Target_SPIDER) * 180. / Pi;
                }
            }
            if (DL_StrVM == 1 && DL_SecVM == 1 && UR_StrVM + DR_StrVM + UL_StrVM + UR_SecVM + DR_SecVM + UL_SecVM < 1) {
                if (TMath::Abs(DL_SecV[0] - DL_StrV[0]) < 0.1) {
                    Si_dE = DL_StrV[0] - piedestal[0];
                    Si_theta  = atan((24. + 0.75 + 24. / 16. * DL_StrVN[0]) / Target_SPIDER) * 180. / Pi;
                }
            }
            GANILData->Fill();
        }
    }
    GANILData->Write();
    outputFile->Close();
}
