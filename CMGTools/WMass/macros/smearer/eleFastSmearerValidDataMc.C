#define eleFastSmearerValidDataMc_cxx
#include "eleFastSmearerValidDataMc.h"
#include <TH2.h>
#include <TStyle.h>
#include <iostream>
#include <TLorentzVector.h>

using namespace std;

TVector3* CSFrame(TLorentzVector dilepton) {

  TVector3* myCSFrameP;
  TVector3  myCSFrame[3];
  myCSFrameP = myCSFrame;

  float pMass = 0.938272;
  float beamE = 4000;
  float hypot = sqrt(pMass*pMass+beamE*beamE);
  float sign = 1.;
  if (dilepton.Z()<0) sign = -1.;
  TLorentzVector proton1 = TLorentzVector(0.,0.,sign*beamE,hypot);  
  TLorentzVector proton2 = TLorentzVector(0.,0.,-sign*beamE,hypot);
  proton1.Boost(-dilepton.BoostVector()); 
  proton2.Boost(-dilepton.BoostVector());

  TVector3 CSAxis = (proton1.Vect().Unit()-proton2.Vect().Unit()).Unit();
  TVector3 yAxis  = (proton1.Vect().Unit()).Cross((proton2.Vect().Unit()));
  yAxis = yAxis.Unit();
  TVector3 xAxis = yAxis.Cross(CSAxis);
  xAxis = xAxis.Unit();
  myCSFrame[0] = xAxis;
  myCSFrame[1] = yAxis;
  myCSFrame[2] = CSAxis;

  return myCSFrameP;
}

float cosThetaCS(TLorentzVector lplus, TLorentzVector lminus) {

  TLorentzVector dilep = lplus + lminus;
  TLorentzVector boostedLep = lminus;
  boostedLep.Boost(-dilep.BoostVector());

  TVector3 *csframe = CSFrame(dilep);
  return cos(boostedLep.Angle(csframe[2]));
}
    
float phiCS(TLorentzVector lplus, TLorentzVector lminus) {

  TLorentzVector dilep = lplus + lminus;
  TLorentzVector boostedLep = lminus;
  boostedLep.Boost(-dilep.BoostVector());

  TVector3 *csframe = CSFrame(dilep);
  float phi = atan2((boostedLep.Vect()*csframe[1]),(boostedLep.Vect()*csframe[0]));
  if(phi<0)
    return phi + 3.1415926536;
  else 
    return phi;
}

void eleFastSmearerValidDataMc::Loop() {

  bool isMC = true;

  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  // Data/MC SF file
  TFile fileSF("electrons_scale_factors.root");
  TH2F *histoSF = (TH2F*)fileSF.Get("electronsDATAMCratio_FO_ID_ISO");

  // xsec weight (chiara: only for DY!)
  float dyXsec = 3503.7;
  float unskimmedN = 30429578.;  // from skimAnalyzerCount/SkimReport.txt 
  float xsecWeight = dyXsec/unskimmedN;

  // Create histograms
  TH1F *Hr_invMass  = new TH1F("Hr_invMass",  "Hr_invMass",  20, 80.,100.);
  TH1F *Hr_pT1      = new TH1F("Hr_pT1",      "Hr_pT1",      55, 25.,80.);
  TH1F *Hr_pT2      = new TH1F("Hr_pT2",      "Hr_pT2",      50, 20.,70.);
  TH1F *Hr_eta1     = new TH1F("Hr_eta1",     "Hr_eta1",     50, -2.5,2.5);
  TH1F *Hr_eta2     = new TH1F("Hr_eta2",     "Hr_eta2",     50, -2.5,2.5);
  TH1F *Hr_pTZ      = new TH1F("Hr_pTZ",      "Hr_pTZ",      80, 0.,80.);
  TH1F *Hr_pTZzoom  = new TH1F("Hr_pTZzoom",  "Hr_pTZzoom",  30, 0.,30.);
  TH1F *Hr_pTZzoom1 = new TH1F("Hr_pTZzoom1", "Hr_pTZzoom1", 20, 0.,20.);
  TH1F *Hr_pTZzoom2 = new TH1F("Hr_pTZzoom2", "Hr_pTZzoom2", 10, 0.,10.);
  TH1F *Hr_phiStarZ      = new TH1F("Hr_phiStarZ",      "Hr_phiStarZ",      30, 0, 3.14);
  TH1F *Hr_cosThetaStarZ = new TH1F("Hr_cosThetaStarZ", "Hr_cosThetaStarZ", 30, -1., 1.);
  TH1F *Hr_yZ      = new TH1F("Hr_yZ",      "Hr_yZ",      50, -3.,3.);

  TH1F *HrEB_pT1   = new TH1F("HrEB_pT1",   "HrEB_pT1",     55, 25.,80.);
  TH1F *HrEB_pT2   = new TH1F("HrEB_pT2",   "HrEB_pT2",     50, 20.,70.);
  TH1F *HrEE_pT1   = new TH1F("HrEE_pT1",   "HrEE_pT1",     55, 25.,80.);
  TH1F *HrEE_pT2   = new TH1F("HrEE_pT2",   "HrEE_pT2",     50, 20.,70.);

  TH1F *HrEBEB_invMass    = new TH1F("HrEBEB_invMass",    "HrEBEB_invMass",    20,80.,100.);
  TH1F *HrNotEBEB_invMass = new TH1F("HrNotEBEB_invMass", "HrNotEBEB_invMass", 20,80.,100.);

  TH1F *HrEBZ_invMass       = new TH1F("HrEBZ_invMass",       "HrEBZ_invMass",       20,80.,100.);
  TH1F *HrEBZ_pTZ           = new TH1F("HrEBZ_pTZ",           "HrEBZ_pTZ",           80, 0.,80.);
  TH1F *HrEBZ_pTZzoom       = new TH1F("HrEBZ_pTZzoom",       "HrEBZ_pTZzoom",       30, 0.,30.);
  TH1F *HrEBZ_pTZzoom1      = new TH1F("HrEBZ_pTZzoom1",      "HrEBZ_pTZzoom1",      20, 0.,20.);
  TH1F *HrEBZ_pTZzoom2      = new TH1F("HrEBZ_pTZzoom2",      "HrEBZ_pTZzoom2",      10, 0.,10.);
  TH1F *HrEBZ_phiStarZ      = new TH1F("HrEBZ_phiStarZ",      "HrEBZ_phiStarZ",      30, 0, 3.14);
  TH1F *HrEBZ_cosThetaStarZ = new TH1F("HrEBZ_cosThetaStarZ", "HrEBZ_cosThetaStarZ", 30, -1., 1.);

  TH1F *HrEEZ_invMass       = new TH1F("HrEEZ_invMass",       "HrEEZ_invMass",       20,80.,100.);
  TH1F *HrEEZ_pTZ           = new TH1F("HrEEZ_pTZ",           "HrEEZ_pTZ",           80, 0.,80.);
  TH1F *HrEEZ_pTZzoom       = new TH1F("HrEEZ_pTZzoom",       "HrEEZ_pTZzoom",       30, 0.,30.);
  TH1F *HrEEZ_pTZzoom1      = new TH1F("HrEEZ_pTZzoom1",      "HrEEZ_pTZzoom1",      20, 0.,20.);
  TH1F *HrEEZ_pTZzoom2      = new TH1F("HrEEZ_pTZzoom2",      "HrEEZ_pTZzoom2",      10, 0.,10.);
  TH1F *HrEEZ_phiStarZ      = new TH1F("HrEEZ_phiStarZ",      "HrEEZ_phiStarZ",      30, 0, 3.14);
  TH1F *HrEEZ_cosThetaStarZ = new TH1F("HrEEZ_cosThetaStarZ", "HrEEZ_cosThetaStarZ", 30, -1., 1.);


  // Event loop 
  cout << "Entries: " << nentries << endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if (ientry % 10000 == 0) cout << "Event " << ientry << endl;    
    if (isMC && (ientry % 2 != 0)) continue;     

    // event weight
    double myWeight = 1.;
    if (isMC) myWeight = puWeight*xsecWeight;

    // HLT bit
    // if (!HLT_DoubleEl) continue;

    // Reco electrons
    int recoFirst = -99;
    int recoSec   = -99;

    // At least 2 leptons
    if (nLepGood<2) continue;

    for (int ilep=0; ilep<nLepGood; ilep++) {

      // electrons only
      if (abs(LepGood_pdgId[ilep])!=11) continue; 

      // trigger match
      // if (LepGood_trgMatch[ilep]<0.5)  continue;

      // eta and pT cut (minimal, probably looser than applied in the analysis)  
      if (fabs(LepGood_eta[ilep])>2.5) continue;
      if (LepGood_pt[ilep]<15) continue;

      // barrel-endcap transition  
      if ( fabs(LepGood_scEta[ilep])>1.442 && fabs(LepGood_scEta[ilep])<1.566 ) continue;  

      // trigger emulation        
      // if (LepGood_eleMVAPreselId[ilep]==0) continue;
      
      // ID + isolation + conv.rejection
      if (LepGood_eleMVAId[ilep]<2)      continue;   
      if (LepGood_relIso04[ilep]>0.15)   continue;  
      if (LepGood_convVetoFull[ilep]==0) continue;   // LepGood_lostHits==0 + LepGood_conv
      
      if (recoFirst<-90)    recoFirst = ilep; 
      else if (recoSec<-90) recoSec = ilep;
      else continue;
    }

    // At least 2 reco electrons...
    if (recoFirst<-90 || recoSec<-90) continue; 


    // ... with mee compatible with mZ
    TLorentzVector recoFirstP(0,0,0,0);     
    TLorentzVector recoSecP(0,0,0,0);     
    recoFirstP.SetPtEtaPhiM(LepGood_pt[recoFirst],LepGood_eta[recoFirst],LepGood_phi[recoFirst],0.);
    recoSecP.SetPtEtaPhiM(LepGood_pt[recoSec],LepGood_eta[recoSec],LepGood_phi[recoSec],0.);
    float invMassReco = (recoFirstP+recoSecP).M();
    if (invMassReco>100 || invMassReco<80) continue;

    // Pt thresholds (can be applied here since I select anyway the two highest pT ones)
    float recoFirstPt  = recoFirstP.Pt();
    float recoSecPt    = recoSecP.Pt();
    if (recoFirstPt<recoSecPt) cout << "problem with ordering" << endl;
    if (recoFirstPt<25) continue;
    if (recoSecPt<20)   continue;

    // Eta
    float recoFirstEta = recoFirstP.Eta();
    float recoSecEta   = recoSecP.Eta();

    // data/MC eff SF
    float sfWFirst = 1.;
    float sfWSec   = 1.;

    if (isMC) {

      int binXSfFirst = 0.;
      int binYSfFirst = 0.;
      int binXSfSec   = 0.;
      int binYSfSec   = 0.;

      if (fabs(recoFirstEta)<=0.8) binXSfFirst = 1;
      else if (fabs(recoFirstEta)<=1.4442 && fabs(recoFirstEta)>0.8)   binXSfFirst = 2;
      else if (fabs(recoFirstEta)<=1.556 && fabs(recoFirstEta)>1.4442) binXSfFirst = 3;  
      else if (fabs(recoFirstEta)<=2. && fabs(recoFirstEta)>1.556)     binXSfFirst = 4;  
      else if (fabs(recoFirstEta)<=2.5 && fabs(recoFirstEta)>2.)       binXSfFirst = 5;    

      if (fabs(recoSecEta)<=0.8) binXSfSec = 1;
      else if (fabs(recoSecEta)<=1.4442 && fabs(recoSecEta)>0.8)   binXSfSec = 2;
      else if (fabs(recoSecEta)<=1.556 && fabs(recoSecEta)>1.4442) binXSfSec = 3;  
      else if (fabs(recoSecEta)<=2. && fabs(recoSecEta)>1.556)     binXSfSec = 4;  
      else if (fabs(recoSecEta)<=2.5 && fabs(recoSecEta)>2.)       binXSfSec = 5;    

      if (recoFirstPt<=15 && recoFirstPt>10) binYSfFirst = 1;
      else if (recoFirstPt<=20 && recoFirstPt>15)  binYSfFirst = 2;
      else if (recoFirstPt<=30 && recoFirstPt>20)  binYSfFirst = 3;
      else if (recoFirstPt<=40 && recoFirstPt>30)  binYSfFirst = 4;
      else if (recoFirstPt<=50 && recoFirstPt>40)  binYSfFirst = 5;
      else if (recoFirstPt>50) binYSfFirst = 6;

      if (recoSecPt<=15 && recoSecPt>10) binYSfSec = 1;
      else if (recoSecPt<=20 && recoSecPt>15)  binYSfSec = 2;
      else if (recoSecPt<=30 && recoSecPt>20)  binYSfSec = 3;
      else if (recoSecPt<=40 && recoSecPt>30)  binYSfSec = 4;
      else if (recoSecPt<=50 && recoSecPt>40)  binYSfSec = 5;
      else if (recoSecPt>50) binYSfSec = 6;

      if (binXSfFirst>0 && binYSfFirst>0 && binYSfFirst>0 && binYSfSec>0) {
	sfWFirst = histoSF->GetBinContent(binXSfFirst,binYSfFirst);
	sfWSec   = histoSF->GetBinContent(binXSfSec,binYSfSec);
      } else {
	cout << recoFirstPt << " " << recoSecPt << " " << recoFirstEta << " " << recoSecEta << endl;
	cout << "problem!" << endl;
      }
    }

    // Fill histos
    Hr_pT1  -> Fill(recoFirstPt, myWeight*sfWFirst*sfWSec);   
    Hr_pT2  -> Fill(recoSecPt,   myWeight*sfWSec*sfWFirst); 
    Hr_eta1 -> Fill(recoFirstEta,myWeight*sfWFirst*sfWSec); 
    Hr_eta2 -> Fill(recoSecEta,  myWeight*sfWSec*sfWFirst); 

    if (fabs(recoFirstEta)<1.5) 
      HrEB_pT1->Fill(recoFirstPt,myWeight*sfWFirst*sfWSec);
    else 
      HrEE_pT1->Fill(recoFirstPt,myWeight*sfWFirst*sfWSec);
    if (fabs(recoSecEta)<1.5) 
      HrEB_pT2->Fill(recoSecPt,myWeight*sfWSec*sfWFirst);
    else 
      HrEE_pT2->Fill(recoSecPt,myWeight*sfWSec*sfWFirst);

    float zPtReco          = (recoFirstP+recoSecP).Perp();   
    float zYReco           = (recoFirstP+recoSecP).Rapidity();   

    recoFirstP.SetPtEtaPhiM(LepGood_pt[recoFirst],LepGood_eta[recoFirst],LepGood_phi[recoFirst],0.);

    float phiStarReco      = 1000.;
    float cosThetaStarReco = 1000.;
    if (LepGood_charge[recoFirst]>0 && LepGood_charge[recoSec]<0) {
      phiStarReco      = phiCS(recoFirstP, recoSecP);
      cosThetaStarReco = cosThetaCS(recoFirstP, recoSecP);
    } else if (LepGood_charge[recoFirst]<0 && LepGood_charge[recoSec]>0) {
      phiStarReco      = phiCS(recoSecP, recoFirstP);
      cosThetaStarReco = cosThetaCS(recoSecP, recoFirstP);
    }

    Hr_invMass       -> Fill(invMassReco,myWeight*sfWFirst*sfWSec);
    Hr_pTZ           -> Fill(zPtReco,myWeight*sfWFirst*sfWSec);
    Hr_pTZzoom       -> Fill(zPtReco,myWeight*sfWFirst*sfWSec);
    Hr_pTZzoom1      -> Fill(zPtReco,myWeight*sfWFirst*sfWSec);
    Hr_pTZzoom2      -> Fill(zPtReco,myWeight*sfWFirst*sfWSec);
    Hr_yZ            -> Fill(zYReco,myWeight*sfWFirst*sfWSec);
    if (phiStarReco<900) {
      Hr_phiStarZ      -> Fill(phiStarReco,myWeight*sfWFirst*sfWSec);
      Hr_cosThetaStarZ -> Fill(cosThetaStarReco,myWeight*sfWFirst*sfWSec);
    }

    bool isEBEB = false;
    if ( fabs(recoFirstEta)<1.5 && fabs(recoSecEta)<1.5 ) isEBEB= true;

    if (isEBEB) 
      HrEBEB_invMass->Fill(invMassReco,myWeight*sfWFirst*sfWSec);
    else 
      HrNotEBEB_invMass->Fill(invMassReco,myWeight*sfWFirst*sfWSec);

    if (fabs(zYReco)<1.5) {
      HrEBZ_invMass->Fill(invMassReco,myWeight*sfWFirst*sfWSec);
      HrEBZ_pTZ->Fill(zPtReco,myWeight*sfWFirst*sfWSec);
      HrEBZ_pTZzoom->Fill(zPtReco,myWeight*sfWFirst*sfWSec);
      HrEBZ_pTZzoom1->Fill(zPtReco,myWeight*sfWFirst*sfWSec);
      HrEBZ_pTZzoom2->Fill(zPtReco,myWeight*sfWFirst*sfWSec);
      if (phiStarReco<900) {
	HrEBZ_phiStarZ->Fill(phiStarReco,myWeight*sfWFirst*sfWSec);
	HrEBZ_cosThetaStarZ->Fill(cosThetaStarReco,myWeight*sfWFirst*sfWSec);
      }
    } else if (fabs(zYReco)>1.5 && fabs(zYReco)<2.5) {
      HrEEZ_invMass->Fill(invMassReco,myWeight*sfWFirst*sfWSec);
      HrEEZ_pTZ->Fill(zPtReco,myWeight*sfWFirst*sfWSec);
      HrEEZ_pTZzoom->Fill(zPtReco,myWeight*sfWFirst*sfWSec);
      HrEEZ_pTZzoom1->Fill(zPtReco,myWeight*sfWFirst*sfWSec);
      HrEEZ_pTZzoom2->Fill(zPtReco,myWeight*sfWFirst*sfWSec);
      if (phiStarReco<900) {
	HrEEZ_phiStarZ->Fill(phiStarReco,myWeight*sfWFirst*sfWSec);
	HrEEZ_cosThetaStarZ->Fill(cosThetaStarReco,myWeight*sfWFirst*sfWSec);
      }
    }


  } // event loop


  // ----------------------------------------
  // cosmetics
  Hr_pT1->SetLineColor(1);
  Hr_pT2->SetLineColor(1);
  Hr_eta1->SetLineColor(1);  
  Hr_eta2->SetLineColor(1);  
  Hr_invMass->SetLineColor(1);
  Hr_pTZ->SetLineColor(1);
  Hr_pTZzoom->SetLineColor(1);
  Hr_pTZzoom1->SetLineColor(1);
  Hr_pTZzoom2->SetLineColor(1);
  Hr_phiStarZ->SetLineColor(1);
  Hr_cosThetaStarZ->SetLineColor(1);
  Hr_yZ->SetLineColor(1);
  HrEB_pT1->SetLineColor(1);
  HrEB_pT2->SetLineColor(1);
  HrEE_pT1->SetLineColor(1);
  HrEE_pT2->SetLineColor(1);
  HrEBEB_invMass->SetLineColor(1);
  HrNotEBEB_invMass->SetLineColor(1);
  HrEBZ_invMass->SetLineColor(1);
  HrEBZ_pTZ->SetLineColor(1);
  HrEBZ_pTZzoom->SetLineColor(1);
  HrEBZ_pTZzoom1->SetLineColor(1);
  HrEBZ_pTZzoom2->SetLineColor(1);
  HrEBZ_phiStarZ->SetLineColor(1);
  HrEBZ_cosThetaStarZ->SetLineColor(1);
  HrEEZ_invMass->SetLineColor(1);
  HrEEZ_pTZ->SetLineColor(1);
  HrEEZ_pTZzoom->SetLineColor(1);
  HrEEZ_pTZzoom1->SetLineColor(1);
  HrEEZ_pTZzoom2->SetLineColor(1);
  HrEEZ_phiStarZ->SetLineColor(1);
  HrEEZ_cosThetaStarZ->SetLineColor(1);

  TFile *fileOut;
  if (!isMC) fileOut = new TFile("validationData.root","RECREATE");
  if (isMC)  fileOut = new TFile("validationMC.root","RECREATE");
  fileOut->cd();
  Hr_pT1->Write();
  Hr_pT2->Write();
  Hr_eta1->Write();
  Hr_eta2->Write();
  Hr_pTZ->Write();
  Hr_pTZzoom->Write();
  Hr_pTZzoom1->Write();
  Hr_pTZzoom2->Write();
  Hr_invMass->Write();
  Hr_phiStarZ->Write();
  Hr_cosThetaStarZ->Write();
  Hr_yZ->Write();
  HrEB_pT1->Write();
  HrEB_pT2->Write();
  HrEE_pT1->Write();
  HrEE_pT2->Write();
  HrEBEB_invMass->Write();
  HrNotEBEB_invMass->Write();
  HrEBZ_invMass->Write();
  HrEBZ_pTZ->Write();
  HrEBZ_pTZzoom->Write();
  HrEBZ_pTZzoom1->Write();
  HrEBZ_pTZzoom2->Write();
  HrEBZ_phiStarZ->Write();
  HrEBZ_cosThetaStarZ->Write();
  HrEEZ_invMass->Write();
  HrEEZ_pTZ->Write();
  HrEEZ_pTZzoom->Write();
  HrEEZ_pTZzoom1->Write();
  HrEEZ_pTZzoom2->Write();
  HrEEZ_phiStarZ->Write();
  HrEEZ_cosThetaStarZ->Write();

} // Loop 
