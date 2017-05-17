#define eleFastSmearer_cxx
#include "eleFastSmearer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>

using namespace std;

// Functions declaration 
void initializeHistograms(vector<vector<TH1F*> > &hists, string name, Int_t NPtBins, Int_t NEtaBins, Int_t nbins, Double_t xmin, Double_t xmax) {

  hists.resize(NPtBins+2);
  for (int i=0; i < NPtBins+2; ++i) {
    hists[i].resize(NEtaBins+2);
    for (int j=0; j < NEtaBins+2; ++j) {
      hists[i][j]= new TH1F(Form("%s_PtBin%d_EtaBin%d",name.c_str(),i,j), "; (Reco Pt - Gen Pt)/Gen Pt; Number of Events", nbins, xmin, xmax);
    }
  }
}

// Bins from 1 to nbins, overflow in nbins+1                            
Int_t FindBin( double value, double bins[], Int_t nbins) {

  Int_t nbinboundaries = nbins+1;
  Int_t bin = 0;
  for (int i=0; i < nbinboundaries; ++i) {
    if (i < nbinboundaries-1) {
      if (value >= bins[i] && value < bins[i+1]) {
        bin = i+1;
        break;
      }
    } else if (i == nbinboundaries-1) {
      if (value >= bins[i]) {
        bin = nbinboundaries;
        break;
      }
    }
  }
  return bin;
}

void eleFastSmearer::Loop() {

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  // Create arrays to store the resolution fit maps
  const Int_t NPtBins  = 12;
  const Int_t NEtaBins = 14;
  double ptBins[NPtBins+1]   = { 15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,120.};
  double etaBins[NEtaBins+1] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4442, 1.566, 1.8, 2.0, 2.2, 2.4, 2.5, 2.6 };

  // Create histograms: resolution maps
  vector<vector<TH1F*> > PtResolution_PtEta;
  initializeHistograms(PtResolution_PtEta, "ElectronsPtResolution", NPtBins, NEtaBins, 100, -0.75, 0.75 );

  // Create histograms: efficiency maps                                                                         
  TH2F *NumEfficiency_PtEta = new TH2F("NumEfficiency_PtEta","NumEfficiency_PtEta",105,15.,120.,52,0.,2.6);
  TH2F *DenEfficiency_PtEta = new TH2F("DenEfficiency_PtEta","DenEfficiency_PtEta",105,15.,120.,52,0.,2.6);

  // Counters                                                                                                       
  int totEvents=0;
  int eleEvents=0;
  int genAcc=0;
  int recoTwo=0;
  int recoOne=0;
  int recoZero=0;

  // Event loop                                                                          
  cout << "Entries: " << nentries << endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // to use different events to determine and validate smearings                       
    if (ientry % 2 == 0) continue;
    if (ientry % 1000001 == 0) cout << "Event " << ientry << endl;
    totEvents++;

    // event weight                                                                                       
    double myWeight = puWeight;

    // gen-level electron selection                            
    int genEleIdx = -99;
    int genPosIdx = -99;
    for (int imc=0; imc<nGenP6StatusThree; imc++) {
      if (genEleIdx>-10 && genPosIdx>-10) continue;
      if (GenP6StatusThree_motherId[imc]==23) {
        if (GenP6StatusThree_pdgId[imc]==-11) genEleIdx = imc;
        if (GenP6StatusThree_pdgId[imc]==+11) genPosIdx = imc;
      }
    }
    if (genEleIdx<-10 && genPosIdx>-10) { cout << "error" << endl; continue; }
    if (genEleIdx>-10 && genPosIdx<-10) { cout << "error" << endl; continue; }
    if (genEleIdx<-10 && genPosIdx<-10) continue;
    eleEvents++;

    TLorentzVector genEle(0,0,0,0);
    TLorentzVector genPos(0,0,0,0);
    genEle.SetPtEtaPhiM(GenP6StatusThree_pt[genEleIdx], GenP6StatusThree_eta[genEleIdx], GenP6StatusThree_phi[genEleIdx], 0.);
    genPos.SetPtEtaPhiM(GenP6StatusThree_pt[genPosIdx], GenP6StatusThree_eta[genPosIdx], GenP6StatusThree_phi[genPosIdx], 0.);
    float genElePt   = GenP6StatusThree_pt[genEleIdx];
    float genPosPt   = GenP6StatusThree_pt[genPosIdx];
    float genEleEta  = GenP6StatusThree_eta[genEleIdx];
    float genPosEta  = GenP6StatusThree_eta[genPosIdx];

    // Remove events with gen leptons outside range 
    if (genElePt<15) continue;
    if (genPosPt<15) continue;
    if (fabs(genEleEta)>2.6) continue;
    if (fabs(genPosEta)>2.6) continue;
    genAcc++;
    
    // Which histo bins 
    int genElePtBin  = FindBin( genElePt, ptBins, NPtBins);
    int genPosPtBin  = FindBin( genPosPt, ptBins, NPtBins);
    int genEleEtaBin = FindBin( fabs(genEleEta) , etaBins, NEtaBins);
    int genPosEtaBin = FindBin( fabs(genPosEta) , etaBins, NEtaBins);

    // MC efficiency: denominator
    DenEfficiency_PtEta->Fill(genElePt, fabs(genEleEta));   // here I fill wo pu-weight => put it at num only                  
    DenEfficiency_PtEta->Fill(genPosPt, fabs(genPosEta));   // here I fill wo pu-weight => put it at num only               

    // trigger bit                                                                              
    // if (!HLT_DoubleEl) continue;                                                                                     

    // Reco electrons and match with genLevel    
    float minDrEle=1000.;
    float minDrPos=1000.;
    int recoEleIdx = -99;
    int recoPosIdx = -99;
    for (int ilep=0; ilep<nLepGood; ilep++) {

      // electrons only                          
      if (abs(LepGood_pdgId[ilep])!=11) continue;

      // trigger match   
      // if (LepGood_trgMatch[ilep]<0.5) continue;  

      // eta and pT cut (minimal, probably looser than applied in the analysis)   
      if (fabs(LepGood_scEta[ilep])>2.5) continue;
      if (LepGood_pt[ilep]<15) continue;

      // barrel-endcap transition    
      if ( fabs(LepGood_scEta[ilep])>1.442 && fabs(LepGood_scEta[ilep])<1.566 ) continue;

      // trigger emulation
      if (LepGood_eleMVAPreselId[ilep]==0) continue;

      // ID + isolation + conv.rejection
      if (LepGood_eleMVAId[ilep]<2)       continue;
      if (LepGood_relIso04[ilep]>0.15)    continue;
      if (LepGood_convVetoFull[ilep]==0)  continue;   // LepGood_lostHits==0 + LepGood_conv

      // for selected electrons
      TLorentzVector thisReco(0,0,0,0);
      thisReco.SetPtEtaPhiM(LepGood_pt[ilep],LepGood_eta[ilep],LepGood_phi[ilep],0.);
      if (thisReco.DeltaR(genEle)<minDrEle) {
        minDrEle = thisReco.DeltaR(genEle);
        recoEleIdx = ilep;
      }
      if (thisReco.DeltaR(genPos)<minDrPos) {
        minDrPos = thisReco.DeltaR(genPos);
        recoPosIdx = ilep;
      }
    }

    if (minDrEle>0.1 && minDrPos>0.1) { recoZero++; continue; }
    if (minDrEle<0.1 && minDrPos<0.1 && recoEleIdx==recoPosIdx) { cout << "problem" << endl; continue; }

    if (minDrEle<0.1) {
      float recoElePt  = LepGood_pt[recoEleIdx];
      NumEfficiency_PtEta->Fill(genElePt, fabs(genEleEta), myWeight);
      PtResolution_PtEta[genElePtBin][genEleEtaBin]->Fill( (recoElePt-genElePt)/genElePt, myWeight);
    }
    if (minDrPos<0.1) {
      float recoPosPt  = LepGood_pt[recoPosIdx];
      NumEfficiency_PtEta->Fill(genPosPt, fabs(genPosEta), myWeight);
      PtResolution_PtEta[genPosPtBin][genPosEtaBin]->Fill( (recoPosPt-genPosPt)/genPosPt, myWeight);
    }
    if ( minDrEle<0.1 && minDrPos<0.1) recoTwo++;
    if ((minDrEle<0.1 && minDrPos>0.1) || (minDrEle>0.1 && minDrPos<0.1)) recoOne++;

  } // event loop                                                                  

  // ----------------------------------------                                                                                                                         
  // Efficiencies              
  TH2F *Efficiency_PtEta = (TH2F*)NumEfficiency_PtEta->Clone("Efficiency_PtEta");
  Efficiency_PtEta->Divide(DenEfficiency_PtEta);

  // Output file 
  TFile *fileOutput = new TFile("PtResolutionAndEffHist.root", "RECREATE");
  for (int i=0; i < NPtBins+2; ++i) {
    for (int j=0; j < NEtaBins+2; ++j) {
      fileOutput->WriteTObject(PtResolution_PtEta[i][j], PtResolution_PtEta[i][j]->GetName(), "WriteDelete");
    }
  }
  fileOutput->WriteTObject(NumEfficiency_PtEta,NumEfficiency_PtEta->GetName(), "WriteDelete");
  fileOutput->WriteTObject(DenEfficiency_PtEta,DenEfficiency_PtEta->GetName(), "WriteDelete");
  fileOutput->WriteTObject(Efficiency_PtEta,Efficiency_PtEta->GetName(), "WriteDelete");

  fileOutput->Close();
  delete fileOutput;


  cout << "Tot number of events = "  << totEvents << endl;
  cout << "Zee events = "            << eleEvents << endl;
  cout << "Gen within acceptance = " << genAcc << endl;
  cout << "ZeroEle reco = " << recoZero << endl;
  cout << "OneEle reco = "  << recoOne  << endl;
  cout << "TwoEle reco = "  << recoTwo  << endl;

}
