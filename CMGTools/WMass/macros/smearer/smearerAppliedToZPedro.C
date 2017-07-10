#define smearerAppliedToZPedro_cxx
#include "smearerAppliedToZPedro.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>

// RooFit headers
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"

using namespace std;

// Functions declaration
UInt_t FindBin( double value, double bins[], UInt_t nbins) {

  UInt_t nbinboundaries = nbins+1;
  UInt_t bin = 0;
  for (uint i=0; i < nbinboundaries; ++i) {
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


void smearerAppliedToZPedro::Loop() {

  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  // Data/MC SF file
  TFile fileSF("electrons_scale_factors.root");
  TH2F *histoSF = (TH2F*)fileSF.Get("electronsDATAMCratio_FO_ID_ISO");

  // Arrays for resolution maps
  const UInt_t NPtBins  = 12;
  const UInt_t NEtaBins = 14;
  double ptBins[NPtBins+1]   = { 15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,120.};
  double etaBins[NEtaBins+1] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4442, 1.566, 1.8, 2.0, 2.2, 2.4, 2.5, 2.6 }; 

  // Read histos with efficiency and resolution maps
  TFile *fileMaps = new TFile("/afs/cern.ch/work/c/crovelli/Wmass/CMSSW_5_3_22_patch1/src/CMGTools/WMass/macros/smearer/PtResolutionAndEffHistOUT.root");
  TH2F *Efficiency_PtEta = (TH2F*)fileMaps->Get("Efficiency_PtEta");
  TH2F *DoubleSidedCBShapeParamArray_mean   = (TH2F*)fileMaps->Get("DoubleSidedCBShapeParamArray_mean");
  TH2F *DoubleSidedCBShapeParamArray_sigma  = (TH2F*)fileMaps->Get("DoubleSidedCBShapeParamArray_sigma");
  TH2F *DoubleSidedCBShapeParamArray_alphaL = (TH2F*)fileMaps->Get("DoubleSidedCBShapeParamArray_alphaL");
  TH2F *DoubleSidedCBShapeParamArray_alphaR = (TH2F*)fileMaps->Get("DoubleSidedCBShapeParamArray_alphaR");
  TH2F *DoubleSidedCBShapeParamArray_nL     = (TH2F*)fileMaps->Get("DoubleSidedCBShapeParamArray_nL");
  TH2F *DoubleSidedCBShapeParamArray_nR     = (TH2F*)fileMaps->Get("DoubleSidedCBShapeParamArray_nR");

  // Create histograms
  TH1F *Hg_invMass = new TH1F("Hg_invMass", "Hg_invMass", 20, 80.,100.);
  TH1F *Hg_pT1     = new TH1F("Hg_pT1",     "Hg_pT1",     55, 25.,80.);
  TH1F *Hg_pT2     = new TH1F("Hg_pT2",     "Hg_pT2",     50, 20.,70.);
  TH1F *Hg_pTZ     = new TH1F("Hg_pTZ",     "Hg_pTZ",     80, 0.,80.);
  TH1F *Hg_pTZzoom = new TH1F("Hg_pTZzoom", "Hg_pTZzoom", 30, 0.,30.);
  TH1F *Hg_yZ      = new TH1F("Hg_yZ",      "Hg_yZ",      50, -3.,3.);
  TH1F *Hg_eta1    = new TH1F("Hg_eta1",    "Hg_eta1",    50, -2.5,2.5);
  TH1F *Hg_eta2    = new TH1F("Hg_eta2",    "Hg_eta2",    50, -2.5,2.5);  

  TH1F *HgEB_pT1   = new TH1F("HgEB_pT1",     "HgEB_pT1",     55, 25.,80.);
  TH1F *HgEB_pT2   = new TH1F("HgEB_pT2",     "HgEB_pT2",     50, 20.,70.);
  TH1F *HgEE_pT1   = new TH1F("HgEE_pT1",     "HgEE_pT1",     55, 25.,80.);
  TH1F *HgEE_pT2   = new TH1F("HgEE_pT2",     "HgEE_pT2",     50, 20.,70.);

  TH1F *HgEBEB_invMass = new TH1F("HgEBEB_invMass", "HgEBEB_invMass", 20,80.,100.);
  TH1F *HgEBEB_pTZ     = new TH1F("HgEBEB_pTZ",     "HgEBEB_pTZ",     80, 0.,80.);
  TH1F *HgEBEB_pTZzoom = new TH1F("HgEBEB_pTZzoom", "HgEBEB_pTZzoom", 30, 0.,30.);

  TH1F *HgNotEBEB_invMass = new TH1F("HgNotEBEB_invMass", "HgNotEBEB_invMass", 20,80.,100.);
  TH1F *HgNotEBEB_pTZ     = new TH1F("HgNotEBEB_pTZ",     "HgNotEBEB_pTZ",     80, 0.,80.);
  TH1F *HgNotEBEB_pTZzoom = new TH1F("HgNotEBEB_pTZzoom", "HgNotEBEB_pTZzoom", 30, 0.,30.);

  // Event loop 
  cout << "Entries: " << nentries << endl;
  for (Long64_t jentry=0; jentry<1000000;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if (ientry % 10000 == 0) cout << "Event " << ientry << endl;    
    if (ientry % 2 != 0) continue;

    // event weight
    double myWeight = 1.;

    // gen-level electron selection
    if (nl!=2) continue;
    //if (abs(pid[0])!=11 || abs(pid[1])!=11) continue;  // chiaraaaaaaaaaaaa!!!
    if (abs(pid[0])!=13 || abs(pid[1])!=13) continue;
    if (charge[0]*charge[1]>0) { cout << "charge error" << endl; continue; }   

    // Before FSR
    TLorentzVector gen0(0,0,0,0); 
    TLorentzVector gen1(0,0,0,0); 
    gen0.SetPtEtaPhiM(dressed_pt[0], dressed_eta[0], dressed_phi[0], dressed_m[0]);
    gen1.SetPtEtaPhiM(dressed_pt[1], dressed_eta[1], dressed_phi[1], dressed_m[1]); 
    float gen0Pt  = dressed_pt[0];
    float gen1Pt  = dressed_pt[1];
    float gen0Eta = dressed_eta[0];
    float gen1Eta = dressed_eta[1];

    // Remove events with gen leptons outside range
    if (gen0Pt<15) continue;
    if (gen1Pt<15) continue;
    if (fabs(gen0Eta)>2.6) continue;
    if (fabs(gen1Eta)>2.6) continue;
    
    // Which histo bins
    int gen0PtBin  = FindBin( gen0Pt, ptBins, NPtBins);    
    int gen1PtBin  = FindBin( gen1Pt, ptBins, NPtBins);
    int gen0EtaBin = FindBin( fabs(gen0Eta) , etaBins, NEtaBins);
    int gen1EtaBin = FindBin( fabs(gen1Eta) , etaBins, NEtaBins);
    if (gen0EtaBin==0 || gen1EtaBin==0 || gen0EtaBin==15 || gen1EtaBin==15) cout << "STRANGE!" << endl;
    if (gen0PtBin==0 || gen1PtBin==0) cout << "STRANGE!" << endl;

    // Smearing
    float theMean0   = DoubleSidedCBShapeParamArray_mean->GetBinContent(gen0PtBin,gen0EtaBin);
    float theMean1   = DoubleSidedCBShapeParamArray_mean->GetBinContent(gen1PtBin,gen1EtaBin);
    float theSigma0  = DoubleSidedCBShapeParamArray_sigma->GetBinContent(gen0PtBin,gen0EtaBin);
    float theSigma1  = DoubleSidedCBShapeParamArray_sigma->GetBinContent(gen1PtBin,gen1EtaBin);
    float theAlphaL0 = DoubleSidedCBShapeParamArray_alphaL->GetBinContent(gen0PtBin,gen0EtaBin);
    float theAlphaL1 = DoubleSidedCBShapeParamArray_alphaL->GetBinContent(gen1PtBin,gen1EtaBin);
    float theAlphaR0 = DoubleSidedCBShapeParamArray_alphaR->GetBinContent(gen0PtBin,gen0EtaBin);
    float theAlphaR1 = DoubleSidedCBShapeParamArray_alphaR->GetBinContent(gen1PtBin,gen1EtaBin);
    float theNL0 = DoubleSidedCBShapeParamArray_nL->GetBinContent(gen0PtBin,gen0EtaBin);
    float theNL1 = DoubleSidedCBShapeParamArray_nL->GetBinContent(gen1PtBin,gen1EtaBin);
    float theNR0 = DoubleSidedCBShapeParamArray_nR->GetBinContent(gen0PtBin,gen0EtaBin);
    float theNR1 = DoubleSidedCBShapeParamArray_nR->GetBinContent(gen1PtBin,gen1EtaBin);

    RooRealVar ptRes0("ptRes0","ptRes0",-0.75,0.75);
    RooRealVar *mean0   = new RooRealVar("mean0","mean0",theMean0);
    RooRealVar *sigma0  = new RooRealVar("sigma0","sigma0",theSigma0);
    RooRealVar *alphaL0 = new RooRealVar("alphaL0","alphaL0",theAlphaL0);
    RooRealVar *alphaR0 = new RooRealVar("alphaR0","alphaR0",theAlphaR0);
    RooRealVar *nL0     = new RooRealVar("nL0","nL0",theNL0);
    RooRealVar *nR0     = new RooRealVar("nR0","nR0",theNR0);
    RooDoubleCB *model0 = new RooDoubleCB("model0","model0",ptRes0,*mean0,*sigma0,*alphaL0,*nL0,*alphaR0,*nR0);

    RooRealVar ptRes1("ptRes1","ptRes1",-0.75,0.75);
    RooRealVar *mean1   = new RooRealVar("mean1","mean1",theMean1);
    RooRealVar *sigma1  = new RooRealVar("sigma1","sigma1",theSigma1);
    RooRealVar *alphaL1 = new RooRealVar("alphaL1","alphaL1",theAlphaL1);
    RooRealVar *alphaR1 = new RooRealVar("alphaR1","alphaR1",theAlphaR1);
    RooRealVar *nL1     = new RooRealVar("nL1","nL1",theNL1);
    RooRealVar *nR1     = new RooRealVar("nR1","nR1",theNR1);
    RooDoubleCB *model1 = new RooDoubleCB("model1","model1",ptRes1,*mean1,*sigma1,*alphaL1,*nL1,*alphaR1,*nR1);

    RooDataSet* data0 = model0->generate(ptRes0,1);
    RooDataSet* data1 = model1->generate(ptRes1,1);
    RooArgSet ras0 = *data0->get(0);
    RooArgSet ras1 = *data1->get(0);
    RooRealVar* var0 = (RooRealVar*)ras0.find("ptRes0");
    RooRealVar* var1 = (RooRealVar*)ras1.find("ptRes1");
    float ptOverPtt0 = 1.+(var0->getVal());
    float ptOverPtt1 = 1.+(var1->getVal());
    delete data0;
    delete data1;

    delete mean0;
    delete sigma0;
    delete alphaL0;
    delete alphaR0;
    delete nL0;
    delete nR0;
    delete model0;
    delete mean1;
    delete sigma1;
    delete alphaL1;
    delete alphaR1;
    delete nL1;
    delete nR1;
    delete model1;

    // smeared quantities
    float smearedPT0  = ptOverPtt0*gen0Pt;
    float smearedPT1  = ptOverPtt1*gen1Pt;
    // not smeared....
    float smearedETA0 = gen0Eta;
    float smearedETA1 = gen1Eta;

    // Efficiency
    int gen0PtBinEff  = int(gen0Pt-15)+1;
    int gen1PtBinEff  = int(gen1Pt-15)+1;  
    int gen0EtaBinEff = int( (fabs(gen0Eta))/0.05 )+1;
    int gen1EtaBinEff = int( (fabs(gen1Eta))/0.05 )+1;
    float theEff0 = Efficiency_PtEta->GetBinContent(gen0PtBinEff,gen0EtaBinEff);
    float theEff1 = Efficiency_PtEta->GetBinContent(gen1PtBinEff,gen1EtaBinEff);

    // Smeared quantities: pT and eta
    float highestPt  = smearedPT0;
    float highestEta = smearedETA0;
    float highestEff = theEff0;
    float lowestPt   = smearedPT1;
    float lowestEta  = smearedETA1;
    float lowestEff  = theEff1;
    if (smearedPT0<smearedPT1) {
      highestPt  = smearedPT1;
      highestEta = smearedETA1;
      highestEff = theEff1;
      lowestPt   = smearedPT0;
      lowestEta  = smearedETA0;
      lowestEff  = theEff0;
    }

    // Smeared quantities: Z
    TLorentzVector gen0Sm(0,0,0,0); 
    TLorentzVector gen1Sm(0,0,0,0); 
    gen0Sm.SetPtEtaPhiM(smearedPT0, dressed_eta[0], dressed_phi[0], dressed_m[0]);
    gen1Sm.SetPtEtaPhiM(smearedPT1, dressed_eta[1], dressed_phi[1], dressed_m[1]);
    float invMassSm = (gen0Sm+gen1Sm).M();
    float zPtSm = (gen0Sm+gen1Sm).Perp();
    float zYSm = (gen0Sm+gen1Sm).Rapidity();
    
    // Further selection to match the one in data
    if (highestPt<25) continue;
    if (lowestPt<20)  continue;
    if (invMassSm<80 || invMassSm>100) continue;

    // Dta/MC efficiency scale factors
    int binXSfFirst = 0.;
    int binYSfFirst = 0.;
    int binXSfSec   = 0.;
    int binYSfSec   = 0.;
    float sfWFirst  = 1.;
    float sfWSec    = 1.;

    if (fabs(highestEta)<=0.8) binXSfFirst = 1;
    else if (fabs(highestEta)<=1.4442 && fabs(highestEta)>0.8)   binXSfFirst = 2;
    else if (fabs(highestEta)<=1.556 && fabs(highestEta)>1.4442) binXSfFirst = 3;  
    else if (fabs(highestEta)<=2. && fabs(highestEta)>1.556)     binXSfFirst = 4;  
    else if (fabs(highestEta)<=2.6 && fabs(highestEta)>2.)       binXSfFirst = 5;    
    
    if (fabs(lowestEta)<=0.8) binXSfSec = 1;
    else if (fabs(lowestEta)<=1.4442 && fabs(lowestEta)>0.8)   binXSfSec = 2;
    else if (fabs(lowestEta)<=1.556 && fabs(lowestEta)>1.4442) binXSfSec = 3;  
    else if (fabs(lowestEta)<=2. && fabs(lowestEta)>1.556)     binXSfSec = 4;  
    else if (fabs(lowestEta)<=2.6 && fabs(lowestEta)>2.)       binXSfSec = 5;    
    
    if (highestPt<=15 && highestPt>10) binYSfFirst = 1;
    else if (highestPt<=20 && highestPt>15)  binYSfFirst = 2;
    else if (highestPt<=30 && highestPt>20)  binYSfFirst = 3;
    else if (highestPt<=40 && highestPt>30)  binYSfFirst = 4;
    else if (highestPt<=50 && highestPt>40)  binYSfFirst = 5;
    else if (highestPt>50) binYSfFirst = 6;
    
    if (lowestPt<=15 && lowestPt>10) binYSfSec = 1;
    else if (lowestPt<=20 && lowestPt>15)  binYSfSec = 2;
    else if (lowestPt<=30 && lowestPt>20)  binYSfSec = 3;
    else if (lowestPt<=40 && lowestPt>30)  binYSfSec = 4;
    else if (lowestPt<=50 && lowestPt>40)  binYSfSec = 5;
    else if (lowestPt>50) binYSfSec = 6;
    
    if (binXSfFirst>0 && binYSfFirst>0 && binYSfFirst>0 && binYSfSec>0) {
      sfWFirst = histoSF->GetBinContent(binXSfFirst,binYSfFirst);
      sfWSec   = histoSF->GetBinContent(binXSfSec,binYSfSec);
    } else {
      cout << highestPt << " " << lowestPt << " " << highestEta << " " << lowestEta << endl;
      cout << "problem!" << endl;
    }
    
    // Plots
    myWeight = myWeight*highestEff*lowestEff*sfWFirst*sfWSec;
    Hg_pT1->Fill(highestPt,myWeight);
    Hg_pT2->Fill(lowestPt,myWeight);
    Hg_eta1->Fill(highestEta,myWeight);
    Hg_eta2->Fill(lowestEta,myWeight);
    Hg_invMass->Fill(invMassSm,myWeight);
    Hg_pTZ->Fill(zPtSm,myWeight);
    Hg_pTZzoom->Fill(zPtSm,myWeight);
    Hg_yZ->Fill(zYSm,myWeight);

    float highestInEB = true;
    float lowestInEB  = true;
    if (fabs(highestEta)>1.5) highestInEB = false;
    if (fabs(lowestEta)>1.5)  lowestInEB  = false; 
    if (highestInEB) HgEB_pT1->Fill(highestPt,myWeight);
    else HgEE_pT1->Fill(highestPt,myWeight);
    if (lowestInEB) HgEB_pT2->Fill(lowestPt,myWeight);
    else HgEE_pT2->Fill(lowestPt,myWeight);
    if (highestInEB && lowestInEB) {
      HgEBEB_invMass->Fill(invMassSm,myWeight);
      HgEBEB_pTZ->Fill(zPtSm,myWeight);
      HgEBEB_pTZzoom->Fill(zPtSm,myWeight);
    } else {
      HgNotEBEB_invMass->Fill(invMassSm,myWeight);
      HgNotEBEB_pTZ->Fill(zPtSm,myWeight);
      HgNotEBEB_pTZzoom->Fill(zPtSm,myWeight);
    }

  } // event loop


  // ----------------------------------------
  TFile fileOut("smearedAppliedToZPedro.root","RECREATE");
  fileOut.cd();
  Hg_pT1->Write();
  Hg_pT2->Write();
  Hg_pTZ->Write();
  Hg_pTZzoom->Write();
  Hg_yZ->Write();
  Hg_invMass->Write();
  Hg_eta1->Write();
  Hg_eta2->Write();
  HgEB_pT1->Write();
  HgEB_pT2->Write();
  HgEE_pT1->Write();
  HgEE_pT2->Write();
  HgEBEB_invMass->Write(); 
  HgEBEB_pTZ->Write(); 
  HgEBEB_pTZzoom->Write();
  HgNotEBEB_invMass->Write(); 
  HgNotEBEB_pTZ->Write(); 
  HgNotEBEB_pTZzoom->Write();

} // Loop 
