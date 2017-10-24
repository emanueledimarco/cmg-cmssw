#define smearerAppliedToZ_cxx
#include "smearerAppliedToZ.h"
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


void smearerAppliedToZ::Loop() {

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
  TH1F *Hg_invMass = new TH1F("Hg_invMass", "Hg_invMass", 20,80.,100.);
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

    // chiara: per avere gli stessi eventi che nella valid
    if (ientry % 2 != 0) continue;

    // event weight - chiara: non metto il peso xsec perche' comunque confronto le shapes
    //                        non metto pu weight che in teoria e' entrato in media nello smearing 
    double myWeight = 1;

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
    // We want 2 electrons at gen level  
    if (genEleIdx<-10 && genPosIdx>-10) { cout << "error" << endl; continue; }
    if (genEleIdx>-10 && genPosIdx<-10) { cout << "error" << endl; continue; }
    if (genEleIdx<-10 || genPosIdx<-10) continue;

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

    // Which histo bins
    int genElePtBin  = FindBin( genElePt, ptBins, NPtBins);    
    int genPosPtBin  = FindBin( genPosPt, ptBins, NPtBins);
    int genEleEtaBin = FindBin( fabs(genEleEta) , etaBins, NEtaBins);
    int genPosEtaBin = FindBin( fabs(genPosEta) , etaBins, NEtaBins);
    if (genEleEtaBin==0 || genPosEtaBin==0 || genEleEtaBin==15 || genPosEtaBin==15) cout << "STRANGE!" << endl;
    if (genElePtBin==0 || genPosPtBin==0) cout << "STRANGE!" << endl;

    // Smearing
    float theMeanEle   = DoubleSidedCBShapeParamArray_mean->GetBinContent(genElePtBin,genEleEtaBin);
    float theMeanPos   = DoubleSidedCBShapeParamArray_mean->GetBinContent(genPosPtBin,genPosEtaBin);
    float theSigmaEle  = DoubleSidedCBShapeParamArray_sigma->GetBinContent(genElePtBin,genEleEtaBin);
    float theSigmaPos  = DoubleSidedCBShapeParamArray_sigma->GetBinContent(genPosPtBin,genPosEtaBin);
    float theAlphaLEle = DoubleSidedCBShapeParamArray_alphaL->GetBinContent(genElePtBin,genEleEtaBin);
    float theAlphaLPos = DoubleSidedCBShapeParamArray_alphaL->GetBinContent(genPosPtBin,genPosEtaBin);
    float theAlphaREle = DoubleSidedCBShapeParamArray_alphaR->GetBinContent(genElePtBin,genEleEtaBin);
    float theAlphaRPos = DoubleSidedCBShapeParamArray_alphaR->GetBinContent(genPosPtBin,genPosEtaBin);
    float theNLEle = DoubleSidedCBShapeParamArray_nL->GetBinContent(genElePtBin,genEleEtaBin);
    float theNLPos = DoubleSidedCBShapeParamArray_nL->GetBinContent(genPosPtBin,genPosEtaBin);
    float theNREle = DoubleSidedCBShapeParamArray_nR->GetBinContent(genElePtBin,genEleEtaBin);
    float theNRPos = DoubleSidedCBShapeParamArray_nR->GetBinContent(genPosPtBin,genPosEtaBin);

    RooRealVar ptResEle("ptResEle","ptResEle",-0.75,0.75);
    RooRealVar *meanEle   = new RooRealVar("meanEle","meanEle",theMeanEle);
    RooRealVar *sigmaEle  = new RooRealVar("sigmaEle","sigmaEle",theSigmaEle);
    RooRealVar *alphaLEle = new RooRealVar("alphaLEle","alphaLEle",theAlphaLEle);
    RooRealVar *alphaREle = new RooRealVar("alphaREle","alphaREle",theAlphaREle);
    RooRealVar *nLEle     = new RooRealVar("nLEle","nLEle",theNLEle);
    RooRealVar *nREle     = new RooRealVar("nREle","nREle",theNREle);
    RooDoubleCB *modelEle = new RooDoubleCB("modelEle","modelEle",ptResEle,*meanEle,*sigmaEle,*alphaLEle,*nLEle,*alphaREle,*nREle);

    RooRealVar ptResPos("ptResPos","ptResPos",-0.75,0.75);
    RooRealVar *meanPos   = new RooRealVar("meanPos","meanPos",theMeanPos);
    RooRealVar *sigmaPos  = new RooRealVar("sigmaPos","sigmaPos",theSigmaPos);
    RooRealVar *alphaLPos = new RooRealVar("alphaLPos","alphaLPos",theAlphaLPos);
    RooRealVar *alphaRPos = new RooRealVar("alphaRPos","alphaRPos",theAlphaRPos);
    RooRealVar *nLPos     = new RooRealVar("nLPos","nLPos",theNLPos);
    RooRealVar *nRPos     = new RooRealVar("nRPos","nRPos",theNRPos);
    RooDoubleCB *modelPos = new RooDoubleCB("modelPos","modelPos",ptResPos,*meanPos,*sigmaPos,*alphaLPos,*nLPos,*alphaRPos,*nRPos);

    RooDataSet* dataEle = modelEle->generate(ptResEle,1);
    RooDataSet* dataPos = modelPos->generate(ptResPos,1);
    RooArgSet rasEle = *dataEle->get(0);
    RooArgSet rasPos = *dataPos->get(0);
    RooRealVar* varEle = (RooRealVar*)rasEle.find("ptResEle");
    RooRealVar* varPos = (RooRealVar*)rasPos.find("ptResPos");
    float ptOverPttEle = 1.+(varEle->getVal());
    float ptOverPttPos = 1.+(varPos->getVal());
    delete dataEle;
    delete dataPos;

    delete meanEle;
    delete sigmaEle;
    delete alphaLEle;
    delete alphaREle;
    delete nLEle;
    delete nREle;
    delete modelEle;
    delete meanPos;
    delete sigmaPos;
    delete alphaLPos;
    delete alphaRPos;
    delete nLPos;
    delete nRPos;
    delete modelPos;

    // smeared quantities
    float smearedPTEle  = ptOverPttEle*genElePt;
    float smearedPTPos  = ptOverPttPos*genPosPt;
    // not smeared....
    float smearedETAEle = genEleEta;
    float smearedETAPos = genPosEta;

    // Efficiency
    int genElePtBinEff  = int(genElePt-15)+1;
    int genPosPtBinEff  = int(genPosPt-15)+1;  
    int genEleEtaBinEff = int( (fabs(genEleEta))/0.05 )+1;
    int genPosEtaBinEff = int( (fabs(genPosEta))/0.05 )+1;
    float theEffEle = Efficiency_PtEta->GetBinContent(genElePtBinEff,genEleEtaBinEff);
    float theEffPos = Efficiency_PtEta->GetBinContent(genPosPtBinEff,genPosEtaBinEff);

    // Smeared quantities: pT and eta
    float highestPt  = smearedPTEle;
    float highestEta = smearedETAEle;
    float highestEff = theEffEle;
    float lowestPt   = smearedPTPos;
    float lowestEta  = smearedETAPos;
    float lowestEff  = theEffPos;
    if (smearedPTEle<smearedPTPos) {
      highestPt  = smearedPTPos;
      highestEta = smearedETAPos;
      highestEff = theEffPos;
      lowestPt   = smearedPTEle;
      lowestEta  = smearedETAEle;
      lowestEff  = theEffEle;
    }

    // Smeared quantities: Z
    TLorentzVector genEleSm(0,0,0,0); 
    TLorentzVector genPosSm(0,0,0,0); 
    genEleSm.SetPtEtaPhiM(smearedPTEle, GenP6StatusThree_eta[genEleIdx], GenP6StatusThree_phi[genEleIdx], 0.);
    genPosSm.SetPtEtaPhiM(smearedPTPos, GenP6StatusThree_eta[genPosIdx], GenP6StatusThree_phi[genPosIdx], 0.);
    float invMassSm = (genEleSm+genPosSm).M();
    float zPtSm     = (genEleSm+genPosSm).Perp();
    float zYSm      = (genEleSm+genPosSm).Rapidity();
    
    // Further selection to match the one in data
    if (highestPt<25) continue;
    if (lowestPt<20)  continue;
    if (invMassSm<80 || invMassSm>100) continue;


    // Data/MC efficiency scale factors
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
  TFile fileOut("smearedAppliedToZ.root","RECREATE");
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
