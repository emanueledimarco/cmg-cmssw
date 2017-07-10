#define eleFastSmearerValid_cxx
#include "eleFastSmearerValid.h"
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


void eleFastSmearerValid::Loop() {

  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  
  // Arrays for resolution maps
  const Int_t NPtBins  = 12;
  const Int_t NEtaBins = 14;
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

  TH1F *Hgne_invMass = new TH1F("Hgne_invMass", "Hgne_invMass", 20, 80.,100.);
  TH1F *Hgne_pT1     = new TH1F("Hgne_pT1",     "Hgne_pT1",     55, 25.,80.);
  TH1F *Hgne_pT2     = new TH1F("Hgne_pT2",     "Hgne_pT2",     50, 20.,70.);
  TH1F *Hgne_pTZ     = new TH1F("Hgne_pTZ",     "Hgne_pTZ",     80, 0.,80.);
  TH1F *Hgne_pTZzoom = new TH1F("Hgne_pTZzoom", "Hgne_pTZzoom", 30, 0.,30.);
  TH1F *Hgne_yZ      = new TH1F("Hgne_yZ",      "Hgne_yZ",      50, -3.,3.);
  TH1F *Hgne_eta1    = new TH1F("Hgne_eta1",    "Hgne_eta1",    50, -2.5,2.5);
  TH1F *Hgne_eta2    = new TH1F("Hgne_eta2",    "Hgne_eta2",    50, -2.5,2.5);

  TH1F *Hr_invMass = new TH1F("Hr_invMass", "Hr_invMass", 20, 80.,100.);
  TH1F *Hr_pT1     = new TH1F("Hr_pT1",     "Hr_pT1",     55, 25.,80.);
  TH1F *Hr_pT2     = new TH1F("Hr_pT2",     "Hr_pT2",     50, 20.,70.);
  TH1F *Hr_pTZ     = new TH1F("Hr_pTZ",     "Hr_pTZ",     80, 0.,80.);
  TH1F *Hr_pTZzoom = new TH1F("Hr_pTZzoom", "Hr_pTZzoom", 30, 0.,30.);
  TH1F *Hr_yZ      = new TH1F("Hr_yZ",      "Hr_yZ",      50, -3.,3.);
  TH1F *Hr_eta1    = new TH1F("Hr_eta1",    "Hr_eta1",    50, -2.5,2.5);
  TH1F *Hr_eta2    = new TH1F("Hr_eta2",    "Hr_eta2",    50, -2.5,2.5);

  // Event loop 
  cout << "Entries: " << nentries << endl;
  for (Long64_t jentry=0; jentry<1000000;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // to use different events to determine and validate smearings
    if (ientry % 2 != 0) continue;
    if (ientry % 10000 == 0) cout << "Event " << ientry << endl;    

    // event weight
    double myWeight   = 1.;                   // I don't use pu weight on smeared quantities (pu should be accounted for in the weights)
    double myWeightPU = myWeight*puWeight;

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

    // Smeared plots
    TLorentzVector genEleSm(0,0,0,0); 
    TLorentzVector genPosSm(0,0,0,0); 
    genEleSm.SetPtEtaPhiM(smearedPTEle, GenP6StatusThree_eta[genEleIdx], GenP6StatusThree_phi[genEleIdx], 0.);
    genPosSm.SetPtEtaPhiM(smearedPTPos, GenP6StatusThree_eta[genPosIdx], GenP6StatusThree_phi[genPosIdx], 0.);
    float invMassSm = (genEleSm+genPosSm).M();
    float zPtSm     = (genEleSm+genPosSm).Perp();
    float zYSm      = (genEleSm+genPosSm).Rapidity();
    float maxSmearPt  = smearedPTEle;
    float minSmearPt  = smearedPTPos;
    float maxSmearEta = smearedETAEle;
    float minSmearEta = smearedETAPos;
    if (smearedPTEle<smearedPTPos) { 
      maxSmearPt  = smearedPTPos;
      minSmearPt  = smearedPTEle;
      maxSmearEta = smearedETAPos;
      minSmearEta = smearedETAEle;
    }
    if ( invMassSm<100 && invMassSm>80 && maxSmearPt>25 && minSmearPt>20 ) {
      Hg_pT1  -> Fill(maxSmearPt,myWeight*theEffEle*theEffPos);     
      Hg_pT2  -> Fill(minSmearPt,myWeight*theEffEle*theEffPos);     
      Hg_eta1 -> Fill(maxSmearEta,myWeight*theEffEle*theEffPos);
      Hg_eta2 -> Fill(minSmearEta,myWeight*theEffEle*theEffPos);
      Hg_invMass->Fill(invMassSm,myWeight*theEffEle*theEffPos);
      Hg_pTZ->Fill(zPtSm,myWeight*theEffEle*theEffPos);
      Hg_pTZzoom->Fill(zPtSm,myWeight*theEffEle*theEffPos);
      Hg_yZ->Fill(zYSm,myWeight*theEffEle*theEffPos);
    }

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
      // if (LepGood_eleMVAPreselId[ilep]==0) continue;
      
      // ID + isolation + conv.rejection
      if (LepGood_eleMVAId[ilep]<2)      continue;   
      if (LepGood_relIso04[ilep]>0.15)   continue;  
      if (LepGood_convVetoFull[ilep]==0) continue;   // LepGood_lostHits==0 + LepGood_conv

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

    // We want 2 electrons at reco level
    if (minDrEle>0.1 || minDrPos>0.1) { continue; }
    if (minDrEle<0.1 && minDrPos<0.1 && recoEleIdx==recoPosIdx) { cout << "problem" << endl; continue; }

    // Filla gen-level
    if ( invMassSm<100 && invMassSm>80 && maxSmearPt>25 && minSmearPt>20 ) {
      Hgne_pT1->Fill(maxSmearPt,myWeight); 
      Hgne_pT2->Fill(minSmearPt,myWeight);    
      Hgne_eta1->Fill(maxSmearEta,myWeight);  
      Hgne_eta2->Fill(minSmearEta,myWeight); 
      Hgne_invMass->Fill(invMassSm,myWeight);
      Hgne_pTZ->Fill(zPtSm,myWeight);
      Hgne_pTZzoom->Fill(zPtSm,myWeight);
      Hgne_yZ->Fill(zYSm,myWeight);
    }

    // Filla reco-level
    float recoElePt  = LepGood_pt[recoEleIdx];
    float recoEleEta = LepGood_eta[recoEleIdx];
    float recoPosPt  = LepGood_pt[recoPosIdx];
    float recoPosEta = LepGood_eta[recoPosIdx];
    float maxRecoPt  = recoElePt;
    float minRecoPt  = recoPosPt;
    float maxRecoEta = recoEleEta;
    float minRecoEta = recoPosEta;
    if (recoPosPt>recoElePt) {
      maxRecoPt  = recoPosPt;
      minRecoPt  = recoElePt;
      maxRecoEta = recoPosEta;
      minRecoEta = recoEleEta;
    }
    TLorentzVector recoEle(0,0,0,0);     
    TLorentzVector recoPos(0,0,0,0);     
    recoEle.SetPtEtaPhiM(LepGood_pt[recoEleIdx],LepGood_eta[recoEleIdx],LepGood_phi[recoEleIdx],0.);
    recoPos.SetPtEtaPhiM(LepGood_pt[recoPosIdx],LepGood_eta[recoPosIdx],LepGood_phi[recoPosIdx],0.);
    float invMassReco = (recoEle+recoPos).M();
    float zPtReco     = (recoEle+recoPos).Perp();   
    float zYReco      = (recoEle+recoPos).Rapidity();   
    if ( invMassReco<100 && invMassReco>80 && maxRecoPt>25 && minRecoPt>20 ) {
      Hr_pT1->Fill(maxRecoPt,myWeightPU);   
      Hr_pT2->Fill(minRecoPt,myWeightPU);   
      Hr_eta1->Fill(maxRecoEta,myWeightPU); 
      Hr_eta2->Fill(minRecoEta,myWeightPU); 
      Hr_invMass->Fill(invMassReco,myWeightPU);
      Hr_pTZ->Fill(zPtReco,myWeightPU);
      Hr_pTZzoom->Fill(zPtReco,myWeightPU);
      Hr_yZ->Fill(zYReco,myWeightPU);
    }

  } // event loop


  // ----------------------------------------
  Hr_pT1->SetLineColor(1);
  Hg_pT1->SetLineColor(2);
  Hgne_pT1->SetLineColor(3);
  Hr_pT2->SetLineColor(1);
  Hg_pT2->SetLineColor(2);
  Hgne_pT2->SetLineColor(3);
  Hr_invMass->SetLineColor(1);
  Hg_invMass->SetLineColor(2);
  Hgne_invMass->SetLineColor(3);
  Hr_pTZ->SetLineColor(1);
  Hg_pTZ->SetLineColor(2);
  Hgne_pTZ->SetLineColor(3);
  Hr_pTZzoom->SetLineColor(1);
  Hg_pTZzoom->SetLineColor(2);
  Hgne_pTZzoom->SetLineColor(3);
  Hr_eta1->SetLineColor(1);  
  Hg_eta1->SetLineColor(2);  
  Hgne_eta1->SetLineColor(3);  
  Hr_eta2->SetLineColor(1);  
  Hg_eta2->SetLineColor(2);  
  Hgne_eta2->SetLineColor(3);  
  Hr_yZ->SetLineColor(1);
  Hg_yZ->SetLineColor(2);
  Hgne_yZ->SetLineColor(3);

  TCanvas cPt1("cPt1","",1);
  Hr_pT1->Draw();
  Hg_pT1->Draw("same");
  cPt1.SaveAs("cPt1.png");

  TCanvas cPt2("cPt2","",1);
  Hr_pT2->Draw();
  Hg_pT2->Draw("same");
  cPt2.SaveAs("cPt2.png");

  TCanvas cMee("cMee","",1);
  Hr_invMass->Draw();
  Hg_invMass->Draw("same");
  cMee.SaveAs("cMee.png");
   
  TCanvas cPtZ("cPtZ","",1);
  Hr_pTZ->Draw();
  Hg_pTZ->Draw("same");
  cPtZ.SaveAs("cPtZ.png");

  TCanvas cPtZzoom("cPtZzoom","",1);
  Hr_pTZzoom->Draw();
  Hg_pTZzoom->Draw("same");
  cPtZzoom.SaveAs("cPtZzoom.png");

  TCanvas cYZ("cYZ","",1);
  Hr_yZ->Draw();
  Hg_yZ->Draw("same");
  cYZ.SaveAs("cYZ.png");

  TCanvas cPt1ne("cPt1ne","",1);
  Hr_pT1->Draw();
  Hgne_pT1->Draw("same");
  cPt1ne.SaveAs("cPt1NoEff.png");

  TCanvas cPt2ne("cPt2ne","",1);
  Hr_pT2->Draw();
  Hgne_pT2->Draw("same");
  cPt2ne.SaveAs("cPt2NoEff.png");

  TCanvas cMeene("cMeene","",1);
  Hr_invMass->Draw();
  Hgne_invMass->Draw("same");
  cMeene.SaveAs("cMeeNoEff.png");
   
  TCanvas cPtZne("cPtZne","",1);
  Hr_pTZ->Draw();
  Hgne_pTZ->Draw("same");
  cPtZne.SaveAs("cPtZNoEff.png");

  TCanvas cPtZneZoom("cPtZneZoom","",1);
  Hr_pTZzoom->Draw();
  Hgne_pTZzoom->Draw("same");
  cPtZneZoom.SaveAs("cPtZNoEffZoom.png");

  TCanvas cYZne("cYZne","",1);
  Hr_yZ->Draw();
  Hgne_yZ->Draw("same");
  cYZne.SaveAs("cYZNoEff.png");

  TCanvas cEta1("cEta1","",1);
  Hr_eta1->Draw();
  Hg_eta1->Draw("same");
  cEta1.SaveAs("cEta1.png");

  TCanvas cEta1ne("cEta1ne","",1);
  Hr_eta1->Draw();
  Hgne_eta1->Draw("same");
  cEta1ne.SaveAs("cEta1NoEff.png");

  TCanvas cEta2("cEta2","",1);
  Hr_eta2->Draw();
  Hg_eta2->Draw("same");
  cEta2.SaveAs("cEta2.png");

  TCanvas cEta2ne("cEta2ne","",1);
  Hr_eta2->Draw();
  Hgne_eta2->Draw("same");
  cEta2ne.SaveAs("cEta2NoEff.png");

  TFile fileOut("validation.root","RECREATE");
  fileOut.cd();
  Hr_pT1->Write();
  Hr_pT2->Write();
  Hr_pTZ->Write();
  Hr_pTZzoom->Write();
  Hr_yZ->Write();
  Hr_invMass->Write();
  Hr_eta1->Write();
  Hr_eta2->Write();
  Hg_pT1->Write();
  Hg_pT2->Write();
  Hg_pTZ->Write();
  Hg_pTZzoom->Write();
  Hg_yZ->Write();
  Hg_invMass->Write();
  Hg_eta1->Write();
  Hg_eta2->Write();
  Hgne_pT1->Write();
  Hgne_pT2->Write();
  Hgne_pTZ->Write();
  Hgne_pTZzoom->Write();
  Hgne_yZ->Write();
  Hgne_invMass->Write();
  Hgne_eta1->Write();
  Hgne_eta2->Write();

  // delete

} // Loop 
