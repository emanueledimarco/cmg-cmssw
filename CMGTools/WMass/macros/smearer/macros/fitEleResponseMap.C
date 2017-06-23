#include <TROOT.h>   
#include <TSystem.h> 
#include <TCanvas.h>      
#include <TFile.h>         
#include <TH1F.h>         
#include <TH2F.h>
#include <vector>           
#include <iostream>         

// RooFit headers
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooFitResult.h"
#include "RooPlot.h"

#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"

void fitEleResponseMap()  {  

  // -------------------------------------------------------
  // Bins
  const int NPtBins  = 12;
  const int NEtaBins = 14;

  // histos to save outputs
  // bin0: underflow; bin NPtBins+1: overflow
  TH2F *DoubleSidedCBShapeParamArray_mean   = new TH2F( "DoubleSidedCBShapeParamArray_mean",   "", NPtBins, 0, NPtBins, NEtaBins, 0, NEtaBins);
  TH2F *DoubleSidedCBShapeParamArray_sigma  = new TH2F( "DoubleSidedCBShapeParamArray_sigma",  "", NPtBins, 0, NPtBins, NEtaBins, 0, NEtaBins);
  TH2F *DoubleSidedCBShapeParamArray_alphaL = new TH2F( "DoubleSidedCBShapeParamArray_alphaL", "", NPtBins, 0, NPtBins, NEtaBins, 0, NEtaBins);
  TH2F *DoubleSidedCBShapeParamArray_nL     = new TH2F( "DoubleSidedCBShapeParamArray_nL",     "", NPtBins, 0, NPtBins, NEtaBins, 0, NEtaBins);
  TH2F *DoubleSidedCBShapeParamArray_alphaR = new TH2F( "DoubleSidedCBShapeParamArray_alphaR", "", NPtBins, 0, NPtBins, NEtaBins, 0, NEtaBins);
  TH2F *DoubleSidedCBShapeParamArray_nR     = new TH2F( "DoubleSidedCBShapeParamArray_nR",     "", NPtBins, 0, NPtBins, NEtaBins, 0, NEtaBins);
  for (uint i=0; i < NPtBins+2; ++i) {        
    for (uint j=0; j < NEtaBins+2; ++j) {
      DoubleSidedCBShapeParamArray_mean->SetBinContent(i,j,1.0);
      DoubleSidedCBShapeParamArray_sigma->SetBinContent(i,j,1.0);
      DoubleSidedCBShapeParamArray_alphaL->SetBinContent(i,j,1.0);
      DoubleSidedCBShapeParamArray_nL->SetBinContent(i,j,1.0);
      DoubleSidedCBShapeParamArray_alphaR->SetBinContent(i,j,1.0);
      DoubleSidedCBShapeParamArray_nR->SetBinContent(i,j,1.0);
    }
  }

  // Read histos with resolution and fit them
  TFile *infile  = new TFile("/afs/cern.ch/work/c/crovelli/Wmass/CMSSW_5_3_22_patch1/src/CMGTools/WMass/macros/smearer/PtResolutionAndEffHist.root"); 

  // Loop over bins
  for (uint i=0; i < NPtBins+2; ++i) {
    for (uint j=0; j < NEtaBins+2; ++j) {
      TH1F *hist = (TH1F*)infile->Get(Form("ElectronsPtResolution_PtBin%d_EtaBin%d",i,j));
      cout << "New round: ptBin = " << i << ", etaBin = " << j << ", check histo: " << hist->GetTitle() << " with " << hist->GetEntries() << " entries " << endl;
      cout << "Numero di entries: " << hist->GetEntries() << ", integral " << hist->Integral() << endl;

      RooRealVar elePtRes("elePtRes","elePtRes",-0.3,0.3);
      elePtRes.setRange("range",-0.3,0.3);
      elePtRes.setBins(100);

      RooDataHist *data=0;
      data = new RooDataHist("data","data",RooArgSet(elePtRes),hist);

      RooRealVar *sigma;
      sigma  = new RooRealVar("sigma","sigma",0.01,0.0001,0.05);            // barrel
      if (j>=7) sigma = new RooRealVar("sigma","sigma",0.02,0.0001,0.07);  // endcap 

      RooRealVar *mean = new RooRealVar("mean","mean",0.05,-0.1,0.1);

      // ad hoc tuning
      if (j==7) {
	mean  = new RooRealVar("mean","mean",0.,-0.1,0.1);
	sigma = new RooRealVar("sigma","sigma",0.03,0.0001,0.07); 
      }
      if (j==9) mean  = new RooRealVar("mean","mean",0.05,-0.1,0.2);
      if (i==4 && (j==11 || j==12) ) {
	sigma = new RooRealVar("sigma","sigma",0.02,0.0001,0.07);
	mean  = new RooRealVar("mean","mean",0.,-0.1,0.1);
      }
      if (i==5 && j==8 ) mean = new RooRealVar("mean","mean",0.,-0.1,0.1);
      if (i==6 && j==12) mean = new RooRealVar("mean","mean",0.,-0.1,0.1);
      if (i==6 && j==2)  mean = new RooRealVar("mean","mean",0.,-0.1,0.1);
      if (i==8 && j==2)  mean = new RooRealVar("mean","mean",0.,-0.1,0.1);
      if (i==10 && j==2) mean = new RooRealVar("mean","mean",0.,-0.1,0.1);
      if (i==12 && j==1) mean = new RooRealVar("mean","mean",0.,-0.1,0.1);
      if (i==12 && j==2) mean = new RooRealVar("mean","mean",0.,-0.1,0.1);
      if (i==12 && j==9) mean = new RooRealVar("mean","mean",0.,-0.1,0.1);
      if (i==13 && j==2) mean = new RooRealVar("mean","mean",0.,-0.1,0.1);
      if (i==13 && j==3) mean = new RooRealVar("mean","mean",0.,-0.1,0.1);
      if (i==13 && j==4) mean = new RooRealVar("mean","mean",0.,-0.1,0.1);
      if (i==8 && j==11) sigma = new RooRealVar("sigma","sigma",0.03,0.0001,0.07);
      if (i==9 && j==12) sigma = new RooRealVar("sigma","sigma",0.03,0.0001,0.07);

      RooRealVar *alphaL = new RooRealVar("alphaL","alphaL",5,0,50);
      RooRealVar *nL     = new RooRealVar("nL","nL",3,0,30);
      if (j==7) alphaL = new RooRealVar("alphaL","alphaL",1,0,10);
      if (j==7) nL = new RooRealVar("nL","nL",5,0,50);

      RooRealVar *alphaR = new RooRealVar("alphaR","alphaR",5,0,50);
      RooRealVar *nR     = new RooRealVar("nR","nR",3,0,30);
      RooDoubleCB *model = new RooDoubleCB("ElePtResModel","ElePtResModel",elePtRes,*mean,*sigma,*alphaL,*nL,*alphaR,*nR);
                                                       
      RooFitResult *fitResult=0;
      fitResult = model->fitTo(*data, RooFit::Strategy(1), RooFit::Save());

      DoubleSidedCBShapeParamArray_mean->SetBinContent(i,j,mean->getVal());
      DoubleSidedCBShapeParamArray_sigma->SetBinContent(i,j,sigma->getVal());
      DoubleSidedCBShapeParamArray_alphaL->SetBinContent(i,j,alphaL->getVal());
      DoubleSidedCBShapeParamArray_nL->SetBinContent(i,j,nL->getVal());
      DoubleSidedCBShapeParamArray_alphaR->SetBinContent(i,j,alphaR->getVal());
      DoubleSidedCBShapeParamArray_nR->SetBinContent(i,j,nR->getVal());

      RooPlot *frame = elePtRes.frame(RooFit::Bins(100));
      data->plotOn(frame,RooFit::MarkerStyle(kFullCircle),RooFit::MarkerSize(0.8),RooFit::DrawOption("ZP"));    
      model->plotOn(frame);

      TCanvas *cv = new TCanvas("cv","cv",800,600);
      frame->Draw();
      cv->SaveAs(Form("ElePtResolution_PtBin%d_EtaBin%d.png",i,j)); 
      
      delete cv;
      delete data;
      delete mean;
      delete sigma;
      delete alphaL;
      delete nL;
      delete alphaR;
      delete nR;
      delete model;

    } // pt bins loop
  }   // eta bins loop


  // Save resolution + efficiency histos
  TH1F *NumEfficiency_PtEta = (TH1F*)infile->Get("NumEfficiency_PtEta");
  TH1F *DenEfficiency_PtEta = (TH1F*)infile->Get("DenEfficiency_PtEta");
  TH1F *Efficiency_PtEta    = (TH1F*)infile->Get("Efficiency_PtEta");

  TFile *outfile = new TFile("PtResolutionAndEffHistOUT.root","RECREATE");
  outfile->cd();
  NumEfficiency_PtEta->Write();
  DenEfficiency_PtEta->Write();
  Efficiency_PtEta->Write();
  outfile->WriteTObject(DoubleSidedCBShapeParamArray_mean,   "DoubleSidedCBShapeParamArray_mean",   "WriteDelete");
  outfile->WriteTObject(DoubleSidedCBShapeParamArray_sigma,  "DoubleSidedCBShapeParamArray_sigma",  "WriteDelete");
  outfile->WriteTObject(DoubleSidedCBShapeParamArray_alphaL, "DoubleSidedCBShapeParamArray_alphaL", "WriteDelete");
  outfile->WriteTObject(DoubleSidedCBShapeParamArray_nL,     "DoubleSidedCBShapeParamArray_nL",     "WriteDelete");
  outfile->WriteTObject(DoubleSidedCBShapeParamArray_alphaR, "DoubleSidedCBShapeParamArray_alphaR", "WriteDelete");
  outfile->WriteTObject(DoubleSidedCBShapeParamArray_nR,     "DoubleSidedCBShapeParamArray_nR",     "WriteDelete");
  delete outfile;

  delete infile;

  return;
}

