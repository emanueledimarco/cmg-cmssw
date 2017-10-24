#include <stdio.h>
#include <cmath>
#include <iostream>

//ROOT header files
#include <TROOT.h>
#include <TAttFill.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TVector3.h>
#include <TVirtualFitter.h>

using namespace std; 

float drawTH1pair(TH1* h1, TH1* h2, 
		 const string& xAxisNameTmp = "", const string& yAxisName = "Events", float lumi=-1, const string& canvasName = "default", 
		 const string& outputDIR = "./", int mycolor=2,
		 const string& legEntry1 = "data", const string& legEntry2 = "MC", const string& ratioPadYaxisName = "data/MC") 
{

  h1->Sumw2();
  h2->Sumw2();

  string xAxisName = "";
  string separator = "::";
  Bool_t setXAxisRangeFromUser = false;
  Double_t xmin = 0;
  Double_t xmax = 0;

  size_t pos = xAxisNameTmp.find(separator);
  if (pos != string::npos) {
    string xrange = "";
    setXAxisRangeFromUser = true;
    xAxisName.assign(xAxisNameTmp, 0, pos); 
    xrange.assign(xAxisNameTmp, pos + separator.size(), string::npos);
    separator = ",";
    pos = xrange.find(separator);
    string numString = ""; 
    numString.assign(xrange,0,pos);
    xmin = std::stod(numString);
    numString.assign(xrange,pos + separator.size(), string::npos);
    xmax = std::stod(numString);
  } else {
    xAxisName = xAxisNameTmp;
  }

  if (yAxisName == "a.u.") {
    h1->Scale(1./h1->Integral());
    h2->Scale(1./h2->Integral());
  }
  else if (lumi>-1) {
    h1->Scale(lumi/h1->Integral());
    h2->Scale(lumi/h2->Integral());
  }

  h1->SetStats(0);
  h2->SetStats(0);

  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);

  TH1* frame =  (TH1*) h1->Clone("frame");
  frame->SetTitle("");
  frame->GetXaxis()->SetLabelSize(0.04);
  frame->SetStats(0);

  h1->SetLineColor(kBlack);
  h1->SetMarkerColor(kBlack);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);

  h1->SetTitle("");
  h1->GetXaxis()->SetLabelSize(0);
  h1->GetYaxis()->SetTitle(yAxisName.c_str());
  h1->GetYaxis()->SetTitleOffset(1.1);
  h1->GetYaxis()->SetTitleSize(0.05);
  h1->GetYaxis()->SetRangeUser(0.0, max(h1->GetMaximum(),h2->GetMaximum()) * 1.2);
  if (setXAxisRangeFromUser) h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->Draw("EP");

  h2->SetTitle("");
  h2->SetLineColor(mycolor);
  h2->SetLineWidth(2);
  h2->Draw("hist E same");

  TLegend leg (0.55,0.7,0.95,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(h1,legEntry1.c_str(),"PLE");
  leg.AddEntry(h2,legEntry2.c_str(),"L");
  leg.Draw("same");
  canvas->RedrawAxis("sameaxis");

  pad2->Draw();
  pad2->cd();

  frame->Reset("ICES");
  frame->GetYaxis()->SetRangeUser(0.8,1.2);
  frame->GetYaxis()->SetNdivisions(5);
  frame->GetYaxis()->SetTitle(ratioPadYaxisName.c_str());
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetTitle(xAxisName.c_str());
  if (setXAxisRangeFromUser) frame->GetXaxis()->SetRangeUser(xmin,xmax);
  frame->GetXaxis()->SetTitleSize(0.05);

  TH1D* ratio = (TH1D*) h1->Clone("ratio");
  TH1D* den_noerr = (TH1D*) h2->Clone("den_noerr");
  TH1D* den = (TH1D*) h2->Clone("den");
  for(int iBin = 1; iBin < den->GetNbinsX()+1; iBin++)
    den_noerr->SetBinError(iBin,0.);

  ratio->Divide(den_noerr);
  den->Divide(den_noerr);
  den->SetFillColor(kGray);
  frame->Draw();
  ratio->SetMarkerSize(0.85);
  ratio->Draw("EPsame");
  den->Draw("E2same");

  TF1* line = new TF1("horiz_line","1",ratio->GetXaxis()->GetBinLowEdge(1),ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()+1));
  line->SetLineColor(mycolor);
  line->SetLineWidth(2);
  line->Draw("Lsame");
  ratio->Draw("EPsame");
  pad2->RedrawAxis("sameaxis");

  // Calculate chi2                                                                                        
  double chi2 = h1->Chi2Test(h2,"CHI2/NDF WW");
  TLegend leg2 (0.14,0.25,0.32,0.28,NULL,"brNDC");
  leg2.SetFillColor(0);
  leg2.SetFillStyle(1);
  leg2.SetBorderSize(0);
  leg2.SetLineColor(0);
  leg2.AddEntry((TObject*)0,Form("#chi^{2}/ndf = %.2f",chi2),"");
  leg2.Draw("same");

  canvas->SaveAs((outputDIR + canvasName + ".png").c_str());
  canvas->SaveAs((outputDIR + canvasName + ".pdf").c_str());
  canvas->SaveAs((outputDIR + canvasName + ".root").c_str());

  if (yAxisName == "a.u.") h1->GetYaxis()->SetRangeUser(max(0.0001,min(h1->GetMinimum(),h2->GetMinimum())*0.8),max(h1->GetMaximum(),h2->GetMaximum())*100);
  else h1->GetYaxis()->SetRangeUser(max(0.001,min(h1->GetMinimum(),h2->GetMinimum())*0.8),max(h1->GetMaximum(),h2->GetMaximum())*100);
  canvas->SetLogy(0);

  delete canvas;
  frame->Reset("ICES");

  return chi2;
}

/*
void validationPlot() {

  // MC for closure test: MC vs smearing
  TFile *fileIN = new TFile("rootFiles/new/validation.root");

  TH1F *Hr_invMass = (TH1F*)fileIN->Get("Hr_invMass");
  TH1F *Hr_pT1 = (TH1F*)fileIN->Get("Hr_pT1");
  TH1F *Hr_pT2 = (TH1F*)fileIN->Get("Hr_pT2");
  TH1F *Hr_pTZ = (TH1F*)fileIN->Get("Hr_pTZ");
  TH1F *Hr_pTZzoom = (TH1F*)fileIN->Get("Hr_pTZzoom");
  TH1F *Hr_yZ   = (TH1F*)fileIN->Get("Hr_yZ");
  TH1F *Hr_eta1 = (TH1F*)fileIN->Get("Hr_eta1");
  TH1F *Hr_eta2 = (TH1F*)fileIN->Get("Hr_eta2");

  TH1F *Hg_invMass = (TH1F*)fileIN->Get("Hg_invMass");
  TH1F *Hg_pT1 = (TH1F*)fileIN->Get("Hg_pT1");
  TH1F *Hg_pT2 = (TH1F*)fileIN->Get("Hg_pT2");
  TH1F *Hg_pTZ = (TH1F*)fileIN->Get("Hg_pTZ");
  TH1F *Hg_pTZzoom = (TH1F*)fileIN->Get("Hg_pTZzoom");
  TH1F *Hg_yZ   = (TH1F*)fileIN->Get("Hg_yZ");
  TH1F *Hg_eta1 = (TH1F*)fileIN->Get("Hg_eta1");
  TH1F *Hg_eta2 = (TH1F*)fileIN->Get("Hg_eta2");

  TH1F *Hgne_invMass = (TH1F*)fileIN->Get("Hgne_invMass");
  TH1F *Hgne_pT1 = (TH1F*)fileIN->Get("Hgne_pT1");
  TH1F *Hgne_pT2 = (TH1F*)fileIN->Get("Hgne_pT2");
  TH1F *Hgne_pTZ = (TH1F*)fileIN->Get("Hgne_pTZ");
  TH1F *Hgne_pTZzoom = (TH1F*)fileIN->Get("Hgne_pTZzoom");
  TH1F *Hgne_yZ   = (TH1F*)fileIN->Get("Hgne_yZ");
  TH1F *Hgne_eta1 = (TH1F*)fileIN->Get("Hgne_eta1");
  TH1F *Hgne_eta2 = (TH1F*)fileIN->Get("Hgne_eta2");

  // Data for data vs MC o vs smearing
  TFile *fileINdata = new TFile("rootFiles/new/validationData.root");
  TH1F *HrData_pT1 = (TH1F*)fileINdata->Get("Hr_pT1");
  TH1F *HrData_pT2 = (TH1F*)fileINdata->Get("Hr_pT2");
  TH1F *HrData_eta1 = (TH1F*)fileINdata->Get("Hr_eta1");
  TH1F *HrData_eta2 = (TH1F*)fileINdata->Get("Hr_eta2");
  TH1F *HrData_invMass = (TH1F*)fileINdata->Get("Hr_invMass");
  TH1F *HrData_pTZ = (TH1F*)fileINdata->Get("Hr_pTZ");
  TH1F *HrData_pTZzoom = (TH1F*)fileINdata->Get("Hr_pTZzoom");
  TH1F *HrData_phiStarZ = (TH1F*)fileINdata->Get("Hr_phiStarZ");
  TH1F *HrData_cosThetaStarZ = (TH1F*)fileINdata->Get("Hr_cosThetaStarZ");
  TH1F *HrData_yZ   = (TH1F*)fileINdata->Get("Hr_yZ");
  TH1F *HrDataEB_pT1 = (TH1F*)fileINdata->Get("HrEB_pT1");
  TH1F *HrDataEB_pT2 = (TH1F*)fileINdata->Get("HrEB_pT2");
  TH1F *HrDataEE_pT1 = (TH1F*)fileINdata->Get("HrEE_pT1");
  TH1F *HrDataEE_pT2 = (TH1F*)fileINdata->Get("HrEE_pT2");
  TH1F *HrDataEBEB_invMass = (TH1F*)fileINdata->Get("HrEBEB_invMass");
  TH1F *HrDataNotEBEB_invMass = (TH1F*)fileINdata->Get("HrNotEBEB_invMass");
  TH1F *HrDataEBZ_invMass = (TH1F*)fileINdata->Get("HrEBZ_invMass");
  TH1F *HrDataEBZ_pTZ = (TH1F*)fileINdata->Get("HrEBZ_pTZ");
  TH1F *HrDataEBZ_pTZzoom = (TH1F*)fileINdata->Get("HrEBZ_pTZzoom");
  TH1F *HrDataEBZ_phiStarZ = (TH1F*)fileINdata->Get("HrEBZ_phiStarZ");
  TH1F *HrDataEBZ_cosThetaStarZ = (TH1F*)fileINdata->Get("HrEBZ_cosThetaStarZ");
  TH1F *HrDataEEZ_invMass = (TH1F*)fileINdata->Get("HrEEZ_invMass");
  TH1F *HrDataEEZ_pTZ = (TH1F*)fileINdata->Get("HrEEZ_pTZ");
  TH1F *HrDataEEZ_pTZzoom = (TH1F*)fileINdata->Get("HrEEZ_pTZzoom");
  TH1F *HrDataEEZ_phiStarZ = (TH1F*)fileINdata->Get("HrEEZ_phiStarZ");
  TH1F *HrDataEEZ_cosThetaStarZ = (TH1F*)fileINdata->Get("HrEEZ_cosThetaStarZ");


  /*
  // Full MC for data vs MC
  TFile *fileINmc = new TFile("rootFiles/new/validationMC.root");
  TH1F *HrMc_pT1 = (TH1F*)fileINmc->Get("Hr_pT1");
  TH1F *HrMc_pT2 = (TH1F*)fileINmc->Get("Hr_pT2");
  TH1F *HrMc_eta1 = (TH1F*)fileINmc->Get("Hr_eta1");
  TH1F *HrMc_eta2 = (TH1F*)fileINmc->Get("Hr_eta2");
  TH1F *HrMc_invMass = (TH1F*)fileINmc->Get("Hr_invMass");
  TH1F *HrMc_pTZ = (TH1F*)fileINmc->Get("Hr_pTZ");
  TH1F *HrMc_pTZzoom = (TH1F*)fileINmc->Get("Hr_pTZzoom");
  TH1F *HrMc_yZ   = (TH1F*)fileINmc->Get("Hr_yZ");
  TH1F *HrMcEB_pT1 = (TH1F*)fileINmc->Get("HrEB_pT1");
  TH1F *HrMcEB_pT2 = (TH1F*)fileINmc->Get("HrEB_pT2");
  TH1F *HrMcEE_pT1 = (TH1F*)fileINmc->Get("HrEE_pT1");
  TH1F *HrMcEE_pT2 = (TH1F*)fileINmc->Get("HrEE_pT2");
  TH1F *HrEBEBMc_invMass = (TH1F*)fileINmc->Get("HrEBEB_invMass");
  TH1F *HrEBEBMc_pTZ = (TH1F*)fileINmc->Get("HrEBEB_pTZ");
  TH1F *HrEBEBMc_pTZzoom = (TH1F*)fileINmc->Get("HrEBEB_pTZzoom");
  TH1F *HrNotEBEBMc_invMass = (TH1F*)fileINmc->Get("HrNotEBEB_invMass");
  TH1F *HrNotEBEBMc_pTZ = (TH1F*)fileINmc->Get("HrNotEBEB_pTZ");
  TH1F *HrNotEBEBMc_pTZzoom = (TH1F*)fileINmc->Get("HrNotEBEB_pTZzoom");
  */

/*
  // Smeared quantities, to be used with whatever generator
  TFile *fileINSm = new TFile("smearedAppliedToZPedro__ZJPtsqmin4.root");
  TH1F *HsmMc_pT1  = (TH1F*)fileINSm->Get("Hg_pT1");
  TH1F *HsmMc_pT2  = (TH1F*)fileINSm->Get("Hg_pT2");
  TH1F *HsmMc_eta1 = (TH1F*)fileINSm->Get("Hg_eta1");
  TH1F *HsmMc_eta2 = (TH1F*)fileINSm->Get("Hg_eta2");
  TH1F *HsmMc_invMass  = (TH1F*)fileINSm->Get("Hg_invMass");
  TH1F *HsmMc_pTZ      = (TH1F*)fileINSm->Get("Hg_pTZ");
  TH1F *HsmMc_pTZzoom  = (TH1F*)fileINSm->Get("Hg_pTZzoom");
  TH1F *HsmMc_phiStarZ = (TH1F*)fileINSm->Get("Hg_phiStarZ");
  TH1F *HsmMc_cosThetaStarZ = (TH1F*)fileINSm->Get("Hg_cosThetaStarZ");
  TH1F *HsmMc_yZ    = (TH1F*)fileINSm->Get("Hg_yZ");
  TH1F *HsmMcEB_pT1 = (TH1F*)fileINSm->Get("HgEB_pT1");
  TH1F *HsmMcEB_pT2 = (TH1F*)fileINSm->Get("HgEB_pT2");
  TH1F *HsmMcEE_pT1 = (TH1F*)fileINSm->Get("HgEE_pT1");
  TH1F *HsmMcEE_pT2 = (TH1F*)fileINSm->Get("HgEE_pT2");
  TH1F *HsmMcEBEB_invMass = (TH1F*)fileINSm->Get("HgEBEB_invMass");
  TH1F *HsmMcNotEBEB_invMass = (TH1F*)fileINSm->Get("HgNotEBEB_invMass");
  TH1F *HsmMcEBZ_invMass = (TH1F*)fileINSm->Get("HgEBZ_invMass");
  TH1F *HsmMcEBZ_pTZ = (TH1F*)fileINSm->Get("HgEBZ_pTZ");
  TH1F *HsmMcEBZ_pTZzoom = (TH1F*)fileINSm->Get("HgEBZ_pTZzoom");
  TH1F *HsmMcEBZ_phiStarZ = (TH1F*)fileINSm->Get("HgEBZ_phiStarZ");
  TH1F *HsmMcEBZ_cosThetaStarZ = (TH1F*)fileINSm->Get("HgEBZ_cosThetaStarZ");
  TH1F *HsmMcEEZ_invMass = (TH1F*)fileINSm->Get("HgEEZ_invMass");
  TH1F *HsmMcEEZ_pTZ = (TH1F*)fileINSm->Get("HgEEZ_pTZ");
  TH1F *HsmMcEEZ_pTZzoom = (TH1F*)fileINSm->Get("HgEEZ_pTZzoom");
  TH1F *HsmMcEEZ_phiStarZ = (TH1F*)fileINSm->Get("HgEEZ_phiStarZ");
  TH1F *HsmMcEEZ_cosThetaStarZ = (TH1F*)fileINSm->Get("HgEEZ_cosThetaStarZ");
*/

  /*
  // closure of method, Mc vs smearing: plots with efficiency
  drawTH1pair(Hr_pT1, Hg_pT1, "leading pT [GeV]::25::80","a.u.",-1,"pt1WithEff","./",2,"MC","smeared",""); 
  drawTH1pair(Hr_pT2, Hg_pT2, "subleading pT [GeV]::20::70","a.u.",-1,"pt2WithEff","./",2,"MC","smeared",""); 
  drawTH1pair(Hr_pTZ, Hg_pTZ, "Z pT [GeV]::0::80","a.u.",-1,"ptZWithEff","./",2,"MC","smeared",""); 
  drawTH1pair(Hr_pTZzoom, Hg_pTZzoom, "Z pT [GeV]::0::30","a.u.",-1,"ptZzoomWithEff","./",2,"MC","smeared",""); 
  drawTH1pair(Hr_invMass, Hg_invMass, "m(ee) [GeV]::80::100","a.u.",-1,"meeWithEff","./",2,"MC","smeared",""); 
  drawTH1pair(Hr_eta1, Hg_eta1, "leading eta::-2.5::2.5","a.u.",-1,"eta1WithEff","./",2,"MC","smeared",""); 
  drawTH1pair(Hr_eta2, Hg_eta2, "subleading eta::-2.5::2.5","a.u.",-1,"eta2WithEff","./",2,"MC","smeared",""); 
  drawTH1pair(Hr_yZ, Hg_yZ, "Z rapidity::-3.::3.","a.u.",-1,"yZWithEff","./",2,"MC","smeared",""); 
  */

  /*
  // closure of method, Mc vs smearing: plots wo efficiency
  drawTH1pair(Hr_pT1, Hgne_pT1, "leading pT [GeV]::25::80","a.u.",-1,"pt1NoEff","./",3,"MC","smeared",""); 
  drawTH1pair(Hr_pT2, Hgne_pT2, "subleading pT [GeV]::20::70","a.u.",-1,"pt2NoEff","./",3,"MC","smeared",""); 
  drawTH1pair(Hr_pTZ, Hgne_pTZ, "Z pT [GeV]::0::80","a.u.",-1,"ptZNoEff","./",3,"MC","smeared",""); 
  drawTH1pair(Hr_pTZzoom, Hgne_pTZzoom, "Z pT [GeV]::0::30","a.u.",-1,"ptZzoomNoEff","./",3,"MC","smeared",""); 
  drawTH1pair(Hr_invMass, Hgne_invMass, "m(ee) [GeV]::80::100","a.u.",-1,"meeNoEff","./",3,"MC","smeared",""); 
  drawTH1pair(Hr_eta1, Hgne_eta1, "leading eta::-2.5::2.5","a.u.",-1,"eta1NoEff","./",3,"MC","smeared",""); 
  drawTH1pair(Hr_eta2, Hgne_eta2, "subleading eta::-2.5::2.5","a.u.",-1,"eta2NoEff","./",3,"MC","smeared",""); 
  drawTH1pair(Hr_yZ, Hgne_yZ, "Z rapidity::-3.::3.","a.u.",-1,"yZNoEff","./",3,"MC","smeared",""); 
  */

/*
  // data/MC 
  drawTH1pair(HrData_pT1, HrMc_pT1, "leading pT [GeV]::25::80","a.u.",19700.,"pt1dataMc","./",4,"Data","MC",""); 
  drawTH1pair(HrData_pT2, HrMc_pT2, "subleading pT [GeV]::20::70","a.u.",19700.,"pt2dataMc","./",4,"Data","MC",""); 
  drawTH1pair(HrData_eta1, HrMc_eta1, "leading eta::-2.5::2.5","a.u.",19700.,"eta1dataMc","./",4,"Data","MC",""); 
  drawTH1pair(HrData_eta2, HrMc_eta2, "subleading eta::-2.5::2.5","a.u.",19700.,"eta2dataMc","./",4,"Data","MC",""); 
  drawTH1pair(HrData_pTZ, HrMc_pTZ, "Z pT [GeV]::0::80","a.u.",19700.,"ptZdataMc","./",4,"Data","MC",""); 
  drawTH1pair(HrData_pTZzoom, HrMc_pTZzoom, "Z pT [GeV]::0::30","a.u.",19700.,"ptZzoomdataMc","./",4,"Data","MC",""); 
  drawTH1pair(HrData_invMass, HrMc_invMass, "m(ee) [GeV]::80::100","a.u.",19700.,"meedataMc","./",4,"Data","MC",""); 
  drawTH1pair(HrData_yZ, HrMc_yZ, "Z rapidity::-3.::3.","a.u.",19700.,"yZdataMc","./",4,"Data","MC",""); 
}
*/

void validationPlot() {

  // Data 
  TFile *fileINdata = new TFile("rootFiles/new/validationData.root");
  TH1F *HrData_pT1 = (TH1F*)fileINdata->Get("Hr_pT1");
  TH1F *HrData_pT2 = (TH1F*)fileINdata->Get("Hr_pT2");
  TH1F *HrData_eta1 = (TH1F*)fileINdata->Get("Hr_eta1");
  TH1F *HrData_eta2 = (TH1F*)fileINdata->Get("Hr_eta2");
  TH1F *HrData_invMass = (TH1F*)fileINdata->Get("Hr_invMass");
  TH1F *HrData_pTZ = (TH1F*)fileINdata->Get("Hr_pTZ");
  TH1F *HrData_pTZzoom  = (TH1F*)fileINdata->Get("Hr_pTZzoom");
  TH1F *HrData_pTZzoom1 = (TH1F*)fileINdata->Get("Hr_pTZzoom1");
  TH1F *HrData_pTZzoom2 = (TH1F*)fileINdata->Get("Hr_pTZzoom2");
  TH1F *HrData_phiStarZ = (TH1F*)fileINdata->Get("Hr_phiStarZ");
  TH1F *HrData_cosThetaStarZ = (TH1F*)fileINdata->Get("Hr_cosThetaStarZ");
  TH1F *HrData_yZ   = (TH1F*)fileINdata->Get("Hr_yZ");
  TH1F *HrDataEB_pT1 = (TH1F*)fileINdata->Get("HrEB_pT1");
  TH1F *HrDataEB_pT2 = (TH1F*)fileINdata->Get("HrEB_pT2");
  TH1F *HrDataEE_pT1 = (TH1F*)fileINdata->Get("HrEE_pT1");
  TH1F *HrDataEE_pT2 = (TH1F*)fileINdata->Get("HrEE_pT2");
  TH1F *HrDataEBEB_invMass = (TH1F*)fileINdata->Get("HrEBEB_invMass");
  TH1F *HrDataNotEBEB_invMass = (TH1F*)fileINdata->Get("HrNotEBEB_invMass");
  TH1F *HrDataEBZ_invMass = (TH1F*)fileINdata->Get("HrEBZ_invMass");
  TH1F *HrDataEBZ_pTZ = (TH1F*)fileINdata->Get("HrEBZ_pTZ");
  TH1F *HrDataEBZ_pTZzoom  = (TH1F*)fileINdata->Get("HrEBZ_pTZzoom");
  TH1F *HrDataEBZ_pTZzoom1 = (TH1F*)fileINdata->Get("HrEBZ_pTZzoom1");
  TH1F *HrDataEBZ_pTZzoom2 = (TH1F*)fileINdata->Get("HrEBZ_pTZzoom2");
  TH1F *HrDataEBZ_phiStarZ = (TH1F*)fileINdata->Get("HrEBZ_phiStarZ");
  TH1F *HrDataEBZ_cosThetaStarZ = (TH1F*)fileINdata->Get("HrEBZ_cosThetaStarZ");
  TH1F *HrDataEEZ_invMass = (TH1F*)fileINdata->Get("HrEEZ_invMass");
  TH1F *HrDataEEZ_pTZ = (TH1F*)fileINdata->Get("HrEEZ_pTZ");
  TH1F *HrDataEEZ_pTZzoom  = (TH1F*)fileINdata->Get("HrEEZ_pTZzoom");
  TH1F *HrDataEEZ_pTZzoom1 = (TH1F*)fileINdata->Get("HrEEZ_pTZzoom1");
  TH1F *HrDataEEZ_pTZzoom2 = (TH1F*)fileINdata->Get("HrEEZ_pTZzoom2");
  TH1F *HrDataEEZ_phiStarZ = (TH1F*)fileINdata->Get("HrEEZ_phiStarZ");
  TH1F *HrDataEEZ_cosThetaStarZ = (TH1F*)fileINdata->Get("HrEEZ_cosThetaStarZ");

  // Smeared quantities, to be used with whatever generator
  //TFile *fileINSm = new TFile("rootFiles/new/smearedAppliedToZPedro__ZJcentral.root");
  //TFile *fileINSm = new TFile("rootFiles/new/smearedAppliedToZPedro__DYToMuMu_M_50_TuneAZ_8TeV_ATLAS.root");
  //TFile *fileINSm = new TFile("rootFiles/new/smearedAppliedToZPedro__ZJPtsqmin4.root");
  TFile *fileINSm = new TFile("rootFiles/new/smearedAppliedToZPedro__ZJcentral__w47.root"); 
  //TFile *fileINSm = new TFile("rootFiles/new/smearedAppliedToZPedro__ZJPtsqmin4__w48.root");
  TH1F *HsmMc_pT1  = (TH1F*)fileINSm->Get("Hg_pT1");
  TH1F *HsmMc_pT2  = (TH1F*)fileINSm->Get("Hg_pT2");
  TH1F *HsmMc_eta1 = (TH1F*)fileINSm->Get("Hg_eta1");
  TH1F *HsmMc_eta2 = (TH1F*)fileINSm->Get("Hg_eta2");
  TH1F *HsmMc_invMass  = (TH1F*)fileINSm->Get("Hg_invMass");
  TH1F *HsmMc_pTZ      = (TH1F*)fileINSm->Get("Hg_pTZ");
  TH1F *HsmMc_pTZzoom  = (TH1F*)fileINSm->Get("Hg_pTZzoom");
  TH1F *HsmMc_pTZzoom1 = (TH1F*)fileINSm->Get("Hg_pTZzoom1");
  TH1F *HsmMc_pTZzoom2 = (TH1F*)fileINSm->Get("Hg_pTZzoom2");
  TH1F *HsmMc_phiStarZ = (TH1F*)fileINSm->Get("Hg_phiStarZ");
  TH1F *HsmMc_cosThetaStarZ = (TH1F*)fileINSm->Get("Hg_cosThetaStarZ");
  TH1F *HsmMc_yZ    = (TH1F*)fileINSm->Get("Hg_yZ");
  TH1F *HsmMcEB_pT1 = (TH1F*)fileINSm->Get("HgEB_pT1");
  TH1F *HsmMcEB_pT2 = (TH1F*)fileINSm->Get("HgEB_pT2");
  TH1F *HsmMcEE_pT1 = (TH1F*)fileINSm->Get("HgEE_pT1");
  TH1F *HsmMcEE_pT2 = (TH1F*)fileINSm->Get("HgEE_pT2");
  TH1F *HsmMcEBEB_invMass = (TH1F*)fileINSm->Get("HgEBEB_invMass");
  TH1F *HsmMcNotEBEB_invMass = (TH1F*)fileINSm->Get("HgNotEBEB_invMass");
  TH1F *HsmMcEBZ_invMass = (TH1F*)fileINSm->Get("HgEBZ_invMass");
  TH1F *HsmMcEBZ_pTZ = (TH1F*)fileINSm->Get("HgEBZ_pTZ");
  TH1F *HsmMcEBZ_pTZzoom  = (TH1F*)fileINSm->Get("HgEBZ_pTZzoom");
  TH1F *HsmMcEBZ_pTZzoom1 = (TH1F*)fileINSm->Get("HgEBZ_pTZzoom1");
  TH1F *HsmMcEBZ_pTZzoom2 = (TH1F*)fileINSm->Get("HgEBZ_pTZzoom2");
  TH1F *HsmMcEBZ_phiStarZ = (TH1F*)fileINSm->Get("HgEBZ_phiStarZ");
  TH1F *HsmMcEBZ_cosThetaStarZ = (TH1F*)fileINSm->Get("HgEBZ_cosThetaStarZ");
  TH1F *HsmMcEEZ_invMass = (TH1F*)fileINSm->Get("HgEEZ_invMass");
  TH1F *HsmMcEEZ_pTZ = (TH1F*)fileINSm->Get("HgEEZ_pTZ");
  TH1F *HsmMcEEZ_pTZzoom  = (TH1F*)fileINSm->Get("HgEEZ_pTZzoom");
  TH1F *HsmMcEEZ_pTZzoom1 = (TH1F*)fileINSm->Get("HgEEZ_pTZzoom1");
  TH1F *HsmMcEEZ_pTZzoom2 = (TH1F*)fileINSm->Get("HgEEZ_pTZzoom2");
  TH1F *HsmMcEEZ_phiStarZ = (TH1F*)fileINSm->Get("HgEEZ_phiStarZ");
  TH1F *HsmMcEEZ_cosThetaStarZ = (TH1F*)fileINSm->Get("HgEEZ_cosThetaStarZ");

  TH1F *HsmMc_pTZzoom_ww[282];
  TH1F *HsmMc_pTZzoom1_ww[282];
  TH1F *HsmMc_pTZzoom2_ww[282];
  for (int ii=0; ii<282; ii++) { 
    TString thisii = Form("%d",ii);
    TString thename  = TString("Hg_pTZzoom_ww[")+thisii+TString("]");
    TString thename1 = TString("Hg_pTZzoom1_ww[")+thisii+TString("]");
    TString thename2 = TString("Hg_pTZzoom2_ww[")+thisii+TString("]");
    HsmMc_pTZzoom_ww[ii]  = (TH1F*)fileINSm->Get(thename);
    HsmMc_pTZzoom1_ww[ii] = (TH1F*)fileINSm->Get(thename1);
    HsmMc_pTZzoom2_ww[ii] = (TH1F*)fileINSm->Get(thename2);
  }

  // chi2 study
  float theweight[282];
  float theweight1[282];
  float theweight2[282];
  float chi2_ptzoom[282];
  float chi2_ptzoom1[282];
  float chi2_ptzoom2[282];
  theweight[0] = 0;

  drawTH1pair(HrData_pT1,  HsmMc_pT1,  "leading pT [GeV]::25::80","a.u.",19700.,"pt1_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrData_pT2,  HsmMc_pT2,  "subleading pT [GeV]::20::70","a.u.",19700.,"pt2_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrData_eta1, HsmMc_eta1, "leading eta::-2.5::2.5","a.u.",19700.,"eta1_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrData_eta2, HsmMc_eta2, "subleading eta::-2.5::2.5","a.u.",19700.,"eta2_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrData_invMass, HsmMc_invMass, "m(ee) [GeV]::80::100","a.u.",19700.,"mee_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrData_yZ, HsmMc_yZ, "Z rapidity::-3.::3.","a.u.",19700.,"yZ_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrData_phiStarZ, HsmMc_phiStarZ, "phi*::0::3.14", "a.u.",19700.,"phiStar_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrData_cosThetaStarZ, HsmMc_cosThetaStarZ, "cos(theta*)::-1::1", "a.u.",19700.,"cosThetaStar_central","./",7,"Data","central, w47",""); 
  //
  drawTH1pair(HrData_pTZ,  HsmMc_pTZ,  "Z pT [GeV]::0::80","a.u.",19700.,"ptZ_central","./",7,"Data","central, w47",""); 
  //
  //drawTH1pair(HrData_pTZzoom,  HsmMc_pTZzoom,  "Z pT [GeV]::0::30","a.u.",19700.,"ptZzoom_central", "./",7,"Data","central, w47",""); 
  //drawTH1pair(HrData_pTZzoom1, HsmMc_pTZzoom1, "Z pT [GeV]::0::20","a.u.",19700.,"ptZzoom1_central","./",7,"Data","central, w47",""); 
  //drawTH1pair(HrData_pTZzoom2, HsmMc_pTZzoom2, "Z pT [GeV]::0::10","a.u.",19700.,"ptZzoom2_central","./",7,"Data","central, w47",""); 
  chi2_ptzoom[0]  = drawTH1pair(HrData_pTZzoom,  HsmMc_pTZzoom,  "Z pT [GeV]::0::30","a.u.",19700.,"ptZzoom_central", "./",7,"Data","central, w47",""); 
  chi2_ptzoom1[0] = drawTH1pair(HrData_pTZzoom1, HsmMc_pTZzoom1, "Z pT [GeV]::0::20","a.u.",19700.,"ptZzoom1_central","./",7,"Data","central, w47",""); 
  chi2_ptzoom2[0] = drawTH1pair(HrData_pTZzoom2, HsmMc_pTZzoom2, "Z pT [GeV]::0::10","a.u.",19700.,"ptZzoom2_central","./",7,"Data","central, w47",""); 
  //
  drawTH1pair(HrDataEB_pT1,  HsmMcEB_pT1,  "EB, leading pT [GeV]::25::80","a.u.",19700.,"pt1EB_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrDataEB_pT2,  HsmMcEB_pT2,  "EB, subleading pT [GeV]::20::70","a.u.",19700.,"pt2EB_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrDataEE_pT1,  HsmMcEE_pT1,  "EE, leading pT [GeV]::25::80","a.u.",19700.,"pt1EE_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrDataEE_pT2,  HsmMcEE_pT2,  "EE, subleading pT [GeV]::20::70","a.u.",19700.,"pt2EE_central","./",7,"Data","central, w47",""); 
  //
  drawTH1pair(HrDataEBEB_invMass, HsmMcEBEB_invMass, "m(ee) [GeV]::80::100","a.u.",19700.,"meeEBEB_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrDataNotEBEB_invMass, HsmMcNotEBEB_invMass, "m(ee) [GeV]::80::100","a.u.",19700.,"meeNotEBEB_central","./",7,"Data","central, w47",""); 
  //
  drawTH1pair(HrDataEBZ_invMass, HsmMcEBZ_invMass, "m(ee) [GeV]::80::100","a.u.",19700.,"meeEBZ_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrDataEEZ_invMass, HsmMcEEZ_invMass, "m(ee) [GeV]::80::100","a.u.",19700.,"meeEEZ_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrDataEBZ_pTZ,     HsmMcEBZ_pTZ,     "Z pT [GeV]::0::80",    "a.u.",19700.,"ptZEBZ_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrDataEEZ_pTZ,     HsmMcEEZ_pTZ,     "Z pT [GeV]::0::80",    "a.u.",19700.,"ptZEEZ_central","./",7,"Data","central, w47",""); 

  drawTH1pair(HrDataEBZ_pTZzoom,  HsmMcEBZ_pTZzoom,  "Z pT [GeV]::0::30", "a.u.",19700.,"ptZzoomEBZ_central", "./",7,"Data","central, w47",""); 
  drawTH1pair(HrDataEEZ_pTZzoom,  HsmMcEEZ_pTZzoom,  "Z pT [GeV]::0::30", "a.u.",19700.,"ptZzoomEEZ_central", "./",7,"Data","central, w47",""); 
  drawTH1pair(HrDataEBZ_pTZzoom1, HsmMcEBZ_pTZzoom1, "Z pT [GeV]::0::20", "a.u.",19700.,"ptZzoom1EBZ_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrDataEEZ_pTZzoom1, HsmMcEEZ_pTZzoom1, "Z pT [GeV]::0::20", "a.u.",19700.,"ptZzoom1EEZ_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrDataEBZ_pTZzoom2, HsmMcEBZ_pTZzoom2, "Z pT [GeV]::0::10", "a.u.",19700.,"ptZzoom2EBZ_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrDataEEZ_pTZzoom2, HsmMcEEZ_pTZzoom2, "Z pT [GeV]::0::10", "a.u.",19700.,"ptZzoom2EEZ_central","./",7,"Data","central, w47",""); 

  drawTH1pair(HrDataEBZ_invMass, HsmMcEBZ_invMass, "m(ee) [GeV]::80::100", "a.u.",19700.,"meeEBZ_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrDataEEZ_invMass, HsmMcEEZ_invMass, "m(ee) [GeV]::80::100", "a.u.",19700.,"meeEEZ_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrDataEBZ_phiStarZ, HsmMcEBZ_phiStarZ, "phi*::0::3.14", "a.u.",19700.,"phiStarEBZ_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrDataEEZ_phiStarZ, HsmMcEEZ_phiStarZ, "phi*::0::3.14", "a.u.",19700.,"phiStarEEZ_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrDataEBZ_cosThetaStarZ, HsmMcEBZ_cosThetaStarZ, "cos(theta*)::-1::1", "a.u.",19700.,"cosThetaStarEBZ_central","./",7,"Data","central, w47",""); 
  drawTH1pair(HrDataEEZ_cosThetaStarZ, HsmMcEEZ_cosThetaStarZ, "cos(theta*)::-1::1", "a.u.",19700.,"cosThetaStarEEZ_central","./",7,"Data","central, w47",""); 

  int okweight  = 1;
  int okweight1 = 1;
  int okweight2 = 1;
  float minchi2  = 1000.;
  float minchi21 = 1000.;
  float minchi22 = 1000.;
  int minchi2weight  = -1;
  int minchi2weight1 = -1;
  int minchi2weight2 = -1;
  for (int ii=1; ii<282; ii++) { 
    TString thisii = Form("%d",ii); 
    TString thename  = TString("ptZzoom_central_ww")+thisii;
    TString thename1 = TString("ptZzoom1_central_ww")+thisii;
    TString thename2 = TString("ptZzoom2_central_ww")+thisii;

    if (HsmMc_pTZzoom_ww[ii]) {
      float mychi2 = drawTH1pair(HrData_pTZzoom, HsmMc_pTZzoom_ww[ii], "Z pT [GeV]::0::30","a.u.",19700.,(string)thename,"./",7,"Data","central","");
      if (fabs(mychi2)<20 && mychi2>0.01) { 
	chi2_ptzoom[okweight] = mychi2; 
	theweight[okweight] = ii;
	okweight++;
	if (mychi2<minchi2) { minchi2 = mychi2; minchi2weight = ii; }
      }
    }

    if (HsmMc_pTZzoom1_ww[ii]) {
      float mychi21 = drawTH1pair(HrData_pTZzoom1, HsmMc_pTZzoom1_ww[ii], "Z pT [GeV]::0::20","a.u.",19700.,(string)thename1,"./",7,"Data","central","");
      if (fabs(mychi21)<20 && mychi21>0.01) { 
	chi2_ptzoom1[okweight1] = mychi21; 
	theweight1[okweight1] = ii;
	okweight1++;
	if (mychi21<minchi21) { minchi21 = mychi21; minchi2weight1 = ii; }
      }
    }

    if (HsmMc_pTZzoom2_ww[ii]) {
      float mychi22 = drawTH1pair(HrData_pTZzoom2, HsmMc_pTZzoom2_ww[ii], "Z pT [GeV]::0::10","a.u.",19700.,(string)thename2,"./",7,"Data","central","");
      if (fabs(mychi22)<20 && mychi22>0.01) { 
	chi2_ptzoom2[okweight2] = mychi22; 
	theweight2[okweight2] = ii;
	okweight2++;
	if (mychi22<minchi22) { minchi22 = mychi22; minchi2weight2 = ii; }
      }
    }
  }

  cout << "pTzoom: best chi2/ndf in 0-30 = " << minchi2  << " corresponding to weight " << minchi2weight  << endl;
  cout << "pTzoom: best chi2/ndf in 0-20 = " << minchi21 << " corresponding to weight " << minchi2weight1 << endl; 
  cout << "pTzoom: best chi2/ndf in 0-10 = " << minchi22 << " corresponding to weight " << minchi2weight2 << endl; 
  
  gStyle->SetOptStat(0);
  TGraph *myGraphPtZoom = new TGraph(okweight,theweight,chi2_ptzoom);
  TCanvas cPtZoom("cPtZoomChi2","",1);
  myGraphPtZoom->SetMarkerSize(1);
  myGraphPtZoom->SetMarkerColor(1);
  myGraphPtZoom->SetMarkerStyle(20);
  myGraphPtZoom->Draw("AP");
  cPtZoom.SaveAs("trendChi2_0-30.png");

  TGraph *myGraphPtZoom1 = new TGraph(okweight1,theweight1,chi2_ptzoom1);
  TCanvas cPtZoom1("cPtZoomChi2","",1);
  myGraphPtZoom1->SetMarkerSize(1);
  myGraphPtZoom1->SetMarkerColor(2);
  myGraphPtZoom1->SetMarkerStyle(21);
  myGraphPtZoom1->Draw("AP");
  cPtZoom1.SaveAs("trendChi2_0-20.png");

  TGraph *myGraphPtZoom2 = new TGraph(okweight2,theweight2,chi2_ptzoom2);
  TCanvas cPtZoom2("cPtZoomChi2","",1);
  myGraphPtZoom2->SetMarkerSize(1);
  myGraphPtZoom2->SetMarkerColor(4);
  myGraphPtZoom2->SetMarkerStyle(22);
  myGraphPtZoom2->Draw("AP");
  cPtZoom2.SaveAs("trendChi2_0-10.png");

  TCanvas cPtZoomAll("cPtZoomChiAll","",1);
  TH2F *myChi2H = new TH2F("myChi2H","myChi2H",100,-1,283,100,0,20);
  myChi2H->GetXaxis()->SetTitle("weight");
  myChi2H->GetYaxis()->SetTitle("#Chi^{2}");
  myChi2H->Draw();
  myGraphPtZoom->Draw("Psame");
  myGraphPtZoom1->Draw("Psame");
  myGraphPtZoom2->Draw("Psame");
  cPtZoomAll.SaveAs("variousTrendChi2.png");
}

