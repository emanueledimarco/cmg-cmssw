//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr  7 14:51:09 2017 by ROOT version 6.02/05
// from TTree data/data
// found on file: /cmsrm/pc28_2/crovelli/data/Wmass/ZJ_ptsqmin4_Scan0.root
//////////////////////////////////////////////////////////

#ifndef smearerAppliedToZPedro_h
#define smearerAppliedToZPedro_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class smearerAppliedToZPedro {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Fixed size dimensions of array or collections stored in the TTree if any.
  
  // Declaration of leaf types
  Int_t           nw;
  Float_t         w[10];   //[nw]
  Int_t           nl;
  Int_t           pid[10];   //[nl]
  Int_t           charge[10];   //[nl]
  Float_t         pt[10];   //[nl]
  Float_t         eta[10];   //[nl]
  Float_t         phi[10];   //[nl]
  Float_t         m[10];   //[nl]
  Float_t         dressed_pt[10];   //[nl]
  Float_t         dressed_eta[10];   //[nl]
  Float_t         dressed_phi[10];   //[nl]
  Float_t         dressed_m[10];   //[nl]
  Float_t         imbalance_pt[10];
  Float_t         imbalance_eta[10];
  Float_t         imbalance_phi[10];
  Float_t         vecbos_pt;
  Float_t         vecbos_eta;
  Float_t         vecbos_phi;
  Float_t         vecbos_m;
  Int_t           id1;
  Int_t           id2;
  Float_t         x1;
  Float_t         x2;
  Float_t         qscale;
  
  // List of branches
  TBranch        *b_nw;   //!
  TBranch        *b_w;   //!
  TBranch        *b_nl;   //!
  TBranch        *b_pid;   //!
  TBranch        *b_charge;   //!
  TBranch        *b_pt;   //!
  TBranch        *b_eta;   //!
  TBranch        *b_phi;   //!
  TBranch        *b_m;   //!
  TBranch        *b_dressed_pt;   //!
  TBranch        *b_dressed_eta;   //!
  TBranch        *b_dressed_phi;   //!
  TBranch        *b_dressed_m;   //!
  TBranch        *b_imbalance_pt;   //!
  TBranch        *b_imbalance_eta;   //!
  TBranch        *b_imbalance_phi;   //!
  TBranch        *b_vecbos_pt;   //!
  TBranch        *b_vecbos_eta;   //!
  TBranch        *b_vecbos_phi;   //!
  TBranch        *b_vecbos_m;   //!
  TBranch        *b_id1;   //!
  TBranch        *b_id2;   //!
  TBranch        *b_x1;   //!
  TBranch        *b_x2;   //!
  TBranch        *b_qscale;   //!
  
  smearerAppliedToZPedro(TTree *tree=0);
  virtual ~smearerAppliedToZPedro();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef smearerAppliedToZPedro_cxx
smearerAppliedToZPedro::smearerAppliedToZPedro(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/cmsrm/pc28_2/crovelli/data/Wmass/Zj_nominal.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("/cmsrm/pc28_2/crovelli/data/Wmass/Zj_nominal.root");
    }
    TDirectory * dir = (TDirectory*)f->Get("/cmsrm/pc28_2/crovelli/data/Wmass/Zj_nominal.root:/analysis");
    dir->GetObject("data",tree);
  }
  Init(tree);
}

smearerAppliedToZPedro::~smearerAppliedToZPedro()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t smearerAppliedToZPedro::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t smearerAppliedToZPedro::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void smearerAppliedToZPedro::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  
  fChain->SetBranchAddress("nw", &nw, &b_nw);
  fChain->SetBranchAddress("w", w, &b_w);
  fChain->SetBranchAddress("nl", &nl, &b_nl);
  fChain->SetBranchAddress("pid", pid, &b_pid);
  fChain->SetBranchAddress("charge", charge, &b_charge);
  fChain->SetBranchAddress("pt", pt, &b_pt);
  fChain->SetBranchAddress("eta", eta, &b_eta);
  fChain->SetBranchAddress("phi", phi, &b_phi);
  fChain->SetBranchAddress("m", m, &b_m);
  fChain->SetBranchAddress("dressed_pt", dressed_pt, &b_dressed_pt);
  fChain->SetBranchAddress("dressed_eta", dressed_eta, &b_dressed_eta);
  fChain->SetBranchAddress("dressed_phi", dressed_phi, &b_dressed_phi);
  fChain->SetBranchAddress("dressed_m", dressed_m, &b_dressed_m);
  fChain->SetBranchAddress("imbalance_pt", imbalance_pt, &b_imbalance_pt);
  fChain->SetBranchAddress("imbalance_eta", imbalance_eta, &b_imbalance_eta);
  fChain->SetBranchAddress("imbalance_phi", imbalance_phi, &b_imbalance_phi);
  fChain->SetBranchAddress("vecbos_pt", &vecbos_pt, &b_vecbos_pt);
  fChain->SetBranchAddress("vecbos_eta", &vecbos_eta, &b_vecbos_eta);
  fChain->SetBranchAddress("vecbos_phi", &vecbos_phi, &b_vecbos_phi);
  fChain->SetBranchAddress("vecbos_m", &vecbos_m, &b_vecbos_m);
  fChain->SetBranchAddress("id1", &id1, &b_id1);
  fChain->SetBranchAddress("id2", &id2, &b_id2);
  fChain->SetBranchAddress("x1", &x1, &b_x1);
  fChain->SetBranchAddress("x2", &x2, &b_x2);
  fChain->SetBranchAddress("qscale", &qscale, &b_qscale);
  Notify();
}

Bool_t smearerAppliedToZPedro::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}

void smearerAppliedToZPedro::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t smearerAppliedToZPedro::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef smearerAppliedToZPedro_cxx
