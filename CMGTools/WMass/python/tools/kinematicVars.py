from CMGTools.TTHAnalysis.treeReAnalyzer import *
from CMGTools.WMass.tools.standaloneElectronCalibrator import ElectronCalibrator
from CMGTools.WMass.tools.PileUpReWeighter import PileUpReWeighter
from CMGTools.WMass.tools.eventVars_wmass import SimpleVBoson
import types, os, math
import numpy as np

class KinematicVars:
    def __init__(self):
        ROOT.gSystem.Load("libFWCoreFWLite.so")
        ROOT.AutoLibraryLoader.enable()
        ROOT.gSystem.Load("libCMGToolsWMass.so")
    def initSample(self,nevents,dataset):
        self.isMC = not any(x in dataset for x in "DoubleMu DoubleEl DoubleEG MuEG MuonEG SingleMu SingleEl".split())
        aifilename = "%s/src/CMGTools/WMass/data/theory/angularcoeffs/Ai_8TeV_atlas.root" % os.environ['CMSSW_BASE']
        tfile = ROOT.TFile.Open(aifilename)
        etabins = ["00","10","20","35"]
        self.Ai = {}
        for i in range(5):
            for e in range(len(etabins)-1):
                aistr = "A%d_8TeV_atlas_eta%s_%s" % (i,etabins[e],etabins[e+1])
                self.Ai[(i,e)] = tfile.Get(aistr)
    def listBranches(self):
        mylist = [ ("nAP","I"), ("AP","F",8,"nAP"),
                   ("z_y","F"), ("costheta_cs","F"), ("phi_cs","F")]
        self.branches = mylist
        print "self.branches = ",self.branches[:]
        return self.branches[:]
    def CSFrame(self,dilepton):
        pMass = 0.938272
        beamE = 8000;
        sign = np.sign(dilepton.Z())
        proton1 = ROOT.TLorentzVector(0.,0.,sign*beamE,hypot(beamE,pMass));  proton2 = ROOT.TLorentzVector(0.,0.,-sign*beamE,hypot(beamE,pMass))
        proton1.Boost(-dilepton.BoostVector()); proton2.Boost(-dilepton.BoostVector())
        CSAxis = (proton1.Vect().Unit()-proton2.Vect().Unit()).Unit()
        yAxis = (proton1.Vect().Unit()).Cross((proton2.Vect().Unit()));
        yAxis = yAxis.Unit();
        xAxis = yAxis.Cross(CSAxis);
        xAxis = xAxis.Unit();
        return (xAxis,yAxis,CSAxis)
    def cosThetaCS(self,lplus,lminus):
        dilep = lplus + lminus
        boostedLep = lminus
        boostedLep.Boost(-dilep.BoostVector())
        csframe = self.CSFrame(dilep)
        return cos(boostedLep.Angle(csframe[2]))
    def phiCS(self,lplus,lminus):
        dilep = lplus + lminus
        boostedLep = lminus
        boostedLep.Boost(-dilep.BoostVector())
        csframe = self.CSFrame(dilep)
        phi = atan2((boostedLep.Vect()*csframe[1]),(boostedLep.Vect()*csframe[0]))
        if(phi<0): return phi + 2*ROOT.TMath.Pi()
        else: return phi
    def P_i(self,theta,phi,i):
        ct = cos(theta); st = sqrt(1-pow(ct,2))
        cp = cos(phi); sp = sqrt(1-pow(cp,2))
        if i==0: return 0.5*(1-3*ct*ct)
        elif i==1: return 2*ct*st*cp
        elif i==2: return 0.5*st*st*(cp*cp-sp*sp)
        elif i==3: return st*cp
        elif i==4: return ct
        elif i==5: return 2*st*st*sp*cp
        elif i==6: return 2*st*ct*sp
        elif i==7: return st*sp
        else: return -999.
    def A_i(self,i,ptz,yz):
        if abs(yz)<1: histo = self.Ai[(i,0)]
        elif abs(yz)<2: histo = self.Ai[(i,1)]
        else: histo = self.Ai[(i,2)]
        if histo.GetEntries() == 0: return 0
        xmin = histo.GetXaxis().GetXmin()+1e-3
        xmax = histo.GetXaxis().GetXmax()-1e-3
        xbin = histo.GetXaxis().FindBin(min(max(ptz,xmin),xmax))
        return histo.GetBinContent(xbin)
    def PtEtaPhiM4V(self,pt,eta,phi,m):
        tlv = ROOT.TLorentzVector()
        tlv.SetPtEtaPhiM(pt,eta,phi,m)
        return tlv
    def __call__(self,event):
        # prepare output
        ret = {}; 
        genp = [p for p in Collection(event,"GenP6StatusThree","nGenP6StatusThree")] if self.isMC else []
        leps = [l for l in Collection(event,"LepGood","nLepGood")]
        leps_4v=[ self.PtEtaPhiM4V(l.pt,l.eta,l.phi,0.511*e-3) for l in leps ]
        
        ret["nAP"] = 5 if len(leps)>=2 else 0
        ret["AP"]=[]
        if len(leps)>=2:
            (lplus,lminus) = (leps_4v[0],leps_4v[1]) if leps[0].charge>0 else (leps_4v[1],leps_4v[0])
            Z = SimpleVBoson(leps_4v[:2])
            costheta_cs = self.cosThetaCS(lplus,lminus); phi_cs = self.phiCS(lplus,lminus);
            ret["costheta_cs"] = costheta_cs; ret["phi_cs"] = phi_cs
            ret["z_y"] = Z.y()
            for i in range(ret["nAP"]):
                ret["AP"].append(self.A_i(i,Z.pt(),Z.y())*self.P_i(acos(costheta_cs),phi_cs,i))
        else:
            ret["costheta_cs"] = ret["phi_cs"] = ret["z_y"] = -999

        ### return
        return ret

if __name__ == '__main__':
    from sys import argv
    file = ROOT.TFile(argv[1])
    tree = file.Get("tree")
    tree.vectorTree = True
    class Tester(Module):
        def __init__(self, name):
            Module.__init__(self,name,None)
            self.sf = EventVarsWmass()
        def analyze(self,ev):
            #if ev.metNoMu_pt < 200: return True
            print "\nrun %6d lumi %4d event %d: metNoMu %d" % (ev.run, ev.lumi, ev.evt, ev.metNoMu_pt)
            print self.sf(ev)
    el = EventLoop([ Tester("tester") ])
    el.loop([tree], maxEvents = 50)

