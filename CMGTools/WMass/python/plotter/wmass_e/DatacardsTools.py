#!/usr/bin/env python
import os.path
import sys,ROOT,os
ROOT.gROOT.SetBatch(True)

class CardsChecker:
    def __init__(self, card_dir, n_masses):
        self.card_dir = card_dir
        self.subdirs = [os.path.basename(x[0]) for x in os.walk(card_dir) if "eta" in x[0]]
        self.datacards = {}
        self.cardinputs = {}
        for subdir in self.subdirs:
            for m in range(int(n_masses)):
                for charge in ["neg","pos"]:
                    key = "wenu_mass%d_%s_%s" % (m,charge,subdir)
                    f_txt = subdir+"/"+key+".card.txt" 
                    f_root = subdir+"/"+key+".input.root"
                    self.datacards[key] = f_txt
                    self.cardinputs[key] = f_root

    def checkCards(self):
        resubcmds = {}
        for key,dc in self.datacards.iteritems():
            if not os.path.exists(self.card_dir+"/"+dc): 
                print "datacard ",dc," is not present in ",self.card_dir
                resubcmds[key] = "bsub -q {queue} -o {dir}/{logfile} {dir}/{srcfile}".format(
                    queue="8nh", dir="/".join([os.getcwd(),self.card_dir,"jobs"]), logfile=key+"_resub.log", srcfile=key+".sh")
        for key,f in self.cardinputs.iteritems():
            f_ok = True
            if not os.path.exists(self.card_dir+"/"+f): 
                print "input root file ",f," is not present in ",self.card_dir
                f_ok = False
            else:
                tfile = ROOT.TFile.Open(self.card_dir+"/"+f)
                if not tfile or tfile.IsZombie():
                    print f, " is Zombie"
                    f_ok = False
            if not f_ok: 
                resubcmds[key] = "bsub -q {queue} -o {dir}/{logfile} {dir}/{srcfile}".format(
                    queue="8nh", dir="/".join([os.getcwd(),self.card_dir,"jobs"]), logfile=key+"_resub.log", srcfile=key+".sh")
        return resubcmds

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] [inputs]")
    parser.add_option("-c", "--check-cards", dest="check", default=False, action="store_true", help="Check if there are all the datacards and ROOT files are fine");
    (options, args) = parser.parse_args()

    if options.check:
        if len(args)<2: print "needed inputs: datacards_dir nmasses"; exit(1)
        cc = CardsChecker(args[0],args[1])
        result = cc.checkCards()
        if len(result)==0: print "All cards are GOOD."
        else: 
            for k,cmd in result.iteritems(): print cmd

