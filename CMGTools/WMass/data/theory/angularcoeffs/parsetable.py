#!/usr/bin/env python
# usage: python parsetable.py A*_8TeV_atlas.txt
import os,re,sys
from array import array
from math import *

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

def parseOneTable(filename,tfile):
    inputfile = open(filename)
    lines = inputfile.readlines()
    for i,line in enumerate(lines):
        l=line.rstrip()
        p = re.compile('(-?0\.\d+\spm\s0\.\d+\spm\s0\.\d+)')
        vals = p.findall(l)
        etastr = "eta%s_%s" % (etaBins[i],etaBins[i+1])
        histoname = basename+"_"+etastr
        histo = ROOT.TH1F(histoname,"",len(ptBins)-1,array('f',ptBins))
        for b,ptb in enumerate(vals):
            ai = [float(f) for f in ptb.split(" pm ")]
            #print "pt bin: [",ptBins[b],",",ptBins[b+1],"] = ",ai
            histo.SetBinContent(b+1,ai[0])
            histo.SetBinError(b+1,hypot(ai[1],ai[2]))
        tfile.WriteTObject(histo)


if __name__ == '__main__':

    args = sys.argv[:]
    if len(args) == 0:
        print "You must provide at least one file to parse"
        exit(0)

    ptBins = [0.0 , 2.5 , 5.0 , 8.0 , 11.4 , 14.9 , 18.5 , 22.0 , 25.5 , 29.0 , 32.6 , 36.4 , 40.4 , 44.9 , 50.2 , 56.4 , 63.9 , 73.4 , 85.4 , 105 , 132 , 173 , 253 , 600 ]
    etaBins = ["00", "10", "20", "35"]

    tfile = ROOT.TFile.Open("Ai_8TeV_atlas.root","recreate")

    for fname in args[1:]:
        print "fname = ", fname
        basename = os.path.basename(fname).split('.')[0]
        print "Creating histograms for ", basename
        parseOneTable(fname,tfile)

    tfile.Close()

    print "DONE."
