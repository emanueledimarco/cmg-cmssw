#!/usr/bin/env python
import os,re,sys

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

args = sys.argv[:]

filename = args[1]
basename = os.path.basename(filename).split('.')[0]
print "Creating histograms for ", basename

ptBins = [0.0 , 2.5 , 5.0 , 8.0 , 11.4 , 14.9 , 18.5 , 22.0 , 25.5 , 29.0 , 32.6 , 36.4 , 40.4 , 44.9 , 50.2 , 56.4 , 63.9 , 73.4 , 85.4 , 105 , 132 , 173 , 253 , 600 ]
etaBins = [0.0, 1.0, 2.0, 3.5]

tfile = ROOT.TFile.Open(basename+".root","recreate")

inputfile = open(filename)
lines = inputfile.readlines()
for i,line in enumerate(lines):
    l=line.rstrip()
    p = re.compile('(-?0\.\d+\spm\s0\.\d+\spm\s0\.\d+)')
    vals = p.findall(l)
    print "Ai for eta bin ",i," = ", vals
    for b,ptb in enumerate(vals):
        ai = [float(f) for f in ptb.split(" pm ")]
        print "pt bin: [",ptBins[b],",",ptBins[b+1],"] = ",ai
