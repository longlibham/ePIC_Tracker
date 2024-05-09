# this is a python script

#=======================================================================
#   Copyright (C) 2024 Univ. of Bham  All rights reserved.
#   
#   		FileName：		FastTrackingPlot.py
#   	 	Author：		LongLI <long.l@cern.ch>
#   		Time：			2024.03.27
#   		Description：
#
#======================================================================

import ROOT
import sys
import os
import argparse
import math
import numpy as np
from array import array


class Bigaus:
    def __call__(self, arr, par):
        rev = par[0]*ROOT.TMath.Gaus(arr[0], par[1], par[2], 0)
        rev += par[3]*ROOT.TMath.Gaus(arr[0], par[4], par[5], 0)
        return rev

class Trigaus:
    def __call__(self, arr, par):
        rev = par[0]*ROOT.TMath.Gaus(arr[0], par[1], par[2], 0)
        rev += par[3]*ROOT.TMath.Gaus(arr[0], par[4], par[5], 0)
        rev += par[6]*ROOT.TMath.Gaus(arr[0], par[7], par[8], 0)
        return rev

def read_tree(Tfile, outfile):
    tree = Tfile.Get('tracks')

    # load the tree
    c1 = ROOT.TCanvas("c1", "c1", 1800, 600)
    h1 = ROOT.TH1D("h1d", "#DeltaP/P;#DeltaP/P;Entries/0.01", 100, 0.5, 0.5)
    h2 = ROOT.TH1D("h2", "DCA;#DCA;Entries/0.001cm", 100, 0.05, 0.05)

    entries = tree.GetEntries()
    np_px = np.zeros(entries)
    np_py = np.zeros(entries)
    np_pz = np.zeros(entries)
    np_pcax = np.zeros(entries)
    np_pcay = np.zeros(entries)
    np_pcaz = np.zeros(entries)
    np_gpx = np.zeros(entries)
    np_gpy = np.zeros(entries)
    np_gpz = np.zeros(entries)
    np_gvx = np.zeros(entries)
    np_gvy = np.zeros(entries)
    np_gvz = np.zeros(entries)
    
    
    for entry in range(entries):
        tree.GetEntry(entry)
        if math.isnan(tree.px):
            continue
        h1.Fill((math.sqrt(tree.px**2+tree.py**2+tree.pz**2)- math.sqrt(tree.gpx**2+tree.gpy**2+tree.gpz**2))/math.sqrt(tree.gpx**2+tree.gpy**2+tree.gpz**2))
        h2.Fill(math.sqrt(tree.pcax**2+tree.pcay**2+tree.pcaz**2) - math.sqrt(tree.gvx**2+tree.gvy**2+tree.gvz**2))

    c1.Divide(2, 1)
    c1.cd(1)
    h1.Fit('gaus','', '', -0.1, 0.1)
    h1.Draw()
    c1.cd(2)
    
    #bigaus = Bigaus()
    trigaus = Trigaus()
    #fit_dca = ROOT.TF1("fit_dca", bigaus, -0.02, 0.02, 6)
    fit_dca = ROOT.TF1("fit_dca", trigaus, -0.02, 0.02, 9)
    fit_dca.SetParameter(0, 500)
    fit_dca.SetParName(0, "A1")
    fit_dca.SetParLimits(0, 0, 10000)
    fit_dca.SetParameter(1, 0)
    fit_dca.SetParName(1, "#mu1")
    fit_dca.SetParLimits(1, -0.01, 0.01)
    fit_dca.SetParameter(2, 0.002)
    fit_dca.SetParName(2, "#sigma1")
    fit_dca.SetParLimits(2, 0.0002, 0.02)
    fit_dca.SetParameter(3, 500)
    fit_dca.SetParName(3, "A2")
    fit_dca.SetParLimits(3, 0, 10000)
    fit_dca.SetParameter(4, 0)
    fit_dca.SetParName(4, "#mu2")
    fit_dca.SetParLimits(4, -0.01, 0.01)
    fit_dca.SetParameter(5, 0.002)
    fit_dca.SetParName(5, "#sigma2")
    fit_dca.SetParLimits(5, 0.0002, 0.02)
    fit_dca.SetParameter(6, 500)
    fit_dca.SetParName(6, "A3")
    fit_dca.SetParLimits(6, 0, 10000)
    fit_dca.SetParameter(7, 0)
    fit_dca.SetParName(7, "#mu3")
    fit_dca.SetParLimits(7, -0.01, 0.01)
    fit_dca.SetParameter(8, 0.002)
    fit_dca.SetParName(8, "#sigma3")
    fit_dca.SetParLimits(8, 0.0002, 0.02)

    h2.Fit("fit_dca", 'S', '', -0.01, 0.01)
    h2.Draw()
    c1.Print(outfile.replace(".root", ".pdf"))

    dp_mean, dp_std, dca_mean, dca_std = h1.GetMean()*100, h1.GetRMS()*100, h2.GetMean()*10000, h2.GetRMS()*10000 # unit %, um

    Tfile.Close()
    return dp_mean, dp_std, dca_mean, dca_std

def process_tree(args):
    pgen = array('d')
    dp_mean = array('d')
    dp_std = array('d')
    dca_mean = array('d')
    dca_std = array('d')
    for i, each in enumerate(args.file_in):
        f = ROOT.TFile(each,"READ")
        p_m, p_s, dca_m, dca_s = read_tree(f, each)
        dp_mean.append(p_m)
        dp_std.append(p_s)
        dca_mean.append(dca_m)
        dca_std.append(dca_s)

        p_num = float(each.split("_")[-1].split("G")[0])
        pgen.append(p_num)


    # prepare data for TGraphErrrors
    num = len(pgen)
    ge1 = ROOT.TGraph(num, pgen, dp_std)
    ge2 = ROOT.TGraphErrors(num, pgen, dca_std)
    
    c2 = ROOT.TCanvas('c2', 'c2', 1800, 600)
    c2.Divide(2,1)
    c2.cd(1)
    ge1.GetXaxis().SetTitle("P [GeV/c]")
    ge1.GetYaxis().SetRangeUser(0, 5)
    ge1.GetYaxis().SetTitle("#DeltaP/P [%]")
    ge1.Draw("AP")
    c2.cd(2)
    ge2.GetXaxis().SetTitle("P [GeV/c]")
    ge2.GetYaxis().SetTitle("DCA [#mum]")
    ge2.GetYaxis().SetRangeUser(0, 50)
    ge2.Draw("AP")
    c2.Print("Summary_plot.pdf")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Read data form FastTracking file and plotting")
    parser.add_argument('file_in', nargs='+', help="input file names")

    args = parser.parse_args()

    process_tree(args)

