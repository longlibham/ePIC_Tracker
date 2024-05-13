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
    
    # fixed pt or fixed p
    fixed = outfile.split('/')[-2].split('_')[-1]

    gp = int(outfile.split('-')[-1].split('G')[0]) 
    bins = gp*200
    # load the tree
    c1 = ROOT.TCanvas("c1", "c1", 1800, 1200)
    ptitle = ''
    if fixed == 'p':
        ptitle = fixed
    elif fixed == 'pT':
        ptitle = 'p_{T}'
    
    title = '#Delta%s/%s;#Delta%s/%s;Entries/%d' %(ptitle, ptitle, ptitle, ptitle, round(0.2/bins, 5))     
    h1 = ROOT.TH1D("h1", title, bins, 0.1, 0.1)
    h3 = ROOT.TH1D("h3", "DCA_{T};DCA_{T};Entries/0.001cm", 100, 0.05, 0.05)
    h4 = ROOT.TH1D("h4", "DCA_{Z};DCA_{Z};Entries/0.001cm", 100, 0.05, 0.05)

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
        denumerator = 0
        if fixed == 'p': 
            denumerator = math.sqrt(tree.gpx**2 + tree.gpy**2 + tree.gpz**2)
        elif fixed == 'pT': 
            denumerator = math.sqrt(tree.gpz**2 + tree.gpy**2)
        else:
            print('please choose the right data fixp/fixpt!') 
            exit(-1)
        h1.Fill((math.sqrt(tree.px**2+tree.py**2)- math.sqrt(tree.gpx**2+tree.gpy**2))/denumerator)
        h3.Fill(tree.dca2d)
        h4.Fill(tree.pcaz - tree.gvz)

    c1.Divide(2, 2)
    c1.cd(1)
    frp_p = h1.Fit('gaus','S', '', h1.GetMean()-h1.GetRMS(), h1.GetMean()+h1.GetRMS())
    h1.GetXaxis().SetRangeUser(-0.1, 0.1)
    h1.Draw()
    c1.cd(3)
   
    ''' 
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
    '''

    #h2.Fit("fit_dca", 'S', '', h2.GetMean()-h2.GetRMS(), h2.GetMean+h2.GetRMS())
    frp_dcat = h3.Fit("gaus", 'S', '', h3.GetMean()-2*h3.GetRMS(), h3.GetMean()+2*h3.GetRMS())
    h3.Draw()
    c1.cd(4)
    frp_dcaz = h4.Fit("gaus", 'S', '', h4.GetMean()- h4.GetRMS(), h4.GetMean()+ h4.GetRMS())
    h4.Draw()
    c1.Print(outfile.replace(".root", ".pdf"))

    dp_mean, dp_std = frp_p.Parameter(1)*100, frp_p.Parameter(2)*100
    dcat_mean, dcat_std, dcaz_mean, dcaz_std = frp_dcat.Parameter(1)*10000, frp_dcat.Parameter(2)*10000, frp_dcaz.Parameter(1)*10000, frp_dcaz.Parameter(2)*10000 # unit %, um

    Tfile.Close()
    return dp_mean, dp_std, dcat_mean, dcat_std, dcaz_mean, dcaz_std

def process_tree(args):
    pgen = array('d')
    dp_mean = array('d')
    dp_std = array('d')
    dcat_mean = array('d')
    dcat_std = array('d')
    dcaz_mean = array('d')
    dcaz_std = array('d')
    # p/pt fixed
    fixed = ''
    for i, each in enumerate(args.file_in):
        f = ROOT.TFile(each,"READ")
        p_m, p_s, dcat_m, dcat_s, dcaz_m, dcaz_s = read_tree(f, each)
        dp_mean.append(p_m)
        dp_std.append(p_s)
        dcat_mean.append(dcat_m)
        dcat_std.append(dcat_s)
        dcaz_mean.append(dcaz_m)
        dcaz_std.append(dcaz_s)

        p_num = float(each.split("_")[-1].split("G")[0].split("-")[0])
        pgen.append(p_num)
        fixed = each.split('/')[-2].split('_')[-1]

    # prepare data for TGraphErrrors
    num = len(pgen)
    print(dcat_std)
    ge1 = ROOT.TGraph(num, pgen, dp_std)
    ge3 = ROOT.TGraph(num, pgen, dcat_std)
    ge4 = ROOT.TGraph(num, pgen, dcaz_std)
    
    c2 = ROOT.TCanvas('c2', 'c2', 1800, 1200)
    c2.Divide(2,2)
    c2.cd(1)


    ptitle = ''
    if fixed == 'p':
        ptitle = 'p'
    elif fixed == 'pT':
        ptitle = 'p_{T}'


    xtitle = '%s [GeV/c]' % ptitle
    ytitle = '#Delta%s/%s ' %(ptitle, ptitle)
    ytitle += '[%]'    
    ge1.GetXaxis().SetTitle(xtitle)
    ge1.GetYaxis().SetRangeUser(0, 2.5)
    ge1.GetYaxis().SetTitle(ytitle)
    ge1.Draw("AP")
    

    c2.cd(3)
    ge3.GetXaxis().SetTitle(xtitle)
    ge3.GetYaxis().SetTitle("DCA_{T} [#mum]")
    ge3.GetYaxis().SetRangeUser(0, 50)
    ge3.Draw("AP")

    c2.cd(4)
    ge4.GetXaxis().SetTitle(xtitle)
    ge4.GetYaxis().SetTitle("DCA_{Z} [#mum]")
    ge4.GetYaxis().SetRangeUser(0, 50)
    ge4.Draw("AP")
    c2.Print(f"plots/Summary_plot_{fixed}.pdf")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Read data form FastTracking file and plotting")
    parser.add_argument('file_in', nargs='+', help="input file names")

    args = parser.parse_args()

    process_tree(args)

