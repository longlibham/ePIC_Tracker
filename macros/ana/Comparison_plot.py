# this is a python script

#=======================================================================
#   Copyright (C) 2024 Univ. of Bham  All rights reserved.
#   
#   		FileName：		Comparison_plot.py
#   	 	Author：		LongLI <long.l@cern.ch>
#   		Time：			2024.06.19
#   		Description：
#
#======================================================================

import os
import sys
import argparse
from ROOT import *


def comparison_plot(args):

    mgp = TMultiGraph()
    mgdcat = TMultiGraph()
    mgdcaz = TMultiGraph()
    
    fixed = ''
    for i, each in enumerate(args.root_input):
        print(each)
        title = each.split('/')[-1].split('_')[-2]
        fixed = each.split('/')[-1].split('_')[-1].split('.')[0]
        f = TFile(each, 'r')
        gep = f.Get('Momentum')
        gedcat = f.Get('DCAT')
        gedcaz = f.Get('DCAZ')


        gep.SetLineColor(i+1)
        gep.SetMarkerSize(1.2)
        gep.SetMarkerColor(i+1)
        gep.SetName(title)
        gep.SetTitle(title)
        mgp.Add(gep, 'AP')

        gedcat.SetLineColor(i+1)
        gedcat.SetMarkerSize(1.2)
        gedcat.SetMarkerColor(i+1)
        gedcat.SetName(title)
        gedcat.SetTitle(title)
        mgdcat.Add(gedcat, 'AP')

        gedcaz.SetLineColor(i+1)
        gedcaz.SetMarkerSize(1.2)
        gedcaz.SetMarkerColor(i+1)
        gedcaz.SetName(title)
        gedcaz.SetName(title)
        mgdcaz.Add(gedcaz, 'AP')


    c1 = TCanvas('c1', 'c1', 1800, 1200)
    c1.Divide(2,2)
    pad1 = c1.cd(1)
    
    xtitle = ''
    if fixed == 'p':
        xtitle = 'p [GeV]'
        ytitle = '#Deltap/p [%]'
    else:
        xtitle = 'p_{T} [GeV]'
        ytitle = '#Deltap_{T}/p_{T} [%]'
    mgp.GetXaxis().SetTitle(xtitle)
    mgp.GetYaxis().SetTitle(ytitle)
    mgp.GetYaxis().SetRangeUser(0, 2.5)
    mgp.Draw("AP")

    lgd = pad1.BuildLegend(0.15, 0.5, 0.4, 0.8)
    lgd.SetTextFont(62)
    lgd.SetFillColorAlpha(0,0)
    lgd.SetBorderSize(0)

    c1.cd(3)

    mgdcat.GetXaxis().SetTitle(xtitle)
    mgdcat.GetYaxis().SetTitle('DCA_{T} [#mum]')
    mgdcat.GetYaxis().SetRangeUser(0, 100)
    mgdcat.Draw("AP")

    c1.cd(4)

    mgdcaz.GetXaxis().SetTitle(xtitle)
    mgdcaz.GetYaxis().SetTitle('DCA_{Z} [#mum]')
    mgdcaz.GetYaxis().SetRangeUser(0, 100)
    mgdcaz.Draw("AP")

    c1.Print(f'plots/Comparison_plots_{fixed}.pdf')





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='ePIC SVT single particle simulation')
    parser.add_argument('root_input', nargs='+', help='input root files')
    
    args = parser.parse_args()
    comparison_plot(args)



