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
import ROOT
from ROOT import TMultiGraph, TCanvas, TFile, TGraph, TGraphErrors, TF1



def comparison_plot(args):

    mgp = TMultiGraph()
    mgdcat = TMultiGraph()
    mgdcaz = TMultiGraph()
    mgeff = TMultiGraph()
    
    fixed = ''
    for i, each in enumerate(args.root_input):
        print(each)
        title = ''
        filename = each.split('/')[-1]
        if('LS' in filename):
            title = filename.split('_')[-3] + '@' + filename.split('_')[-2]
        else:
            title = each.split('/')[-1].split('_')[-2]
        
        fixed = each.split('/')[-1].split('_')[-1].split('.')[0]
        f = TFile(each, 'r')
        gep = f.Get('Momentum')
        gedcat = f.Get('DCAT')
        gedcaz = f.Get('DCAZ')
        geeff = f.Get('Eff')


        gep.SetLineColor(i+1)
        gep.SetMarkerSize(1.2)
        gep.SetMarkerColor(i+1)
        gep.SetName(title)
        gep.SetTitle(title)

        if('LS' in filename):
            if('CurvedOB' in filename):
                gep.SetMarkerStyle(4)
                gep.SetMarkerColor(1)
                if('2' in filename):
                    gep.SetMarkerColor(4)
                elif('3' in filename):
                    gep.SetMarkerColor(2)
                if('Simple'in filename):
                    gep.SetMarkerColor(4)

            elif('FlatOB' in filename):
                gep.SetMarkerStyle(8)
                gep.SetMarkerColor(4)
                if('2' in filename):
                    gep.SetMarkerColor(1)
                elif('3' in filename):
                    gep.SetMarkerColor(2)

        mgp.Add(gep, 'AP')

        geeff.SetLineColor(i+1)
        geeff.SetMarkerSize(1.2)
        geeff.SetMarkerColor(i+1)
        geeff.SetName(title)
        geeff.SetTitle(title)

        if('LS' in filename):
            if('CurvedOB' in filename):
                geeff.SetMarkerStyle(4)
                geeff.SetMarkerColor(1)
                if('2' in filename):
                    geeff.SetMarkerColor(4)
                elif('3' in filename):
                    geeff.SetMarkerColor(2)
                if('Simple' in filename):
                    geeff.SetMarkerColor(4)

            elif('FlatOB' in filename):
                geeff.SetMarkerStyle(8)
                geeff.SetMarkerColor(4)
                if('2' in filename):
                    geeff.SetMarkerColor(4)
                elif('3' in filename):
                    geeff.SetMarkerColor(2)
        mgeff.Add(geeff, 'AP')

        gedcat.SetLineColor(i+1)
        gedcat.SetMarkerSize(1.2)
        gedcat.SetMarkerColor(i+1)
        gedcat.SetName(title)
        gedcat.SetTitle(title)

        if('LS' in filename):
            if('CurvedOB' in filename):
                gedcat.SetMarkerStyle(4)
                gedcat.SetMarkerColor(1)
                if('2' in filename):
                    gedcat.SetMarkerColor(4)
                elif('3' in filename):
                    gedcat.SetMarkerColor(2)

                if('Simple' in filename):
                    gedcat.SetMarkerColor(4)

            elif('FlatOB' in filename):
                gedcat.SetMarkerStyle(8)
                gedcat.SetMarkerColor(1)
                if('2' in filename):
                    gedcat.SetMarkerColor(4)
                elif('3' in filename):
                    gedcat.SetMarkerColor(2)    
        mgdcat.Add(gedcat, 'AP')

        gedcaz.SetLineColor(i+1)
        gedcaz.SetMarkerSize(1.2)
        gedcaz.SetMarkerColor(i+1)
        gedcaz.SetName(title)
        gedcaz.SetName(title)

        if('LS' in filename):
            if('CurvedOB' in filename):
                gedcaz.SetMarkerStyle(4)
                gedcaz.SetMarkerColor(2)
                if('2' in filename):
                    gedcaz.SetMarkerColor(4)
                elif('3' in filename):
                    gedcaz.SetMarkerColor(2)
                if('Simple' in filename):
                    gedcaz.SetMarkerColor(4)


            elif('FlatOB' in filename):
                gedcaz.SetMarkerStyle(8)
                gedcaz.SetMarkerColor(4)
                if('2' in filename):
                    gedcaz.SetMarkerColor(4)
                elif('3' in filename):
                    gedcaz.SetMarkerColor(2)
        mgdcaz.Add(gedcaz, 'AP')

    title_size = 0.08
    label_size = 0.08
    lmargin = 0.2
    bmargin = 0.2
    c1 = TCanvas('c1', 'c1', 1800, 1200)
    c1.Divide(2,2)
    pad1 = c1.cd(1)
    pad1.SetGridy(1)
    pad1.SetLeftMargin(lmargin)
    pad1.SetBottomMargin(bmargin)

    xtitle = ''
    if fixed == 'p':
        xtitle = 'p [GeV]'
        ytitle = '#sigma_{p}/p [%]'
    else:
        xtitle = 'p_{T} [GeV]'
        ytitle = '#sigma_{p_{T}}/p_{T} [%]'
    mgp.GetXaxis().SetTitle(xtitle)
    mgp.GetYaxis().SetTitle(ytitle)
    mgp.GetXaxis().SetTitleSize(title_size)
    mgp.GetXaxis().SetLabelSize(label_size)
    mgp.GetYaxis().SetTitleSize(title_size)
    mgp.GetYaxis().SetLabelSize(label_size)
    mgp.GetYaxis().SetRangeUser(0, 1.5)
    mgp.Draw("AP")

    # Draw YP requirement

    fdcat = TF1("yp", "sqrt(0.5*0.5 + 0.05*0.05*x*x)", 0.1, 20)
    fdcat.SetTitle("YP requirement")
    fdcat.SetLineColor(4)
    fdcat.SetLineStyle(2)
    fdcat.Draw("SAME")



    lgd = pad1.BuildLegend(0.22, 0.48, 0.55, 0.88, '', 'P')
    lgd.SetTextFont(62)
    lgd.SetFillColorAlpha(0,0)
    lgd.SetBorderSize(0)

    pad2 = c1.cd(2)
    pad2.SetGridy(1)
    pad2.SetLeftMargin(lmargin)
    pad2.SetBottomMargin(bmargin)
    mgeff.GetXaxis().SetTitle(xtitle)
    mgeff.GetYaxis().SetTitle('Eff. [%]')
    mgeff.GetXaxis().SetTitleSize(title_size)
    mgeff.GetXaxis().SetLabelSize(label_size)
    mgeff.GetYaxis().SetTitleSize(title_size)
    mgeff.GetYaxis().SetLabelSize(label_size)
    mgeff.GetYaxis().SetRangeUser(0, 105)
    mgeff.Draw("AP")

    pad3 = c1.cd(3)
    pad3.SetGridy(1)
    pad3.SetLeftMargin(lmargin)
    pad3.SetBottomMargin(bmargin)
    mgdcat.GetXaxis().SetTitle(xtitle)
    mgdcat.GetYaxis().SetTitle('DCA_{T} [#mum]')
    mgdcat.GetXaxis().SetTitleSize(title_size)
    mgdcat.GetXaxis().SetLabelSize(label_size)
    mgdcat.GetYaxis().SetTitleSize(title_size)
    mgdcat.GetYaxis().SetLabelSize(label_size)
    mgdcat.GetYaxis().SetRangeUser(0, 100)
    mgdcat.Draw("AP")

    f = TF1("yp", "sqrt(5*5 + 20*20/(x*x))", 0.1, 20)
    f.SetTitle("YP requirement")
    f.SetLineColor(4)
    f.SetLineStyle(2)
    f.Draw("SAME")


    pad4 = c1.cd(4)
    pad4.SetGridy(1)
    pad4.SetLeftMargin(lmargin)
    pad4.SetBottomMargin(bmargin)
    mgdcaz.GetXaxis().SetTitle(xtitle)
    mgdcaz.GetYaxis().SetTitle('DCA_{Z} [#mum]')
    mgdcaz.GetXaxis().SetTitleSize(title_size)
    mgdcaz.GetXaxis().SetLabelSize(label_size)
    mgdcaz.GetYaxis().SetTitleSize(title_size)
    mgdcaz.GetYaxis().SetLabelSize(label_size)
    mgdcaz.GetYaxis().SetRangeUser(0, 100)
    mgdcaz.Draw("AP")

    fdcaz = TF1("yp", "sqrt(5*5 + 20*20/(x*x))", 0.1, 20)
    fdcaz.SetTitle("YP requirement")
    fdcaz.SetLineColor(4)
    fdcaz.SetLineStyle(2)
    fdcaz.Draw("SAME")

    c1.Print(f'plots/Comparison_plots_{fixed}.pdf')





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='ePIC SVT single particle simulation')
    parser.add_argument('root_input', nargs='+', help='input root files')
    
    args = parser.parse_args()
    comparison_plot(args)



