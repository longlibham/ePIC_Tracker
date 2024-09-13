# this is a python script

#=======================================================================
#   Copyright (C) 2024 Univ. of Bham  All rights reserved.
#   
#   		FileName：		autorun.py
#   	 	Author：		LongLI <long.l@cern.ch>
#   		Time：			2024.05.10
#   		Description：
#
#======================================================================
import os
import sys
import argparse
import math
import ROOT


def autorun(args):
    save_path = ''
    if args.use_pt:
        save_path = os.path.join(args.outroot, 'fix_pT/')
    else:
        save_path = os.path.join(args.outroot, 'fix_p/')
    
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    
    print(save_path)

    for p in args.prange:
        use_pt = str(args.use_pt).lower()
        cmd = f'root -b -q Fun4All_ePIC_Tracker.C\({args.evts}\,\\"{args.particlename}\\"\,{use_pt}\,{p}\,{p}\,{args.emin}\,{args.emax}\,{args.phimin}\,{args.phimax}\,\\"{save_path}\\"\)'
        print(cmd)
        os.system(cmd)






if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Automatic run for Fun4All_ePIC_SVT_OB')
    parser.add_argument('--evts', type=int, default=10000, help='events for each run')
    parser.add_argument('--outroot', type=str, default='../data/', help="Destination of the output root files")
    parser.add_argument('--particlename', type=str, default='pi-', help='Particle name used in the simulation')
    parser.add_argument('--use_pt', type=bool, default=False, help='use P/PT for particle generation P by default')
    parser.add_argument('--prange', nargs='+', type=float, default=[0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0], help="Momentum range of the generated particles")
    parser.add_argument('--emin', type=float, default=-0.8, help='start of the eta range')
    parser.add_argument('--emax', type=float, default=0.8, help='stop of the eta range')
    parser.add_argument('--phimin', type=str, default='-M_PI/2.', help='start of the phi range')
    parser.add_argument('--phimax', type=str, default='M_PI/2.', help='stop of the phi range')
    args = parser.parse_args()

    autorun(args)

        

