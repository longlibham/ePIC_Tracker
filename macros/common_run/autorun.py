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

def autorun(args):
    prange = range(args.pmin, args.pmax)
    for p in prange:
        use_pt = str(args.use_pt).lower()
        cmd = f"root -b -q Fun4All_ePIC_Tracker.C\({args.evts}\,{use_pt}\,{p}\,{p}\,{args.emin}\,{args.emax}\)"
        print(cmd)
        os.system(cmd)





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Automatic run for Fun4All_ePIC_SVT_OB')
    parser.add_argument('--pmin', type=int, default=1, help='start of the momentum range')
    parser.add_argument('--pmax', type=int, default=1, help='stop of the momentum range')
    parser.add_argument('--emin', type=float, default=-0.1, help='start of the eta range')
    parser.add_argument('--emax', type=float, default=0.1, help='stop of the eta range')
    parser.add_argument('--evts', type=int, default=10000, help='events for each run')
    parser.add_argument('--use_pt', type=bool, default=False, help='use P/PT for particle generation P by default')
    print("works here")
    args = parser.parse_args()

    autorun(args)

        

