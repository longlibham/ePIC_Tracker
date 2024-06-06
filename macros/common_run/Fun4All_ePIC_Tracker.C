/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		Fun4All_ePIC_Tracker.cxx
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.06.05
*   	Description：	
*
=========================================================================*/

#ifndef MACRO_FUN4ALLEPICTRACKER_C
#define MACRO_FUN4ALLEPICTRACKER_C

#include <GlobalVariables.C>

#include <DisplayOn.C>
#include "Fun4All_ePIC_Tracker_Setup.C"
#include "../detectors/Fun4All_particle_generator.C"
#include <G4_Input.C>
#include <TROOT.h>

#include <fun4all/Fun4AllServer.h>
#include <phool/recoConsts.h>
#include <RooUnblindPrecision.h>

#include <iostream>
#include <fstream>
#include <sstream>

R__LOAD_LIBRARY(libfun4all.so)
using namespace std;

int Fun4All_ePIC_Tracker(
    const int nEvents = 10000,
    bool use_pt = false, 
    const double pmin= 10., 
    const double pmax = 10., 
    const double emin = -0.1, 
    const double emax = 0.1, 
    const double phimin = -M_PI/2.,
    const double phimax = M_PI/2.,
    const int skip = 0,
    const string& outfile = ""
){
    ////////////////////////////////////
	//Make the server
	///////////////////////////////////
	Fun4AllServer* se = Fun4AllServer::instance();
	se->Verbosity(0);
	recoConsts* rc = recoConsts::instance();
	// if you want to fix the random seed to reproduce results
	// set this flag
	// nail th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is down so I know what the first event looks like...
	// 	rc->Set_IntFlag("RANDOMSEED", 12345);
	

	//
	// particle generator
	//
	
	//the PHG4ParticleGenerator makes cones using phi and eta 
//	PHG4ParticleGenerator* gen = new PHG4ParticleGenerator("PGENERATOR");
//	gen->set_name("pi-");
//	gen->set_vtx(0, 0, 0);
//	gen->set_eta_range(emin, emax);
//	gen->set_mom_range(pmin, pmax);
	//gen->set_z_range(0, 0);
//	gen->set_phi_range(0, M_PI/2.);
//	se->registerSubsystem(gen);

	// particle generator according to Fun4All_particle_generator.C    --Long LI 2024-05-14
	if (use_pt)
		Enable::USE_PT = true;
	else
		Enable::USE_P = true;
	Input::SIMPLE = true;
	Fun4All_particle_generator(nEvents, "pi-", pmin, pmax, emin, emax, phimin, phimax);
	
	//
	//Geant4 setup
	//
	
	PHG4Reco* g4Reco = new PHG4Reco();
	g4Reco->set_field(1.7);  // Telsa

    //enable beam pipe

    //enable SVT_IB
    Enable::ePIC_SVTIB = true;
    Enable::ePIC_SVTIB_OVERLAPCHECK = true;

    //enable SVT_OB

    //enable blackhole
    Enable::BLACKHOLE = true;
    BlackHoleGeometry::visible = true;

    //enable tracking
    Enable::TRACKING = true;
    Enable::TRACKING_EVAL = Enable::TRACKING && true;
    G4TRACKING::DISPLACED_VERTEX = true; // this option exclude vertex in the track fitting and use RAVE to reconstruct primary and 2
                                        //projections to calorimeters
    


    //Initialsize the selected sybsystems
    ePICInit();

    if (!Input::READHITS){
        ePICSetup();
    }

    //
    // Tracking and PID
    //

    if(Enable::TRACKING) Tracking_Reco();
    ostringstream oss;
    oss.str("");
	if(Enable::USE_PT)
		oss<<"../data/fixed_pT/FastTrackingEval_"<<pmin<<"-"<<pmax<<"GeV.root";
	else if(Enable::USE_P)
		oss<<"../data/fixed_p/FastTrackingEval_"<<pmin<<"-"<<pmax<<"GeV.root";

    if(Enable::TRACKING_EVAL) Tracking_Eval(oss.str().c_str());

    if(!outfile.empty() && Enable::DSTOUT){
        Fun4AllOutputManager* out = new Fun4AllDstOutputManager("DSTOUT", outfile);
		se->registerOutputManager(out);
    }

    //
    // Set up Input Managers
    InputManagers();

    //
    // Event processing
    //

    if(Enable::DISPLAY){
        DisplayOn();
        gROOT->ProcessLine("Fun4AllServer* se = Fun4AllServer::instance();");
        gROOT->ProcessLine("PHG4Reco* g4 = (PHG4Reco*) se->getSubsysReco(\"PHG4RECO\");");

        cout<<"------------------------------------------------------------------------"<<endl;
        cout<<"You are in event display mode. Run one event with"<<endl;
        cout<<"se->run(1)"<<endl;
        cout<<"Run Geant4 command with following examples"<<endl;
        gROOT->ProcessLine("displaycmd()");
        return 0;

    }

    if(nEvents < 0){
        return 0;
    }

    // if we run any of the particle generators and use 0 it'll run forever
    if(nEvents == 0 && !Input::READHITS && !Input::HEPMC && !Input::READEIC){

        cout<<"using 0 for number of events is a bad idea when using particle generators"<<endl;
        cout<<"it will run forever, so I just return without running anything"<<endl;
        return 0;
    }
    se->skip(skip);
    se->run(nEvents);

    //
    // Exit
    //

    se->End();
    cout<<"All Done!"<<endl;

    delete se;

    gSystem->Exit(0);
    return 0;
}
#endif