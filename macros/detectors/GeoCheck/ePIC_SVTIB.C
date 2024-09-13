/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_SVTIB.cxx
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.07.19
*   	Description：	
*
=========================================================================*/


#ifndef MACRO_FUN4ALLEPICSVTIB_C
#define MACRO_FUN4ALLEPICSVTIB_C

#include <iostream>
#include <fstream>
#include <sstream>

#include <epic_svtib/ePIC_SVTIB_Subsystem.h>
#include <g4detectors/PHG4DetectorSubsystem.h>
#include <g4histos/G4HitNtuple.h>

#include <g4main/PHG4ParticleGenerator.h>
#include <g4main/PHG4ParticleGeneratorBase.h>
#include <g4main/PHG4Reco.h>
#include <g4main/PHG4SimpleEventGenerator.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllDummyInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>

#include <phool/recoConsts.h>

#include <g4detectors/PHG4CylinderSubsystem.h>

#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <g4trackfastsim/PHG4TrackFastSimEval.h>

#include <g4main/PHG4TruthSubsystem.h>


// user defined
#include "../Fun4All_particle_generator.C"
#include <GlobalVariables.C>
#include <G4_Input.C>
#include <G4_Global.C>

using namespace std;

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libePIC_SVTIB_Detector.so)
R__LOAD_LIBRARY(libg4histos.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)

int ePIC_SVTIB(const int nEvents = 10000, bool use_pt = false, const double pmin= 10., const double pmax = 10., const double emin = -0.1, const double emax = 0.1, const string& outfile = ""){
	if (pmin > pmax){
		cout<<"pmin should be less than pmax, please choose the correct momentum range..."<<endl;
		gSystem->Exit(-1);
	}

	////////////////////////////////////
	//Make the server
	///////////////////////////////////
	Fun4AllServer* se = Fun4AllServer::instance();
	se->Verbosity(0);
	recoConsts* rc = recoConsts::instance();
	
	// particle generator according to Fun4All_particle_generator.C    --Long LI 2024-05-14
	if (use_pt)
		Enable::USE_PT = true;
	else
		Enable::USE_P = true;
	Input::SIMPLE = true;
	Fun4All_particle_generator(nEvents, "pi-", pmin, pmax, emin, emax, -M_PI/2., M_PI/2.);
	
	//
	//Geant4 setup
	//
	
	PHG4Reco* g4Reco = new PHG4Reco();
	g4Reco->set_field(1.7);  // Telsa
		
	//
	// build my SVT IB layers
	//
	int ib_layers = 3;
    double si_radius[3] = {3.6, 4.8, 12.0};
    double si_thickness = 50e-4;

	double matrix_length = 0.3571;
	double switch_length = 0.002;
	double backbone_length = 0.006;
	int ntile = 25;

    double lec_length = 0.45;
    double rec_length = 0.15;
    double periphery_width = 0.0525;
    double tile_width = 0.9782;
	
	ePIC_SVTIB_Subsystem* svt_ib;
	for(int i = 0; i<ib_layers; i++){	
		
		svt_ib = new ePIC_SVTIB_Subsystem("SVTIB", i);
		
		// set the parameters
		svt_ib->set_double_param("radius", si_radius[i]);
        svt_ib->set_double_param("si_thickness", si_thickness);
		
        svt_ib->set_double_param("matrix_length", matrix_length);
		svt_ib->set_double_param("switch_length", switch_length);
		svt_ib->set_double_param("backbone_length", backbone_length);
		svt_ib->set_int_param("ntile", ntile);
		
        svt_ib->set_double_param("lec_length", lec_length);
        svt_ib->set_double_param("rec_length", rec_length);
        svt_ib->set_double_param("periphery_width", periphery_width);
        svt_ib->set_double_param("tile_width", tile_width);


		
		svt_ib->SetActive();
        svt_ib->OverlapCheck(1);

		g4Reco->registerSubsystem(svt_ib);
	
	}
	


	// Black hole swallows everything - prevent loopers from returning
	// to inner detectors, length is given by default eta = +- 1.1 range
	PHG4CylinderSubsystem* cyl;
	cyl = new PHG4CylinderSubsystem("BlackHole", 0);
	cyl->set_double_param("radius", 80);
	cyl->set_double_param("thickness", 0.1);
	cyl->SetActive();
	cyl->BlackHole();
	g4Reco->registerSubsystem(cyl);


	PHG4TruthSubsystem* truth = new PHG4TruthSubsystem();
	g4Reco->registerSubsystem(truth);

	se->registerSubsystem(g4Reco);

	//
	//fast pattern recognition and full Kalman filter
	//output evaluation file for truth track and reco tracks are PHG4TruthInfoContainer
	//
	
	// SVT IB
	PHG4TrackFastSim* kalman = new PHG4TrackFastSim("PHG4TrackFastSim");
	kalman->Verbosity(0);
	kalman->set_use_vertex_in_fitting(false);
	
	kalman->set_sub_top_node_name("SVTX");
	kalman->set_trackmap_out_name("SvtxTrackMap");
	
	
	// SVT IB
	ostringstream oss;
	for(int i = 0; i < ib_layers; i++){
		oss.str("");
		oss<<"G4HIT_SVTIB_"<<i;
		kalman->add_phg4hits(
				oss.str().c_str(),					//const std::string& phg4hitsNames,
				PHG4TrackFastSim::Cylinder, 	// const DETECTOR_TYPE phg4dettype,	Vertical_Plane/Cylinder
				50e-4,							//radial-resolution [cm]
				10./10000./sqrt(12.), 							//azimuthal-resolution [cm]
				10./10000./sqrt(12.), 							//z-resolution [cm]
				1, 								//efficiency
				0 								//noise hits
				);
	}
	
	se->registerSubsystem(kalman);
	se->Print("NODETREE");	
	
	///////////////////////////////////////
	//IO managers
	//////////////////////////////////////
	PHG4TrackFastSimEval* fast_sim_eval = new PHG4TrackFastSimEval("FastTrackingEval");
	oss.str("");
	if(Enable::USE_PT)
		oss<<"../data/fixed_pT/FastTrackingEval_"<<pmin<<"-"<<pmax<<"GeV.root";
	else if(Enable::USE_P)
		oss<<"../data/fixed_p/FastTrackingEval_"<<pmin<<"-"<<pmax<<"GeV.root";
	
	string evalfile = oss.str();
	fast_sim_eval->set_filename(evalfile);
	se->registerSubsystem(fast_sim_eval);

	//
	//output DST file for further offline analysis
	//
	if(!outfile.empty()){
		Fun4AllOutputManager* out = new Fun4AllDstOutputManager("DSTOUT", outfile);
		se->registerOutputManager(out);
	}
	Fun4AllInputManager* in = new Fun4AllDummyInputManager("JADE");
	se->registerInputManager(in);

	// events = 0 => run forever
	if (nEvents <= 0){
		return 0;
	}

	se->run(nEvents);
	//svt_ob->Print();
	se->End();
	cout<< "All Done!!!" <<endl;
	delete se;
	gSystem->Exit(0);

	return 0;

}

PHG4ParticleGenerator* get_gen(const char* name="PGENERATOR"){
	Fun4AllServer* se = Fun4AllServer::instance();
	PHG4ParticleGenerator* pgun = (PHG4ParticleGenerator*) se->getSubsysReco(name);
	return pgun;
}

#endif