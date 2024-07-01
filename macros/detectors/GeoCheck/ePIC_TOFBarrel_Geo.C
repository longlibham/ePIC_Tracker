/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_TOFBarrel_Geo.C
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.06.21
*   	Description：	
*
=========================================================================*/
#ifndef FUN4ALL_EPICTOFBARRELGEO_C
#define FUN4ALL_EPICTOFBARRELGEO_C
#include <iostream>
#include <fstream>
#include <sstream>

#include <epic_tofbarrel/ePIC_TOFBarrel_Subsystem.h>
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


R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libePIC_TOFBarrel_Detector.so)
R__LOAD_LIBRARY(libg4histos.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)

int ePIC_TOFBarrel_Geo(const int nEvents = 10000, bool use_pt = false, const double pmin= 10., const double pmax = 10., const double emin = -0.1, const double emax = 0.1, const string& outfile = ""){
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

	PHG4Reco* g4Reco = new PHG4Reco();
	g4Reco->set_field(1.7);  // Telsa

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

	//
	// build my MPGD layers
	//
	int nlayers = 1;
	double zmin = -112.5;
	double zmax = 174.0;
	double rtof = 64.6;
	double sensor_width = 3.2;
	double module_width = 5.6;
	double coolingtube_width = 0.75;
	double sensor_t = 0.03;
	double hybrid_t = 2*0.008125;
	double CFskin_t = 2*0.0075;
	double CoolingTube_t = 0.08;
	double Coolant_t = 0.08;
	double CFoam_t = 2*0.29;
	double CHoneycomb_t = 2*0.29;

    int stave_number = 144;

	
	ePIC_TOFBarrel_Subsystem* TOFBarrel;
	for(int ilayer = 0; ilayer<nlayers; ilayer++){	
		
		TOFBarrel = new ePIC_TOFBarrel_Subsystem("TOFBarrel", ilayer);
		TOFBarrel->set_double_param("radius", rtof);
		TOFBarrel->set_double_param("zmin", zmin);
		TOFBarrel->set_double_param("zmax", zmax);
		TOFBarrel->set_double_param("sensor_width", sensor_width);
		TOFBarrel->set_double_param("module_width", module_width);
		TOFBarrel->set_double_param("coolingtube_width", coolingtube_width);
		TOFBarrel->set_double_param("sensor_thickness", sensor_t);
		TOFBarrel->set_double_param("hybrid_thickness", hybrid_t);
		TOFBarrel->set_double_param("CFskin_thickness", CFskin_t);
		TOFBarrel->set_double_param("CoolingTube_thickness", CoolingTube_t);
		TOFBarrel->set_double_param("Coolant_thickness", Coolant_t);
		TOFBarrel->set_double_param("CFoam_thickness", CFoam_t);
		TOFBarrel->set_double_param("CHoneycomb_thickness", CHoneycomb_t);

		TOFBarrel->set_int_param("stave_number", stave_number);
		
		TOFBarrel->SetActive();
		TOFBarrel->OverlapCheck(1);
//		TOFBarrel->SuperDetector("SVTOB");
		g4Reco->registerSubsystem(TOFBarrel);
	}

	// Black hole swallows everything - prevent loopers from returning
	// to inner detectors, length is given by default eta = +- 1.1 range
	PHG4CylinderSubsystem* cyl;
	cyl = new PHG4CylinderSubsystem("BlackHole", 0);
	cyl->set_double_param("radius", rtof + module_width + 0.1);
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
	
	kalman->set_sub_top_node_name("TOFBarrel");
	kalman->set_trackmap_out_name("TOFBarrelTrackMap");
	
	
	kalman->add_phg4hits(
			"G4HIT_TOFBarrel",					//const std::string& phg4hitsNames,
			PHG4TrackFastSim::Cylinder, 	// const DETECTOR_TYPE phg4dettype,	Vertical_Plane/Cylinder
			0.3,							//radial-resolution [cm]
			0.015, 							//azimuthal-resolution [cm]
			0.015, 							//z-resolution [cm]
			1, 								//efficiency
			0 								//noise hits
			);
	
	
	se->registerSubsystem(kalman);
	se->Print("NODETREE");	
	
	///////////////////////////////////////
	//IO managers
	//////////////////////////////////////
	PHG4TrackFastSimEval* fast_sim_eval = new PHG4TrackFastSimEval("FastTrackingEval");
	ostringstream oss;
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
