/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_TOFBarrel_Geo.C
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.06.21
*   	Description：	
*
=========================================================================*/
#ifndef FUN4ALL_EPICOUTERMPGDGEO_C
#define FUN4ALL_EPICOUTERMPGDGEO_C
#include <iostream>
#include <fstream>
#include <sstream>

#include <epic_outermpgd/ePIC_OuterMPGD_Subsystem.h>
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
R__LOAD_LIBRARY(libePIC_OuterMPGD_Detector.so)
R__LOAD_LIBRARY(libg4histos.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)

int ePIC_OuterMPGD_Geo(const int nEvents = 10000, bool use_pt = false, const double pmin= 10., const double pmax = 10., const double emin = -0.1, const double emax = 0.1, const string& outfile = ""){
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
	double r_outermpgd = 72.5;
	double module_width = 36.;
	double zmin = -164.5;
	double zmax = 174.5;
	double window_t = 50./10000.;
	double windowgap_t = 2./10.;
	double driftgap_t = 3./10.;
	double foilcu_t = 5./10000.;
	double readoutelectrode_t = 10./10000.;
	double foilkapton_t = 50./10000.;
	double readoutnomex_t = 50./10000.;
	double readoutkapton_t = 50./10000.;
	double pcb_t = 2.8/10.; 
	double carbonfoam_t = 3.5/10.;
	int stave_z = 2;
	int stave_phi = 12;

	
	ePIC_OuterMPGD_Subsystem* OuterMPGD;
	for(int ilayer = 0; ilayer<nlayers; ilayer++){	
		
		OuterMPGD = new ePIC_OuterMPGD_Subsystem("OuterMPGD", ilayer);
		OuterMPGD->set_double_param("radius", r_outermpgd);
		OuterMPGD->set_double_param("module_width", module_width);
		OuterMPGD->set_double_param("zmin", zmin);
		OuterMPGD->set_double_param("zmax", zmax);
		OuterMPGD->set_double_param("window_thickness", window_t);
		OuterMPGD->set_double_param("windowgap_thickness", windowgap_t);
		OuterMPGD->set_double_param("driftgap_thickness", driftgap_t);
		OuterMPGD->set_double_param("foilcu_thickness", foilcu_t);
		OuterMPGD->set_double_param("readoutelectrode_thickness", readoutelectrode_t);
		OuterMPGD->set_double_param("foilkapton_thickness", foilkapton_t);
		OuterMPGD->set_double_param("readoutnomex_thickness", readoutnomex_t);
		OuterMPGD->set_double_param("readoutkapton_thickness", readoutkapton_t);
		OuterMPGD->set_double_param("pcb_thickness", pcb_t);
		OuterMPGD->set_double_param("carbonfoam_thickness", carbonfoam_t);

		OuterMPGD->set_int_param("stave_phi", stave_phi);
		OuterMPGD->set_int_param("stave_z", stave_z);

		
		OuterMPGD->SetActive();
		OuterMPGD->OverlapCheck(1);
//		TOFBarrel->SuperDetector("SVTOB");
		g4Reco->registerSubsystem(OuterMPGD);
	}

	// Black hole swallows everything - prevent loopers from returning
	// to inner detectors, length is given by default eta = +- 1.1 range
	PHG4CylinderSubsystem* cyl;
	cyl = new PHG4CylinderSubsystem("BlackHole", 0);
	cyl->set_double_param("radius", r_outermpgd + 2.5);
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
	
	kalman->set_sub_top_node_name("OuterMPGD");
	kalman->set_trackmap_out_name("OuterMPGDTrackMap");
	
	
	kalman->add_phg4hits(
			"G4HIT_OuterMPGD",					//const std::string& phg4hitsNames,
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
