/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		Fun4All_ePIC_SVT_OB.C
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.04.28
*   	Description：	
*
=========================================================================*/

#ifndef MACRO_FUN4ALLEPICSVTOB_C
#define MACRO_FUN4ALLEPICSVTOB_C

#include <iostream>
#include <fstream>
#include <sstream>

#include <epic_svt_ob/ePIC_SVT_OB_Subsystem.h>
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

using namespace std;

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libePIC_SVT_OB_Detector.so)
R__LOAD_LIBRARY(libg4histos.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)

int Fun4All_ePIC_SVT_OB(const int nEvents = 10000, const string& evalfile="FastTrackingEval.root", const string& outfile = ""){
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
	PHG4ParticleGenerator* gen = new PHG4ParticleGenerator("PGENERATOR");
	gen->set_name("pi-");
	gen->set_vtx(0, 0, 0);
	gen->set_eta_range(-0.8, 0.8);
	gen->set_mom_range(10.0, 10.0);
	//gen->set_z_range(0, 0);
	gen->set_phi_range(0, M_PI/2.);
	se->registerSubsystem(gen);
	
	//
	//Geant4 setup
	//
	
	PHG4Reco* g4Reco = new PHG4Reco();
	g4Reco->set_field(1.7);  // Telsa
	
	//g4Reco->save_DST_geometry(false);
	// try non default physics lists
	// g4Reco->SetPhysicsList("FTFP_BERT_HP")

	//
	//build SVT IB layers
	//
	const int ib_layers = 3;
	double si_thickness = 0.005;  // cm by default
	double svxrad[ib_layers] = {3.6, 4.8, 12.0};
	double length[ib_layers] = {27., 27., 27.};  // -1 use eta coverage to determine length
	PHG4CylinderSubsystem *cyl;
	// here is our silicon:
	for (int ilayer = 0; ilayer < ib_layers; ilayer++)
	{
		cyl = new PHG4CylinderSubsystem("SVTX", ilayer);
		cyl->set_double_param("radius", svxrad[ilayer]);
		cyl->set_string_param("material", "G4_Si");  // Silicon (G4 definition)
		cyl->set_double_param("thickness", si_thickness);
		cyl->SetActive();
		cyl->SuperDetector("SVTX");
		if (length[ilayer] > 0)
		{
		  cyl->set_double_param("length", length[ilayer]);
		}
		g4Reco->registerSubsystem(cyl);
	  }
		
	//
	// build my SVT OB layers
	//
	// L3/L4 OB
	double r_inner[2] = {27.1, 41.8};
	double r_outer[2] = {27.7, 42.4};
	double carbon_thickness = 0.5*0.03;
   	double carbon_length[2] = {54.31, 83.75};
	double carbon_width = 3.92;
	//double si_thickness = 0.005;
	si_thickness = 0.005;
	double si_length[2] = {15.0529, 12.8880};
	double si_width = 3.3128;
	double si_carbon_gap = 0.5;
	int n_silicon_z[2] = {4, 8};
	int n_stave_phi[2] = {46, 70};
	double las_overlap[2] = {2.55, 2.98};
	const int ob_layers = 2;
	
	ePIC_SVT_OB_Subsystem* svt_ob;
	for(int i = 0; i<ob_layers; i++){	
		
		svt_ob = new ePIC_SVT_OB_Subsystem("SVTOB", i);
		
		// set the parameters
		svt_ob->set_double_param("r_inner", r_inner[i]);
		svt_ob->set_double_param("r_outer", r_outer[i]);
		// carbon support 
		svt_ob->set_double_param("carbon_thickness", carbon_thickness);
		svt_ob->set_double_param("carbon_length", carbon_length[i]);
		svt_ob->set_double_param("carbon_width", carbon_width);
		// LAS geo 
		svt_ob->set_double_param("si_thickness", si_thickness);
		svt_ob->set_double_param("si_length", si_length[i]);
		svt_ob->set_double_param("si_width", si_width);
		svt_ob->set_double_param("si_carbon_gap", si_carbon_gap);
		svt_ob->set_int_param("n_silicon_z", n_silicon_z[i]);
		svt_ob->set_int_param("n_stave_phi", n_stave_phi[i]);

		// overlaps of LAS
		
		svt_ob->set_double_param("las_overlap", las_overlap[i]);
		
		svt_ob->SetActive();
//		svt_ob->SuperDetector("SVTOB");
		
		g4Reco->registerSubsystem(svt_ob);
	
	}
	


	// Black hole swallows everything - prevent loopers from returning
	// to inner detectors, length is given by default eta = +- 1.1 range
//	PHG4CylinderSubsystem* cyl;
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
	
	
	kalman->add_phg4hits(
			"G4HIT_SVTX",					//const std::string& phg4hitsNames,
			PHG4TrackFastSim::Cylinder, 	// const DETECTOR_TYPE phg4dettype,	Vertical_Plane/Cylinder
			50e-4,							//radial-resolution [cm]
			5e-4, 							//azimuthal-resolution [cm]
			5e-4, 							//z-resolution [cm]
			1, 								//efficiency
			0 								//noise hits
			);
	
	// SVT OB
	ostringstream oss;
	// add si tracker
	//PHG4TrackFastSim* kalman; //= new PHG4TrackFastSim("PHG4TrackFastSim");
	for(int i = 0; i < ob_layers; i++){
		//kalman = new PHG4TrackFastSim("PHG4TrackFastSim");
		//kalman->set_use_vertex_in_fitting(false);
		//oss.str("");
		//oss<<"SVTOB_"<<i;
		//kalman->set_sub_top_node_name(oss.str().c_str());
		//oss.str("");
		//oss<<"SvtobTrackMap_"<<i;
		//kalman->set_trackmap_out_name(oss.str().c_str());

		oss.str("");
		oss<<"G4HIT_SVTOB_"<<i;
		kalman->add_phg4hits(
				oss.str().c_str(),					//const std::string& phg4hitsNames,
				PHG4TrackFastSim::Cylinder, 	// const DETECTOR_TYPE phg4dettype,	Vertical_Plane/Cylinder
				50e-4,							//radial-resolution [cm]
				5e-4, 							//azimuthal-resolution [cm]
				5e-4, 							//z-resolution [cm]
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


	// this (dummy) input manager just drives the event loop
//	Fun4AllInputManager* in = new Fun4AllDummyInputManager("Dummy");
//	se->registerInputManager(in);
	
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
