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


// user defined
#include "../Fun4All_particle_generator.C"
#include <GlobalVariables.C>
#include <G4_Input.C>
#include <G4_Global.C>

using namespace std;

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libePIC_SVT_OB_Detector.so)
R__LOAD_LIBRARY(libg4histos.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)

int Fun4All_ePIC_SVT_OB(const int nEvents = 10000, bool use_pt = false, const double pmin= 10., const double pmax = 10., const double emin = -0.1, const double emax = 0.1, const string& outfile = ""){
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
	Fun4All_particle_generator(nEvents, "pi-", pmin, pmax, emin, emax, -M_PI/2., M_PI/2.);
	
	//
	//Geant4 setup
	//
	
	PHG4Reco* g4Reco = new PHG4Reco();
	g4Reco->set_field(1.7);  // Telsa
	
	//g4Reco->save_DST_geometry(false);
	// try non default physics lists
	// g4Reco->SetPhysicsList("FTFP_BERT_HP")

		
	//
	// build my SVT OB layers
	//
	// L3/L4 OB
	
	double r_inner[2] = {26.475, 41.675};//{22., 43.8};//{26.5, 39.5};
	double r_outer[2] = {27.925, 43.125};//{23.2, 45.};//{27.7, 40.7};

	double carbon_thickness = 100./10000.; //0.2/100*26.356;
	double carbon_length[2] = {53.2, 89.064};//{44.532, 89.064};//{53.2, 79.8};
	double carbon_width = 3.92;
	double cf_leftendcap = 1.2;
	double cf_rightendcap = 1.2;
	double cf_margin = 0.1; //carbon_thickness;
	double cf_curve_radius = 9.0115;
	double cf_center_height = 0.85;
	double cf_edge_height = 0.351;
	double oring_radius = 0.25;
	double oring_spacing = 0.3;
	double cf_csupport_width = 0.5;
	double k9_center_height = 0.35;
	double si_thickness = 0.05/100*9.37;
	
	double si_width = 3.9128;
	double matrix_length = 0.3571;
	double switch_length = 0.002;
	double backbone_length = 0.006;
	int nmatrix = 3;
	int ntile[2] = {12, 10};//{12, 12};

	double n_silicon_z[2] = {4, 8};//{4, 6};
	double n_stave_phi[2] = {46, 70};
	double lec_length = 0.45;
	double rec_length = 0.15;
	double anc_length = 1.0;
	double anc_thickness = 300./10000.;
	double las_airspace = 0.6;
	double peri_width = 0.0525;
	double kapton_thickness = 50./10000.;
	const int ob_layers = 1;
	
	ePIC_SVT_OB_Subsystem* svt_ob;
	for(int i = 0; i<ob_layers; i++){	
		
		svt_ob = new ePIC_SVT_OB_Subsystem("SVTOB", i);
		
		// set the parameters
		svt_ob->set_double_param("r_inner", r_inner[i]);
		svt_ob->set_double_param("r_outer", r_outer[i]);
		svt_ob->set_double_param("carbon_thickness", carbon_thickness);
		svt_ob->set_double_param("carbon_length", carbon_length[i]);
		svt_ob->set_double_param("carbon_width", carbon_width);
		svt_ob->set_double_param("cf_leftendcap", cf_leftendcap);
		svt_ob->set_double_param("cf_rightendcap", cf_rightendcap);
		svt_ob->set_double_param("cf_margin", cf_margin);
		svt_ob->set_double_param("cf_curve_radius", cf_curve_radius);
		svt_ob->set_double_param("cf_center_height", cf_center_height);
		svt_ob->set_double_param("cf_edge_height", cf_edge_height);
		svt_ob->set_double_param("oring_radius", oring_radius);
		svt_ob->set_double_param("oring_spacing", oring_spacing);
		svt_ob->set_double_param("cf_csupport_width", cf_csupport_width);
		svt_ob->set_double_param("k9_center_height", k9_center_height);
		svt_ob->set_double_param("si_thickness", si_thickness);

		svt_ob->set_double_param("si_width", si_width);
		svt_ob->set_double_param("matrix_length", matrix_length);
		svt_ob->set_double_param("switch_length", switch_length);
		svt_ob->set_double_param("backbone_length", backbone_length);
		svt_ob->set_int_param("nmatrix", nmatrix);
		svt_ob->set_int_param("ntile", ntile[i]);


		svt_ob->set_int_param("n_silicon_z", n_silicon_z[i]);
		svt_ob->set_int_param("n_stave_phi", n_stave_phi[i]);
		svt_ob->set_double_param("lec_length", lec_length);
		svt_ob->set_double_param("rec_length", rec_length);
		svt_ob->set_double_param("anc_length", anc_length);
		svt_ob->set_double_param("anc_thickness", anc_thickness);
		svt_ob->set_double_param("las_airspace", las_airspace);
		svt_ob->set_double_param("periphery_width", peri_width);
		svt_ob->set_double_param("kapton_thickness", kapton_thickness);

		
		svt_ob->SetActive();
		svt_ob->OverlapCheck(1);
//		svt_ob->SuperDetector("SVTOB");
		
		g4Reco->registerSubsystem(svt_ob);
	
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

	// SVT OB
	ostringstream oss;
	for(int i = 0; i < ob_layers; i++){
		oss.str("");
		oss<<"G4HIT_SVTOB_"<<i;
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
