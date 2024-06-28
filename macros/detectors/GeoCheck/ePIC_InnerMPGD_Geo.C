/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_MPGD.C
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.06.21
*   	Description：	
*
=========================================================================*/
#ifndef FUN4ALL_EPICINNERMPGDSINGLE_C
#define FUN4ALL_EPICINNERMPGDSINGLE_C
#include <iostream>
#include <fstream>
#include <sstream>

#include <epic_innermpgd/ePIC_InnerMPGD_Subsystem.h>
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
R__LOAD_LIBRARY(libePIC_InnerMPGD_Detector.so)
R__LOAD_LIBRARY(libg4histos.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)

int ePIC_InnerMPGD_Geo(const int nEvents = 10000, bool use_pt = false, const double pmin= 10., const double pmax = 10., const double emin = -0.1, const double emax = 0.1, const string& outfile = ""){
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
	double rmpgd = 55.;
    double coat_t = 5.;
    double z_min = -105.;
    double z_max = 143.;
    double KaptonOverlay_t = 50./10000.;
    double CuGround_t = 1.58/10000.;
    double pcb_t = 100./10000.;
    double CuStrip_t = 12./10000.;
    double KaptonStrip_t = 75./10000.;
    double ResistiveStrip_t = 128./10000.;
    double gas_t = 20./10000.;
    double mesh_t = 18./10000.;
    double GasGap_t = 3000./10000.;
    double DriftCuElectrode_t = 5./10000.;
    double DriftKapton_t = 250./10000.;
    double DriftCuGround_t = 0.41/10000.;
    double fudge_t = 570./10000.;

    int stave_number = 128;

	
	ePIC_InnerMPGD_Subsystem* innerMPGD;
	for(int ilayer = 0; ilayer<nlayers; ilayer++){	
		
		innerMPGD = new ePIC_InnerMPGD_Subsystem("InnerMPGD", ilayer);
		innerMPGD->set_double_param("radius", rmpgd);
        innerMPGD->set_double_param("coat_thickness", coat_t);
        innerMPGD->set_double_param("z_min", z_min);
        innerMPGD->set_double_param("z_max", z_max);
        innerMPGD->set_double_param("KaptonOverlay_thickness", KaptonOverlay_t);
        innerMPGD->set_double_param("CuGround_thickness", CuGround_t);
        innerMPGD->set_double_param("PCB_thickness", pcb_t);
        innerMPGD->set_double_param("CuStrip_thickness", CuStrip_t);
        innerMPGD->set_double_param("KaptonStrip_thickness", KaptonStrip_t);
        innerMPGD->set_double_param("ResistiveStrip_thickness", ResistiveStrip_t);
        innerMPGD->set_double_param("gas_thickness", gas_t);
        innerMPGD->set_double_param("mesh_thickness", mesh_t);
        innerMPGD->set_double_param("GasGap_thickness", GasGap_t);
        innerMPGD->set_double_param("DriftCuElectrode_thickness", DriftCuElectrode_t);
        innerMPGD->set_double_param("DriftKapton_thickness", DriftKapton_t);
        innerMPGD->set_double_param("DriftCuGround_thickness", DriftCuGround_t);
        innerMPGD->set_double_param("fudge_thickness", fudge_t);
		innerMPGD->set_int_param("stave_number", stave_number);
		
		innerMPGD->SetActive();
//		innerMPGD->SuperDetector("SVTOB");
		g4Reco->registerSubsystem(innerMPGD);
	}

	// Black hole swallows everything - prevent loopers from returning
	// to inner detectors, length is given by default eta = +- 1.1 range
	PHG4CylinderSubsystem* cyl;
	cyl = new PHG4CylinderSubsystem("BlackHole", 0);
	cyl->set_double_param("radius", rmpgd + coat_t/2. + 0.1);
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
	
	kalman->set_sub_top_node_name("InnerMPGD");
	kalman->set_trackmap_out_name("InnerMPGDTrackMap");
	
	
	kalman->add_phg4hits(
			"G4HIT_InnerMPGD",					//const std::string& phg4hitsNames,
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
