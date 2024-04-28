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

using namespace std;

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libePIC_SVT_OB_Detector.so)
R__LOAD_LIBRARY(libg4histos.so)

void Fun4All_ePIC_SVT_OB(int nEvents = 1){
	////////////////////////////////////
	//Make the server
	///////////////////////////////////
	Fun4AllServer* se = Fun4AllServer::instance();
	recoConsts* rc = recoConsts::instance();
	// if you want to fix the random seed to reproduce results
	// set this flag
	// nail this down so I know what the first event looks like...
	// 	rc->Set_IntFlag("RANDOMSEED", 12345);
	

	//
	// particle generator
	//
	
	//the PHG4ParticleGenerator makes cones using phi and eta 
	PHG4ParticleGenerator* gen = new PHG4ParticleGenerator();
	gen->set_name("pi-");
	gen->set_vtx(0, 0, 0);
	gen->set_eta_range(-0.8, 0.8);
	gen->set_mom_range(10.0, 10.0);
	//gen->set_z_range(0, 0);
	gen->set_phi_range(0, TMath::Pi());
	se->registerSubsystem(gen);
	
	//
	//Geant4 setup
	//
	
	PHG4Reco* g4Reco = new PHG4Reco();
	g4Reco->Set_field(1.7);  // Telsa
	g4Reco->save_DST_geometry(false);
	// try non default physics lists
	// g4Reco->SetPhysicsList("FTFP_BERT_HP")
	
	//
	// build my SVT OB layers
	//
	
	ePIC_SVT_OB_Subsystem* svt_ob = new ePIC_SVT_OB_Subsystem("SVT_OB_0", 0);
	g4Reco->registerSubsystem(svt_ob);

	se->registerSubsystem(g4Reco);


	/////////////////////////////////////////
	//Fun4All module
	////////////////////////////////////////
	
	G4HitNtuple* hits = new G4HitNtuple("hits");
	hits->AddNode("SVT_OB", 0);
	se->registerSubsystem(hits);

	///////////////////////////////////////
	//IO managers
	//////////////////////////////////////
	

	// this (dummy) input manager just drives the event loop
	Fun4AllInputManager* in = new Fun4AllDummyInputManager("Dummy");
	se->registerInputManager(in);

	// events = 0 => run forever
	if (nEvents <= 0){
		return 0;
	}

	se->run(nEvents);
	svt_ob->Print();
	se->End();
	cout<< "All Done!!!" <<endl;
	delete se;
	gSystem->Exit(0);

}

#endif
