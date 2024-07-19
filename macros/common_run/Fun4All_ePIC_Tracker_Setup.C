/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		Fun4All_ePIC_Tracker_Setup.C
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.04.28
*   	Description：	
*
=========================================================================*/

#ifndef MACRO_FUN4ALLEPICTRACKERSETUP_C
#define MACRO_FUN4ALLEPICTRACKERSETUP_C

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <epic_svtib/ePIC_SVTIB_Subsystem.h>
#include <epic_svt_ob/ePIC_SVT_OB_Subsystem.h>
#include <epic_innermpgd/ePIC_InnerMPGD_Subsystem.h>
#include <epic_tofbarrel/ePIC_TOFBarrel_Subsystem.h>
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
#include "../detectors/Fun4All_particle_generator.C"
#include "../detectors/ePIC_Tracking.C"
#include "../detectors/ePIC_SVT_Pipe.C"
#include "../detectors/ePIC_SVT_IB.C"
#include "../detectors/ePIC_SVT_OB.C"
#include "../detectors/ePIC_InnerMPGD.C"
#include "../detectors/ePIC_TOFBarrel.C"
#include "../detectors/ePIC_OuterMPGD.C"
#include <GlobalVariables.C>
#include "../detectors/ePIC_BlackHole.C"
#include "../detectors/ePIC_World.C"
#include <G4_Input.C>
#include <G4_Global.C>

using namespace std;

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libePIC_SVTIB_Detector.so)
R__LOAD_LIBRARY(libePIC_SVT_OB_Detector.so)
R__LOAD_LIBRARY(libePIC_InnerMPGD_Detector.so)
R__LOAD_LIBRARY(libePIC_TOFBarrel_Detector.so)
R__LOAD_LIBRARY(libePIC_OuterMPGD_Detector.so)
R__LOAD_LIBRARY(libg4histos.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)

void ePICInit(){
		//particle generator check 
		if(Enable::IP6 and Enable::IP8){
			cout<<"Can not enable Enable::IP6 and Enable::IP8 at the same time!"<<endl;
			gSystem->Exit(1);
		}

		if(Enable::TRACKING) TrackingInit();
		//init the beam pipe
		if(Enable::PIPE) PipeInit();
		// init the SVT inner barrel
		if(Enable::ePIC_SVTIB) ePIC_SVTIB_Init();
		//init the SVT outer barrel
		if(Enable::ePIC_SVTOB) ePIC_SVTOB_Init();
		//init the Inner MPGD barrel
		if(Enable::ePIC_InnerMPGD) ePIC_InnerMPGD_Init();
		//init the TOFBarrel
		if(Enable::ePIC_TOFBarrel) ePIC_TOFBarrel_Init();
		//init the OuterMPGD
		if(Enable::ePIC_OuterMPGD) ePIC_OuterMPGD_Init();

}

void ePICSetup(){
	
	Fun4AllServer* se = Fun4AllServer::instance();
	PHG4Reco* g4Reco = new PHG4Reco();
	g4Reco->set_field(1.7); //Telsa
	WorldInit(g4Reco);

	//g4Reco->set_rapidity_coverage(1.1);
	if(G4P6DECAYER::decayType != EDecayType::kAll){
		g4Reco->set_force_decay(G4P6DECAYER::decayType);
	}
	
	double radius = 0.;
	// beam pipe
	if(Enable::PIPE) radius = Pipe(g4Reco, radius);
	//tracking service

	// SVT IB layers
	if(Enable::ePIC_SVTIB) radius = ePIC_SVT_IB(g4Reco, 3, radius);

	//SVT OB layers
	if(Enable::ePIC_SVTOB) radius = ePIC_SVT_OB(g4Reco, 2, radius);

	//Inner MPGD layers
	if(Enable::ePIC_InnerMPGD) radius = ePIC_InnerMPGD(g4Reco, 1, radius);

	// TOF barrel layers
	if(Enable::ePIC_TOFBarrel) radius = ePIC_TOFBarrel(g4Reco, 1, radius);

	//Outer MPGD layer
	if(Enable::ePIC_OuterMPGD) radius = ePIC_OuterMPGD(g4Reco, 1, radius);
	
	//BlackHole if enabled, needs infor from all previus sub detectors for dimensions
	if(Enable::BLACKHOLE) BlackHole(g4Reco, radius);

	PHG4TruthSubsystem* truth = new PHG4TruthSubsystem();
	g4Reco->registerSubsystem(truth);
	//finally adjust the world size in case the default is too small
	WorldSize(g4Reco, radius);

	se->registerSubsystem(g4Reco);
	return ;
}



#endif
