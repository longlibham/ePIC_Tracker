/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_MPGD.C
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.06.21
*   	Description：	
*
=========================================================================*/
#ifndef FUN4ALL_EPICOUTERMPGD_C
#define FUN4ALL_EPICOUTERMPGD_C

#include <g4detectors/PHG4CylinderSubsystem.h>

#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <g4trackfastsim/PHG4TrackFastSimEval.h>

#include <g4main/PHG4Reco.h>
#include <g4main/PHG4TruthSubsystem.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllDummyInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>

#include "ePIC_Tracking.C"
#include "GlobalVariables.C"

#include <phool/recoConsts.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libePIC_OuterMPGD_Detector.so)

using namespace std;


namespace Enable{
	bool ePIC_OuterMPGD = false;
	bool ePIC_OuterMPGD_OVERLAPCHECK = false;
}  // name space Enable

void ePIC_OuterMPGD_Init(){

}



void OuterMPGDFastKalmanFilterConfig(PHG4TrackFastSim* kalman_filter, int ilay, double radius, bool addproj){

	kalman_filter->add_phg4hits(
			Form("G4HIT_OuterMPGD_%d", ilay),
			PHG4TrackFastSim::Cylinder,
			0.3/sqrt(12), //radial-resolution [cm]
			0.015, // azimuthal-resolution [cm]
			0.015, // z-resolution [cm]
			1,   //efficiency
			0    // noise
			);
	kalman_filter->add_cylinder_state(Form("OuterMPGD_%d", ilay), radius);
	if(addproj)TRACKING::ProjectionNames.insert(Form("OuterMPGD_%d", ilay));
}


double ePIC_OuterMPGD(PHG4Reco* g4Reco, const int nlayers = 1, double radius = 0){

	if (nlayers > 1){
		cout<<"OuterMPGD layer should not exceed 1 for the moment!"<<endl;
		exit(1);
	}
	
	bool OverlapCheck = Enable::OVERLAPCHECK || Enable::ePIC_OuterMPGD_OVERLAPCHECK;

    //
	// build my MPGD layers
	//
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
	int stave_z = 2;
	int stave_phi = 12;

	if(r_outermpgd <= radius){
		cout<<"Geometry overlap happens, please check the radius of each layer!"<<endl;
		exit(-1);
	}

	
	
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

		OuterMPGD->set_int_param("stave_phi", stave_phi);
		OuterMPGD->set_int_param("stave_z", stave_z);
		
		OuterMPGD->SetActive();
//		innerMPGD->SuperDetector("SVTOB");

		OuterMPGD->OverlapCheck(OverlapCheck);
		g4Reco->registerSubsystem(OuterMPGD);

		OuterMPGDFastKalmanFilterConfig(TRACKING::FastKalmanFilter, ilayer, r_outermpgd, false);//true);
		
	}

	// update the BlackHole geometry 
	BlackHoleGeometry::max_radius = r_outermpgd + 5.;  // redundancy for black hole 5. cm
	BlackHoleGeometry::min_z = zmin ;
	BlackHoleGeometry::max_z = zmax ;
	BlackHoleGeometry::gap = no_overlapp;

	return r_outermpgd + 2.5;
}

#endif
