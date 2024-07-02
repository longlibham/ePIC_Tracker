/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_TOFBarrel.cxx
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.07.02
*   	Description：	
*
=========================================================================*/
#ifndef FUN4ALL_EPICTOFBARREL_C
#define FUN4ALL_EPICTOFBARREL_C

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
R__LOAD_LIBRARY(libePIC_TOFBarrel_Detector.so)

using namespace std;

namespace Enable{
	bool ePIC_TOFBarrel = false;
	bool ePIC_TOFBarrel_OVERLAPCHECK = false;
}  // name space Enable

void ePIC_TOFBarrel_Init(){

}



void TOFBarrelFastKalmanFilterConfig(PHG4TrackFastSim* kalman_filter, int ilay, double radius, bool addproj){

	kalman_filter->add_phg4hits(
			Form("G4HIT_TOFBarrel_%d", ilay),
			PHG4TrackFastSim::Cylinder,
			0.03/sqrt(12), //radial-resolution [cm]
			0.01/sqrt(12), // azimuthal-resolution [cm]
			1./sqrt(12), // z-resolution [cm]
			1,   //efficiency
			0    // noise
			);
	kalman_filter->add_cylinder_state(Form("TOFBarrel_%d", ilay), radius);
	if(addproj)TRACKING::ProjectionNames.insert(Form("TOFBarrel_%d", ilay));
}


double ePIC_TOFBarrel(PHG4Reco* g4Reco, const int nlayers = 1, double radius = 0){

	if (nlayers > 1){
		cout<<"TOFBarrel layer should not exceed 1 for the moment!"<<endl;
		exit(1);
	}
	
	bool OverlapCheck = Enable::OVERLAPCHECK || Enable::ePIC_TOFBarrel_OVERLAPCHECK;

    //
	// build TOF layers layers
	//
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
	double CFoam_t = 2*0.265;
	double CHoneycomb_t = 2*0.265;

    int stave_number = 144;


	if(rtof <= radius){
		cout<<"Geometry overlap happens, please check the radius of each layer!"<<endl;
		exit(-1);
	}

	
	
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
//		innerMPGD->SuperDetector("SVTOB");

		TOFBarrel->OverlapCheck(OverlapCheck);
		g4Reco->registerSubsystem(TOFBarrel);

		TOFBarrelFastKalmanFilterConfig(TRACKING::FastKalmanFilter, ilayer, rtof, false);//true);
		
	}

	// update the BlackHole geometry 
    double r_boundary = rtof + sensor_t + 2*hybrid_t + 2*CFskin_t + CoolingTube_t + Coolant_t + CFoam_t + 1.; // add 1 cm due to tilt
	BlackHoleGeometry::max_radius = r_boundary;
	BlackHoleGeometry::min_z = zmin;
	BlackHoleGeometry::max_z = zmax;
	BlackHoleGeometry::gap = no_overlapp;

	return r_boundary;
}

#endif