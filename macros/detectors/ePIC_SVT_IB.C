/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_SVT_IB.C
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.05.28
*   	Description：	
*
=========================================================================*/
#ifndef FUN4ALL_EPICSVTIB_C
#define FUN4ALL_EPICSVTIB_C

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
R__LOAD_LIBRARY(libePIC_SVTIB_Detector.so)

using namespace std;

namespace SVTIB{
	double si_radius[3] = {3.6, 4.8, 12.0};
    double si_mat = 0.05;

	double matrix_length = 0.3571;
	double switch_length = 0.002;
	double backbone_length = 0.006;
	int ntile = 25;

    double lec_length = 0.45;
    double rec_length = 0.15;
    double peri_width = 0.0525; //0;
    double tile_width = 0.9782;

}
namespace Enable{
	bool ePIC_SVTIB = false;
	bool ePIC_SVTIB_OVERLAPCHECK = false;
}  // name space Enable

void ePIC_SVTIB_Init(){

}



void SVTIBFastKalmanFilterConfigSVTX(PHG4TrackFastSim* kalman_filter, int ilay, double radius, bool addproj){

	kalman_filter->add_phg4hits(
			Form("G4HIT_SVTXIB_%d", ilay),
			PHG4TrackFastSim::Cylinder,
			50./10000./sqrt(12.), //radial-resolution [cm]
			20./10000./sqrt(12.), // azimuthal-resolution [cm]
			20./10000./sqrt(12.), // z-resolution [cm]
			1,   //efficiency
			0    // noise
			);
	kalman_filter->add_cylinder_state(Form("SVTXIB_%d", ilay), radius);
	if(addproj)TRACKING::ProjectionNames.insert(Form("SVTXIB_%d", ilay));

}


double ePIC_SVT_IB(PHG4Reco* g4Reco, const int nlayers = 3, double radius = 0){

	if (nlayers > 3){
		cout<<"Silicon VTX layers should not exceed 3!"<<endl;
		exit(1);
	}
	
	bool OverlapCheck = Enable::OVERLAPCHECK || Enable::ePIC_SVTIB_OVERLAPCHECK;

	if(SVTIB::si_radius[0] <= radius){
		cout<<"Geometry overlap happens, please check the radius of each layer!"<<endl;
		exit(-1);
	}

	ePIC_SVTIB_Subsystem* svt_ib;
	for (int ilayer = 0; ilayer < nlayers; ilayer++){
		svt_ib = new ePIC_SVTIB_Subsystem("SVTXIB", ilayer);
		svt_ib->set_double_param("radius", SVTIB::si_radius[ilayer]);
		svt_ib->set_double_param("si_thickness", SVTIB::si_mat/100.*9.37);

        svt_ib->set_double_param("matrix_length", SVTIB::matrix_length);
		svt_ib->set_double_param("switch_length", SVTIB::switch_length);
		svt_ib->set_double_param("backbone_length", SVTIB::backbone_length);
		svt_ib->set_int_param("ntile", SVTIB::ntile);

		svt_ib->set_double_param("lec_length", SVTIB::lec_length);
		svt_ib->set_double_param("rec_length", SVTIB::rec_length);
		svt_ib->set_double_param("periphery_width", SVTIB::peri_width);
		svt_ib->set_double_param("tile_width", SVTIB::tile_width);
    	svt_ib->SetActive();
    //	svt_ib->SuperDetector("SVTX");

		svt_ib->OverlapCheck(OverlapCheck);
    	g4Reco->registerSubsystem(svt_ib);

		SVTIBFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilter, ilayer, SVTIB::si_radius[ilayer], false);//true);
		SVTIBFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilterInnerTrack, ilayer, SVTIB::si_radius[ilayer], false);
		SVTIBFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilterSiliconTrack, ilayer, SVTIB::si_radius[ilayer], false);
	
	}
	
	// update the BlackHole geometry 
	double ib_length = SVTIB::lec_length + SVTIB::ntile*(SVTIB::backbone_length + 3*(SVTIB::matrix_length + SVTIB::switch_length)) + SVTIB::rec_length;
	BlackHoleGeometry::max_radius = SVTIB::si_radius[nlayers-1];
	BlackHoleGeometry::min_z = -ib_length/2.;
	BlackHoleGeometry::max_z = ib_length/2.;
	BlackHoleGeometry::gap = no_overlapp;

	return SVTIB::si_radius[nlayers-1];
}

#endif
