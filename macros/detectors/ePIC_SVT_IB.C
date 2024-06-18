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
#include "ePIC_SVT_IB.C"
#include "GlobalVariables.C"

#include <phool/recoConsts.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

R__LOAD_LIBRARY(libg4detectors.so)

using namespace std;

namespace Enable{
	bool ePIC_SVTIB = false;
	bool ePIC_SVTIB_OVERLAPCHECK = false;
}  // name space Enable

void ePIC_SVTIB_Init(){

}



void SVTIBFastKalmanFilterConfigSVTX(PHG4TrackFastSim* kalman_filter, int ilay, double radius, bool addproj){

	kalman_filter->add_phg4hits(
			Form("G4HIT_SVTX_%d", ilay),
			PHG4TrackFastSim::Cylinder,
			50./10000., //radial-resolution [cm]
			10./10000./sqrt(12.), // azimuthal-resolution [cm]
			10./10000./sqrt(12.), // z-resolution [cm]
			1,   //efficiency
			0    // noise
			);
	kalman_filter->add_cylinder_state(Form("SVTX_%d", ilay), radius);
	if(addproj)TRACKING::ProjectionNames.insert(Form("SVTX_%d", ilay));

}


double ePIC_SVT_IB(PHG4Reco* g4Reco, const int nlayers = 3, double radius = 0){

	if (nlayers > 5){
		cout<<"Silicon VTX layers should not exceed 5!"<<endl;
		exit(1);
	}
	
	bool OverlapCheck = Enable::OVERLAPCHECK || Enable::ePIC_SVTIB_OVERLAPCHECK;

	PHG4CylinderSubsystem* cyl(nullptr);

	// for ePIC-SVT 5 layers of silicon
  	double si_mat[5] = {0.05, 0.05, 0.05, 0.25, 0.55}; 
  	double svxrad[5] = {3.6, 4.8, 12.0, 27.0, 42.0};
  	double length[5] = {27., 27., 27., 54., 84.};  // -1 use eta coverage to determine length
	if(svxrad[0] <= radius){
		cout<<"Geometry overlap happens, please check the radius of each layer!"<<endl;
		exit(-1);
	}

	double gap = 0.1;	

	for (int ilayer = 0; ilayer < nlayers; ilayer++){
    	cyl = new PHG4CylinderSubsystem("SVTX", ilayer);
    	cyl->set_double_param("radius", svxrad[ilayer]);
    	cyl->set_string_param("material", "G4_Si");  // Silicon (G4 definition)
    	cyl->set_double_param("thickness", si_mat[ilayer]/100.*9.37);
		cyl->set_double_param("place_z", 0.);
    	cyl->SetActive();
    //	cyl->SuperDetector("SVTX");
    	if (length[ilayer] > 0){
      		cyl->set_double_param("length", length[ilayer]);
    	}
		
		cyl->OverlapCheck(OverlapCheck);
    	g4Reco->registerSubsystem(cyl);

		SVTIBFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilter, ilayer, svxrad[ilayer], false);//true);
		SVTIBFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilterInnerTrack, ilayer, svxrad[ilayer], false);
		SVTIBFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilterSiliconTrack, ilayer, svxrad[ilayer], false);
	
	}
	
	// update the BlackHole geometry 
	BlackHoleGeometry::max_radius = svxrad[nlayers-1];
	BlackHoleGeometry::min_z = -length[nlayers-1]/2.;
	BlackHoleGeometry::max_z = length[nlayers-1]/2.;
	BlackHoleGeometry::gap = gap;

	return svxrad[nlayers-1];
}

#endif
