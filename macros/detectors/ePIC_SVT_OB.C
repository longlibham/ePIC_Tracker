/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_SVT_OB.C
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.06.19
*   	Description：	
*
=========================================================================*/
#ifndef FUN4ALL_EPICSVTOB_C
#define FUN4ALL_EPICSVTOB_C

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
R__LOAD_LIBRARY(libePIC_SVT_OB_Detector.so)

using namespace std;

namespace Enable{
	bool ePIC_SVTOB = false;
	bool ePIC_SVTOB_OVERLAPCHECK = false;
}  // name space Enable

void ePIC_SVTOB_Init(){

}



void SVTOBFastKalmanFilterConfigSVTX(PHG4TrackFastSim* kalman_filter, int ilay, double radius, bool addproj){

	kalman_filter->add_phg4hits(
			Form("G4HIT_SVTXOB_%d", ilay),
			PHG4TrackFastSim::Cylinder,
			50./10000./sqrt(12.), //radial-resolution [cm]
			20./10000./sqrt(12.), // azimuthal-resolution [cm]
			20./10000./sqrt(12.), // z-resolution [cm]
			1,   //efficiency
			0    // noise
			);
	kalman_filter->add_cylinder_state(Form("SVTXOB_%d", ilay), radius);
	if(addproj)TRACKING::ProjectionNames.insert(Form("SVTXOB_%d", ilay));
}


double ePIC_SVT_OB(PHG4Reco* g4Reco, const int nlayers = 2, double radius = 0){

	if (nlayers > 2){
		cout<<"Silicon VTX OB layers should not exceed 2!"<<endl;
		exit(1);
	}
	
	bool OverlapCheck = Enable::OVERLAPCHECK || Enable::ePIC_SVTOB_OVERLAPCHECK;

    //
	// build my SVT OB layers
	//
	// L3/L4 OB
	double r_inner[2] = {27.1, 41.8};
	double r_outer[2] = {27.7, 42.4};

	if(r_inner[0] <= radius){
		cout<<"Geometry overlap happens, please check the radius of each layer!"<<endl;
		exit(-1);
	}

	double carbon_thickness = 0.5*0.03;//0.5*0.03;
   	double carbon_length[2] = {54.31, 83.75};
	double carbon_width = 3.92;
	double si_thickness = 0.05/100.*9.37; //0.005;
	double si_length[2] = {15.0529, 12.8880};
	double si_width = 3.3128;
	double si_carbon_gap = 0.5;
	int n_silicon_z[2] = {4, 8};
	int n_stave_phi[2] = {46, 70};
	double las_overlap[2] = {2.55, 2.98};
	double gap = 0.1;
	
	ePIC_SVT_OB_Subsystem* svt_ob;
	for(int ilayer = 0; ilayer<nlayers; ilayer++){	
		
		svt_ob = new ePIC_SVT_OB_Subsystem("SVTXOB", ilayer);
		
		// set the parameters
		svt_ob->set_double_param("r_inner", r_inner[ilayer]);
		svt_ob->set_double_param("r_outer", r_outer[ilayer]);
		// carbon support 
		svt_ob->set_double_param("carbon_thickness", carbon_thickness);
		svt_ob->set_double_param("carbon_length", carbon_length[ilayer]);
		svt_ob->set_double_param("carbon_width", carbon_width);
		// LAS geo 
		svt_ob->set_double_param("si_thickness", si_thickness);
		svt_ob->set_double_param("si_length", si_length[ilayer]);
		svt_ob->set_double_param("si_width", si_width);
		svt_ob->set_double_param("si_carbon_gap", si_carbon_gap);
		svt_ob->set_int_param("n_silicon_z", n_silicon_z[ilayer]);
		svt_ob->set_int_param("n_stave_phi", n_stave_phi[ilayer]);

		// overlaps of LAS
		
		svt_ob->set_double_param("las_overlap", las_overlap[ilayer]);
		
		svt_ob->SetActive();
//		svt_ob->SuperDetector("SVTOB");

		svt_ob->OverlapCheck(OverlapCheck);
		g4Reco->registerSubsystem(svt_ob);

		SVTOBFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilter, ilayer, (r_inner[ilayer]+r_outer[ilayer])/2., false);//true);
		SVTOBFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilterInnerTrack, ilayer, (r_inner[ilayer]+r_outer[ilayer])/2., false);
		SVTOBFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilterSiliconTrack, ilayer, (r_inner[ilayer]+r_outer[ilayer])/2., false);
	
	
	}

	// update the BlackHole geometry 
	BlackHoleGeometry::max_radius = r_outer[nlayers-1];
	BlackHoleGeometry::min_z = -carbon_length[nlayers-1]/2.;
	BlackHoleGeometry::max_z = carbon_length[nlayers-1]/2.;
	BlackHoleGeometry::gap = gap;

	return r_outer[nlayers-1];
}

#endif