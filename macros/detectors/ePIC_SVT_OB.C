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

namespace ePIC_SVTOB{
	double r_inner[2] = {27.0, 41.7};
	double r_outer[2] = {28.2, 42.9};

	double carbon_thickness = 0.2/100*26.356;
	double carbon_length[2] = {53.2, 89.04};
	double carbon_width = 3.92;
	double cf_leftendcap = 1.2;
	double cf_rightendcap = 1.2;
	double cf_margin = 0.1;
	double cf_curve_radius = 9.0115;
	double cf_center_height = 0.85;
	double cf_edge_height = 0.351;
	double oring_radius = 0.25;
	double oring_spacing = 0.3;
	double cf_csupport_width = 0.5;
	double k9_center_height = 0.5;
	double si_thickness = 0.05/100*9.37;
	double si_length[2] = {13.0, 10.83};
	double si_width = 3.9128;
	double n_silicon_z[2] = {4, 8};
	double n_stave_phi[2] = {46, 70};
	double lec_length = 0.45;
	double rec_length = 0.15;
	double anc_length = 1.0;
	double anc_thickness = 300./10000.;
	double las_airspace = 0.6;
	double peri_width = 0.0525;
	double kapton_thickness = 100./10000.;


}

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

	if(ePIC_SVTOB::r_inner[0] <= radius){
		cout<<"Geometry overlap happens, please check the radius of each layer!"<<endl;
		exit(-1);
	}


	
	ePIC_SVT_OB_Subsystem* svt_ob;
	for(int ilayer = 0; ilayer<nlayers; ilayer++){	
		
		svt_ob = new ePIC_SVT_OB_Subsystem("SVTXOB", ilayer);
		
		// set the parameter
		svt_ob->set_double_param("r_inner", ePIC_SVTOB::r_inner[ilayer]);
		svt_ob->set_double_param("r_outer", ePIC_SVTOB::r_outer[ilayer]);
		svt_ob->set_double_param("carbon_thickness", ePIC_SVTOB::carbon_thickness);
		svt_ob->set_double_param("carbon_length", ePIC_SVTOB::carbon_length[ilayer]);
		svt_ob->set_double_param("carbon_width", ePIC_SVTOB::carbon_width);
		svt_ob->set_double_param("cf_leftendcap", ePIC_SVTOB::cf_leftendcap);
		svt_ob->set_double_param("cf_rightendcap", ePIC_SVTOB::cf_rightendcap);
		svt_ob->set_double_param("cf_margin", ePIC_SVTOB::cf_margin);
		svt_ob->set_double_param("cf_curve_radius", ePIC_SVTOB::cf_curve_radius);
		svt_ob->set_double_param("cf_center_height", ePIC_SVTOB::cf_center_height);
		svt_ob->set_double_param("cf_edge_height", ePIC_SVTOB::cf_edge_height);
		svt_ob->set_double_param("oring_radius", ePIC_SVTOB::oring_radius);
		svt_ob->set_double_param("oring_spacing", ePIC_SVTOB::oring_spacing);
		svt_ob->set_double_param("cf_csupport_width", ePIC_SVTOB::cf_csupport_width);
		svt_ob->set_double_param("k9_center_height", ePIC_SVTOB::k9_center_height);
		svt_ob->set_double_param("si_thickness", ePIC_SVTOB::si_thickness);
		svt_ob->set_double_param("si_length", ePIC_SVTOB::si_length[ilayer]);
		svt_ob->set_double_param("si_width", ePIC_SVTOB::si_width);
		svt_ob->set_int_param("n_silicon_z", ePIC_SVTOB::n_silicon_z[ilayer]);
		svt_ob->set_int_param("n_stave_phi", ePIC_SVTOB::n_stave_phi[ilayer]);
		svt_ob->set_double_param("lec_length", ePIC_SVTOB::lec_length);
		svt_ob->set_double_param("rec_length", ePIC_SVTOB::rec_length);
		svt_ob->set_double_param("anc_length", ePIC_SVTOB::anc_length);
		svt_ob->set_double_param("anc_thickness", ePIC_SVTOB::anc_thickness);
		svt_ob->set_double_param("las_airspace", ePIC_SVTOB::las_airspace);
		svt_ob->set_double_param("periphery_width", ePIC_SVTOB::peri_width);
		svt_ob->set_double_param("kapton_thickness", ePIC_SVTOB::kapton_thickness);

		
		svt_ob->SetActive();
//		svt_ob->SuperDetector("SVTOB");

		svt_ob->OverlapCheck(OverlapCheck);
		g4Reco->registerSubsystem(svt_ob);

		SVTOBFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilter, ilayer, (ePIC_SVTOB::r_inner[ilayer]+ePIC_SVTOB::r_outer[ilayer])/2., false);//true);
		SVTOBFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilterInnerTrack, ilayer, (ePIC_SVTOB::r_inner[ilayer]+ ePIC_SVTOB::r_outer[ilayer])/2., false);
		SVTOBFastKalmanFilterConfigSVTX(TRACKING::FastKalmanFilterSiliconTrack, ilayer, (ePIC_SVTOB::r_inner[ilayer]+ ePIC_SVTOB::r_outer[ilayer])/2., false);
	
	
	}

	// update the BlackHole geometry 
	BlackHoleGeometry::max_radius = ePIC_SVTOB::r_outer[nlayers-1];
	BlackHoleGeometry::min_z = -ePIC_SVTOB::carbon_length[nlayers-1]/2.;
	BlackHoleGeometry::max_z = ePIC_SVTOB::carbon_length[nlayers-1]/2.;
	BlackHoleGeometry::gap = no_overlapp;

	return ePIC_SVTOB::r_outer[nlayers-1];
}

#endif