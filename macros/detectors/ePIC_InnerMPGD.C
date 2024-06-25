/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_MPGD.C
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.06.21
*   	Description：	
*
=========================================================================*/
#ifndef FUN4ALL_EPICINNERMPGD_C
#define FUN4ALL_EPICINNERMPGD_C

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
R__LOAD_LIBRARY(libePIC_InnerMPGD_Detector.so)

using namespace std;


namespace Enable{
	bool ePIC_InnerMPGD = false;
	bool ePIC_InnerMPGD_OVERLAPCHECK = false;
}  // name space Enable

void ePIC_InnerMPGD_Init(){

}



void InnerMPGDFastKalmanFilterConfig(PHG4TrackFastSim* kalman_filter, int ilay, double radius, bool addproj){

	kalman_filter->add_phg4hits(
			Form("G4HIT_InnerMPGD_%d", ilay),
			PHG4TrackFastSim::Cylinder,
			0.3, //radial-resolution [cm]
			0.015, // azimuthal-resolution [cm]
			0.015, // z-resolution [cm]
			1,   //efficiency
			0    // noise
			);
	kalman_filter->add_cylinder_state(Form("InnerMPGD_%d", ilay), radius);
	if(addproj)TRACKING::ProjectionNames.insert(Form("InnerMPGD_%d", ilay));
}


double ePIC_InnerMPGD(PHG4Reco* g4Reco, const int nlayers = 1, double radius = 0){

	if (nlayers > 1){
		cout<<"InnerMPGD layer should not exceed 1 for the moment!"<<endl;
		exit(1);
	}
	
	bool OverlapCheck = Enable::OVERLAPCHECK || Enable::ePIC_InnerMPGD_OVERLAPCHECK;

    //
	// build my MPGD layers
	//
	double rmpgd = 55.;
    double coat_t = 5.;
    double z_min = -105.;
    double z_max = 143.;
    double KaptonOverlay_t = 50./10000.;
    double CuGround_t = 1.58/10000.;
    double pcb_t = 100./10000.;
    double CuStrip_t = 12./10000.;
    double KaptonStrip_t = 75./10000.;
    double ResistiveStrip_t = 128./10000;
    double gas_t = 20./10000.;
    double mesh_t = 18./10000.;
    double GasGap_t = 3000./10000.;
    double DriftCuElectrode_t = 5./10000.;
    double DriftKapton_t = 250./10000.;
    double DriftCuGround_t = 0.41/10000.;
    double fudge_t = 570./10000.;

    int stave_number = 128;


	if(rmpgd <= radius){
		cout<<"Geometry overlap happens, please check the radius of each layer!"<<endl;
		exit(-1);
	}

	
	
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
        innerMPGD->set_double_param("ResitiveStrip_thickness", ResistiveStrip_t);
        innerMPGD->set_double_param("gas_thickness", gas_t);
        innerMPGD->set_double_param("mesh_thickness", mesh_t);
        innerMPGD->set_double_param("GasGap_thickness", GasGap_t);
        innerMPGD->set_double_param("DriftCuElectrode_thickness", DriftCuElectrode_t);
        innerMPGD->set_double_param("DriftKapton_thickness", DriftKapton_t);
        innerMPGD->set_double_param("DriftCuGround_thickness", DriftCuGround_t);
        innerMPGD->set_double_param("fudge_thickness", fudge_t);
		
		innerMPGD->SetActive();
//		innerMPGD->SuperDetector("SVTOB");

		innerMPGD->OverlapCheck(OverlapCheck);
		g4Reco->registerSubsystem(innerMPGD);

        double r_sensitive = rmpgd + KaptonOverlay_t + CuGround_t + pcb_t + CuStrip_t + KaptonStrip_t + ResistiveStrip_t + gas_t + mesh_t;
		InnerMPGDFastKalmanFilterConfig(TRACKING::FastKalmanFilter, ilayer, r_sensitive, false);//true);
		
	
	
	}

	// update the BlackHole geometry 
	BlackHoleGeometry::max_radius = rmpgd + coat_t/2.;
	BlackHoleGeometry::min_z = min_z;
	BlackHoleGeometry::max_z = max_z;
	BlackHoleGeometry::gap = no_overlapp;

	return rmpgd + coat_t/2.;
}

#endif
