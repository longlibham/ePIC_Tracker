/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_tracking.C
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.05.28
*   	Description：	
*
=========================================================================*/
#ifndef MACRO_EPICTRACKING_C
#define MACRO_EPICTRACKING_C

#include <GlobalVariables.C>
#include <trackreco/PHRaveVertexing.h>
#include <g4trackfastsim/PHG4TrackFastSimEval.h>
#include <fun4all/Fun4AllServer.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

R__LOAD_LIBRARY(libtrack_reco.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)

using namespace std;

namespace Enable{
	bool TRACKING = false;
	bool TRACKING_EVAL = false;
	int TRACKING_VERBOSITY = 0;
} // namespace Enable

namespace TRACKING{
	//string TrackNodeName = "TrackMap";
	string TrackTopName = "SVTX";
	PHG4TrackFastSim* FastKalmanFilter(nullptr);
	PHG4TrackFastSim* FastKalmanFilterSiliconTrack(nullptr);
	PHG4TrackFastSim* FastKalmanFilterInnerTrack(nullptr);
	std::set<std::string> ProjectionNames;
}

namespace G4TRACKING
{
  bool DISPLACED_VERTEX = false;
  bool PROJECTION_EEMC = false;
  bool PROJECTION_EHCAL = false;
  bool PROJECTION_CEMC = false;
  bool PROJECTION_BECAL = false;
  bool PROJECTION_HCALIN = false;
  bool PROJECTION_HCALOUT = false;
  bool PROJECTION_FEMC = false;
  bool PROJECTION_FHCAL = false;
  bool PROJECTION_LFHCAL = false;
}  // namespace G4TRACKING



void TrackingInit(){
	TRACKING::FastKalmanFilter = new PHG4TrackFastSim("PHG4TrackFastSim");
	TRACKING::FastKalmanFilterSiliconTrack = new PHG4TrackFastSim("FastKalmanFilterSiliconTrack");
	TRACKING::FastKalmanFilterInnerTrack = new PHG4TrackFastSim("FastKalmanFilterInnerTrack");
}

void InitFastKalmanFilter(PHG4TrackFastSim* kalman_filter){
	// kalman_filter->Smearing(false);
	if(G4TRACKING::DISPLACED_VERTEX){
		// do not use truth vertex in the track fitting,
		// which would lead to worse momentum resolution for prompt tracks
		// but this allows displaced track analysis including CCA and vertex finding
		kalman_filter->set_use_vertex_in_fitting(false);
		kalman_filter->set_vertex_xy_resolution(0);
		kalman_filter->set_vertex_z_resolution(0);
		kalman_filter->enable_vertexing(true);
	}
	else{
		kalman_filter->set_use_vertex_in_fitting(true);
		kalman_filter->set_vertex_xy_resolution(5e-4);
		kalman_filter->set_vertex_z_resolution(5e-4);
	}
	kalman_filter->set_sub_top_node_name("TRACKS");
}

void Tracking_Reco(){
	int verbosity = std::max(Enable::VERBOSITY, Enable::TRACKING_VERBOSITY);

	// fun4all server
	//
	Fun4AllServer* se = Fun4AllServer::instance();

	if(TRACKING::FastKalmanFilter == nullptr){
		cout<<__PRETTY_FUNCTION__<<": missing the expected initialization for TRACKING::FastKalmanFilter."<<endl;
		exit(1);
	}

	InitFastKalmanFilter(TRACKING::FastKalmanFilter);
	TRACKING::FastKalmanFilter->Verbosity(verbosity);
	TRACKING::FastKalmanFilter->set_sub_top_node_name(TRACKING::TrackTopName);
	TRACKING::FastKalmanFilter->set_trackmap_out_name(TRACKING::TrackNodeName);

	se->registerSubsystem(TRACKING::FastKalmanFilter);

	//next, tracks with partial usage of the tracker stack
	if(TRACKING::FastKalmanFilterInnerTrack == nullptr){
		cout<<__PRETTY_FUNCTION__<<": missing the expected initialization for TRACKING::FastKalmanFilterInnerTrack."<<endl;
		exit(1);
	}

	InitFastKalmanFilter(TRACKING::FastKalmanFilterInnerTrack);
	TRACKING::FastKalmanFilterInnerTrack->Verbosity(verbosity);
	TRACKING::FastKalmanFilterInnerTrack->set_trackmap_out_name("InnerTrackMap");
	TRACKING::FastKalmanFilterInnerTrack->enable_vertexing(false);
	se->registerSubsystem(TRACKING::FastKalmanFilterInnerTrack);

	if(TRACKING::FastKalmanFilterSiliconTrack == nullptr){
		cout<<__PRETTY_FUNCTION__<<": missing the expected initialization for TRACKING::FastKalmanFilterSiliconTrack."<<endl;
		exit(1);
	}

	InitFastKalmanFilter(TRACKING::FastKalmanFilterSiliconTrack);
	TRACKING::FastKalmanFilterSiliconTrack->Verbosity(verbosity);
	TRACKING::FastKalmanFilterSiliconTrack->set_trackmap_out_name("SiliconTrackMap");
	TRACKING::FastKalmanFilterSiliconTrack->enable_vertexing(false);
	se->registerSubsystem(TRACKING::FastKalmanFilterSiliconTrack);

	return;

}

void Tracking_Eval(const std::string &outputfile){
	int verbosity = std::max(Enable::VERBOSITY, Enable::TRACKING_VERBOSITY);
	//
	// Fun4All server
	//
	Fun4AllServer* se = Fun4AllServer::instance();

	//
	// Fast Tracking evaluation
	//
	
	PHG4TrackFastSimEval* fast_sim_eval = new PHG4TrackFastSimEval("FastTrackingEval");
	fast_sim_eval->set_trackmapname(TRACKING::TrackNodeName);
	fast_sim_eval->set_filename(outputfile);
	fast_sim_eval->Verbosity(verbosity);


	cout<<"Tracking_Eval(): configuration of track projections in PHG4TrackFastSimEval"<<endl;
	cout<<"/*std::set<std::string>*/ TRACKING::ProjectionNames = {";
	bool first = true;
	for(const string& proj: TRACKING::ProjectionNames){
		if(first)
			first = false;
		else
			cout<<", ";
		cout<<"\""<<proj<<"\"";

		fast_sim_eval->AddProjection(proj);
	}
	cout<<"};"<<endl;

	se->registerSubsystem(fast_sim_eval);

	// now partial track fits
	fast_sim_eval = new PHG4TrackFastSimEval("FastTrackingEval_InnerTrackMap");
	fast_sim_eval->set_trackmapname("InnerTrackMap");
	fast_sim_eval->set_filename(outputfile + ".InnerTrackMap.root");
	fast_sim_eval->Verbosity(verbosity);
	se->registerSubsystem(fast_sim_eval);

	fast_sim_eval = new PHG4TrackFastSimEval("FastTrackingEval_SiliconTrackMap");
	fast_sim_eval->set_trackmapname("SiliconTrackMap");
	fast_sim_eval->set_filename(outputfile + ".SiliconTrackMap.root");
	fast_sim_eval->Verbosity(verbosity);
	se->registerSubsystem(fast_sim_eval);


}


#endif




