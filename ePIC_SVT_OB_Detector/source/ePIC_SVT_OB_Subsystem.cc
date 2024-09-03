/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_SVT_OB_Subsystem.cc
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.04.18
*   	Description：	
*
=========================================================================*/

#include <iostream>
#include "ePIC_SVT_OB_Subsystem.h"
#include "ePIC_SVT_OB_Detector.h"
#include "ePIC_SVT_OB_SteppingAction.h"
#include <fstream>
#include <sstream>

#include <phparameter/PHParameters.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>
#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <cmath>

using namespace std;

ePIC_SVT_OB_Subsystem::ePIC_SVT_OB_Subsystem(const std::string& name, const int lyr)
	: PHG4DetectorSubsystem(name, lyr)
	  , m_Detector(nullptr)
	  , m_SteppingAction(nullptr){
	  
		  InitializeParameters();
}

int ePIC_SVT_OB_Subsystem::InitRunSubsystem(PHCompositeNode* topNode){
	
	// use world material if material was not set so far
	if(GetParams()->get_string_param("material") == "WorldMaterial"){
		recoConsts* rc = recoConsts::instance();
		GetParams()->set_string_param("material", rc->get_StringFlag("WorldMaterial"));
	}



	PHNodeIterator iter(topNode);
	PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode* >(iter.findFirst("PHCompositeNode", "DST"));
	PHNodeIterator dstIter(dstNode);
	PHCompositeNode* DetNode = dynamic_cast<PHCompositeNode* > (dstIter.findFirst("PHCompositeNode", Name()));
	
	if(!DetNode){
		DetNode = new PHCompositeNode(Name());
		dstNode->addNode(DetNode);
	}

	string g4hitnodename = "G4HIT_" + Name();
	PHG4HitContainer* g4_hits = findNode::getClass<PHG4HitContainer>(DetNode, g4hitnodename);

	if(!g4_hits){
		g4_hits = new PHG4HitContainer(g4hitnodename);
		DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, g4hitnodename, "PHObject"));
	}

	//create detector
	m_Detector = new ePIC_SVT_OB_Detector(this, topNode, GetParams(), Name(), GetLayer());
	m_Detector -> OverlapCheck(CheckOverlap());

	// create stepping action
	m_SteppingAction = new ePIC_SVT_OB_SteppingAction(m_Detector);

	return 0;

}

int ePIC_SVT_OB_Subsystem::process_event(PHCompositeNode* topNode){
	// pass top node to stepping action so that it gets 
	// relevant nodes need internally
	if(m_SteppingAction){
		m_SteppingAction->SetInterfacePointers(topNode);
	}

	return 0;
}

void ePIC_SVT_OB_Subsystem::Print(const string& what) const{
	
	//
	if (m_Detector){
		m_Detector->Print(what);
	}
	return;
}

PHG4Detector* ePIC_SVT_OB_Subsystem::GetDetector(void) const {
	
	return m_Detector;
}

void ePIC_SVT_OB_Subsystem::SetDefaultParameters(){
	
	set_default_double_param("r_inner", NAN);
	set_default_double_param("r_outer", NAN);
	set_default_double_param("place_x", 0.);
	set_default_double_param("place_y", 0.);
	set_default_double_param("place_z", 0.);
	set_default_double_param("carbon_thickness", 0.);
	set_default_double_param("carbon_length", 0.);
	set_default_double_param("carbon_width", 0.);
	set_default_double_param("cf_leftendcap", 0.);
	set_default_double_param("cf_rightendcap", 0.);
	set_default_double_param("cf_margin", 0.);
	set_default_double_param("cf_curve_radius", 0.);
	set_default_double_param("cf_center_height", 0.);
	set_default_double_param("cf_edge_height", 0.);

	set_default_double_param("oring_radius", 0.);
	set_default_double_param("oring_spacing", 0.);

	set_default_double_param("cf_csupport_width", 0.);

	set_default_double_param("k9_center_height", 0.);

	set_default_double_param("si_thickness", 0.);
	set_default_double_param("si_length", 0.);
	set_default_double_param("si_width", 0.);
	set_default_int_param("n_silicon_z", 0);
	set_default_int_param("n_stave_phi", 0);
	set_default_double_param("lec_length", 0.);
	set_default_double_param("rec_length", 0.);
	set_default_double_param("anc_length", 0.);
	set_default_double_param("anc_thickness", 0.);
	set_default_double_param("las_airspace", 0.);
	set_default_double_param("periphery_width", 0.);
	set_default_double_param("kapton_thickness", 0.);


	// place holder, will be replaced by world material if not set by other means(macro)
	set_default_string_param("material", "WorldMaterial");
}
