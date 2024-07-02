/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_OuterMPGD_Subsystem.cc
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.07.02
*   	Description：	
*
=========================================================================*/

#include <iostream>
#include "ePIC_OuterMPGD_Subsystem.h"
#include "ePIC_OuterMPGD_Detector.h"
#include "ePIC_OuterMPGD_SteppingAction.h"
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


ePIC_OuterMPGD_Subsystem::ePIC_OuterMPGD_Subsystem(const std::string& name, const int lyr)
	: PHG4DetectorSubsystem(name, lyr)
	  , m_Detector(nullptr)
	  , m_SteppingAction(nullptr){
	  
		  InitializeParameters();
}

int ePIC_OuterMPGD_Subsystem::InitRunSubsystem(PHCompositeNode* topNode){
	
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
	m_Detector = new ePIC_OuterMPGD_Detector(this, topNode, GetParams(), Name(), GetLayer());
	m_Detector -> OverlapCheck(CheckOverlap());

	// create stepping action
	m_SteppingAction = new ePIC_OuterMPGD_SteppingAction(m_Detector);

	return 0;

}

int ePIC_OuterMPGD_Subsystem::process_event(PHCompositeNode* topNode){
	// pass top node to stepping action so that it gets 
	// relevant nodes need internally
	if(m_SteppingAction){
		m_SteppingAction->SetInterfacePointers(topNode);
	}

	return 0;
}

void ePIC_OuterMPGD_Subsystem::Print(const string& what) const{
	
	//
	if (m_Detector){
		m_Detector->Print(what);
	}
	return;
}

PHG4Detector* ePIC_OuterMPGD_Subsystem::GetDetector(void) const {
	
	return m_Detector;
}

void ePIC_OuterMPGD_Subsystem::SetDefaultParameters(){
	
	set_default_double_param("radius", NAN);
    set_default_double_param("module_width", NAN);
    set_default_double_param("zmin", NAN);
    set_default_double_param("zmax", NAN);
    set_default_double_param("window_thickness", 0);
    set_default_double_param("windowgap_thickness", 0);
    set_default_double_param("driftgap_thickness", 0);
    set_default_double_param("foilcu_thickness", 0);
    set_default_double_param("readoutelectrode_thickness", 0);
    set_default_double_param("foilkapton_thickness", 0);
    set_default_double_param("readoutnomex_thickness", 0);
    set_default_double_param("readoutkapton_thickness", 0);
    set_default_double_param("pcb_thickness", 0);
    set_default_int_param("stave_phi", 0);
    set_default_int_param("stave_z", 0);


	// place holder, will be replaced by world material if not set by other means(macro)
	set_default_string_param("material", "WorldMaterial");
}
