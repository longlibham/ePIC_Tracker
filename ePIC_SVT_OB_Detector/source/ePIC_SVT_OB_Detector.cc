/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_SVT_OB_Detector.cc
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.04.08
*   	Description：	
*
=========================================================================*/

#include <iostream>
#include "ePIC_SVT_OB_Detector.h"
#include <fstream>
#include <sstream>

#include <g4main/PHG4Detector.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Color.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4Transform3D.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <iostream>

class G4VSolid;
class CompositeNode;

using namespace std;

ePIC_SVT_OB_Detector::ePIC_SVT_OB_Detector(PHG4Subsystem* subsys, PHCompositeNode* Node, const std::string &dname):PHG4Detector(subsys, Node, dname){

}

int ePIC_SVT_OB_Detector::IsInDetector(G4PhysicalVolume* volume) const{
	set<PHG4PhysicalVolume*>::const_iterator iter = m_PhysicalVolumesSet.find(volume);
	if(iter != m_PhysicalVolumesSet.end()){
		return 1;
	}
	return 0;
}

void ePIC_SVT_OB_Detector::ConstructMe(G4LogicalVolume* logicWorld){
	// L3 barrel 
	double l3_length = 43.3333 * cm;
	double l3_width = 3.9127 * cm;
	double carbon_t = 0.05 * cm;
	double si_t = 0.005 * cm;
	double si_c_gap = 0.3 * cm;
	double l3_gap = 0.5 * cm;

	int las_per_strip = 2;
	int l3_nstrips = 36;

	double l3_r = 22.2 *cm;
	dphi = twopi/l3_nstrips;

	// silicon layer
	G4VSolid* solid_si_l3 = new G4Box("L3_Si_box", l3_length/2., l3_width/2., si_t/2.);
	G4LogicalVolume* si_logical_l3 = new G4LogicalVolume(solid_si_l3, G4Material::GetMaterial("G4_Si"), "L3_Si_logical");
	G4VisAttributes* sivis = new G4VisAttributes(G4Color(G4Colour::Yellow()));
	sivis->SetForceSolid(true);
	si_logical_l3->SetVisAttributes(sivis);

	// Carbon layer
	G4VSolid* solid_carbon_l3 = new G4Box("L3_carbon_box", l3_length/2., l3_width/2., carbon_t/2.);
	G4LogicalVolume* carbon_logical_l3 = new G4LogicalVolume(solid_carbon_l3, G4Material::GetMaterial("G4_C"), "L3_carbon_logical");
	G4VisAttributes* cvis = new G4VisAttributes(G4Color(G4Colour::Black()));
	cvis->SetForceSolid(true);
	carbon_logical_l3->SetVisAttributes(cvis);
	

	// placement     ***** most important *****


}

