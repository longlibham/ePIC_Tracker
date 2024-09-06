/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_SVTIB_Detector.cc
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.07.18
*   	Description：	
*
=========================================================================*/

#include <iostream>
#include "ePIC_SVTIB_Detector.h"
#include <fstream>
#include <sstream>

#include <g4main/PHG4Detector.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4Color.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4Transform3D.hh>
#include <Geant4/G4VisAttributes.hh>
#include <phparameter/PHParameters.h>
#include <Geant4/G4NistManager.hh>

#include <cmath>
#include <iostream>

#include <TSystem.h>

class G4VSolid;
class PHCompositeNode;

using namespace std;

ePIC_SVTIB_Detector::ePIC_SVTIB_Detector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string &dname, const int lyr)
	:PHG4Detector(subsys, Node, dname)
	, m_Params(parameters)
	, m_Layer(lyr){

}

int ePIC_SVTIB_Detector::IsInDetector(G4VPhysicalVolume* volume) const{
	set<G4VPhysicalVolume*>::const_iterator iter = m_PhysicalVolumesSet.find(volume);
	if(iter != m_PhysicalVolumesSet.end()){
		return 1;
	}
	return 0;
}

void ePIC_SVTIB_Detector::ConstructMe(G4LogicalVolume* logicWorld){
	// get essential parameter passing from PHParameters
    double radius = m_Params->get_double_param("radius") * cm;
    double si_thickness = m_Params->get_double_param("si_thickness") * cm;
    double length = m_Params->get_double_param("length") * cm;
    double lec_length = m_Params->get_double_param("lec_length") * cm;
    double rec_length = m_Params->get_double_param("rec_length") * cm;
    double peri_width = m_Params->get_double_param("periphery_width") * cm;
    double tile_width = m_Params->get_double_param("tile_width") * cm;

    int nphi = 2*M_PI*radius/tile_width;
   

    // construct the original cylinder
    if(!isfinite(radius) || !isfinite(length) || !isfinite(lec_length)
    || !isfinite(rec_length) || !isfinite(peri_width) || !isfinite(tile_width)){
        cout << "PHWHERE: " << "Bad parameter for " << GetName() << endl;
        cout << "radius: " << radius << " length: "<< length <<endl;
        cout << "lec_length: " << lec_length << " rec_length: " << rec_length << endl;
        cout << "peri_width: " << peri_width << " tile_width: " << tile_width << endl;
        gSystem->Exit(-1);
    }

    G4VSolid* si_cyl = new G4Tubs("si_cylinder", radius, radius + si_thickness, length/2., 0., 2*M_PI);
    G4VSolid* si_sub = NULL;
    if (peri_width > 0) si_sub = new G4Box("si_substract", peri_width/2., (si_thickness+m_redundancy)/2., (length + m_redundancy)/2.);

    // subtract the insensitive periphery area.
    G4VSolid* si_sensitive_solid = nullptr;
    ostringstream oss;

    for(int i = 0; i < nphi; i++){
        
        if(peri_width <= 0) continue;  // for no dead area.

        double phi = i * 2*M_PI/nphi;
        G4RotationMatrix rotm = G4RotationMatrix();
        rotm.rotateZ(M_PI/2. + phi);
        double radius_sub = radius + si_thickness/2.;
        G4ThreeVector pos = G4ThreeVector(
            radius_sub*cos(phi),
            radius_sub*sin(phi),
            0
            );

        G4Transform3D transform = G4Transform3D(rotm, pos); 
        oss.str("");
        oss<<"si_sensitive_"<<i;
        if(i == 0) si_sensitive_solid = new G4SubtractionSolid(oss.str().c_str(), si_cyl, si_sub, transform);
        else si_sensitive_solid = new G4SubtractionSolid(oss.str().c_str(), si_sensitive_solid, si_sub, transform);
    }

    // build the LogicalVolume for sensitive area.
    G4Colour col_si_sensitive = G4Colour(251./255., 248./255., 0., 1.);
    G4Colour col_si_insensitive = G4Colour(185./255., 99./255., 83./255., 1.);

    G4NistManager* nist = G4NistManager::Instance();
    G4Material* mat_si = nist->FindOrBuildMaterial("G4_Si");

    G4LogicalVolume* si_sen_logic = NULL;
    if(si_sensitive_solid)
        si_sen_logic = new G4LogicalVolume(si_sensitive_solid, mat_si, "SiSensorLogic");
    else
        si_sen_logic = new G4LogicalVolume(si_cyl, mat_si, "SiSensorLogic");

    G4VisAttributes* sensor_vis = new G4VisAttributes(col_si_sensitive);
    sensor_vis->SetForceSolid(true);
    si_sen_logic->SetVisAttributes(sensor_vis);

    G4VPhysicalVolume* sensor_phy =
    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0., 0., 0),
        si_sen_logic,
        "sensor_layer",
        logicWorld,
        0,
        0,
        OverlapCheck()
    );

    m_PhysicalVolumesSet.insert(sensor_phy);

    if(peri_width > 0){
         // for dead area
        G4VSolid* insen_solid = new G4Box("Insensitive_si", peri_width/2., si_thickness/2., length/2.);
        G4LogicalVolume* insen_logic = new G4LogicalVolume(insen_solid, mat_si, "SiSubLogic");
        G4VisAttributes* insen_vis = new G4VisAttributes(col_si_insensitive);
        insen_vis->SetForceSolid(true);
        insen_logic->SetVisAttributes(insen_vis);
        
        for(int i = 0; i < nphi; i++){
            double phi = i * 2*M_PI/nphi;
            G4RotationMatrix rotm = G4RotationMatrix();
            rotm.rotateZ(M_PI/2. + phi);
            G4ThreeVector pos = G4ThreeVector(
                radius*cos(phi),
                radius*sin(phi),
                0
                );

            G4Transform3D transform = G4Transform3D(rotm, pos); 
            new G4PVPlacement(
                transform,
                insen_logic,
                "insen_si_layer",
                logicWorld,
                0,
                i,
                OverlapCheck()
            );
        }
    
        // Left/right endcap
        G4VSolid* lec_solid = new G4Tubs("lec_cylinder", radius, radius+si_thickness, lec_length/2., 0, 2*M_PI);
        G4LogicalVolume* lec_logic = new G4LogicalVolume(lec_solid, mat_si, "LECLogic");
        lec_logic->SetVisAttributes(insen_vis);

        new G4PVPlacement(
            nullptr,
            G4ThreeVector(0, 0, -(length + lec_length)/2. ),
            lec_logic,
            "LeftEndcap_cyl",
            logicWorld,
            0,
            0,
            OverlapCheck()
        );

        G4VSolid* rec_solid = new G4Tubs("rec_cylinder", radius, radius+si_thickness, rec_length/2., 0, 2*M_PI);
        G4LogicalVolume* rec_logic = new G4LogicalVolume(rec_solid, mat_si, "RECLogic");
        rec_logic->SetVisAttributes(insen_vis);

        new G4PVPlacement(
            nullptr,
            G4ThreeVector(0, 0, (length + rec_length)/2. + 2*m_nooverlap),
            rec_logic,
            "RightEndcap_cyl",
            logicWorld,
            0,
            0,
            OverlapCheck()
        );
    
    }
    
    return;


}

void ePIC_SVTIB_Detector::Print(const std::string &what) const{
	std::cout<<"ePIC SVT IB Detector: "<<std::endl;
	if(what == "ALL" || what == "VOLUME"){
		std::cout<<"Version 0.1"<<std::endl;
	}
	return;
}
