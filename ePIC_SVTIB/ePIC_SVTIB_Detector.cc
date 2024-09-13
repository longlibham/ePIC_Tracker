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

    double matrix_length = m_Params->get_double_param("matrix_length") * cm;
    double switch_length = m_Params->get_double_param("switch_length") * cm;
    double backbone_length = m_Params->get_double_param("backbone_length") * cm;
    double ntile = m_Params->get_int_param("ntile");
    
    double lec_length = m_Params->get_double_param("lec_length") * cm;
    double rec_length = m_Params->get_double_param("rec_length") * cm;
    
    double peri_width = m_Params->get_double_param("periphery_width") * cm;
    double tile_width = m_Params->get_double_param("tile_width") * cm;

    
    int nphi = 2*M_PI*radius/tile_width;


    // construct the original cylinder
    if(!isfinite(radius) ||!isfinite(si_thickness) ||!isfinite(matrix_length) 
    || !isfinite(switch_length) || !isfinite(backbone_length) || !isfinite(ntile) 
    ||!isfinite(lec_length) || !isfinite(rec_length) || !isfinite(peri_width) 
    || !isfinite(tile_width)){
        cout << "PHWHERE: " << "Bad parameter for " << GetName() << endl;
        cout << "radius: " << radius << " matrix_length: "<< matrix_length <<endl;
        cout<<"switch_length: "<< switch_length<<" backbone_length: "<<backbone_length<<endl;
        cout<<"ntile: " << ntile <<endl;
        cout << "lec_length: " << lec_length << " rec_length: " << rec_length << endl;
        cout << "peri_width: " << peri_width << " tile_width: " << tile_width << endl;
        gSystem->Exit(-1);
    }

    
    // color and material setting
    G4Colour col_si_sensitive = G4Colour(251./255., 248./255., 0., 1.);
    G4Colour col_si_insensitive = G4Colour(185./255., 99./255., 83./255., 1.);

    G4NistManager* nist = G4NistManager::Instance();
    G4Material* mat_si = nist->FindOrBuildMaterial("G4_Si");
    G4Material* mat_air = nist->FindOrBuildMaterial("G4_AIR");

    // VisAttribute
    G4VisAttributes* sensor_vis = new G4VisAttributes(col_si_sensitive);
    sensor_vis->SetForceSolid(true);
    G4VisAttributes* insen_si_vis = new G4VisAttributes(col_si_insensitive);
    insen_si_vis->SetForceSolid(true);


    // BUILD air container for the tile
    double tileContainer_length = 3*(matrix_length + switch_length) + backbone_length;
    double ib_length = ntile*tileContainer_length;
    G4VSolid* tileContainer_solid = new G4Tubs("tileContainer", radius - m_nooverlap, radius + si_thickness + m_nooverlap, tileContainer_length/2., 0, 2*M_PI);
    G4LogicalVolume* tileContainerLogic = new G4LogicalVolume(tileContainer_solid, mat_air, "TileContainerLogic");

    // build the air container for the matrix
    G4VSolid* mtxContainer_solid = new G4Tubs("mtxContainer_solid", radius - m_nooverlap, radius + si_thickness + m_nooverlap, matrix_length/2., 0., 2*M_PI);
    G4LogicalVolume* mtxContainerLogic = new G4LogicalVolume(mtxContainer_solid, mat_air, "mtxContainerLogic");
    

    // place the sensitive/insensitive si to the mtxContainer

    // build the sensor matrix 
    G4VSolid* matrix_cyl = new G4Tubs("matrix_cylinder", radius, radius + si_thickness, matrix_length/2., 0., 2*M_PI);
    G4VSolid* si_sub = new G4Tubs("si_substract", radius - m_redundancy/2., radius + si_thickness + m_redundancy/2., matrix_length/2. + m_redundancy, 0., peri_width/radius);

    ostringstream oss;
    G4RotationMatrix rotm = G4RotationMatrix();
    G4Transform3D transform;
    G4ThreeVector pos = G4ThreeVector(0., 0., 0.);
    double dphi = 2*M_PI/nphi;

    for(int i = 0; i < nphi; i++){
        rotm.rotateZ(dphi);       
        transform = G4Transform3D(rotm, pos); 

        oss.str("");
        oss<<"si_sensitive_"<<i;
        matrix_cyl = new G4SubtractionSolid(oss.str().c_str(), matrix_cyl, si_sub, transform);
       

    }


    // build the LogicalVolume for sensitive area.
    G4LogicalVolume* si_sen_logic = new G4LogicalVolume(matrix_cyl, mat_si, "SiSensorLogic");
    si_sen_logic->SetVisAttributes(sensor_vis);

    // put the sensor in the mtxContainer
    G4VPhysicalVolume* sensor_phy = 
    new G4PVPlacement(
        nullptr,
        pos,
        si_sen_logic,
        "sensor_layer",
        mtxContainerLogic,
        0,
        0,
        OverlapCheck()
    );
    m_PhysicalVolumesSet.insert(sensor_phy);
    
    // build the insensitive periphery and fill in the gap
    G4VSolid* peri_solid = new G4Tubs("PeripherySolid", radius, radius + si_thickness, matrix_length/2., 0., peri_width/radius);
    G4LogicalVolume* periLogic = new G4LogicalVolume(peri_solid, mat_si, "PeripheryLogic");
    periLogic->SetVisAttributes(insen_si_vis);


    rotm = G4RotationMatrix();
    for(int i =0; i<nphi; i++){
        rotm.rotateZ(dphi);
        transform = G4Transform3D(rotm, pos);
        oss.str("");
        oss<<"periphery_"<<i;

        new G4PVPlacement(
            transform,
            periLogic,
            oss.str().c_str(),
            mtxContainerLogic,
            0,
            i,
            OverlapCheck()
        );
    }

    // build the backbone insensitive tub
    G4VSolid* backbone_solid = new G4Tubs("BackBone_tub", radius, radius + si_thickness, backbone_length/2., 0., 2*M_PI);
    G4LogicalVolume* backboneLogic = new G4LogicalVolume(backbone_solid, mat_si, "BackboneLogic");
    backboneLogic->SetVisAttributes(insen_si_vis);

    // build the power switch insensitive tub
    G4VSolid* switch_solid = new G4Tubs("PowerSwitchSolid", radius, radius + si_thickness, switch_length/2., 0., 2*M_PI);
    G4LogicalVolume* switchLogic = new G4LogicalVolume(switch_solid, mat_si, "SwitchLogic");
    switchLogic->SetVisAttributes(insen_si_vis);

    //place the sen/insen silicon into the tileContainer
    // 1. backbone
    double current_z = -tileContainer_length/2. + backbone_length/2.;
    pos = G4ThreeVector(0., 0., current_z);
    new G4PVPlacement(
        nullptr,
        pos,
        backboneLogic,
        "BackboneLayer",
        tileContainerLogic,
        0,
        0,
        OverlapCheck()
    );

    // 2. mtxContainer & power switch
    for(int i =0; i<3; i++){
        if(i == 0) current_z += (backbone_length + matrix_length)/2.;
        else current_z += (switch_length + matrix_length)/2.;
        pos = G4ThreeVector(0, 0, current_z);
        new G4PVPlacement(
            nullptr,
            pos,
            mtxContainerLogic,
            "mtxContainerLayer",
            tileContainerLogic,
            0,
            i,
            OverlapCheck()
        );

        current_z += (matrix_length + switch_length)/2.;
        pos = G4ThreeVector(0, 0, current_z);
        new G4PVPlacement(
            nullptr,
            pos,
            switchLogic,
            "PowerSwitchLogic",
            tileContainerLogic,
            0,
            i,
            OverlapCheck()
        );

    }

    // put the tileContainer in the World
    
    current_z = 0.;
    for(int i = 0; i<ntile; i++){
        pos = G4ThreeVector(0., 0., -ib_length/2. + (2*i +1)*tileContainer_length/2.);
        new G4PVPlacement(
            nullptr,
            pos,
            tileContainerLogic,
            "tileContainerLayer",
            logicWorld,
            0,
            i,
            OverlapCheck()
        );

    }

    // build the LEC/REC
    G4VSolid* lec_solid = new G4Tubs("LECSolid", radius, radius + si_thickness, lec_length/2., 0, 2*M_PI);
    G4LogicalVolume* lecLogic = new G4LogicalVolume(lec_solid, mat_si, "LECLogic");
    lecLogic->SetVisAttributes(insen_si_vis);

    current_z = -ib_length/2. - lec_length/2.;
    pos = G4ThreeVector(0., 0., current_z);
    new G4PVPlacement(
        nullptr,
        pos,
        lecLogic,
        "LECLayer",
        logicWorld,
        0,
        0,
        OverlapCheck()
    );

    G4VSolid* rec_solid = new G4Tubs("RECSolid", radius, radius + si_thickness, rec_length/2., 0, 2*M_PI);
    G4LogicalVolume* recLogic = new G4LogicalVolume(rec_solid, mat_si, "RECLogic");
    recLogic->SetVisAttributes(insen_si_vis);

    current_z = ib_length/2. + rec_length/2.;
    pos = G4ThreeVector(0., 0., current_z);
    new G4PVPlacement(
        nullptr,
        pos,
        recLogic,
        "RECLayer",
        logicWorld,
        0,
        0,
        OverlapCheck()
    );



    return;

    // if(peri_width > 0){
    //      // for dead area
    //     G4VSolid* insen_solid = new G4Box("Insensitive_si", peri_width/2., si_thickness/2., length/2.);
    //     G4LogicalVolume* insen_logic = new G4LogicalVolume(insen_solid, mat_si, "SiSubLogic");
    //     G4VisAttributes* insen_vis = new G4VisAttributes(col_si_insensitive);
    //     insen_vis->SetForceSolid(true);
    //     insen_logic->SetVisAttributes(insen_vis);
        
    //     for(int i = 0; i < nphi; i++){
    //         double phi = i * 2*M_PI/nphi;
    //         G4RotationMatrix rotm = G4RotationMatrix();
    //         rotm.rotateZ(M_PI/2. + phi);
    //         G4ThreeVector pos = G4ThreeVector(
    //             radius*cos(phi),
    //             radius*sin(phi),
    //             0
    //             );

    //         G4Transform3D transform = G4Transform3D(rotm, pos); 
    //         new G4PVPlacement(
    //             transform,
    //             insen_logic,
    //             "insen_si_layer",
    //             logicWorld,
    //             0,
    //             i,
    //             OverlapCheck()
    //         );
    //     }
    
    //     // Left/right endcap
    //     G4VSolid* lec_solid = new G4Tubs("lec_cylinder", radius, radius+si_thickness, lec_length/2., 0, 2*M_PI);
    //     G4LogicalVolume* lec_logic = new G4LogicalVolume(lec_solid, mat_si, "LECLogic");
    //     lec_logic->SetVisAttributes(insen_vis);

    //     new G4PVPlacement(
    //         nullptr,
    //         G4ThreeVector(0, 0, -(length + lec_length)/2. ),
    //         lec_logic,
    //         "LeftEndcap_cyl",
    //         logicWorld,
    //         0,
    //         0,
    //         OverlapCheck()
    //     );

    //     G4VSolid* rec_solid = new G4Tubs("rec_cylinder", radius, radius+si_thickness, rec_length/2., 0, 2*M_PI);
    //     G4LogicalVolume* rec_logic = new G4LogicalVolume(rec_solid, mat_si, "RECLogic");
    //     rec_logic->SetVisAttributes(insen_vis);

    //     new G4PVPlacement(
    //         nullptr,
    //         G4ThreeVector(0, 0, (length + rec_length)/2. + 2*m_nooverlap),
    //         rec_logic,
    //         "RightEndcap_cyl",
    //         logicWorld,
    //         0,
    //         0,
    //         OverlapCheck()
    //     );
    
    // }
    

}

void ePIC_SVTIB_Detector::Print(const std::string &what) const{
	std::cout<<"ePIC SVT IB Detector: "<<std::endl;
	if(what == "ALL" || what == "VOLUME"){
		std::cout<<"Version 0.1"<<std::endl;
	}
	return;
}
