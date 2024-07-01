/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_TOFBarrel_Detector.cc
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.07.01
*   	Description：	
*
=========================================================================*/

#include <iostream>
#include "ePIC_TOFBarrel_Detector.h"
#include <fstream>
#include <sstream>

#include <g4main/PHG4Detector.h>

#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Color.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4Transform3D.hh>
#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4NistManager.hh>
#include <Geant4/G4Element.hh>

#include <phparameter/PHParameters.h>

#include <cmath>
#include <iostream>

#include <TSystem.h>

class G4VSolid;
class PHCompositeNode;

using namespace std;

ePIC_TOFBarrel_Detector::ePIC_TOFBarrel_Detector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string &dname, const int lyr)
	:PHG4Detector(subsys, Node, dname)
	, m_Params(parameters)
	, m_Layer(lyr){

        DefineMaterials();

}

int ePIC_TOFBarrel_Detector::IsInDetector(G4VPhysicalVolume* volume) const{
	set<G4VPhysicalVolume*>::const_iterator iter = m_PhysicalVolumesSet.find(volume);
	if(iter != m_PhysicalVolumesSet.end()){
		return 1;
	}
	return 0;
}


void ePIC_TOFBarrel_Detector::DefineMaterials(){
    G4NistManager* man = G4NistManager::Instance();
    bool isotopes = false;

    // build NOVEC_7200
    G4Element* elC = man->FindOrBuildElement("C", isotopes);
    G4Element* elH = man->FindOrBuildElement("H", isotopes);
    G4Element* elO = man->FindOrBuildElement("O", isotopes);
    G4Element* elF = man->FindOrBuildElement("F", isotopes);

    double density = 1.4811 * g/cm3;
    G4Material* mat_NOVEC7200 = new G4Material("NOVEC_7200", density, 4);
    mat_NOVEC7200->AddElement(elC, 6);
    mat_NOVEC7200->AddElement(elH, 5);
    mat_NOVEC7200->AddElement(elO, 1);
    mat_NOVEC7200->AddElement(elF, 9);

    cout<<*(G4Material::GetMaterialTable())<<endl;

    
}

void ePIC_TOFBarrel_Detector::ConstructMe(G4LogicalVolume* logicWorld){

    double zmax = m_Params->get_double_param("zmax") * cm;
    double zmin = m_Params->get_double_param("zmin") * cm;
    double radius = m_Params->get_double_param("radius") * cm;
    double sensor_width = m_Params->get_double_param("sensor_width") * cm;
    double module_width = m_Params->get_double_param("module_width") * cm;
    double coolingtube_width = m_Params->get_double_param("coolingtube_width") * cm;
    double sensor_t = m_Params->get_double_param("sensor_thickness") * cm;
    double hybrid_t = m_Params->get_double_param("hybrid_thickness") * cm;
    double CFskin_t = m_Params->get_double_param("CFskin_thickness") * cm;
    double CoolingTube_t = m_Params->get_double_param("CoolingTube_thickness") * cm;
    double Coolant_t = m_Params->get_double_param("Coolant_thickness") * cm;
    double CFoam_t = m_Params->get_double_param("CFoam_thickness") * cm;
    double CHoneycomb_t = m_Params->get_double_param("CHoneycomb_thickness") * cm;

 	//int stave_number = m_Params->get_int_param("stave_number"); 

    if(!std::isfinite(zmax) || !std::isfinite(zmin) || !std::isfinite(radius) ||
        !std::isfinite(sensor_width) || !std::isfinite(sensor_t) || !std::isfinite(hybrid_t) ||
        !std::isfinite(CFskin_t) || !std::isfinite(CoolingTube_t) || !std::isfinite(Coolant_t) ||
        !std::isfinite(CFoam_t) || !std::isfinite(CHoneycomb_t) || !std::isfinite(CFskin_t) ||
        !std::isfinite(CFskin_t)){
            cout<<"PHWHERE "<< ": Bad Parameters for "<< GetName() << endl;
            cout<<"zmax: "<<zmax<<" zmin: "<<zmin<<" radius: "<<radius<<endl;
            cout<<"sensor_width: "<<sensor_width<<" sensor_t: "<<sensor_t<<" hybrid_t: "<<hybrid_t<<endl;
            cout<<"CFskin_t: "<<CFskin_t<<" CoolingTube_t: "<<CoolingTube_t<<" Coolant_t: "<<Coolant_t<<endl;
            cout<<"CFoam_t: "<<CFoam_t<<" CHoneycomb_t: "<<CHoneycomb_t<<endl;
        }


    double stave_length = zmax - zmin;
    //double z_offset = (zmax + zmin)/2.;
    
    // construct the gas cylinder to contain the inner AC-LGAD TOF detector
    G4Colour col_kapton = G4Colour(146./255., 12./255., 201./255., 1.);
    G4Colour col_al = G4Colour(235./255., 230./255., 237./255.);
    G4Colour col_si = G4Colour(251./255., 248./255., 0., 1.);
    G4Colour col_cf = G4Colour(0., 0., 0., 1.);
    G4Colour col_novec7200 = G4Colour(0.45, 0.25, 0.0, 0.4);

    G4NistManager* nist = G4NistManager::Instance();

    G4Material* mat_si = nist->FindOrBuildMaterial("G4_Si");
    G4VSolid* sensor_solid = new G4Box("sensor_plane", sensor_width/2., sensor_t/2., stave_length/2.);
    G4LogicalVolume* sensor_logic = new G4LogicalVolume(sensor_solid, mat_si, "SensorLogic");
    G4VisAttributes* sensor_vis = new G4VisAttributes(col_si); 
    sensor_vis->SetForceSolid(true);
    sensor_logic->SetVisAttributes(sensor_vis);


    // build a container for all the insensitive components
    G4Material* mat_con = nist->FindOrBuildMaterial("mRICH_Air");
    double cont_t = 2*hybrid_t + 2*CFskin_t + Coolant_t + CoolingTube_t+ CFoam_t;
    G4VSolid* container_solid = new G4Box("contain_plane", module_width/2., cont_t/2., stave_length/2.);
    G4LogicalVolume* container_logic = new G4LogicalVolume(container_solid, mat_con, "ContainerLogic");
    
    G4Material* mat_kapton = nist->FindOrBuildMaterial("G4_KAPTON");
    G4VSolid* hybrid_solid = new G4Box("hybrid_plane", module_width/2., hybrid_t/2., stave_length/2.);
    G4LogicalVolume* hybrid_logic = new G4LogicalVolume(hybrid_solid, mat_kapton, "HybridLogic");
    G4VisAttributes* hybrid_vis = new G4VisAttributes(col_kapton);
    hybrid_vis->SetForceSolid(true);
    hybrid_logic->SetVisAttributes(hybrid_vis);

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0, (cont_t - hybrid_t)/2., 0),
        hybrid_logic,
        "hybrid_top",
        container_logic,
        0,
        0,
        OverlapCheck()

    );

    G4Material* mat_cfskin = nist->FindOrBuildMaterial("CFRP_INTT");   // it should be CFRPMix2 in eic-shell, need to be modifed
    G4VSolid* cfskin_solid = new G4Box("cfskin_plane", module_width/2., CFskin_t/2., stave_length/2.);
    G4LogicalVolume* cfskin_logic = new G4LogicalVolume(cfskin_solid, mat_cfskin, "CFskinLogic");
    G4VisAttributes* cfskin_vis = new G4VisAttributes(col_cf);
    cfskin_vis->SetForceSolid(true);
    cfskin_logic->SetVisAttributes(cfskin_vis);

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0, (cont_t - CFskin_t)/2. - hybrid_t, 0),
        cfskin_logic,
        "cfskin_top",
        container_logic,
        0,
        0,
        OverlapCheck()
    );
    
    G4Material* mat_al = nist->FindOrBuildMaterial("G4_Al");
    G4VSolid* CoolingTube_solid = new G4Box("CoolingTube_plane", coolingtube_width/2., CoolingTube_t/2., stave_length/2.);
    G4LogicalVolume* CoolingTube_logic = new G4LogicalVolume(CoolingTube_solid, mat_al, "CoolTubeLogic");
    G4VisAttributes* CoolingTube_vis = new G4VisAttributes(col_al);
    CoolingTube_vis->SetForceSolid(true);
    CoolingTube_logic->SetVisAttributes(CoolingTube_vis);

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0, (cont_t - CoolingTube_t)/2. - hybrid_t - CFskin_t, 0),
        CoolingTube_logic,
        "CoolingTube_stave",
        container_logic,
        0,
        0,
        OverlapCheck()
    );



    G4Material* mat_NOVEC7200 = nist->FindOrBuildMaterial("NOVEC_7200");
    G4VSolid* Coolant_solid = new G4Box("Coolant_plane", coolingtube_width/2., Coolant_t/2., stave_length/2.);
    G4LogicalVolume* Coolant_logic = new G4LogicalVolume(Coolant_solid, mat_NOVEC7200, "CoolantLogic");
    G4VisAttributes* Coolant_vis = new G4VisAttributes(col_novec7200);
    Coolant_vis->SetForceSolid(true);
    Coolant_logic->SetVisAttributes(Coolant_vis);

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0, (cont_t - Coolant_t)/2. - hybrid_t - CFskin_t - CoolingTube_t, 0),
        Coolant_logic,
        "Coolant_stave",
        container_logic,
        0,
        0,
        OverlapCheck()
    );


    double cfoam_width = sensor_width - 5.*mm;
    G4Material* mat_cfoam = nist->FindOrBuildMaterial("CF");
    G4VSolid* CFoam_solid = new G4Box("CF_plane", cfoam_width/2., CFoam_t/2., stave_length/2.);
    G4LogicalVolume* CFoam_logic = new G4LogicalVolume(CFoam_solid, mat_cfoam, "CarbonFoamLogic");
    G4VisAttributes* CFoam_vis = new G4VisAttributes(col_cf);
    CFoam_vis->SetForceSolid(true);
    CFoam_logic->SetVisAttributes(CFoam_vis);

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(-cfoam_width/2., -(cont_t - CFoam_t)/2. + hybrid_t + CFskin_t, 0),
        CFoam_logic,
        "CFoam_stave",
        container_logic,
        0,
        0,
        OverlapCheck()
    );

    double choneycomb_width = module_width - cfoam_width -0.1*cm;
    G4VSolid* CHoneycomb_solid = new G4Box("CHoneycomb_plane", choneycomb_width/2., CHoneycomb_t/2., stave_length/2.);
    G4LogicalVolume* CHoneycomb_logic = new G4LogicalVolume(CHoneycomb_solid, mat_cfskin, "CHoneycombLogic");
    G4VisAttributes* CHoneycomb_vis = new G4VisAttributes(col_al);
    CHoneycomb_vis->SetForceSolid(true);
    CHoneycomb_logic->SetVisAttributes(CHoneycomb_vis);

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(choneycomb_width/2., -(cont_t - CFoam_t)/2. + hybrid_t + CFskin_t, 0),
        CHoneycomb_logic,
        "CHoneycomb_stave",
        container_logic,
        0,
        0,
        OverlapCheck()
    );

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0, -(cont_t - CFskin_t)/2. + hybrid_t, 0),
        cfskin_logic,
        "cfskin_bot",
        container_logic,
        0,
        0,
        OverlapCheck()
    );

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0,-(cont_t - hybrid_t)/2., 0),
        hybrid_logic,
        "hybrid_bot",
        container_logic,
        0,
        0,
        OverlapCheck()
    );
    
    // place the container
    G4VPhysicalVolume* sensor_phy = 
    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0, (cont_t + sensor_t)/2., 0),
        sensor_logic,
        "sensor_stave",
        logicWorld,
        0,
        0,
        OverlapCheck()
    );

    m_PhysicalVolumesSet.insert(sensor_phy);

    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0, 0, 0),
        container_logic,
        "container_stave",
        logicWorld,
        0,
        0,
        OverlapCheck()

    );

    // //placement of the detector layers 
    // double no_overlap = 0.008*cm;
    // for(int istave =0; istave<stave_number; istave++){
    //     //Drift Cu ground
        
    //     double phi = istave * 2*M_PI/stave_number;
    //     G4RotationMatrix rotm = G4RotationMatrix();
    //     double tilt = 20./280.*M_PI;
    //     rotm.rotateZ(M_PI/2. + phi + tilt);

    //     G4ThreeVector position = G4ThreeVector(radius*std::cos(phi), radius*std::sin(phi), z_offset);
    //     G4Transform3D transform = G4Transform3D(rotm, position);
    //     G4VPhysicalVolume* sensor_phy = 
    //     new G4PVPlacement(
    //         transform, //rotation, positon
    //         sensor_logic,
    //         "sensor_stave",
    //         logicWorld,
    //         0,
    //         istave,
    //         OverlapCheck()
    //     );

    //     m_PhysicalVolumesSet.insert(sensor_phy);
       
    // //     no_overlap = 0.004*cm;
    // //     double radius_hybrid = radius + sensor_t/2. + hybrid_t/2. + no_overlap;
    // //     position = G4ThreeVector(
    // //         radius_hybrid*std::cos(phi), 
    // //         radius_hybrid*std::sin(phi), 
    // //         z_offset
    // //         );
        
    // //     transform = G4Transform3D(rotm, position);
    // //     new G4PVPlacement(
    // //         transform,
    // //         hybrid_logic,
    // //         "hybrid_stave",
    // //         logicWorld,
    // //         0,
    // //         istave,
    // //         OverlapCheck()
    // //     );

    // //     no_overlap = 0.006*cm;
    // //     double radius_cft = radius_hybrid + hybrid_t/2. + CFskin_t/2. + no_overlap;
    // //     position = G4ThreeVector(
    // //         radius_cft*std::cos(phi),
    // //         radius_cft*std::sin(phi),
    // //         z_offset
    // //         );
    // //     transform = G4Transform3D(rotm, position);
    // //     new G4PVPlacement(
    // //         transform,
    // //         cfskin_logic,
    // //         "CFskin_stave",
    // //         logicWorld,
    // //         0,
    // //         istave,
    // //         OverlapCheck()
    // //     );
        
    // //     double radius_ct = radius_cft + CFskin_t/2. + CoolingTube_t/2. + no_overlap;
    // //     position = G4ThreeVector(
    // //         radius_ct*std::cos(phi),
    // //         radius_ct*std::sin(phi),    
    // //         z_offset
    // //         );
    // //     transform = G4Transform3D(rotm, position);
    // //     new G4PVPlacement(
    // //         transform,
    // //         CoolingTube_logic,
    // //         "CoolingTube_stave",
    // //         logicWorld,
    // //         0,
    // //         istave,
    // //         OverlapCheck()
    // //     );

    // //     no_overlap = 0.007*cm;
    // //     double radius_coolant = radius_ct + CoolingTube_t/2. + Coolant_t/2. + no_overlap;
    // //     position = G4ThreeVector(
    // //         radius_coolant*std::cos(phi),
    // //         radius_coolant*std::sin(phi),
    // //         z_offset
    // //    );
    // //     transform = G4Transform3D(rotm, position);
    // //     new G4PVPlacement(
    // //         transform,
    // //         Coolant_logic,
    // //         "Coolant_stave",
    // //         logicWorld,
    // //         0,
    // //         istave,
    // //         OverlapCheck()
    // //     );

    // //     no_overlap = 0.001*cm;
    // //     double radius_cf = radius_coolant + Coolant_t/2. + CFoam_t/2. + no_overlap;
    // //     position = G4ThreeVector(
    // //         radius_cf*std::cos(phi),
    // //         radius_cf*std::sin(phi),
    // //         z_offset
    // //     );

    // //     transform = G4Transform3D(rotm, position);
    // //     new G4PVPlacement(
    // //         transform,
    // //         CFoam_logic,
    // //         "CFoam_stave",
    // //         logicWorld,
    // //         0,
    // //         istave,
    // //         OverlapCheck()
    // //     );

    // //    no_overlap = 0.25*cm;
    //     double radius_cc = radius + sensor_t/2. + CHoneycomb_t/2. + no_overlap;
    //     position = G4ThreeVector(
    //         radius_cc*std::cos(phi),
    //         radius_cc*std::sin(phi),
    //         z_offset
    //     );

    //     transform = G4Transform3D(rotm, position);
    //     new G4PVPlacement(
    //         transform,
    //         CHoneycomb_logic,
    //         "CHoneycomb_stave",
    //         logicWorld,
    //         0,
    //         istave,
    //         OverlapCheck()
    //     );


    // }

    return ;
}

void ePIC_TOFBarrel_Detector::Print(const string& what)const{
    std::cout<<"ePIC Inner MPGD Detector: "<<std::endl;
	if(what == "ALL" || what == "VOLUME"){
		std::cout<<"Version 0.1"<<std::endl;
	}
	return;
}