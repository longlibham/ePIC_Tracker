/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_InnerMPGD_Detector.cc
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.06.21
*   	Description：	
*
=========================================================================*/

#include <iostream>
#include "ePIC_InnerMPGD_Detector.h"
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

ePIC_InnerMPGD_Detector::ePIC_InnerMPGD_Detector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string &dname, const int lyr)
	:PHG4Detector(subsys, Node, dname)
	, m_Params(parameters)
	, m_Layer(lyr){

        DefineMaterials();

}

int ePIC_InnerMPGD_Detector::IsInDetector(G4VPhysicalVolume* volume) const{
	set<G4VPhysicalVolume*>::const_iterator iter = m_PhysicalVolumesSet.find(volume);
	if(iter != m_PhysicalVolumesSet.end()){
		return 1;
	}
	return 0;
}

void ePIC_InnerMPGD_Detector::DefineMaterials(){
    G4NistManager* man = G4NistManager::Instance();
    bool isotopes = false;
    
    // build the MPGD gas
    G4Material* Butane = man->FindOrBuildMaterial("G4_BUTANE");
    G4Material* Ar = man->FindOrBuildMaterial("G4_Ar");
    double density = 1.714e-3*g/cm3;
    G4Material* MPGDGAS = new G4Material("MPGDGAS", density, ncomponents = 2);
    MPGDGAS->AddMaterial(Butane, 30*perCent);
    MPGDGAS->AddMaterial(Ar, 70*perCent);

    G4Element* elC = man->FindOrBuildElement("C", isotopes);
    G4Element* elH = man->FindOrBuildElement("H", isotopes);
    G4Element* elO = man->FindOrBuildElement("O", isotopes);

    desnsity = 1.86*g/cm3
    G4Material* mat_fr4 = new G4Material("Fr-4", density, ncomponents = 3);
    mat_fr4->AddElement(elC, 11./26.*100*perCent);
    mat_fr4->AddElement(elH, 12./26.*100*perCent);
    mat_fr4->AddElement(elO, 3./26.*100*perCent);

    
}

void ePIC_InnerMPGD_Detector::ConstructMe(G4LogicalVolume* logicWorld){

    double radius = m_Params->get_double_param("radius");
	double coat_t = m_params->get_double_param("coat_thickness") * cm;
    double z_min = m_Params->get_double_param("z_min") * cm;
    double z_max = m_Params->get_double_param("z_max") * cm;
    double KaptonOverlay_t = m_Params->get_double_param("KaptonOverlay_thickness") * cm;
    double CuGround_t = m_Params->get_double_param("CuGround_thickness") * cm;
    double pcb_t = m_Params->get_double_param("PCB_thickness") * cm;
    double CuStrip_t = m_Params->get_double_param("CuStrip_thickness") * cm;
    double KaptonStrip_t = m_Params->get_double_param("KaptonStrip_thickness") * cm;
    double ResistiveStrip_t m_Params->get_double_param("ResistiveStrip_thickness") * cm;
    double gas_t = m_Params->get_double_param("gas_thickness") * cm;
    double mesh_t = m_Params->get_double_param("mesh_thickness") * cm;
    double GasGap_t = m_Params->get_double_param("GasGap_thickness") * cm;
    double DriftCuElectrode_t = m_Params->get_double_param("DriftCuElectrode_thickness") * cm;
    double DriftKapton_t = m_Params->get_double_param("DriftKapton_thickness") * cm;
    double DriftCuGround_t = m_Params->get_double_param("DriftCuGround_thickness") * cm;
    double fudge_t = m_Params->get_double_param("fudge_thickness") * cm;

	int stave_number = m_Params->get_int_param("stave_number") * cm;

    if(!std::isfinite(radius) || !std::isfinite(mpgd_t) || !std::isfinite(z_min) || !std::isfinite(z_max) || 
        !std::isfinite(KaptonOverlay_t) || !std::isfinite(CuGround_t) || !std::isfinite(pcb_t) || !std::isfinite(CuStrip_t) ||
        !std::isfinite(KaptonStrip_t) || !std::isfinite(ResistiveStrip_t) || !std::isfinite(gas_t) || 
        !std::isfinite(mesh_t) || !std::isfinite(GasGap_t) || !std::isfinite(DriftCuElectrode_t) || 
        !std::isfinite(DriftKapton_t) || !std::isfinite(DriftCuGround_t) || !std::isfinite(fudge_t)){
            cout<<"PHWHERE "<< ": Bad Parameters for "<< GetName() << endl;
            cout<<"radius: "<< radius << " mpgd_t: "<<mpgd_t<<" z_min: "<<z_min<<" z_max: "<<endl;
            cout<<"KaptonOverlay_thickness: "<<KaptonOverlay_t<<" CuGround_thickness: "<<CuGround_t<<" pcb_t: "<<pcb_t<<endl;
            cout<<"CuStrip_thickness: "<<CuStrip_t<<" KaptonStrip_thickness: "<<KaptonStrip_t<<" ResistiveStrip_thickness: "<<ResistiveStrip_t<<endl;
            cout<<"gas_thickness: "<<gas_t<<" mesh_thickness: "<<mesh_t<<" GasGap_thickness: "<<GasGap_t<<" DriftCuElectrode_thickness: "<<DriftCuElectrode_t<<endl;
            cout<<"DriftKapton_thickness: "<<DriftKapton_t<<" DriftCuGround_thickness: "<<DriftCuGround_t<<" fudge_thickness: "<<fudge_t<<endl;
        }

    // construct the gas cylinder to contain the inner MPGD detector
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* mat_mpgd_gas = nist->FindOrBuildMaterial("MPGD_GAS");
    G4Colour* col_mpgd_gas = new G4Colour(0.45, 0.25, 0.0, 0.4);
    double dangle = 2*P_PI/stave_number;
    double offset_z = z_max + z_min; 

    G4VSolid* coat_solid = new G4Tubs("ePIC_MPGD_coat", radius - MPGD_t/2., radius + MPGD_t/2., (z_max - z_min)/2., 0., 2*M_PI);
    G4LogicalVolume* coat_logic = new G4LogicalVolume(coat_solid, mat_mpgd_gas, "coatLogic");
    G4VisAttributes* coat_vis = new  G4VisAttributes(G4Color(col_mpgd_gas)); // set the coat color to be brown and alpha 0.4
    coat_vis->SetForceSolid(true);
    coat_logic->SetVisAttributes(coat_vis);

    // placement of the Ar-butane coat
    
    G4RotationMatrix rotm = new G4RotationMatrix();
    rotm.rotateZ(M_PI/2. + phi);

    new G4PVPlacement(rotm,  //no rotation
                    G4ThreeVector(0., 0., offset_z),    //at (0, 0, 0)
                    coat_logic,         //logical volume
                    "GAS_COAT",         //name
                    logicWorld,         // mother volume
                    false,              //no boolean operation
                    0,                  // copy number
                    fCheckOverlaps
    );


    // Detector construction
    //drift copper ground
    G4Material* mat_cu = nist->FindOrBuildMaterial("G4_Cu");
    G4Colour* col_cu = new G4Colour(221./255., 125./255., 4./255., 1.);
    double width = radius*std::tan(2*M_PI/stave_number);
    G4VSolid* DriftCuGround_solid = new G4Box("DriftCuGround_plane", width/2., DriftCuGround_t/2., (z_max - z_min)/2.);
    G4LogicalVolume* DriftCuGround_logic = new G4LogicalVolume(DriftCuGround_solid, mat_cu, "DriftCuGroundLogic");
    G4VisAttributes* DriftCuGround_vis = new G4VisAttributes(G4Color(col_cu)); 
    DriftCuGround_vis->SetForceSolid(true);
    DriftCuGround_logic->SetVisAttributes(DriftCuGround_vis);

    // drift kapton
    G4Material* mat_kapton = nist->FindOrBuildMaterial("G4_KAPTON");
    G4Colour* col_kapton = new G4Colour(146./255., 12./255., 201./255., 1.);

    G4VSolid* DriftKapton_solid = new G4Box("DriftKapton_plane", width/2., DriftKapton_t/2., (z_max - z_min)/2.);
    G4LogicalVolume* DriftKapton_logic = new G4LogicalVolume(DriftKapton_solid, mat_kapton, "DriftKaptonLogic");
    G4VisAttributes* DrfitKapton_vis = new G4VisAttributes(G4Color(col_kapton));
    DriftKapton_vis->SetForceSolid(true);
    DriftKapton_logic->SetVisAttributes(DriftKapton_logic);

    //Drift copper electrode
    G4VSolid* DriftCuElectrode_solid = new G4Box("DriftCuElectrode_plane", width/2., DriftCuElectrode_t/2., (z_max-z_min)/2.);
    G4LogicalVolume* DriftCuElectrode_logic = new G4LogicalVolume(DriftCuElectrode_solid, mat_cu, "DriftCuElectrodeLogic");
    G4VisAttributes* DriftCuElectrode_vis = new G4VisAttributes(G4Color(col_cu));
    DriftCuElectrode_vis->SetForceSolid(true);
    DriftCuElectrode_logic->SetVisAttributes(DriftCuElectrode_vis);

    //gas gap; sensitive area
    G4VSolid* gasgap_solid = new G4Box("GasGap_plane", width/2., GasGap_t/2., (z_max-z_min)/2.);
    G4LogicalVolume* gasgap_logic = new G4LogicalVolume(gasgap_solid, mat_mpgd_gas, "GasGapLogic");
    G4VisAttributes* gasgap_vis = new G4VisAttributes(G4Color(col_mpgd_gas));
    gasgap_vis->SetForceSolid(true);
    gasgap_logic->SetVisAttributes(gasgap_vis);

    // Mesh
    G4Material* mat_mesh = nist->FindOrBuildMaterial("G4_Fe");
    G4Colour* col_inox = G4Colour(235./255., 230./255., 237./255.)
    G4VSolid* mesh_solid = new G4Box("mesh_plane", width/2., mesh_t/2., (z_max-zmin)/2.);
    G4LogicalVolume* mesh_logic = new G4LogicalVolume(mesh_solid, mat_mesh, "MeshLogic");
    G4VisAttributes* mesh_vis = new G4VisAttributes(G4Color(col_inox));
    mesh_vis->SetForceSolid(true);
    mesh_logic->SetVisAttributes(mesh_vis);
    
    //fudge
    G4VSolid* fudge_solid = new G4Box("fudge_plane", width/2., fudge_t/2., (z_max-z_min)/2.);
    G4LogicalVolume* fudge_logic = new G4LogicalVolume(fudge_solid, mat_kapton, "FudgeLogic");
    G4VisAttributes* fudge_vis = new G4VisAttributes(G4Color(col_kapton));
    fudge_vis->SetForceSolid(true);
    fudge_logic->SetVisAttributes(fudge_vis);

    //gas
    G4VSolid* gas_solid = new G4Box("Gas_plane", width/2., gas_t/2., (z_max-z_min)/2.);
    G4LogicalVolume* gas_logic = new G4LogicalVolume(gas_solid, mat_mpgd_gas, "GasLogic");
    G4VisAttributes* gas_vis = new G4VisAttributes(G4Color(col_mpgd_gas));
    gas_vis->SetForceSolid(true);
    gas_logic->SetVisAttributes(gas_vis);

    //ResistiveStrips
    G4Material* mat_rstrip = nist->FindOrBuildMaterial("G4_C");
    G4Colour* col_carbon = new G4Colour(0., 0., 0., 1.);
    G4VSolid* ResistiveStrip_solid = new G4Box("ResistiveStrip_plane", width/2., ResistiveStrip_t/2., (z_max-z_min)/2.);
    G4LogicalVolume* ResistiveStrip_logic = new G4LogicalVolume(ResistiveStrip_solid, mat_rstrip, "ResistiveStripLogic");
    G4VisAttributes* ResistiveStrip_vis = new G4VisAttributes(G4Color(col_carbon));
    ResistiveStrip_vis->SetForceSolid(true);
    ResistiveStrip_logic->SetVisAttributes(ResistiveStrip_vis);


    //Kapton strip
    G4VSolid* KaptonStrip_solid = new G4Box("KaptonStrip_plane", width/2., KaptonStrip_t/2., (z_max-z_min)/2.);
    G4LogicalVolume*  KaptonStrip_logic = new G4LogicalVolume(KaptonStrip_solid, mat_kapton, "KaptonStripLogic");
    G4VisAttributes* KaptonStrip_vis = new G4VisAttributes(G4Color(col_kapton));
    KaptonStrip_vis->SetForceSolid(true);
    KaptonStrip_logic->SetVisAttributes(KaptonStrip_vis);

    //copper strip
    G4VSolid* CuStrip_solid = new G4Box("CuStrip_plane", width/2., CuStrip_t/2., (z_max-z_min)/2.);
    G4LogicalVolume* CuStrip_logic = new G4LogicalVolume(CuStrip_solid, mat_cu, "CuStripLogic");
    G4VisAttributes* CuStrip_vis = new G4VisAttributes(G4Color(col_cu));
    CuStrip_vis->SetForceSolid(true);
    CuStrip_logic->SetVisAttributes(CuStrip_vis);

    //PCB
    G4Material* mat_fr4 = nist->FindOrBuildMaterial("Fr-4");
    G4Colour* col_fr4 = new G4Colour(0., 67., 0., 1.);
    G4VSolid* pcb_solid = new G4Box("pcb_plane", width/2., pcb_t/2., (z_max-z_min)/2.);
    G4LogicalVolume* pcb_logic = new G4LogicalVolume(pcb_solid, mat_fr4, "PCBLogic");
    G4VisAttributes* pcb_vis = new G4VisAttributes(G4Color(col_fr4));
    pcb_vis->SetForceSolid(true);
    pcb_logic->SetVisAttributes(pcb_vis);

    //coppr ground 
    G4VSolid* CuGround_solid = new G4Box("CuGround_plane", width/2., CuGround_t/2., (z_max-z_min)/2.);
    G4LogicalVolume* CuGround_logic = new G4LogicalVolume(CuGround_solid, mat_cu, "CuGroundLogic");
    G4VisAttributes* CuGround_vis = new G4VisAttributes(G4Color(col_cu));
    CuGround_vis->SetForceSolid(true);
    CuGround_logic->SetVisAttributes(CuGround_vis);

    //Kapton overlay
    G4VSolid* KaptonOverlay_solid = new G4Box("KaptonOverlay_plane", width/2., KaptonOverlay_t/2., (z_max-z_min)/2.);
    G4LogicalVolume* KaptonOverlay_logic = new G4LogicalVolume(KaptonOverlay_solid, mat_kapton, "KaptonOverlayLogic");
    G4VisAttributes* KaptonOverlay_vis = new G4VisAttributes(G4Color(col_kapton));
    KaptonOverlay_vis->SetForceSolid(true);
    KaptonOverlay_logic->SetVisAttributes(KaptonOverlay_vis);

    //placement of the detector layers inside the Ar-Butane coat
    for(int istave =0; istave<stave_number; istave++){
        //Drift Cu ground
        double phi = istave * 2*M_PI/stave_number;
        G4ThreeVector position = G4ThreeVector(radius*std::cos(phi), radius*std::sin(phi), offset_z);
        G4Transform3D transform = G4Transform3D(rotm, position);
        new G4PVPlacement(
            transform, //rotation, positon
            DriftCuGround_logic,
            "DriftCuGround_stave",
            coat_logic,
            0,
            istave,
            OverlapCheck()
        );

        //Drift Kapton
        radius_dk = radius + DriftCuGround_t;
        position = G4ThreeVector(radius_dk*std::cos(phi), radius_dk*std::sin(phi), offset_z);
        transform = G4Transform3D(rotm, position);
        new G4PVPlacement(
            transform,
            DriftKapton_logic,
            "DriftKapton_stave",
            coat_logic,
            0,
            istave,
            OverlapCheck()
        );

        //Drift Cu Electrode
        radius_dce = radius_dk + DriftKapton_t;
        position = G4ThreeVector(radius_dce*std::cos(phi), radius_dce*std::sin(phi), offset_z);
        transform = G4Transform3D(rotm, position);
        new G4PVPlacement(
            transform,
            DriftCuElectrode_logic,
            "DriftCuElectrode_stave",
            coat_logic,
            0,
            istave,
            OverlapCheck()
        );

        //gas gap
        radius_gg = radius_dce + DriftCuElectrode_t;
        position = G4ThreeVector(radius_gg*std::cos(phi), radius_gg*std::sin(phi), offset_z);
        transform = G4Transform3D(rotm, position);
        
        G4VPhysicalVolume* gasgap_phy =
        new G4PVPlacement(
            transform,
            gasgap_logic,
            "gasgap_stave",
            coat_logic,
            0,
            istave,
            OverlapCheck()
        );

        // add sensitive detector
        m_PhysicalVolumesSet.insert(gasgap_phy);


        //Mesh
        radius_mesh = radius_gg + GasGap_t;
        position = G4ThreeVector(radius_mesh*std::cos(phi), radius_mesh*std::sin(phi), offset_z);
        transform = G4Transform3D(rotm, position);
        new G4PVPlacement(
            transform,
            mesh_logic,
            "mesh_stave",
            coat_logic,
            0,
            istave,
            OverlapCheck()
        );

        //fudge
        radius_fudge = radius_mesh + mesh_t;
        position = G4ThreeVector(radius_fudge*std::cos(phi), radius_fudge*std::sin(phi), offset_z);
        transform = G4Transform3D(rotm, position);
        new G4PVPlacement(
            transform,
            fudge_logic,
            "fudge_stave",
            coat_logic,
            0,
            istave,
            OverlapCheck()
        );

        //gas
        radius_g = radius_fudge + fudge_t;
        position = G4ThreeVector(radius_g*std::cos(phi), radius_g*std::sin(phi), offset_z);
        transform = G4Transform3D(rotm, position);
        new G4PVPlacement(
            transform,
            gas_logic,
            "gas_stave",
            coat_logic,
            0,
            istave,
            OverlapCheck()
        );

        //Resistive strip
        radius_rs = radius_g + gas_t;
        position = G4ThreeVector(radius_rs*std::cos(phi), radius_rs*std::sin(phi), offset_z);
        transform = G4Transform3D(rotm, position);
        new G4PVPlacement(
            transform,
            ResistiveStrip_logic,
            "ResistiveStrip_stave",
            coat_logic,
            0,
            istave,
            OverlapCheck()
        );

        //Kapton strip
        radius_ks = radius_rs + ResistiveStrip_t;
        position = G4ThreeVector(radius_ks*std::cos(phi), radius_ks*std::sin(phi), offset_z);
        transform = G4Transform3D(rotm, position);
        new G4PVPlacement(
            transform,
            KaptonStrip_logic,
            "KaptonStrip_stave",
            coat_logic,
            0, 
            istave,
            OverlapCheck()
        );

        //Cu strip
        radius_cus = radius_ks + KaptonStrip_t;
        position = G4ThreeVector(radius_cus*std::cos(phi), radius_cus*std::sin(phi), offset_z);
        transform = G4Transform3D(rotm, position);
        new G4PVPlacement(
            transform,
            CuStrip_logic,
            "CuStrip_stave",
            coat_logic,
            0,
            istave,
            OverlapCheck()
        );

        //PCB
        radius_pcb = radius_cus + CuStrip_t;
        position = G4ThreeVector(radius_pcb*std::cos(phi), radius_pcb*std::sin(phi), offset_z);
        transform = G4Transform3D(rotm, position);
        new G4PVPlacement(
            transform,
            pcb_logic,
            "PCB_stave",
            coat_logic,
            0,
            istave,
            OverlapCheck()
        );

        //Cu ground
        radius_cug = radius_pcb + pcb_t;
        position = G4ThreeVector(radius_cug*std::cos(phi), radius_cug*std::sin(phi), offset_z);
        transform = G4Transform3D(rotm, position);
        new G4PVPlacement(
            transform,
            CuGround_logic,
            "CuGround_stave",
            coat_logic,
            0,
            istave,
            OverlapCheck()
        );

        //Kapton Overlay
        radius_ko = radius_cug + CuGround_t;
        position = G4ThreeVector(radius_ko*std::cos(phi), radius_ko*std::sin(phi), offset_z);
        transform = G4Transform3D(rotm, position);
        new G4PVPlacement(
            transform,
            KaptonOverlay_logic,
            "KaptonOverlay_stave",
            coat_logic,
            0,
            istave,
            OverlapCheck()
        );

    }

    return ;
}

void ePIC_InnerMPGD_Detector::Print(const string& what){
    std::cout<<"ePIC Inner MPGD Detector: "<<std::endl;
	if(what == "ALL" || what == "VOLUME"){
		std::cout<<"Version 0.1"<<std::endl;
	}
	return;
}