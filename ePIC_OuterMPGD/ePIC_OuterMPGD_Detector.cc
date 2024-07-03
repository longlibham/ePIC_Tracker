/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_OuterMPGD_Detector.cc
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.07.02
*   	Description：	
*
=========================================================================*/

#include <iostream>
#include "ePIC_OuterMPGD_Detector.h"
#include <fstream>
#include <sstream>

#include <g4main/PHG4Detector.h>

#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Box.hh>
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
#include <Geant4/G4NistManager.hh>
#include <Geant4/G4Element.hh>

#include <phparameter/PHParameters.h>

#include <cmath>
#include <iostream>

#include <TSystem.h>

class G4VSolid;
class PHCompositeNode;

using namespace std;

ePIC_OuterMPGD_Detector::ePIC_OuterMPGD_Detector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string &dname, const int lyr)
	:PHG4Detector(subsys, Node, dname)
	, m_Params(parameters)
	, m_Layer(lyr){

        DefineMaterials();

}

int ePIC_OuterMPGD_Detector::IsInDetector(G4VPhysicalVolume* volume) const{
	set<G4VPhysicalVolume*>::const_iterator iter = m_PhysicalVolumesSet.find(volume);
	if(iter != m_PhysicalVolumesSet.end()){
		return 1;
	}
	return 0;
}

void ePIC_OuterMPGD_Detector::DefineMaterials(){
    // G4NistManager* man = G4NistManager::Instance();
    // bool isotopes = false;
    
    // // build the MPGD gas
    // G4Material* Butane = man->FindOrBuildMaterial("G4_BUTANE");
    // G4Material* Ar = man->FindOrBuildMaterial("G4_Ar");
    // double density = 1.714e-3*g/cm3;
    // G4Material* MPGDGAS = new G4Material("MPGDGAS", density, 2);
    // MPGDGAS->AddMaterial(Butane, 30*perCent);
    // MPGDGAS->AddMaterial(Ar, 70*perCent);

    // G4Element* elC = man->FindOrBuildElement("C", isotopes);
    // G4Element* elH = man->FindOrBuildElement("H", isotopes);
    // G4Element* elO = man->FindOrBuildElement("O", isotopes);

    // density = 1.86*g/cm3;
    // G4Material* mat_fr4 = new G4Material("Fr-4", density, 3);
    // mat_fr4->AddElement(elC, 11./26.*100*perCent);
    // mat_fr4->AddElement(elH, 12./26.*100*perCent);
    // mat_fr4->AddElement(elO, 3./26.*100*perCent);

    // cout<<*(G4Material::GetMaterialTable())<<endl;

    
}

void ePIC_OuterMPGD_Detector::ConstructMe(G4LogicalVolume* logicWorld){

    double radius = m_Params->get_double_param("radius") * cm;
	double module_width = m_Params->get_double_param("module_width") * cm;
    double zmin = m_Params->get_double_param("zmin") * cm;
    double zmax = m_Params->get_double_param("zmax") * cm;
    double window_t = m_Params->get_double_param("window_thickness") * cm;
    double windowgap_t = m_Params->get_double_param("windowgap_thickness") * cm;
    double driftgap_t = m_Params->get_double_param("driftgap_thickness") * cm;
    double foilcu_t = m_Params->get_double_param("foilcu_thickness") * cm;
    double readoutelectrode_t = m_Params->get_double_param("readoutelectrode_thickness") * cm;
    double foilkapton_t = m_Params->get_double_param("foilkapton_thickness") * cm;
    double readoutnomex_t = m_Params->get_double_param("readoutnomex_thickness") * cm;
    double readoutkapton_t = m_Params->get_double_param("readoutkapton_thickness") * cm;
    double pcb_t = m_Params->get_double_param("pcb_thickness") * cm;
    double cf_t = m_Params->get_double_param("carbonfoam_thickness") * cm;


 	int stave_phi = m_Params->get_int_param("stave_phi");
    int stave_z = m_Params->get_int_param("stave_z"); 

    if(!std::isfinite(radius) || !std::isfinite(module_width) || !std::isfinite(zmin) ||
        !std::isfinite(zmax) || !std::isfinite(window_t) || !std::isfinite(windowgap_t) ||
        !std::isfinite(driftgap_t) || !std::isfinite(foilcu_t) || !std::isfinite(readoutelectrode_t) ||
        !std::isfinite(foilkapton_t) || !std::isfinite(readoutnomex_t) || !std::isfinite(readoutkapton_t) ||
        !std::isfinite(pcb_t) || std::isfinite(stave_phi) || std::isfinite(stave_z)
    
    ){
            cout<<"PHWHERE "<< ": Bad Parameters for "<< GetName() << endl;
            cout<< "radius: " << radius << " module_width: " << module_width <<endl;
            cout<< "zmin: " << zmin << " zmax: " << zmax <<endl;
            cout<< "window_t: " << window_t << " windowgap_t: " << windowgap_t <<endl;
            cout<< "driftgap_t: " << driftgap_t << " foilcu_t: " << foilcu_t <<endl;
            cout<< "readoutelectrode_t: " << readoutelectrode_t << " foilkapton_t: " << foilkapton_t <<endl;
            cout<< "readoutnomex_t: " << readoutnomex_t << " readoutkapton: " << readoutkapton_t <<endl;
            cout<< "pcb_t: " << pcb_t << " stave_phi: " << stave_phi << " stave_z: " << stave_z <<endl;
           
        }

//     // construct the gas cylinder to contain the Outer MPGD detector
    // drift gap
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* mat_Ar = nist->FindOrBuildMaterial("G4_Ar");
    G4Colour col_Ar = G4Colour(0.45, 0.25, 0.0, 0.4);
    G4Colour col_kapton = G4Colour(214./255., 116./255., 7./255., 0.7);
    G4Colour col_cu = G4Colour(246./255., 167./255., 1./255., 1.);
    G4Colour col_nomex = G4Colour(1., 1., 1., 0.7);
    G4Colour Col_fr4 = G4Colour(26./255., 117./255., 1./255., 1.);
    G4Colour col_cf = G4Color(0., 0., 0., 1.);

    double no_overlap_z = 0.1 * cm; 
    double offset_z = (zmax + zmin)/2.; 
    double z_length = (zmax - zmin)/2. - no_overlap_z/2.; 
    G4VSolid* driftgap_solid = new G4Box("driftgap_plane", module_width/2., driftgap_t/2., z_length/2.);
    G4LogicalVolume* driftgap_logic = new G4LogicalVolume(driftgap_solid, mat_Ar, "DriftGapLogic");
    G4VisAttributes* driftgap_vis = new  G4VisAttributes(col_Ar); // set the coat color to be brown and alpha 0.4
    driftgap_vis->SetForceSolid(true);
    driftgap_logic->SetVisAttributes(driftgap_vis);

    // window gasgap
    G4VSolid* windowgap_solid = new G4Box("windowgasgap_plane", module_width/2., windowgap_t/2., z_length/2.);
    G4LogicalVolume* windowgap_logic = new G4LogicalVolume(windowgap_solid, mat_Ar, "WindowGasGapLogic");
    G4VisAttributes* windowgap_vis = new G4VisAttributes(col_Ar);
    windowgap_vis->SetForceSolid(true);
    windowgap_logic->SetVisAttributes(windowgap_vis);

    // window
    G4Material* mat_kapton = nist->FindOrBuildMaterial("G4_KAPTON");
    G4VSolid* window_solid = new G4Box("window_palne", module_width/2., window_t/2., z_length/2.);
    G4LogicalVolume* window_logic = new G4LogicalVolume(window_solid, mat_kapton, "WinsdowLogic");
    G4VisAttributes* window_vis = new G4VisAttributes(col_kapton);
    window_vis->SetForceSolid(true);
    window_logic->SetVisAttributes(window_vis);

    // cathode kapton
    G4VSolid* foilkapton_solid = new G4Box("folikapton_plane", module_width/2., foilkapton_t/2., z_length/2.);
    G4LogicalVolume* foilkapton_logic = new G4LogicalVolume(foilkapton_solid, mat_kapton, "FoilKaptonLogic");
    G4VisAttributes* foilkapton_vis = new G4VisAttributes(col_kapton);
    foilkapton_vis->SetForceSolid(true);
    foilkapton_logic->SetVisAttributes(foilkapton_vis);

    //cathode cu
    G4Material* mat_cu = nist->FindOrBuildMaterial("G4_Cu");
    G4VSolid* foilcu_solid = new G4Box("foilcu_plane", module_width/2., foilcu_t/2., z_length/2.);
    G4LogicalVolume* foilcu_logic = new G4LogicalVolume(foilcu_solid, mat_cu, "FoilCuLogic");
    G4VisAttributes* foilcu_vis = new G4VisAttributes(col_cu);
    foilcu_vis->SetForceSolid(true);
    foilcu_logic->SetVisAttributes(foilcu_vis);

    //RWELL cu
    G4VSolid* rwellcu_solid = new G4Box("rwellcu_plane", module_width/2., foilcu_t/2., z_length/2.);
    G4LogicalVolume* rwellcu_logic = new G4LogicalVolume(rwellcu_solid, mat_cu, "RWELLCuLogic");
    G4VisAttributes* rwellcu_vis = new G4VisAttributes(col_cu);
    rwellcu_vis->SetForceSolid(true);
    rwellcu_logic->SetVisAttributes(rwellcu_vis);

    //RWELL kapton
    G4VSolid* rwellkapton_solid = new G4Box("rwellkapton_plane", module_width/2., foilkapton_t/2., z_length/2.);
    G4LogicalVolume* rwellkapton_logic = new G4LogicalVolume(rwellkapton_solid, mat_kapton, "RWELLKaptonLogic");
    G4VisAttributes* rwellkapton_vis = new G4VisAttributes(col_kapton);
    rwellkapton_vis->SetForceSolid(true);
    rwellkapton_logic->SetVisAttributes(rwellkapton_vis);

    //Nomex
    G4Material* mat_nomex = nist->FindOrBuildMaterial("NOMEX");
    G4VSolid* nomex_solid = new G4Box("nomex_plane", module_width/2., readoutnomex_t/2., z_length/2.);
    G4LogicalVolume* nomex_logic = new G4LogicalVolume(nomex_solid, mat_nomex, "NomexLogic");
    G4VisAttributes* nomex_vis = new G4VisAttributes(col_nomex);
    nomex_vis->SetForceSolid(true);
    nomex_logic->SetVisAttributes(nomex_vis);

    //readout electrodes
    G4VSolid* readouteletrode_solid = new G4Box("readoutelectrode_plane", module_width/2., readoutelectrode_t/2., z_length/2.);
    G4LogicalVolume* readoutelectrode_logic = new G4LogicalVolume(readouteletrode_solid, mat_cu, "ReadoutElectrodeLogic");
    G4VisAttributes* readoutelectrode_vis = new G4VisAttributes(col_cu);
    readoutelectrode_vis->SetForceSolid(true);
    readoutelectrode_logic->SetVisAttributes(readoutelectrode_vis);

    //readoutkapton
    G4VSolid* readoutkapton_solid = new G4Box("readoutkapton_plane", module_width/2., readoutkapton_t/2., z_length/2.);
    G4LogicalVolume* readoutkapton_logic = new G4LogicalVolume(readoutkapton_solid, mat_kapton, "ReadoutKaptonLogic");
    G4VisAttributes* readoutkapton_vis = new G4VisAttributes(col_kapton);
    readoutkapton_vis->SetForceSolid(true);
    readoutkapton_logic->SetVisAttributes(readoutkapton_vis);

    //PCB
    G4Material* mat_fr4 = nist->FindOrBuildMaterial("FR4");
    G4VSolid* pcb_solid = new G4Box("pcb_plane", module_width/2., pcb_t/2., z_length/2.);
    G4LogicalVolume* pcb_logic = new G4LogicalVolume(pcb_solid, mat_fr4, "PCBLogic");
    G4VisAttributes* pcb_vis = new G4VisAttributes(Col_fr4);
    pcb_vis->SetForceSolid(true);
    pcb_logic->SetVisAttributes(pcb_vis);

    // carbon foam service
    G4Material* mat_cf = nist->FindOrBuildMaterial("CF");
    G4VSolid* cf_solid = new G4Box("cf_plane", module_width/2., cf_t/2., z_length/2.);
    G4VSolid* aircut_solid = new G4Box("aircut_plane", (module_width - 5.)/2., (cf_t + driftgap_t)/2., (z_length-5.)/2.);
    G4VSolid* carbonfoamservice_solid = new G4SubtractionSolid("carbonfoam_plane", cf_solid, aircut_solid);
    G4LogicalVolume* carbonfoamservice_logic = new G4LogicalVolume(carbonfoamservice_solid, mat_cf, "CarbonFoamLogic");
    G4VisAttributes* carbonfoam_vis = new G4VisAttributes(col_cf);
    carbonfoam_vis->SetForceSolid(true);
    carbonfoamservice_logic->SetVisAttributes(carbonfoam_vis);

    // displacement of the modules
    for(int nz = 0; nz < stave_z; nz++){
        for(int nphi = 0; nphi < stave_phi; nphi++){
            double phi = nphi * 2 * M_PI/stave_phi;
            G4RotationMatrix rotm = G4RotationMatrix();
            rotm.rotateZ(M_PI/2. + phi);

            double posz = offset_z + (2*nz-1)*(no_overlap_z/2. + z_length/2.);
            G4ThreeVector position = G4ThreeVector(
                radius*std::cos(phi), 
                radius*std::sin(phi), 
                posz
                );

            int cpno = nz*stave_phi+nphi;
            G4Transform3D transform = G4Transform3D(rotm, position);
            G4VPhysicalVolume* gasgap_phy = 
            new G4PVPlacement(
                transform,
                driftgap_logic,
                "driftgap_stave",
                logicWorld,
                0,
                cpno,
                OverlapCheck()
            );

            m_PhysicalVolumesSet.insert(gasgap_phy);

            // carbon foam service 
            double radius_cf = radius - driftgap_t/2. - cf_t/2.;
            position = G4ThreeVector(
                radius_cf*std::cos(phi),
                radius_cf*std::sin(phi),
                posz
            );

            transform = G4Transform3D(rotm, position);
            new G4PVPlacement(
                transform,
                carbonfoamservice_logic,
                "carbonfoam_stave",
                logicWorld,
                0,
                cpno,
                OverlapCheck()
            );


            double radius_wgg = radius + driftgap_t/2. + windowgap_t/2.;
            position = G4ThreeVector(
                radius_wgg*std::cos(phi),
                radius_wgg*std::sin(phi),
                posz
            );
            
            transform = G4Transform3D(rotm, position);
            new G4PVPlacement(
                transform,
                windowgap_logic,
                "windowgap_stave",
                logicWorld,
                0,
                cpno,
                OverlapCheck()
            );

            double radius_w = radius_wgg + windowgap_t/2. + window_t/2.;
            position = G4ThreeVector(
                radius_w*std::cos(phi),
                radius_w*std::sin(phi),
                posz
            );

            transform = G4Transform3D(rotm, position);
            new G4PVPlacement(
                transform,
                window_logic,
                "window_stave",
                logicWorld,
                0,
                cpno,
                OverlapCheck()
            );

            double radius_fk = radius_w + window_t/2. + foilkapton_t/2.;
            position = G4ThreeVector(
                radius_fk*std::cos(phi),
                radius_fk*std::sin(phi),
                posz
            );

            transform = G4Transform3D(rotm, position);
            new G4PVPlacement(
                transform,
                foilkapton_logic,
                "foilkapton_stave",
                logicWorld,
                0,
                cpno,
                OverlapCheck()
            );

            double radius_fcu = radius_fk + foilkapton_t/2. + foilcu_t/2.;
            position = G4ThreeVector(
                radius_fcu*std::cos(phi),
                radius_fcu*std::sin(phi),
                posz
            );

            transform = G4Transform3D(rotm, position);
            new G4PVPlacement(
                transform,
                foilcu_logic,
                "foilcu_stave",
                logicWorld,
                0,
                cpno,
                OverlapCheck()
            );

            double radius_rwcu = radius_fcu + foilcu_t/2. + foilcu_t/2.;
            position = G4ThreeVector(
                radius_rwcu*std::cos(phi),
                radius_rwcu*std::sin(phi),
                posz
            );

            transform = G4Transform3D(rotm, position);
            new G4PVPlacement(
                transform,
                rwellcu_logic,
                "rwellcu_stave",
                logicWorld,
                0,
                cpno,
                OverlapCheck()
            );

            double radius_rwk = radius_rwcu + foilcu_t/2. + foilkapton_t/2.;
            position = G4ThreeVector(
                radius_rwk*std::cos(phi),
                radius_rwk*std::sin(phi),
                posz
            );

            transform = G4Transform3D(rotm, position);
            new G4PVPlacement(
                transform,
                rwellkapton_logic,
                "rwellkapton_stave",
                logicWorld,
                0,
                cpno,
                OverlapCheck()
            );

            double radius_nomex = radius_rwk + foilkapton_t/2. + readoutnomex_t/2.;
            position = G4ThreeVector(
                radius_nomex*std::cos(phi),
                radius_nomex*std::sin(phi),
                posz
            );

            transform = G4Transform3D(rotm, position);
            new G4PVPlacement(
                transform,
                nomex_logic,
                "nomex_stave",
                logicWorld,
                0,
                cpno,
                OverlapCheck()
            );

            double radius_roe = radius_nomex + readoutnomex_t/2. + readoutelectrode_t/2.;
            position = G4ThreeVector(
                radius_roe*std::cos(phi),
                radius_roe*std::sin(phi),
                posz
            );
            
            transform = G4Transform3D(rotm, position);
            new G4PVPlacement(
                transform,
                readoutelectrode_logic,
                "readoutelectrode_stave",
                logicWorld,
                0,
                cpno,
                OverlapCheck()
            );

            double radius_rok = radius_roe + readoutelectrode_t/2. + readoutkapton_t/2.;
            position = G4ThreeVector(
                radius_rok*std::cos(phi),
                radius_rok*std::sin(phi),
                posz
            );

            transform = G4Transform3D(rotm, position);
            new G4PVPlacement(
                transform,
                readoutkapton_logic,
                "readoutkapton_stave",
                logicWorld,
                0,
                cpno,
                OverlapCheck()
            );

            double radius_pcb = radius_rok + readoutkapton_t/2. + pcb_t/2.;
            position = G4ThreeVector(
                radius_pcb*std::cos(phi),
                radius_pcb*std::sin(phi),
                posz
            );

            transform = G4Transform3D(rotm, position);
            new G4PVPlacement(
                transform,
                pcb_logic,
                "pcb_stave",
                logicWorld,
                0,
                cpno,
                OverlapCheck()
            );

        }
    }


    return ;
}

void ePIC_OuterMPGD_Detector::Print(const string& what)const{
    std::cout<<"ePIC Outer MPGD Detector: "<<std::endl;
	if(what == "ALL" || what == "VOLUME"){
		std::cout<<"Version 0.1"<<std::endl;
	}
	return;
}