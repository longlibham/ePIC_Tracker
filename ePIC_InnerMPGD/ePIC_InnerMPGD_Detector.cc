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
    G4Material* CO2 = man->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
    G4Material* Ar = man->FindOrBuildMaterial("G4_Ar");
    double density = 1.714e-3*g/cm3;
    G4Material* MPGDGAS = new G4Material("MPGDGAS", density, ncomponents = 2);
    MPGDGAS->AddMaterial(CO2, 30*perCent);
    MPGDGAS->AddMaterial(Ar, 70*perCent);


    
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

    G4VSolid* coat_solid = new G4Tubs("ePIC_MPGD_coat", radius - MPGD_t/2., radius + MPGD_t/2., (z_max - z_min)/2., 0., 2*M_PI);
    G4LogicalVolume* coat_logic = new G4LogicalVolume(coat_solid, mat_mpgd_gas, "coatLogic");

    new G4PVPlacement(nullptr,  //no rotation
                    G4ThreeVector(),    //at (0, 0, 0)
                    coat_logic,         //logical volume
                    "GAS_COAT",         //name
                    logicWorld,         // mother volume
                    false,              //no boolean operation
                    0,                  // copy number
                    fCheckOverlaps
    );

    // Detector construction
    G4VSolid* DriftCuGround_solid = new G4Box("DriftCuGround_plane", );

}