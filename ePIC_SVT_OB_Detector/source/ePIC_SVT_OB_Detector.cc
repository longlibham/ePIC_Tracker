
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
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4Transform3D.hh>
#include <Geant4/G4VisAttributes.hh>
#include <phparameter/PHParameters.h>

#include <cmath>
#include <iostream>

#include <TSystem.h>

class G4VSolid;
class PHCompositeNode;

using namespace std;

ePIC_SVT_OB_Detector::ePIC_SVT_OB_Detector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string &dname, const int lyr)
	:PHG4Detector(subsys, Node, dname)
	, m_Params(parameters)
	, m_Layer(lyr){

}

int ePIC_SVT_OB_Detector::IsInDetector(G4VPhysicalVolume* volume) const{
	set<G4VPhysicalVolume*>::const_iterator iter = m_PhysicalVolumesSet.find(volume);
	if(iter != m_PhysicalVolumesSet.end()){
		return 1;
	}
	return 0;
}

void ePIC_SVT_OB_Detector::ConstructMe(G4LogicalVolume* logicWorld){
	// get essential parameter passing from PHParameters
	double r_inner = m_Params->get_double_param("r_inner") * cm;
	double r_outer = m_Params->get_double_param("r_outer") * cm;
	// carbon support
	double c_t = m_Params->get_double_param("carbon_thickness") * cm;
	double c_l = m_Params->get_double_param("carbon_length") * cm;
	double c_w = m_Params->get_double_param("carbon_width") * cm;
	// silicon geo parameters
	double si_t = m_Params->get_double_param("si_thickness") * cm;
	double si_l = m_Params->get_double_param("si_length") * cm;
	double si_w = m_Params->get_double_param("si_width") * cm;
	double si_c_gap = m_Params->get_double_param("si_carbon_gap") * cm;
	int  nz = m_Params->get_int_param("n_silicon_z");
	int  nphi = m_Params->get_int_param("n_stave_phi");
	
	// overlaps of LAS
	double las_ol = m_Params->get_double_param("las_overlap") * cm;

	double las_tot = si_l*nz - las_ol*(nz-1);
	double z_align = (c_l - las_tot)/2.;

	if(!std::isfinite(r_inner) || !std::isfinite(r_outer) || !std::isfinite(c_t)||
			!std::isfinite(c_l) || !std::isfinite(c_w) || !std::isfinite(si_t) || 
			!std::isfinite(si_l) || !std::isfinite(si_w) || !std::isfinite(nz) || 
			!std::isfinite(nphi)){
		std::cout<<"PHWHERE "<< ": Bad Parameters for "<< GetName() << std::endl;
		std::cout<<"r_inner: "<< r_inner <<" r_outer: "<< r_outer << std::endl;
		std::cout<<"carbon_thickness: "<< c_t <<" length: "<< c_l <<" width: "<< c_w << std::endl;
		std::cout<<"silicon_thickness: "<< si_t <<" length: "<< si_l <<" width: "<< si_w << std::endl;
		std::cout<<"number of LAS in Z: "<< nz << " number of stave in phi: "<< nphi << std::endl;
		gSystem->Exit(1);
	}

	// build the carbon stave
	G4VSolid* c_stave = new G4Box("SVT_OB_stave", c_w/2.0, c_t/2.0, c_l/2.0);
	G4LogicalVolume* stave_logical = new G4LogicalVolume(c_stave, G4Material::GetMaterial("G4_C"), "StaveLogical");
	G4VisAttributes* c_vis = new G4VisAttributes(G4Color(G4Colour::Grey()));
	c_vis->SetForceSolid(true);
	stave_logical->SetVisAttributes(c_vis);
	
	// build LAS
	G4VSolid* las_box = new G4Box("LAS_box", si_w/2.0, si_t/2.0, si_l/2.0);
	G4LogicalVolume* las_logical = new G4LogicalVolume(las_box, G4Material::GetMaterial("G4_Si"), "LASLogial");
	G4VisAttributes* las_vis = new G4VisAttributes(G4Color(G4Colour::Yellow()));
	las_vis->SetForceSolid(true);
	las_logical->SetVisAttributes(las_vis);

	// placement of staves
	if (r_inner + si_t + si_c_gap + c_t >= r_outer){
		std::cout<<"Overlaps will happen between inner and outer layers, please check parameters blew: "<<std::endl;
		std::cout<<"r_inner: "<< r_inner <<std::endl;
		std::cout<<"silicon thickness: "<< si_t << std::endl;
		std::cout<<"Gap between si and carbon: "<< si_c_gap << std::endl;
		std::cout<<"Carbon thickness: "<< c_t <<std::endl;
		gSystem->Exit(1);
	}

	for(int i =0; i < nz; i++){ // loop LAS' in z direction
		if(i == 0){   // place the carbon stave first
			for(int j = 0; j < nphi; j++){
				double phi = j*2*M_PI/nphi;
				G4RotationMatrix rotm = G4RotationMatrix();
				rotm.rotateZ(M_PI/2. + phi);
				G4ThreeVector uz = G4ThreeVector(std::cos(phi), std::sin(phi), 0.);
				G4ThreeVector position = G4ThreeVector(0., 0., 0.);
				if(j%2 == 0){
					position = r_inner*uz;
				}
				else{
					position = r_outer*uz;
				}
				G4Transform3D transform = G4Transform3D(rotm, position);

				
//				G4VPhysicalVolume* c_phy = 
				new G4PVPlacement(
					transform, // rotation, position
					stave_logical,  // logical volume
					"Carbon_stave",
					logicWorld,
					0,   // no boolean operation
					j,   // copy number
					OverlapCheck()
					
					);

			}
		}
		// place the LAS 
		for(int j =0; j < nphi; j++){
			double phi = j*2*M_PI/nphi;
			G4RotationMatrix rotm = G4RotationMatrix();
			rotm.rotateZ(M_PI/2. + phi);
			double pz = -c_l/2 + z_align + (0.5+i) * si_l - i * las_ol;
			double r_eff = 0.;
			if(i%2 == 0 && j%2 == 0){
				r_eff = r_inner + c_t + si_c_gap;
			} 
			else if(i%2 == 0 && j%2 == 1){
				r_eff = r_outer + c_t + si_c_gap;
			}
			else if(i%2 == 1 && j%2 == 0){
				r_eff = r_inner - si_c_gap;
			}
			else{
				r_eff = r_outer - si_c_gap;
			}
			double px = r_eff * std::cos(phi);
			double py = r_eff * std::sin(phi);

			G4ThreeVector position = G4ThreeVector(px, py, pz);
			G4Transform3D transform = G4Transform3D(rotm, position);

			G4VPhysicalVolume* si_phy = 
			new G4PVPlacement(
						transform,
						las_logical,
						"LAS",
						logicWorld,
						0,
						i*nphi +j,
						OverlapCheck()
					);
			m_PhysicalVolumesSet.insert(si_phy);

			cout<<"Process done copy NO."<<i*nphi+j<<endl;
		}

	}

	return;

}

void ePIC_SVT_OB_Detector::Print(const std::string &what) const{
	std::cout<<"ePIC SVT OB Detector: "<<std::endl;
	if(what == "ALL" || what == "VOLUME"){
		std::cout<<"Version 0.1"<<std::endl;
	}
	return;
}
