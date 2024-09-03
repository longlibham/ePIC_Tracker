
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
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4NistManager.hh>
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

	double lec_length = m_Params->get_double_param("lec_length") * cm;
	double rec_length = m_Params->get_double_param("rec_length") * cm;
	double anc_length = m_Params->get_double_param("anc_length") * cm;
	double anc_thickness = m_Params->get_double_param("anc_thickness") * cm;
	double las_airspace = m_Params->get_double_param("las_airspace") * cm;
	double kapton_thickness = m_Params->get_double_param("kapton_thickness") * cm;

	double peri_width = m_Params->get_double_param("periphery_width") * cm;
	

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
	// color settings 
	G4Colour col_si_insensitive = G4Colour(185./255., 99./255., 83./255., 1.);
	G4Colour col_air = G4Colour(126./255., 134./255., 134./255., 0.);
	G4Colour col_kapton = G4Colour(45./255., 179./255., 7./255., .7);

	G4NistManager* nist = G4NistManager::Instance();
	G4Material* mat_cf = nist->FindOrBuildMaterial("CF");
	G4VSolid* c_stave = new G4Box("SVT_OB_stave", c_w/2.0, c_t/2.0, c_l/2.0);
	G4LogicalVolume* stave_logical = new G4LogicalVolume(c_stave, mat_cf, "StaveLogical");
	G4VisAttributes* c_vis = new G4VisAttributes(G4Color(G4Colour::Grey()));
	c_vis->SetForceSolid(true);
	stave_logical->SetVisAttributes(c_vis);

	// build the air box to contain the las
	G4Material* mat_air = nist->FindOrBuildMaterial("G4_AIR");
	double airbox_length = anc_length + las_airspace + lec_length + si_l + rec_length;
	double airbox_thickness = anc_thickness + kapton_thickness;
	G4VSolid* airbox = new G4Box("airbox", si_w/2., airbox_thickness/2., airbox_length/2.);
	G4LogicalVolume* airbox_logic = new G4LogicalVolume(airbox, mat_air, "AirBoxLogic");
	G4VisAttributes* air_vis = new G4VisAttributes(col_air);
	air_vis->SetForceSolid(1);
	airbox_logic->SetVisAttributes(air_vis);
	
	// build LAS
	G4Material* mat_si = nist->FindOrBuildMaterial("G4_Si");
	G4VSolid* las_box = new G4Box("LAS_box", si_w/2.0, si_t/2.0, si_l/2.0);
	G4VSolid* sub_box = new G4Box("sub_box", (peri_width + m_nooverlap*cm)/2., 10., (si_l + m_nooverlap*cm)/2. );

	G4VSolid* sen_box = new G4SubtractionSolid("LAS_sub0", las_box, sub_box, 0,  
	G4ThreeVector(-(si_w - peri_width)/2., 0., 0.));
	sen_box = new G4SubtractionSolid("LAS_sub1", sen_box, sub_box, 0,
	G4ThreeVector(-peri_width/2., 0., 0.));
	sen_box = new G4SubtractionSolid("LAS_sub2", sen_box, sub_box, 0,  
	G4ThreeVector(peri_width/2., 0., 0.));
	sen_box = new G4SubtractionSolid("LAS_sensitive", sen_box, sub_box, 0, 
	G4ThreeVector((si_w - peri_width)/2., 0., 0.));

	G4LogicalVolume* las_logical = new G4LogicalVolume(sen_box, mat_si, "LASLogial");
	G4VisAttributes* las_vis = new G4VisAttributes(G4Color(G4Colour::Yellow()));
	las_vis->SetForceSolid(true);
	las_logical->SetVisAttributes(las_vis);

	G4LogicalVolume* las_sensor_logic = new G4LogicalVolume(las_box, mat_si, "LASensorLogic");
	las_sensor_logic->SetVisAttributes(las_vis);

	G4VPhysicalVolume* las_phys = 
	new G4PVPlacement(
		nullptr,
		G4ThreeVector(0, 0, -(airbox_length/2. - anc_length - las_airspace - lec_length - si_l/2.)),
		las_logical,
		"las_sensor",
		airbox_logic,
		0,
		0,
		OverlapCheck()
	);
	m_PhysicalVolumesSet.insert(las_phys);

	G4VisAttributes* insen_vis = new G4VisAttributes(col_si_insensitive);
	insen_vis->SetForceSolid(1);

	// lec 
	G4VSolid* lec_box = new G4Box("lec_box", si_w/2., si_t/2., lec_length/2.);
	G4LogicalVolume* lec_logic = new G4LogicalVolume(lec_box, mat_si, "LECLogic");
	lec_logic->SetVisAttributes(insen_vis);

	new G4PVPlacement(
		nullptr,
		G4ThreeVector(0, 0, -(airbox_length/2. - anc_length - las_airspace - lec_length/2.)),
		lec_logic,
		"lec_si",
		airbox_logic,
		0,
		0,
		OverlapCheck()
	);

	
	G4VSolid* rec_box = new G4Box("rec_box", si_w/2., si_t/2., rec_length/2.);
	G4LogicalVolume* rec_logic = new G4LogicalVolume(rec_box, mat_si, "RECLogic");
	rec_logic->SetVisAttributes(insen_vis);

	new G4PVPlacement(
		nullptr,
		G4ThreeVector(0, 0, (airbox_length - rec_length)/2.),
		rec_logic,
		"rec_si",
		airbox_logic,
		0,
		0,
		OverlapCheck()
	);

	// ancASIC
	G4VSolid* anc_box = new G4Box("anc_box", anc_length/2., anc_thickness/2., anc_length/2.);
	G4LogicalVolume* anc_logic = new G4LogicalVolume(anc_box, mat_si, "AncLogic");
	anc_logic->SetVisAttributes(insen_vis);

	new G4PVPlacement(
		nullptr,
		G4ThreeVector(-si_w/4., -(airbox_thickness - kapton_thickness)/2., -(airbox_length/2. - anc_length/2.)),
		anc_logic,
		"ancASIC",
		airbox_logic,
		0,
		0,
		OverlapCheck()
	);

	new G4PVPlacement(
		nullptr,
		G4ThreeVector(si_w/4, -(airbox_thickness - kapton_thickness)/2., -(airbox_length/2. - anc_length/2.)),
		anc_logic,
		"ancASIC",
		airbox_logic,
		0,
		1,
		OverlapCheck()
	);
	
	G4VSolid* periphery_box = new G4Box("periphery_box", peri_width/2., si_t/2., si_l/2. );
	G4LogicalVolume* periphery_logic = new G4LogicalVolume(periphery_box, mat_si, "Periphery_logic");
	periphery_logic->SetVisAttributes(insen_vis);

	new G4PVPlacement(
		nullptr,
		G4ThreeVector(-(si_w - peri_width)/2., 0, -(airbox_length/2. - anc_length - las_airspace - lec_length - si_l/2.)),
		periphery_logic,
		"periphery",
		airbox_logic,
		0,
		0,
		OverlapCheck()
	);

	new G4PVPlacement(
		nullptr,
		G4ThreeVector(-peri_width/2., 0, -(airbox_length/2. - anc_length - las_airspace - lec_length - si_l/2.)),
		periphery_logic,
		"periphery",
		airbox_logic,
		0,
		1,
		OverlapCheck()
	);

	new G4PVPlacement(
		nullptr,
		G4ThreeVector(peri_width/2., 0, -(airbox_length/2. - anc_length - las_airspace - lec_length - si_l/2.)),
		periphery_logic,
		"periphery",
		airbox_logic,
		0,
		2,
		OverlapCheck()
	);

	new G4PVPlacement(
		nullptr,
		G4ThreeVector((si_w - peri_width)/2., 0, -(airbox_length/2. - anc_length - las_airspace - lec_length - si_l/2.)),
		periphery_logic,
		"periphery",
		airbox_logic,
		0,
		3,
		OverlapCheck()
	);

	//kapton
	G4Material* mat_kapton = nist->FindOrBuildMaterial("G4_KAPTON");
	G4VSolid* kapton_box = new G4Box("kapton_layer", si_w/2., si_t/2., (las_airspace + anc_length)/2.);
	G4LogicalVolume* kapton_logic = new G4LogicalVolume(kapton_box, mat_kapton, "KaptonLogic");
	G4VisAttributes* kapton_vis = new G4VisAttributes(col_kapton);
	kapton_vis->SetForceSolid(1);
	kapton_logic->SetVisAttributes(kapton_vis);

	new G4PVPlacement(
		nullptr,
		G4ThreeVector(0, (airbox_thickness - kapton_thickness)/2., -(airbox_length - las_airspace - anc_length)/2.),
		kapton_logic,
		"kapton_layer",
		airbox_logic,
		0,
		0,
		OverlapCheck()
	);

	// placement of staves
	if (r_inner + si_t + c_t >= r_outer){
		std::cout<<"Overlaps will happen between inner and outer layers, please check parameters blew: "<<std::endl;
		std::cout<<"r_inner: "<< r_inner <<std::endl;
		std::cout<<"silicon thickness: "<< si_t << std::endl;
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
				r_eff = r_inner + c_t/2. + airbox_thickness/2.;
			} 
			else if(i%2 == 0 && j%2 == 1){
				r_eff = r_outer + c_t/2. + airbox_thickness/2.;
			}
			else if(i%2 == 1 && j%2 == 0){
				r_eff = r_inner - c_t/2. - airbox_thickness/2.;
			}
			else{
				r_eff = r_outer - c_t/2. - airbox_thickness/2.;
			}
			double px = r_eff * std::cos(phi);
			double py = r_eff * std::sin(phi);

			G4ThreeVector position = G4ThreeVector(px, py, pz);
			G4Transform3D transform = G4Transform3D(rotm, position);

			G4VPhysicalVolume* las_sensor_phys = 
			new G4PVPlacement(
						transform,
						airbox_logic,
						"LAS_stave",
						logicWorld,
						0,
						i*nphi +j,
						OverlapCheck()
					);

			m_PhysicalVolumesSet.insert(las_sensor_phys);

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
