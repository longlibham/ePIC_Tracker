
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
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Color.hh>
#include <Geant4/G4UnionSolid.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4IntersectionSolid.hh>
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
		DefineMaterials();
}

int ePIC_SVT_OB_Detector::IsInDetector(G4VPhysicalVolume* volume) const{
	set<G4VPhysicalVolume*>::const_iterator iter = m_PhysicalVolumesSet.find(volume);
	if(iter != m_PhysicalVolumesSet.end()){
		return 1;
	}
	return 0;
}

void ePIC_SVT_OB_Detector::DefineMaterials(){
	// build the k9 foam material
	G4NistManager* nist = G4NistManager::Instance();
	G4Material* mat_air = nist->FindOrBuildMaterial("G4_AIR");
	G4Material* mat_c = nist->FindOrBuildMaterial("G4_C");

	double density = 0.11*g/cm3;
	G4Material* mat_k9 = new G4Material("K9Foam", density, 2);
	mat_k9->AddMaterial(mat_air, 97*perCent);
	mat_k9->AddMaterial(mat_c, 3*perCent);
	cout<<*(G4Material::GetMaterialTable())<<endl;

}

void ePIC_SVT_OB_Detector::ConstructMe(G4LogicalVolume* logicWorld){
	// get essential parameter passing from PHParameters
	double r_inner = m_Params->get_double_param("r_inner") * cm;
	double r_outer = m_Params->get_double_param("r_outer") * cm;
	// carbon support
	double c_t = m_Params->get_double_param("carbon_thickness") * cm;
	double c_l = m_Params->get_double_param("carbon_length") * cm;
	double c_w = m_Params->get_double_param("carbon_width") * cm;
	double cf_curve_radius = m_Params->get_double_param("cf_curve_radius") * cm;
	double cf_center_height = m_Params->get_double_param("cf_center_height") * cm;
	double cf_edge_height = m_Params->get_double_param("cf_edge_height") * cm;
	double cf_margin = m_Params->get_double_param("cf_margin") * cm;
	double cf_lec = m_Params->get_double_param("cf_leftendcap") * cm;
	double cf_rec = m_Params->get_double_param("cf_rightendcap") * cm;

	// o-rings for cooling
	double oring_r = m_Params->get_double_param("oring_radius") * cm;
	double oring_spacing = m_Params->get_double_param("oring_spacing")*cm;
	
	// c-support cf
	double cf_csupport_width = m_Params->get_double_param("cf_csupport_width")*cm;

	//for K9 foam
	double k9_center_height = m_Params->get_double_param("k9_center_height")*cm; 

	// silicon geo parameters
	double si_t = m_Params->get_double_param("si_thickness") * cm;
	double si_w = m_Params->get_double_param("si_width") * cm;
	double matrix_length = m_Params->get_double_param("matrix_length") * cm;
	double switch_length = m_Params->get_double_param("switch_length") * cm;
	double backbone_length = m_Params->get_double_param("backbone_length") * cm;
	int ntile = m_Params->get_int_param("ntile");
	int nmatrix = m_Params->get_int_param("nmatrix");


	double lec_length = m_Params->get_double_param("lec_length") * cm;
	double rec_length = m_Params->get_double_param("rec_length") * cm;
	double anc_length = m_Params->get_double_param("anc_length") * cm;
	double anc_thickness = m_Params->get_double_param("anc_thickness") * cm;
	double las_airspace = m_Params->get_double_param("las_airspace") * cm;
	double kapton_thickness = m_Params->get_double_param("kapton_thickness") * cm;

	double peri_width = m_Params->get_double_param("periphery_width") * cm;
	

	int  nz = m_Params->get_int_param("n_silicon_z");
	int  nphi = m_Params->get_int_param("n_stave_phi");
	
	double las_active_length = ntile*(backbone_length + nmatrix*(matrix_length + switch_length));


	if(!std::isfinite(r_inner) || !std::isfinite(r_outer) || !std::isfinite(c_t)||
			!std::isfinite(c_l) || !std::isfinite(c_w) || !std::isfinite(si_t) || 
			!std::isfinite(matrix_length) || !std::isfinite(switch_length) || !std::isfinite(backbone_length) ||
			!std::isfinite(nmatrix) || !std::isfinite(ntile) || !std::isfinite(si_w) || !std::isfinite(nz) || 
			!std::isfinite(nphi)){
		std::cout<<"PHWHERE "<< ": Bad Parameters for "<< GetName() << std::endl;
		std::cout<<"r_inner: "<< r_inner <<" r_outer: "<< r_outer << std::endl;
		std::cout<<"carbon_thickness: "<< c_t <<" length: "<< c_l <<" width: "<< c_w << std::endl;
		std::cout<<"silicon_thickness: "<< si_t <<" length: "<< matrix_length << " switch length: "<<switch_length<<endl;
		std::cout<<"backbone length: "<<backbone_length<<" nmatrix: "<<nmatrix<<" ntile: "<<ntile<<" width: "<< si_w << std::endl;
		std::cout<<"number of LAS in Z: "<< nz << " number of stave in phi: "<< nphi << std::endl;
		gSystem->Exit(1);
	}

	// build the carbon stave
	// color settings 
	G4Colour col_si_insensitive = G4Colour(185./255., 99./255., 83./255., 1.);
	G4Colour col_air = G4Colour(126./255., 134./255., 134./255., 0.);
	G4Colour col_kapton = G4Colour(45./255., 179./255., 7./255., .7);
	G4Colour col_k9foam = G4Colour(137./255., 248./255., 173./255., 1.);

	G4NistManager* nist = G4NistManager::Instance();
	G4Material* mat_cf = nist->FindOrBuildMaterial("CF");
	G4Material* mat_air = nist->FindOrBuildMaterial("G4_AIR");
	G4Material* mat_kapton = nist->FindOrBuildMaterial("G4_KAPTON");
	G4Material* mat_si = nist->FindOrBuildMaterial("G4_Si");
	G4Material* mat_k9 = nist->FindOrBuildMaterial("K9Foam");


	// G4VisAttributes
	G4VisAttributes* c_vis = new G4VisAttributes(G4Color(G4Colour::Grey()));
	c_vis->SetForceSolid(true);
	G4VisAttributes* air_vis = new G4VisAttributes(col_air);
	air_vis->SetForceSolid(1);
	G4VisAttributes* kapton_vis = new G4VisAttributes(col_kapton);
	kapton_vis->SetForceSolid(1);
	G4VisAttributes* k9foam_vis = new G4VisAttributes(col_k9foam);
	k9foam_vis->SetForceSolid(1);
	G4VisAttributes* las_vis = new G4VisAttributes(G4Color(G4Colour::Yellow()));
	las_vis->SetForceSolid(true);
	G4VisAttributes* insen_vis = new G4VisAttributes(col_si_insensitive);
	insen_vis->SetForceSolid(1);



	double tub_rinner = 0.;
	double tub_router = cf_curve_radius;
	double tub_length = c_l + cf_lec + cf_rec;
	double stave_length = tub_length;
	double sta_phi = 0;
	double sto_phi = 2*M_PI;

	G4VSolid* cf_tub1 = new G4Tubs("carbon_tub1", tub_rinner, tub_router, tub_length/2., sta_phi, sto_phi);
	G4VSolid* cf_tub2 = new G4Tubs("carbon_tub2", tub_rinner, tub_router, tub_length/2., sta_phi, sto_phi);

	G4VSolid* cf_stave = nullptr;
	cf_stave = new G4IntersectionSolid("cf_stave_curved", cf_tub1, cf_tub2, nullptr, G4ThreeVector(0, 2*cf_curve_radius - cf_center_height, 0));

	// polish the edge
	double intersection_chord_length = 2*sqrt(
		cf_curve_radius*cf_curve_radius -
		(cf_curve_radius - cf_center_height/2.)*(cf_curve_radius - cf_center_height/2.)
	);
	double strip_x = cf_edge_height/2.*(cf_curve_radius - cf_center_height/2.)*2/intersection_chord_length;
	strip_x = intersection_chord_length/2. - si_w/2. - cf_margin/2.;


	double box_x = 2*strip_x;
	double box_y = 100.;
	double box_z = stave_length + 10*m_nooverlap;
	G4VSolid* cf_strip_box = new G4Box("strip_box", box_x/2., box_y/2., box_z/2.);
	cf_stave = new G4SubtractionSolid("cf_stave_strip_left", cf_stave, cf_strip_box, nullptr, 
		G4ThreeVector(-intersection_chord_length/2., cf_curve_radius - cf_center_height/2., 0));
	cf_stave = new G4SubtractionSolid("cf_stave_strip_right", cf_stave, cf_strip_box, nullptr,
		G4ThreeVector(intersection_chord_length/2., cf_curve_radius - cf_center_height/2., 0));
		

	// stave container
	tub_length += m_nooverlap;
	cf_tub1 = new G4Tubs("carbon_tub1", tub_rinner, tub_router, tub_length/2., sta_phi, sto_phi);
	cf_tub2 = new G4Tubs("carbon_tub2", tub_rinner, tub_router, tub_length/2., sta_phi, sto_phi);
	G4VSolid* cf_stave_container = new G4IntersectionSolid("cf_stave_container", cf_tub1, cf_tub2, nullptr,
		G4ThreeVector(0, 2*cf_curve_radius - cf_center_height - 2*anc_thickness - 2*kapton_thickness - 2*c_t, 0.));
	

	G4LogicalVolume* cf_stave_container_logic = new G4LogicalVolume(cf_stave_container, mat_air, "StaveContainerLogic");
	cf_stave_container_logic->SetVisAttributes(air_vis);


	// // strip the center
	box_x = si_w;//intersection_chord_length - 2*strip_x - cf_margin;
	box_y = 1. * cm;
	box_z = stave_length - 2*c_t;

	cf_strip_box = new G4Box("strip_box_center", box_x/2., box_y/2., box_z/2.);
	cf_stave = new G4SubtractionSolid("cf_stave_strip_center", cf_stave, cf_strip_box, nullptr, 
		G4ThreeVector(0, cf_curve_radius - cf_center_height/2., 0.));



	// o-rings
	G4VSolid* cf_strip_oring = new G4Tubs("strip_oring", 0., oring_r, c_l + 2*c_t + m_nooverlap, 0., 2*M_PI);
	cf_stave = new G4SubtractionSolid("cf_stave_strip_oring_left", cf_stave, cf_strip_oring, nullptr,
		G4ThreeVector(-(oring_spacing+oring_r), cf_curve_radius - cf_center_height/2., 0.));//-stave_length/2.));
	cf_stave = new G4SubtractionSolid("cf_stave_strip_oring_right", cf_stave, cf_strip_oring, nullptr,
		G4ThreeVector((oring_spacing+oring_r), cf_curve_radius - cf_center_height/2., 0.));//-stave_length/2.));
	
	//C-shape support structure
	box_x = cf_csupport_width;
	box_y = k9_center_height + 2*c_t;
	G4VSolid* cf_csupport_box = new G4Box("cf_csupport_box", box_x/2., box_y/2., box_z/2.);
	
	box_y = box_y - 2*c_t;
	box_z = box_z + 2*c_t;
	G4VSolid* cf_csupport_sub = new G4Box("cf_csupport_sub", box_x/2., box_y/2., box_z/2.);
	
	G4VSolid* cf_csupport_solid = nullptr;
	cf_csupport_solid = new G4SubtractionSolid("cf_csupport_subleft", cf_csupport_box, cf_csupport_sub, nullptr,
		G4ThreeVector(-(c_t+box_x)/2., 0., 0.));
	cf_csupport_solid = new G4SubtractionSolid("cf_csupport_subright", cf_csupport_solid, cf_csupport_sub, nullptr,
		G4ThreeVector((c_t+box_x)/2., 0., 0.));

	G4LogicalVolume* cf_csupport_logic = new G4LogicalVolume(cf_csupport_solid, mat_cf, "CSupportLogic");
	cf_csupport_logic->SetVisAttributes(c_vis);

	
	
	//shell surface
	tub_rinner = cf_curve_radius;
	tub_router = cf_curve_radius + c_t;
	tub_length = stave_length;
	sta_phi = 0;
	sto_phi =si_w/cf_curve_radius;

	G4VSolid* cf_shell_top = new G4Tubs("cf_shell_top", tub_rinner, tub_router, tub_length/2., sta_phi, sto_phi);
	G4VSolid* cf_shell_bot = new G4Tubs("cf_shell_bot", tub_rinner, tub_router, tub_length/2., sta_phi, sto_phi);

	box_x = 10.*cm;
	box_y = si_w*2.+ m_nooverlap;
	box_z = lec_length + las_active_length + rec_length;

	cf_strip_box = new G4Box("strip_box", box_x/2., box_y/2., box_z/2.);
	double posx = 0.;
	double posy = 0.;
	double posz = 0.;

	G4RotationMatrix rotm = G4RotationMatrix();
	
	G4Transform3D trans; //G4Transform3D(rotm, G4ThreeVector(0, cf_center_height, 0));


	for(int j = 0; j<nz/2; j++){  // n in z-axis

		posy = si_w/2.;
		
		posz = -stave_length/2. + cf_lec + j*(box_z + las_active_length) + box_z/2.;
		posx = cf_curve_radius;

		cf_shell_top = new G4SubtractionSolid("cf_shell_top", cf_shell_top, cf_strip_box, nullptr,
			G4ThreeVector(posx, posy, posz));

		posz = -stave_length/2. + cf_lec + (j+1)*las_active_length + j*box_z + box_z/2.;
		cf_shell_bot = new G4SubtractionSolid("cf_shell_bot", cf_shell_bot, cf_strip_box, nullptr,
			G4ThreeVector(posx, posy, posz)); 
	}


	G4LogicalVolume* cf_stave_logic = new G4LogicalVolume(cf_stave, mat_cf, "CFLogical");
	cf_stave_logic->SetVisAttributes(c_vis);

	new G4PVPlacement(
		nullptr,
		G4ThreeVector(0, -anc_thickness - kapton_thickness - c_t, 0),
		cf_stave_logic,
		"CF_stave",
		cf_stave_container_logic,
		0,
		0,
		OverlapCheck()
	);

	G4LogicalVolume* cf_shell_top_logic = new G4LogicalVolume(cf_shell_top, mat_cf, "CF_ShellTOPLogic");
	cf_shell_top_logic->SetVisAttributes(c_vis);

	rotm.rotateZ(M_PI/2. - si_w/(2*cf_curve_radius));
	trans = G4Transform3D(rotm, G4ThreeVector(0., -anc_thickness -kapton_thickness, 0.));
	new G4PVPlacement(
		trans,
		cf_shell_top_logic,
		"CF_ShellTOP",
		cf_stave_container_logic,
		0,
		0,
		0
	);

	G4LogicalVolume* cf_shell_bot_logic = new G4LogicalVolume(cf_shell_bot, mat_cf, "CF_ShellBOTLogic");
	cf_shell_bot_logic->SetVisAttributes(c_vis);

	rotm.rotateZ(M_PI);
	double cf_stave_container_ymin = cf_curve_radius - cf_center_height -2*anc_thickness -2*kapton_thickness; 
	trans = G4Transform3D(rotm, G4ThreeVector(0., cf_stave_container_ymin + cf_curve_radius + anc_thickness + kapton_thickness, 0.));
	new G4PVPlacement(
		trans,
		cf_shell_bot_logic,
		"CF_ShellBOT",
		cf_stave_container_logic,
		0,
		0,
		0
	);

	//placement of the c-shape support
	new G4PVPlacement(
		nullptr,
		G4ThreeVector(0., cf_curve_radius -anc_thickness -kapton_thickness - cf_center_height/2., (cf_lec - cf_rec)/2.),
		cf_csupport_logic,
		"CSupport_stave",
		cf_stave_container_logic,
		0,
		0,
		OverlapCheck()
	);

	// for the K9 foam
	tub_rinner = 0.;
	tub_router = cf_curve_radius;
	tub_length = lec_length + rec_length;
	sta_phi = 0.;
	sto_phi = 2*M_PI;
	G4VSolid* k9tub1 = new G4Tubs("k9tub1", tub_rinner, tub_router, tub_length/2., sta_phi, sto_phi);
	G4VSolid* k9tub2 = new G4Tubs("k9tub2", tub_rinner, tub_router, tub_length/2., sta_phi, sto_phi);

	G4VSolid* k9_curved_solid = new G4IntersectionSolid("k9_curved_solid", k9tub1, k9tub2, nullptr,
		G4ThreeVector(0, 2*tub_router - k9_center_height, 0.));
	
	box_x = 10.*cm + m_nooverlap;
	box_y = 10.*cm;
	box_z = lec_length+rec_length+m_nooverlap;

	cf_strip_box = new G4Box("strip_box", box_x/2., box_y/2., box_z/2.);
	k9_curved_solid = new G4SubtractionSolid("k9_curved_strip_left", k9_curved_solid, cf_strip_box, nullptr,
		G4ThreeVector(-(si_w + box_x - cf_margin - m_nooverlap)/2., cf_curve_radius-k9_center_height/2., 0.));
	k9_curved_solid = new G4SubtractionSolid("k9_curved_strip_right", k9_curved_solid, cf_strip_box, nullptr,
		G4ThreeVector((si_w + box_x - cf_margin - m_nooverlap)/2., cf_curve_radius-k9_center_height/2., 0.));

	box_x = 2*c_t + m_nooverlap;
	cf_strip_box = new G4Box("strip_box", box_x/2., box_y/2., box_z/2.);
	k9_curved_solid = new G4SubtractionSolid("k9_curved_strip_center", k9_curved_solid, cf_strip_box, nullptr,
		G4ThreeVector(0., cf_curve_radius-k9_center_height/2., 0.));

	G4LogicalVolume* k9foam_logic = new G4LogicalVolume(k9_curved_solid, mat_k9, "K9FoamLogic");	
	k9foam_logic->SetVisAttributes(k9foam_vis);

	for(int i = 0; i<nz/2; i++){
		posz = -stave_length/2. + cf_lec + i*(las_active_length + lec_length);
		posy = -anc_thickness -kapton_thickness -(cf_center_height - k9_center_height)/2.;
		new G4PVPlacement(
			nullptr,
			G4ThreeVector(0., posy, posz),
			k9foam_logic,
			"k9_foam_phys",
			cf_stave_container_logic,
			0,
			2*i,
			OverlapCheck()
		);

		posz = stave_length/2. - cf_rec - las_active_length - i*(las_active_length + lec_length);
		new G4PVPlacement(
			nullptr,
			G4ThreeVector(0., posy, posz),
			k9foam_logic,
			"k9_foam_phys",
			cf_stave_container_logic,
			0,
			2*i+1,
			OverlapCheck()

		);

	}


	// air box for LAS
	tub_rinner = cf_curve_radius - c_t/2.;
	tub_router = cf_curve_radius + c_t/2.;
	tub_length = lec_length + las_active_length + rec_length;
	double las_container_length = tub_length;
	sta_phi = 0.;
	sto_phi = si_w/cf_curve_radius;
	G4VSolid* las_container = new G4Tubs("LAS_container", tub_rinner, tub_router, tub_length/2., sta_phi, sto_phi);
	G4LogicalVolume* las_container_logic = new G4LogicalVolume(las_container, mat_air, "LASContainerLogic");
	las_container_logic->SetVisAttributes(air_vis);

	
	// build the LAS
	// 1. the tileContainer
	tub_rinner = cf_curve_radius - si_t/2.;
	tub_router = cf_curve_radius + si_t/2.;
	double tile_length = backbone_length + nmatrix*(matrix_length + switch_length);
	G4VSolid* tileContainer_solid = new G4Tubs("tileContainer_solid", tub_rinner, tub_router, tile_length/2., sta_phi, sto_phi);
	G4LogicalVolume* tileContainerLogic = new G4LogicalVolume(tileContainer_solid, mat_air, "tileContainerLogic");


	// 2. the mtxContainer
	tub_length = matrix_length;

	G4VSolid* mtxContainer_solid = new G4Tubs("mtxContainer_solid", tub_rinner, tub_router, tub_length/2., sta_phi, sto_phi);
	G4LogicalVolume* mtxContainerLogic = new G4LogicalVolume(mtxContainer_solid, mat_air, "mtxContainerLogic");

	//3. the pixel matrix
	G4VSolid* pixelmatrix_solid = new G4Tubs("pixelmatrix_solid", tub_rinner, tub_router, tub_length/2., sta_phi, sto_phi);

	
	tub_rinner = cf_curve_radius - 10*si_t;
	tub_router = cf_curve_radius + 10*si_t;
	tub_length = matrix_length + m_nooverlap;
	sto_phi = peri_width/cf_curve_radius;
	G4VSolid* strip_tub = new G4Tubs("strip_tub", tub_rinner, tub_router, tub_length/2., sta_phi, sto_phi);
	
	// for periphery 

	tub_rinner = cf_curve_radius - si_t/2.;
	tub_router = cf_curve_radius + si_t/2.;
	tub_length = matrix_length;
	G4VSolid* las_peri = new G4Tubs("LAS_periphery", tub_rinner, tub_router, tub_length/2., sta_phi, sto_phi);
	G4LogicalVolume* las_peri_logic = new G4LogicalVolume(las_peri, mat_si, "LASPeripheryLogic");
	las_peri_logic->SetVisAttributes(insen_vis);

	rotm = G4RotationMatrix();
	double rot_angle = 0.;
	rotm.rotateZ(rot_angle);
	G4ThreeVector pos = G4ThreeVector(0., 0., 0.);
	trans = G4Transform3D(rotm, pos);
	pixelmatrix_solid = new G4SubtractionSolid("LAS_Curve_strip_periphery1", pixelmatrix_solid, strip_tub, trans);
	
	new G4PVPlacement(
		trans,
		las_peri_logic,
		"LAS_Periphery1",
		mtxContainerLogic,
		0,
		0,
		OverlapCheck()
	);

	rot_angle = (si_w/2. - peri_width)/cf_curve_radius;
	rotm.rotateZ(rot_angle);
	trans = G4Transform3D(rotm, pos);
	pixelmatrix_solid = new G4SubtractionSolid("LAS_Curve_strip_periphery2", pixelmatrix_solid, strip_tub, trans);

	new G4PVPlacement(
		trans,
		las_peri_logic,
		"LAS_Periphery2",
		mtxContainerLogic,
		0,
		1,
		OverlapCheck()
	);

	rot_angle = peri_width/cf_curve_radius;
	rotm.rotateZ(rot_angle);
	trans = G4Transform3D(rotm,pos);
	pixelmatrix_solid = new G4SubtractionSolid("LAS_Curve_strip_periphery3", pixelmatrix_solid, strip_tub, trans);

	new G4PVPlacement(
		trans,
		las_peri_logic,
		"LAS_periphery3",
		mtxContainerLogic,
		0,
		2,
		OverlapCheck()
	);

	rot_angle = (si_w/2. - peri_width)/cf_curve_radius;
	rotm.rotateZ(rot_angle);
	trans = G4Transform3D(rotm, pos);
	pixelmatrix_solid = new G4SubtractionSolid("LAS_Curve_strip_periphery4", pixelmatrix_solid, strip_tub, trans);

	new G4PVPlacement(
		trans,
		las_peri_logic,
		"LAS_periphery4",
		mtxContainerLogic,
		0,
		3,
		OverlapCheck()
	);	

	

	G4LogicalVolume* pixelmatrixLogic = new G4LogicalVolume(pixelmatrix_solid, mat_si, "pixelmatrixLogic");
	pixelmatrixLogic->SetVisAttributes(las_vis);

	G4VPhysicalVolume* las_phys = 
	new G4PVPlacement(
		nullptr,
		pos,
		pixelmatrixLogic,
		"LAS_Stave",
		mtxContainerLogic,
		0,
		0,
		OverlapCheck()
	);

	m_PhysicalVolumesSet.insert(las_phys);


	// 4. the backbone 
	tub_length = backbone_length;
	sto_phi = si_w/cf_curve_radius;
	G4VSolid* backbone_solid = new G4Tubs("backbone_solid", tub_rinner, tub_router, tub_length/2., sta_phi, sto_phi);
	G4LogicalVolume* backboneLogic = new G4LogicalVolume(backbone_solid, mat_si, "backboneLogic");
	backboneLogic->SetVisAttributes(insen_vis);

	// 5. the power switch
	tub_length = switch_length;
	G4VSolid* switch_solid = new G4Tubs("switch_solid", tub_rinner, tub_router, tub_length/2., sta_phi, sto_phi);
	G4LogicalVolume* switchLogic = new G4LogicalVolume(switch_solid, mat_si, "switchLogic");
	switchLogic->SetVisAttributes(insen_vis);

	// place the backbone matrix and the power switch in the tileContainer
	double current_z = -tile_length/2. + backbone_length/2.;
	pos = G4ThreeVector(0., 0., current_z);
	new G4PVPlacement(
		nullptr,
		pos,
		backboneLogic,
		"backboneLayer",
		tileContainerLogic,
		0,
		0,
		OverlapCheck()
	);

	for(int i =0; i<nmatrix; i++){
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
			0,
			OverlapCheck()
		);

		current_z += (matrix_length + switch_length)/2.;
		pos = G4ThreeVector(0, 0, current_z);
		new G4PVPlacement(
			nullptr,
			pos,
			switchLogic,
			"switchLayer",
			tileContainerLogic,
			0,
			0,
			OverlapCheck()
		);

	}

	// place the tile in to the las box
	for(int i = 0; i<ntile; i++){
		current_z  = -las_container_length/2. + lec_length + (2*i+1)*tile_length/2.;
		pos = G4ThreeVector(0, 0, current_z);
		new G4PVPlacement(
			nullptr,
			pos,
			tileContainerLogic,
			"tileContainerLayer",
			las_container_logic,
			0,
			i,
			OverlapCheck()
		);

	}


	tub_length = lec_length;

	G4VSolid* lec_solid = new G4Tubs("lec_solid", tub_rinner, tub_router, tub_length/2., sta_phi, sto_phi);
	tub_length = rec_length;
	G4VSolid* rec_solid = new G4Tubs("rec_solid", tub_rinner, tub_router, tub_length/2., sta_phi, sto_phi);
	
	G4LogicalVolume* lec_logic =new G4LogicalVolume(lec_solid, mat_si, "LECLogic");
	lec_logic->SetVisAttributes(insen_vis);
	G4LogicalVolume* rec_logic = new G4LogicalVolume(rec_solid, mat_si, "RECLogic");
	rec_logic->SetVisAttributes(insen_vis);

	pos = G4ThreeVector(0., 0., -las_container_length/2. + lec_length/2.);
	new G4PVPlacement(
		nullptr,
		pos,
		lec_logic,
		"LAS_LEC",
		las_container_logic,
		0,
		0,
		OverlapCheck()
	);

	pos = G4ThreeVector(0., 0., las_container_length/2. - rec_length/2.);
	new G4PVPlacement(
		nullptr,
		pos,
		rec_logic,
		"LAS_REC",
		las_container_logic,
		0,
		0,
		OverlapCheck()
	);

	// kapton 
	tub_rinner = cf_curve_radius;
	tub_router = cf_curve_radius + kapton_thickness;
	double kapton_width = anc_length + m_nooverlap;
	sta_phi = 0;
	sto_phi = si_w/cf_curve_radius;

	G4VSolid* kapton_solid = new G4Tubs("kapton_solid", tub_rinner, tub_router, kapton_width/2., sta_phi, sto_phi);
	G4LogicalVolume* kapton_logic = new G4LogicalVolume(kapton_solid, mat_kapton, "KaptonLogic");
	kapton_logic->SetVisAttributes(kapton_vis);


	// Anc-ASIC 
	box_x = anc_length;
	box_y = anc_thickness;
	box_z = anc_length;
	G4VSolid* anc_solid = new G4Box("anc_solid", box_x/2., box_y/2., box_z/2.);
	G4LogicalVolume* anc_logic = new G4LogicalVolume(anc_solid
	, mat_si, "AncASICLogic");
	anc_logic->SetVisAttributes(insen_vis);


	//place the las container
	ostringstream oss;
	G4ThreeVector pos_kapton;
	G4ThreeVector pos_ancleft;
	G4ThreeVector pos_ancright;

	G4RotationMatrix rotanc_left;
	rotanc_left.rotateZ(si_w/(4*cf_curve_radius));
	G4RotationMatrix rotanc_right;
	rotanc_right.rotateZ(-si_w/(4*cf_curve_radius));
	double dt = cf_curve_radius*(1 - cos(si_w/(4*cf_curve_radius)));

	
	G4Transform3D trans_ancleft;
	G4Transform3D trans_ancright;


	for(int i =0; i<nz; i++){
		
		rotm = G4RotationMatrix();
		rotm.rotateZ(M_PI/2. - si_w/(2*cf_curve_radius));
		if(i >= nz/2) rotm.rotateY(M_PI);
		
		if(i%2 == 0){
			pos = G4ThreeVector(0., -anc_thickness-kapton_thickness, -stave_length/2. + cf_lec + i/2*(las_active_length + las_container_length) + las_container_length/2.);
			if( i < nz/2){   // top left
				pos_kapton = G4ThreeVector(0., -anc_thickness, -stave_length/2. + cf_lec + i/2*(las_active_length + las_container_length) - las_airspace/2. - kapton_width/2.);
				// anc asic
				posx = -si_w/4.;//- (cf_curve_radius + c_t + kapton_thickness + m_nooverlap + m_nooverlap + anc_thickness/2.)*sin(si_w/(2*cf_curve_radius));
				posy = cf_curve_radius - dt - anc_thickness;//(cf_curve_radius + c_t + kapton_thickness + m_nooverlap + m_nooverlap + anc_thickness/2.)*cos(si_w/(2*cf_curve_radius));
				posz = -stave_length/2. + cf_lec + i/2*(las_active_length + las_container_length) - las_airspace/2. - kapton_width/2.;
				pos_ancleft = G4ThreeVector(posx, posy, posz);
				trans_ancleft = G4Transform3D(rotanc_left, pos_ancleft);
				
				posx = si_w/4.; //(cf_curve_radius + c_t + kapton_thickness + m_nooverlap + m_nooverlap + anc_thickness/2.)*sin(si_w/(2*cf_curve_radius));
				pos_ancright = G4ThreeVector(posx, posy, posz);
				trans_ancright = G4Transform3D(rotanc_right, pos_ancright);
			}

			else{	// top right
				pos_kapton = G4ThreeVector(0., c_t + kapton_thickness + m_nooverlap, stave_length/2. - cf_rec - (nz - 1 - i)/2*(las_active_length + las_container_length) - las_active_length + las_airspace/2. + kapton_width/2.);
				// anc asic
				posx = -si_w/4.; //- (cf_curve_radius + c_t + kapton_thickness + m_nooverlap + m_nooverlap + anc_thickness/2.)*sin(si_w/(2*cf_curve_radius));
				posy = cf_curve_radius - dt - anc_thickness;//(cf_curve_radius + c_t + kapton_thickness + m_nooverlap + m_nooverlap + anc_thickness/2.)*cos(si_w/(2*cf_curve_radius));
				posz = stave_length/2. - cf_rec - (nz - 1 - i)/2*(las_active_length + las_container_length) - las_active_length + las_airspace/2. + kapton_width/2.;
				pos_ancleft = G4ThreeVector(posx, posy, posz);
				trans_ancleft = G4Transform3D(rotanc_left, pos_ancleft);
				
				posx = si_w/4.;//(cf_curve_radius + c_t + kapton_thickness + m_nooverlap + m_nooverlap + anc_thickness/2.)*sin(si_w/(2*cf_curve_radius));
				pos_ancright = G4ThreeVector(posx, posy, posz);
				trans_ancright = G4Transform3D(rotanc_right, pos_ancright);
			}
				
		}
		else{  
			rotm.rotateZ(M_PI);
			pos = G4ThreeVector(0., cf_curve_radius + cf_stave_container_ymin + anc_thickness + kapton_thickness, -stave_length/2. + cf_lec + (i/2+1)*las_active_length + i/2*las_container_length + las_container_length/2.);
			if(i < nz/2){  // bottom left
				pos_kapton = G4ThreeVector(0., cf_curve_radius + cf_stave_container_ymin + anc_thickness, -stave_length/2. + cf_lec + (i/2 +1)*las_active_length + i/2*las_container_length -las_airspace/2. - kapton_width/2.);
				
				//anc asic
				posx = -si_w/4.;//-(cf_curve_radius - cf_center_height - 2*c_t - m_nooverlap - m_nooverlap - anc_thickness/2.)*sin(si_w/(2*cf_curve_radius));
				posy = cf_stave_container_ymin + dt + anc_thickness;// (cf_curve_radius - cf_center_height - 2*c_t - m_nooverlap - m_nooverlap - anc_thickness/2.)*cos(si_w/(2*cf_curve_radius));
				posz = -stave_length/2. + cf_lec + (i/2 +1)*las_active_length + i/2*las_container_length -las_airspace/2. - kapton_width/2.;
				pos_ancleft = G4ThreeVector(posx, posy, posz);
				trans_ancleft = G4Transform3D(rotanc_right, pos_ancleft);

				posx = si_w/4.; //(cf_curve_radius - cf_center_height - 2*c_t - m_nooverlap - m_nooverlap - anc_thickness/2.)*sin(si_w/(2*cf_curve_radius));
				pos_ancright = G4ThreeVector(posx, posy, posz);
				trans_ancright = G4Transform3D(rotanc_left, pos_ancright);
			}

			else{ // bottom right
				pos_kapton = G4ThreeVector(0., cf_curve_radius + cf_stave_container_ymin + anc_thickness + kapton_thickness, stave_length/2. - cf_rec - (nz - 1 - i)/2*(las_active_length + las_container_length) + las_airspace/2. + kapton_width/2.);

				//anc asic
				posx = -si_w/4.; //-(cf_curve_radius - cf_center_height - 2*c_t - m_nooverlap - m_nooverlap - anc_thickness/2.)*sin(si_w/(2*cf_curve_radius));
				posy = cf_stave_container_ymin + dt + anc_thickness; //(cf_curve_radius - cf_center_height - 2*c_t - m_nooverlap - m_nooverlap - anc_thickness/2.)*cos(si_w/(2*cf_curve_radius));
				posz = stave_length/2. - cf_rec - (nz - 1 - i)/2*(las_active_length + las_container_length) + las_airspace/2. + kapton_width/2.;
				pos_ancleft = G4ThreeVector(posx, posy, posz);
				trans_ancleft = G4Transform3D(rotanc_right, pos_ancleft);

				posx = si_w/4.; //(cf_curve_radius - cf_center_height - 2*c_t - m_nooverlap - m_nooverlap - anc_thickness/2.)*sin(si_w/(2*cf_curve_radius));
				pos_ancright = G4ThreeVector(posx, posy, posz);
				trans_ancright = G4Transform3D(rotanc_left, pos_ancright);
			}
		}
		
		trans = G4Transform3D(rotm, pos);
		oss.str("");
		if(i%2 == 0) oss<<"LAS_TOP"<<i<<endl;
		else oss<<"LAS_BOT"<<i<<endl;
		new G4PVPlacement(
			trans,
			las_container_logic,
			oss.str().c_str(),
			cf_stave_container_logic,
			0,
			i,
			OverlapCheck()
		);

		//placement of the kapton foil
		if(i%2 == 0) oss<<"kapton_TOP"<<i<<endl;
		else oss<<"kapton_BOT"<<i<<endl;
		trans = G4Transform3D(rotm, pos_kapton);
		new G4PVPlacement(
			trans,
			kapton_logic,
			oss.str().c_str(),
			cf_stave_container_logic,
			0,
			i,
			0 //OverlapCheck()
		);
		
		// left anc-asic
		new G4PVPlacement(
			trans_ancleft,
			anc_logic,
			"ancASIC_left",
			cf_stave_container_logic,
			0,
			i,
			0 //OverlapCheck()
		);

		// right one
		new G4PVPlacement(
			trans_ancright,
			anc_logic,
			"ancASIC_right",
			cf_stave_container_logic,
			0,
			i,
			0 //OverlapCheck()
		);
	}

	// kapton foil on the stave edge
	box_x = kapton_thickness;
	box_y = cf_edge_height - 2*m_nooverlap;
	box_z = stave_length;

	G4VSolid* edge_kapton_solid = new G4Box("edge_kapton", box_x/2., box_y/2., box_z/2.);
	G4LogicalVolume* edge_kapton_logic = new G4LogicalVolume(edge_kapton_solid, mat_kapton, "EdgeKaptonLogic");
	edge_kapton_logic->SetVisAttributes(kapton_vis);
	//placement
	posx = strip_x - intersection_chord_length/2. - box_x;
	posy = cf_curve_radius - anc_thickness - kapton_thickness - cf_center_height/2.;
	posz = 0.;
	new G4PVPlacement(
		nullptr,
		G4ThreeVector(posx, posy, posz),
		edge_kapton_logic,
		"edge_kapton_left",
		cf_stave_container_logic,
		0,
		0,
		OverlapCheck()
	);
	posx = intersection_chord_length/2. - strip_x + box_x;
	new G4PVPlacement(
		nullptr,
		G4ThreeVector(posx, posy, posz),
		edge_kapton_logic,
		"edge_kapton_right",
		cf_stave_container_logic,
		0,
		0,
		OverlapCheck()
	);

	// place the curved stave
	rotm = G4RotationMatrix();
	for(int i = 0; i < nphi; i++){

		double dphi = 2*M_PI/nphi;
		double cur_phi = M_PI/2. + i*2*M_PI/nphi;
		double r_eff = 0;
		if(i%2 == 0) r_eff = r_inner - cf_curve_radius;
		else r_eff = r_outer - cf_curve_radius;

		posx = r_eff*cos(cur_phi);
		posy = r_eff*sin(cur_phi);
		posz = 0.;
		
		rotm.rotateZ(dphi);
		pos = G4ThreeVector(posx, posy, posz);
		trans = G4Transform3D(rotm, pos);

		new G4PVPlacement(
			trans,
			cf_stave_container_logic,
			"CFStaveContainer",
			logicWorld,
			0,
			i,
			OverlapCheck()
		);

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
