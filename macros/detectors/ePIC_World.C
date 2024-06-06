/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_World.cxx
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.06.03
*   	Description：	
*
=========================================================================*/
#ifndef MACRO_EPICWORLD_C
#define MACRO_EPICWORLD_C

#include <iostream>
#include <fstream>
#include <sstream>

#include <GlobalVariables.C>
#include <g4main/PHG4Reco.h>

using namespace std;

R__LOAD_LIBRARY(libg4testbench.so)

namespace ePICWORLD{
    double AddSpace = 100.;
    string WorldMaterial = "G4_AIR";
    string PhysicsList = "FTFP_BERT";
} // namespace ePIC_WORLD


void WorldInit(PHG4Reco* g4Reco){
    g4Reco->SetWorldMaterial(ePICWORLD::WorldMaterial);
    g4Reco->SetPhysicsList(ePICWORLD::PhysicsList);
}

void WorldSize(PHG4Reco * g4Reco, double radius){
    double world_radius = std::max(BlackHoleGeometry::max_radius + BlackHoleGeometry::gap, radius);
    g4Reco->SetWorldSizeY(std::max(g4Reco->GetWorldSizeY(), world_radius + ePICWORLD::AddSpace));
    // ePIC world is a symetric cylinder so the center is at (0, 0, 0), pick the largest of abs(min_z) || abs(max_z)
    double min_zval = std::min((BlackHoleGeometry::min_z - BlackHoleGeometry::gap), -((g4Reco->GetWorldSizeZ() - 100)/2.));
    double max_zval = std::max((BlackHoleGeometry::max_z + BlackHoleGeometry::gap), (g4Reco->GetWorldSizeZ() - 100)/2.);
    double final_zval = std::max(fabs(min_zval), fabs(max_zval)) + ePICWORLD::AddSpace;
    g4Reco->SetWorldSizeZ(std::max(g4Reco->GetWorldSizeZ(), 2*(final_zval)));
    return ;

}

#endif

