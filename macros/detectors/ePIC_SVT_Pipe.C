/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_SVT_Pipe.C
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.06.20
*   	Description：	
*
=========================================================================*/
#ifndef MACRO_EPICSVTPIPE_C
#define MACRO_EPICSVTPIPE_C

#include <GlobalVariables.C>
#include <g4detectors/PHG4ConeSubsystem.h>
#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4detectors/PHG4GDMLSubsystem.h>
#include <g4main/PHG4Reco.h>

#include <iostream>
#include <fstream>
#include <sstream>

R__LOAD_LIBRARY(libg4detectors.so)

using namespace std;

// This creates the Enable Flag to be used in the main steering macro
namespace Enable
{
  bool PIPE = false;
  bool PIPE_ABSORBER = false;
  bool PIPE_OVERLAPCHECK = false;
  int PIPE_VERBOSITY = 0;
}  // namespace Enable

namespace ePICPIPE
{
  // Central pipe dimension
  // Extracted via mechanical model: Detector chamber 3-20-20
  // directly implimenting the central Be section in G4 cylinder for max speed simulation in the detector region.
  // The jointer lip structure of the pipe R = 3.2cm x L=5mm is ignored here
  double be_pipe_radius = 3.1000;
  double Au_coating_thickness = 2e-4; // 2um Au coating
  double be_pipe_thickness = 3.1762 - be_pipe_radius;  // 760 um for sPHENIX
  double be_pipe_length_plus = 66.8;                   // +z beam pipe extend.
  double be_pipe_length_neg = -79.8;                   // -z beam pipe extend.
  bool use_forward_pipes = true;
}  // namespace ePICPIPE

void PipeInit()
{
  if (ePICPIPE::use_forward_pipes)
  {
    BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 23.);
    BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, 500.);
    BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -463.);
  }
  else
  {
    BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, ePICPIPE::be_pipe_radius + ePICPIPE::be_pipe_thickness);
    BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, ePICPIPE::be_pipe_length_plus);
    BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, ePICPIPE::be_pipe_length_neg);
  }
}

//! construct beam pipe
double Pipe(PHG4Reco* g4Reco,
            double radius)
{
  bool AbsorberActive = Enable::ABSORBER || Enable::PIPE_ABSORBER;
  bool OverlapCheck = Enable::PIPE_OVERLAPCHECK or Enable::OVERLAPCHECK; // ||  to suppress GDML validation messages
  int verbosity = std::max(Enable::VERBOSITY, Enable::PIPE_VERBOSITY);
  // process pipe extentions?
  const bool do_pipe_hadron_forward_extension = ePICPIPE::use_forward_pipes && true;
  const bool do_pipe_electron_forward_extension = ePICPIPE::use_forward_pipes && true;

  const double be_pipe_length = ePICPIPE::be_pipe_length_plus - ePICPIPE::be_pipe_length_neg;  // pipe length
  const double be_pipe_center = 0.5 * (ePICPIPE::be_pipe_length_plus + ePICPIPE::be_pipe_length_neg);

  if (radius > ePICPIPE::be_pipe_radius)
  {
    cout << "inconsistency: radius: " << radius
         << " larger than pipe inner radius: " << ePICPIPE::be_pipe_radius << endl;
    gSystem->Exit(-1);
  }

  // mid-rapidity beryillium pipe
  PHG4CylinderSubsystem* cyl = new PHG4CylinderSubsystem("VAC_BE_PIPE", 0);
  cyl->set_double_param("radius", 0.0);
  cyl->set_int_param("lengthviarapidity", 0);
  cyl->set_double_param("length", be_pipe_length);
  cyl->set_double_param("place_z", be_pipe_center);
  cyl->set_string_param("material", "G4_Galactic");
  cyl->set_double_param("thickness", ePICPIPE::be_pipe_radius- ePICPIPE::Au_coating_thickness);
  cyl->SuperDetector("PIPE");
  cyl->OverlapCheck(OverlapCheck);
  g4Reco->registerSubsystem(cyl);

  // 2um Au coating for X-ray absorption
  cyl = new PHG4CylinderSubsystem("Au_PIPE", 2);
  cyl->set_double_param("radius", ePICPIPE::be_pipe_radius - ePICPIPE::Au_coating_thickness);
  cyl->set_int_param("lengthviarapidity", 0);
  cyl->set_double_param("length", be_pipe_length);
  cyl->set_double_param("place_z", be_pipe_center);
  cyl->set_string_param("material", "G4_Au");
  cyl->set_double_param("thickness", ePICPIPE::Au_coating_thickness);
  cyl->SuperDetector("PIPE");
  cyl->OverlapCheck(OverlapCheck);
  if (AbsorberActive) cyl->SetActive();
  g4Reco->registerSubsystem(cyl);

  cyl = new PHG4CylinderSubsystem("BE_PIPE", 1);
  cyl->set_double_param("radius", ePICPIPE::be_pipe_radius);
  cyl->set_int_param("lengthviarapidity", 0);
  cyl->set_double_param("length", be_pipe_length);
  cyl->set_double_param("place_z", be_pipe_center);
  cyl->set_string_param("material", "G4_Be");
  cyl->set_double_param("thickness", ePICPIPE::be_pipe_thickness);
  cyl->SuperDetector("PIPE");
  cyl->OverlapCheck(OverlapCheck);
  if (AbsorberActive) cyl->SetActive();
  g4Reco->registerSubsystem(cyl);

  radius = ePICPIPE::be_pipe_radius + ePICPIPE::be_pipe_thickness;

  radius += no_overlapp;

  if (do_pipe_electron_forward_extension)
  {
    if (Enable::IP6)
    {
      PHG4GDMLSubsystem* gdml = new PHG4GDMLSubsystem("ElectronForwardChamber");
      gdml->set_string_param("GDMPath", string(getenv("CALIBRATIONROOT")) + "/Beam/ConstructSimplifiedBeamChamber.gdml");
      gdml->set_string_param("TopVolName", "ElectronForwardChamber");
      gdml->OverlapCheck(OverlapCheck);
      g4Reco->registerSubsystem(gdml);
    }
    if (Enable::IP8)
    {
      cout <<__PRETTY_FUNCTION__<<" IP8 beam chamber not defined yet! Consider disable ePICPIPE::use_forward_pipes"<<endl;
      exit(1);
    }
  }

  if (do_pipe_hadron_forward_extension)
  {
    if (Enable::IP6)
    {
      PHG4GDMLSubsystem* gdml = new PHG4GDMLSubsystem("HadronForwardChamber");
      gdml->set_string_param("GDMPath", string(getenv("CALIBRATIONROOT")) + "/Beam/ConstructSimplifiedBeamChamber.gdml");
      gdml->set_string_param("TopVolName", "HadronForwardChamber");
      gdml->OverlapCheck(OverlapCheck);
      g4Reco->registerSubsystem(gdml);
    }
    else if (Enable::IP8)
    {
      cout <<__PRETTY_FUNCTION__<<" IP8 beam chamber not defined yet! Consider disable ePICPIPE::use_forward_pipes"<<endl;
      exit(1);
    }
  }

  return radius;
}
#endif

