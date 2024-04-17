/*---------------------------------------------------------------------*
 * Barrel tracker designed by LANL EIC team                            *
 * See technical notes for details: arXiv:2009.02888                   *
 * Contact Ping and Xuan @LANL for questions:                          *
 *   Xuan: xuanli@lanl.gov                                             *
 *   Ping: cpwong@lanl.gov                                             *
 *---------------------------------------------------------------------*/

#ifndef MACRO_G4FSTEIC_C
#define MACRO_G4FSTEIC_C

#include "GlobalVariables.C"

#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4detectors/PHG4SectorSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>

#include <g4main/PHG4Reco.h>

// #include <disk1/Disk1Subsystem.h>
// #include <disk3p/Disk3pSubsystem.h>
// #include <disk4p/Disk4pSubsystem.h>
// #include <disk5p/Disk5pSubsystem.h>
// #include <disk4n/Disk4nSubsystem.h>
// #include <disk5n/Disk5nSubsystem.h>
#include "/home/scm/Fun4All/myinstall/include/disk1/Disk1Subsystem.h"
#include "/home/scm/Fun4All/myinstall/include/disk3p/Disk3pSubsystem.h"
#include "/home/scm/Fun4All/myinstall/include/disk4p/Disk4pSubsystem.h"
#include "/home/scm/Fun4All/myinstall/include/disk5p/Disk5pSubsystem.h"
#include "/home/scm/Fun4All/myinstall/include/disk4n/Disk4nSubsystem.h"
#include "/home/scm/Fun4All/myinstall/include/disk5n/Disk5nSubsystem.h"

#include <string>

R__LOAD_LIBRARY(libg4detectors.so)

R__LOAD_LIBRARY(libDisk1.so)
R__LOAD_LIBRARY(libDisk3p.so)
R__LOAD_LIBRARY(libDisk4p.so)
R__LOAD_LIBRARY(libDisk5p.so)
R__LOAD_LIBRARY(libDisk4n.so)
R__LOAD_LIBRARY(libDisk5n.so)

int make_LANL_FST_station(string name, PHG4Reco *g4Reco, double zpos, double Rmin,
                          double Rmax, double offset, double pitch);
int make_supportCyl(string name, PHG4Reco *g4Reco,
                    double r, double t, double length);
//-----------------------------------------------------------------------------------//
namespace Enable
{
  static bool FST = false;
  bool FST_OVERLAPCHECK = false;
}  // namespace Enable

namespace G4FST
{
  namespace SETTING
  {
    bool FST_TPC = false;
    bool SUPPORTCYL = false;
  }  // namespace SETTING
}  // namespace G4FST

//-----------------------------------------------------------------------------------//
void FST_Init()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 48.);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, 127.);
  if (G4FST::SETTING::SUPPORTCYL)
  {
    BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -127.);
  }
}
//-----------------------------------------------------------------------------------//

/* 

The idea is to calculate all the rmax rmin based on z and theta1 and theta2 in a python script (obj_fun.py)
Once calculated just pass the disk dimensions into this script

*/

void FSTSetup(PHG4Reco *g4Reco)
{
  gSystem->Load("libDisk1.so");
  gSystem->Load("libDisk3p.so");
  gSystem->Load("libDisk4p.so");
  gSystem->Load("libDisk5p.so");
  gSystem->Load("libDisk4n.so");
  gSystem->Load("libDisk5n.so");

  const double cm = PHG4Sector::Sector_Geometry::Unit_cm();
  const double mm = .1 * cm;
  const double um = 1e-3 * mm;
  //const double max_radius = 50.;



  const double bkwd_z[] = {25, 45, 70, 100, 200};
  const double bkwd_offset[] = {0., 0., 0., 0.2, 0.7};
  //double bkwd_rmin[] = {3.6, 3.6, 3.6, 3.9, 4.5};
  double bkwd_rmin[] = {4.54, 4.54, 4.54, 5.33, 5.65};
  //double bkwd_rmax[] = {19., 43., 43., 43., 59.};
  double bkwd_rmax[] = {25., 42.5, 43.12, 43.12, 100.};
  const int n_bkwd_disk = sizeof(bkwd_z) / sizeof(*bkwd_z);
  for (unsigned int i = (n_bkwd_disk -1); i < n_bkwd_disk; i++)
  {
/*
    // Below was made to auto calculate the min and max Radius
    if(bkwd_z[i] < uRwell1_e_length) bkwd_rmax[i] = std::min(max_radius, e_slope1*bkwd_z[i] + e_intercept1) - 0.5;
    else if (bkwd_z[i] >= uRwell1_e_length && bkwd_z[i] <= (uRwell1_e_length + uRwell_plateau_length)){bkwd_rmax[i] = uRwell1_radius - 1.5;}
    else if(bkwd_z[i] > (uRwell1_e_length + uRwell_plateau_length)){bkwd_rmax[i] = std::min(max_radius, e_slope2*bkwd_z[i] + e_intercept2) - 0.5;}
    else {cout << "Cannot calculate the RMax exiting" << endl; gSystem->Exit(0); }

    if(bkwd_z[i]>79.8 && bkwd_z[i]>0) bkwd_rmin[i] = (0.0521*bkwd_z[i] + 1.0);
    else bkwd_rmin[i] = 3.3;
*/  
    make_LANL_FST_station(Form("EST_%i", i), g4Reco, -1*bkwd_z[i], bkwd_rmax[i]-0.0001, bkwd_rmax[i], bkwd_offset[i], 10e-4);  //cm
  }


  const double fwd_z[] = {25, 45, 70, 100, 135};
  const double fwd_offset[] = {0., 0., 0., -0.8, -1.7};
  //double fwd_rmin[] = {3.6, 3.6, 3.6, 4.5, 5.4};
  double fwd_rmin[] = {4.54, 4.54, 4.63, 5.61, 6.78};
  //double fwd_rmax[] = {19., 43., 43., 43., 53.};
  double fwd_rmax[] = {25., 42.5, 43.12, 43.12, 43.12};
  const int n_fwd_disk = sizeof(fwd_z) / sizeof(*fwd_z);
  // for (unsigned int i = 0; i < n_fwd_disk; i++)
  // {

/*
    if(fwd_z[i] < uRwell1_h_length) fwd_rmax[i] = std::min(max_radius, h_slope1*fwd_z[i] + h_intercept1) - 0.5;
    else if (fwd_z[i] >= uRwell1_h_length && fwd_z[i] <= (uRwell1_h_length + uRwell_plateau_length)){fwd_rmax[i] = uRwell1_radius - 1.5;}
    //else if(fwd_z[i] > (uRwell1_h_length + uRwell_plateau_length)){fwd_rmax[i] = std::min(max_radius, h_slope2*fwd_z[i] + h_intercept2) - 0.5;}
    else if(fwd_z[i] > (uRwell1_h_length + uRwell_plateau_length) && fwd_z[i] <= 130.){ fwd_rmax[i] = std::min(max_radius, h_slope2*fwd_z[i] + h_intercept2) - 0.5;}
    //else if(fwd_z[i] > 113.){fwd_rmax[i] = h_slope2*(fwd_z[i] + 5.) + h_intercept2 - 0.25;}
    else {cout << "Cannot calculate the RMax exiting" << endl; gSystem->Exit(0); }

    if(fwd_z[i]>66.8 && fwd_z[i]>0) fwd_rmin[i] = (0.0521*fwd_z[i] + 1.0);
    else fwd_rmin[i] = 3.3;
*/
    // make_LANL_FST_station(Form("FST_%i", i), g4Reco, fwd_z[i], fwd_rmin[i], fwd_rmax[i], fwd_offset[i], 10e-4);  //cm
  // }

// Positive Disks
  Disk1Subsystem *disk1;
  Disk3pSubsystem *disk3p;
  Disk4pSubsystem *disk4p;
  Disk5pSubsystem *disk5p;
  Disk4nSubsystem *disk4n;
  Disk5nSubsystem *disk5n;

  std::string diskName = "FBVS_" + std::to_string(0);
  disk1 = new Disk1Subsystem(diskName);
  disk1->set_string_param("material", "G4_Si");
  disk1->set_double_param("r_max", 23.);
  disk1->set_double_param("place_z", 25.);
  disk1->set_double_param("rot_z", 180.);
  disk1->SetActive();
  disk1->SuperDetector(diskName);
  
  g4Reco->registerSubsystem(disk1);
  
  
  diskName = "FBVS_" + std::to_string(1);
  disk1 = new Disk1Subsystem(diskName);
  disk1->set_string_param("material", "G4_Si");
  disk1->set_double_param("r_max", 43.);
  disk1->set_double_param("place_z", 45.);
  disk1->set_double_param("rot_z", 180.);
  disk1->SetActive();
  disk1->SuperDetector(diskName);
  
  g4Reco->registerSubsystem(disk1);
  
  
  diskName = "FBVS_" + std::to_string(2);
  disk3p = new Disk3pSubsystem(diskName);
  disk3p->set_string_param("material", "G4_Si");
  disk3p->set_double_param("r_max", 43.);
  disk3p->set_double_param("place_z", 70.);
  disk3p->set_double_param("rot_z", 180.);
  disk3p->SetActive();
  disk3p->SuperDetector(diskName);
   
  g4Reco->registerSubsystem(disk3p);
  
  
  diskName = "FBVS_" + std::to_string(3);
  disk4p = new Disk4pSubsystem(diskName);
  disk4p->set_string_param("material", "G4_Si");
  disk4p->set_double_param("r_max", 43.);
  disk4p->set_double_param("place_z", 100.);
  disk4p->set_double_param("rot_z", 180.);
  disk4p->SetActive();
  disk4p->SuperDetector(diskName);
   
  g4Reco->registerSubsystem(disk4p);
  
  
  diskName = "FBVS_" + std::to_string(4);
  disk5p = new Disk5pSubsystem(diskName);
  disk5p->set_string_param("material", "G4_Si");
  disk5p->set_double_param("r_max", 43.);
  disk5p->set_double_param("place_z", 135.);
  disk5p->set_double_param("rot_z", 180.);
  disk5p->SetActive();
  disk5p->SuperDetector(diskName);
   
  g4Reco->registerSubsystem(disk5p);


  // Negative Disks
  diskName = "FBVS_" + std::to_string(5);
  disk1 = new Disk1Subsystem(diskName);
  disk1->set_string_param("material", "G4_Si");
  disk1->set_double_param("r_max", 23.);
  disk1->set_double_param("place_z", -25.);
  disk1->set_double_param("rot_z", 180.);
  disk1->SetActive();
  disk1->SuperDetector(diskName);
  
  g4Reco->registerSubsystem(disk1);
  
  
  diskName = "FBVS_" + std::to_string(6);
  disk1 = new Disk1Subsystem(diskName);
  disk1->set_string_param("material", "G4_Si");
  disk1->set_double_param("r_max", 43.);
  disk1->set_double_param("place_z", -45.);
  disk1->set_double_param("rot_z", 180.);
  disk1->SetActive();
  disk1->SuperDetector(diskName);
  
  g4Reco->registerSubsystem(disk1);
  
  
  diskName = "FBVS_" + std::to_string(7);
  disk1 = new Disk1Subsystem(diskName);
  disk1->set_string_param("material", "G4_Si");
  disk1->set_double_param("r_max", 43.);
  disk1->set_double_param("place_z", -70.);
  disk1->set_double_param("rot_z", 180.);
  disk1->SetActive();
  disk1->SuperDetector(diskName);
  
  g4Reco->registerSubsystem(disk1);
  
  
  diskName = "FBVS_" + std::to_string(8);
  disk4n = new Disk4nSubsystem(diskName);
  disk4n->set_string_param("material", "G4_Si");
  disk4n->set_double_param("r_max", 43.);
  disk4n->set_double_param("place_z", -100);
  disk4n->set_double_param("rot_z", 180.);
  disk4n->SetActive();
  disk4n->SuperDetector(diskName);
  
  g4Reco->registerSubsystem(disk4n);
  
  
  diskName = "FBVS_" + std::to_string(9);
  disk5n = new Disk5nSubsystem(diskName);
  disk5n->set_string_param("material", "G4_Si");
  disk5n->set_double_param("r_max", 43.);
  disk5n->set_double_param("place_z", -135.);
  disk5n->set_double_param("rot_z", 180.);
  disk5n->SetActive();
  disk5n->SuperDetector(diskName);
  
  g4Reco->registerSubsystem(disk5n);
  

  if (TRACKING::FastKalmanFilter)
    {
      for (int i{0}; i<10; i++){
	std::string diskName = "FBVS_" + std::to_string(i) + string("_0");
	TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + diskName,           //      const std::string& phg4hitsNames,
						 PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
						 10.e-4 / sqrt(12.),                 //      const float radres,
						 10.e-4 / sqrt(12.),                 //      const float phires,
						 50e-4 / sqrt(12.),                 //      const float lonres, *ignored in plane detector*
						 1.,                                 //      const float eff,
						 0);                                //      const float noise
	
      }
    }







  if (G4FST::SETTING::SUPPORTCYL)
  {
    double gap = 8;                                                               //cm
    double tSupport = 0.2;                                                        //cm
    make_supportCyl("FSTSupportCyl", g4Reco, 50.1 + gap, tSupport, 125.0 * 2.0);  //cm
  }
}
//-----------------------------------------------------------------------------------//
int make_LANL_FST_station(string name, PHG4Reco *g4Reco,
                          double zpos, double Rmin, double Rmax, double offset, double pitch)  //silicon thickness
{
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::FST_OVERLAPCHECK;
  
  double min_polar_angle = atan2(Rmin, zpos);
  double max_polar_angle = atan2(Rmax, zpos);
  
  PHG4CylinderSubsystem *fst;
  fst = new PHG4CylinderSubsystem(name);
  
  fst->set_string_param("material", "G4_Si");
  fst->set_double_param("radius", Rmin);
  fst->set_double_param("thickness", Rmax - Rmin);
  fst->set_double_param("place_x", offset);
  fst->set_double_param("place_z", zpos);
  fst->set_double_param("length", 0.24/100.*9.37);
  // fst->SetActive();
  fst->OverlapCheck(OverlapCheck);
  fst->SuperDetector(name);
  
  
  g4Reco->registerSubsystem(fst);
  // 
  // if (TRACKING::FastKalmanFilter)
  // {
  // TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + name,           //      const std::string& phg4hitsNames,
  // PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
  // pitch / sqrt(12.),                 //      const float radres,
  // pitch / sqrt(12.),                 //      const float phires,
  // 50e-4 / sqrt(12.),                 //      const float lonres, *ignored in plane detector*
  // 0.95,                                 //      const float eff,
  // 0);                                //      const float noise
  // }
  return 0;
}
//-----------------------------------------------------------------------------------//
int make_supportCyl(string name, PHG4Reco *g4Reco, double r, double t, double length)
{
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::FST_OVERLAPCHECK;

  PHG4CylinderSubsystem *cyl = new PHG4CylinderSubsystem(name, 5);
  cyl->set_double_param("radius", r);
  cyl->set_double_param("length", length);
  cyl->set_string_param("material", "CFRP_INTT");  // borrow carbon fiber reinforced polymer used in sPHENIX silicon tracker support
  cyl->set_double_param("thickness", t);
  cyl->set_double_param("place_x", 0.);
  cyl->set_double_param("place_y", 0.);
  cyl->set_double_param("place_z", 0);
  cyl->OverlapCheck(OverlapCheck);  //true);//overlapcheck);
  cyl->SetActive(0);
  //cyl->SuperDetector("");
  cyl->OverlapCheck(Enable::FST_OVERLAPCHECK);  //OverlapCheck);

  g4Reco->registerSubsystem(cyl);
  return 0;
}
#endif

