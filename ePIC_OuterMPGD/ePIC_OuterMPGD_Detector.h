/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_OuterMPGD_Detector.h
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.07.02
*   	Description：	
*
=========================================================================*/

#ifndef EPIC_OUTERMPGD_DETECTOR_H_
#define EPIC_OUTERMPGD_DETECTOR_H_
#include <iostream>
#include <fstream>
#include <sstream>

#include <g4main/PHG4Detector.h>

#include <set>
#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;
class PHParameters;

using namespace std;

class ePIC_OuterMPGD_Detector: public PHG4Detector{

    public:
        //constructor
        ePIC_OuterMPGD_Detector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string &dname, const int layer = 0);
		// destructor
		virtual ~ePIC_OuterMPGD_Detector(){}

		// geometry builder
		virtual void ConstructMe(G4LogicalVolume* world);
		virtual void Print(const std::string& what= "ALL") const;
		
		int IsInDetector(G4VPhysicalVolume* ) const;

		void SuperDetector(const std::string& name){m_SuperDetector = name;}
		const std::string SuperDetector(){return m_SuperDetector;}
		int get_Layer() const {return m_Layer;}

	private:
		//define material for MPGD
		void DefineMaterials();
		PHParameters* m_Params;
		std::set<G4VPhysicalVolume*> m_PhysicalVolumesSet;
		std::string m_SuperDetector;
		int m_Layer;
};

#endif
