/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_SVT_OB_Detector.h
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.04.08
*   	Description：	
*
=========================================================================*/

#ifndef EPIC_SVT_OB_DETECTOR_H_
#define EPIC_SVT_OB_DETECTOR_H_
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

class ePIC_SVT_OB_Detector : public PHG4Detector{
	public:
		// constructor
		ePIC_SVT_OB_Detector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string &dname);
		// destructor
		virtual ~ePIC_SVT_OB_Detector(){}

		// geometry builder
		virtual void ConstructMe(G4LogicalVolume* world);
		virtual void Print(const std::string& what= "ALL") const;
		
		int IsInDetector(G4VPhysicalVolume* ) const;

		void SuperDetector(const std::string& name){m_SuperDetector = name;}
		const std::string SuperDetector(){return m_SuperDetector;}

	private:
		std::set<G4VPhysicalVolume*> m_PhysicalVolumesSet;
		std::string m_SuperDetector;
};


#endif
