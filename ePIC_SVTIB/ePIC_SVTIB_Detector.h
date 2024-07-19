/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_SVTIB_Detector.h
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.07.18
*   	Description：	
*
=========================================================================*/

#ifndef EPIC_SVTIB_DETECTOR_H_
#define EPIC_SVTIB_DETECTOR_H_
#include <iostream>
#include <fstream>
#include <sstream>

#include <g4main/PHG4Detector.h>
#include <Geant4/G4SystemOfUnits.hh>

#include <set>
#include <string>


class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;
class PHParameters;

class ePIC_SVTIB_Detector : public PHG4Detector{
	public:
		// constructor
		ePIC_SVTIB_Detector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string &dname, const int layer = 0);
		// destructor
		virtual ~ePIC_SVTIB_Detector(){}

		// geometry builder
		virtual void ConstructMe(G4LogicalVolume* world);
		virtual void Print(const std::string& what= "ALL") const;
		
		int IsInDetector(G4VPhysicalVolume* ) const;

		void SuperDetector(const std::string& name){m_SuperDetector = name;}
		const std::string SuperDetector(){return m_SuperDetector;}
		int get_Layer() const {return m_Layer;}

	private:
		PHParameters* m_Params;
		double m_redundancy = 0.1 * cm;
        double m_nooverlap = 0.00001 * cm;
		std::set<G4VPhysicalVolume*> m_PhysicalVolumesSet;
		std::string m_SuperDetector;
		int m_Layer;
};


#endif
