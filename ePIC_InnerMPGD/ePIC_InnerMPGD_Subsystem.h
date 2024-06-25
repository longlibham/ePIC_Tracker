/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_InnerMPGD_Subsystem.h
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.06.25
*   	Description：	
*
=========================================================================*/

#ifndef EPIC_INNERMPGD_SUBSYSTEM_H_
#define EPIC_INNERMPGD_SUBSYSTEM_H_
#include <iostream>
#include <fstream>
#include <sstream>

#include <g4main/PHG4Subsystem.h>
#include <g4detectors/PHG4DetectorSubsystem.h>
#include <string>

class PHCompositeNode;
class PHG4Detector;
class ePIC_InnerMPGD_Detector;
class PHG4SteppingAction;

using namespace std;


class ePIC_InnerMPGD_Subsystem : public PHG4DetectorSubsystem{
	public:
		// constructor
		ePIC_InnerMPGD_Subsystem(const std::string& name = "InnerMPGD", const int layer = 0);

		// destructor
		virtual ~ePIC_InnerMPGD_Subsystem(){}

		// creates the m_Detector object and place it on the node tree, under "DETECTORS" node (or whatever) 
		// creates the stepping action and place it on the node tree, under "ACTIONS" node
		// creates relevant hit nodes that will be populated by the stepping action and stored in the output DST

		int InitRunSubsystem(PHCompositeNode*) override;
		virtual PHG4Detector* GetDetector() const;
		virtual PHG4SteppingAction* GetSteppingAction() const {return m_SteppingAction;}

		// Print info (from SubsysReco)
		virtual void Print(const std::string& what = "ALL") const;

		int process_event(PHCompositeNode*) override;

	private:

		// detector geo
		// define from PHG4Detector
		ePIC_InnerMPGD_Detector* m_Detector;

		void SetDefaultParameters() override;

		// particle tracking 'stepping' action
		// derives from PHG4SteppingActions
		PHG4SteppingAction* m_SteppingAction;
};


#endif
