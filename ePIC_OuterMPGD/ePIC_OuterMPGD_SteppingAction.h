/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_OuterMPGD_SteppingAction.h
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.07.02
*   	Description：	
*
=========================================================================*/

#ifndef EPIC_OUTERMPGD_STEPPINGACTION_H_
#define EPIC_OUTERMPGD_STEPPINGACTION_H_
#include <iostream>
#include <fstream>
#include <sstream>

#include <g4main/PHG4SteppingAction.h>
#include <string>

using namespace std;

class G4Step;
class G4VPhysicalVolume;
class PHComporiteNode;
class PHG4Hit;
class PHG4HitContainer;
class ePIC_OuterMPGD_Detector;

class ePIC_OuterMPGD_SteppingAction : public PHG4SteppingAction{
	public:
		// constructor
		ePIC_OuterMPGD_SteppingAction(ePIC_OuterMPGD_Detector*);

		// destructor
		virtual ~ePIC_OuterMPGD_SteppingAction();
		
		// stepping action
		virtual bool UserSteppingAction(const G4Step*, bool);

		// reimplemented from base class
		virtual void SetInterfacePointers(PHCompositeNode*);


	private:
		// pointer to the detector
		ePIC_OuterMPGD_Detector* m_Detector;
		//pointer to hit container
		PHG4HitContainer* m_HitContainer;
		PHG4Hit* m_Hit;
		PHG4HitContainer* m_SaveHitContainer;

		G4VPhysicalVolume* m_SaveVolPre;
		G4VPhysicalVolume* m_SaveVolPost;

		int m_SaveTrackId;
		int m_SavePreStepStatus;
		int m_SavePostStepStatus;
		double m_EdepSum;
		double m_EionSum;
};

#endif