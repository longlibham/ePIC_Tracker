/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ePIC_TOFBarrel_SteppingAction.h
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.07.01
*   	Description：	
*
=========================================================================*/

#ifndef EPIC_TOFBARREL_STEPPINGACTION_H_
#define EPIC_TOFBARREL_STEPPINGACTION_H_
#include <iostream>
#include <fstream>
#include <sstream>

#include <g4main/PHG4SteppingAction.h>
#include <string>

class G4Step;
class G4VPhysicalVolume;
class PHComporiteNode;
class PHG4Hit;
class PHG4HitContainer;
class ePIC_TOFBarrel_Detector;

using namespace std;

class ePIC_TOFBarrel_SteppingAction : public PHG4SteppingAction{
	public:
		// constructor
		ePIC_TOFBarrel_SteppingAction(ePIC_TOFBarrel_Detector*);

		// destructor
		virtual ~ePIC_TOFBarrel_SteppingAction();
		
		// stepping action
		virtual bool UserSteppingAction(const G4Step*, bool);

		// reimplemented from base class
		virtual void SetInterfacePointers(PHCompositeNode*);


	private:
		// pointer to the detector
		ePIC_TOFBarrel_Detector* m_Detector;
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
