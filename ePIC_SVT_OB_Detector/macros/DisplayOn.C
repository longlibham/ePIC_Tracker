/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		DisplayOn.C
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.04.24
*   	Description：	
*
=========================================================================*/
#ifndef MACRO_DISPLAYON_C
#define MACRO_DISPLAYON_C

#include <g4main/PHG4Reco.h>
#include <fun4all/Fun4AllServer.h>

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

// this starts the QT based G4 Gui which takes control
// when x'ed out it will return a pointor to PHG4Reco so
// the gui can be started again
//
PHG4Reco* QTGui(){
	Fun4AllServer* se = Fun4AllServer::instance();
	PHG4Reco* g4 = (PHG4Reco*) se->getSubsysReco("PHG4Reco");
	g4->InitRun(se->topNode());
	g4->ApplyDisplayAction();
	g4->StartGui();
	return g4;
}

//stupid macro to turn on the geant4 display 
//we ask a Fun4All for a pointer to PHG4Reco 
//using the ApplyCommand will start up the G4
//cmd and interpreter and graphics system the 
//vis.mac contains the necessary command to start
//up the visualisation, the next event will be 
//displayed. Do not excute the macro before the 
//PHG4Reco is registered with Fun4All
//
PHG4Reco* DisplayOn(const char* mac = "vis.mac"){
	char cmd[100];
	Fun4AllServer* se = Fun4AllServer::instance();
	PHG4Reco* g4 = (PHG4Reco*) se->getSubsysReco("PHG4Reco");
	g4->InitRun(se->topNode());
	g4->ApplyDisplayAction();
	sprintf(cmd, "/control/excute/ %s", mac);
	g4->ApplyCommand(cmd);
	return g4;
}

// print out the command I always forget
//
void displaycmd(){
	cout << "draw axis: " << endl;
	cout << "g4->ApplyCommand(\"/vis/scene/add/axes 0 0 50 cm\")" << endl;
	cout << "zoom" <<endl;
	cout << "g4->ApplyCommand(\"/vis/viewer/zoom 1\")" << endl;
	cout << "viewpoint: " << endl;
	cout << "g4->ApplyCommand(\"/vis/viewer/set/viewpointThetaPhi 0 0 \")" << endl;
	cout << "panTo: " << endl;
	cout << "g4->ApplyCommand(\"/vis/viewer/panTo 0 0 cm\")" << endl;
	cout << "print to eps: " << endl;
	cout << "g4->ApplyCommand(\"/vis/ogl/printEPS \")" << endl;
	cout << "set background color: " << endl;
	cout << "g4->ApplyCommand(\"vis/viewer/set/background white\")" << endl;
}

#endif
