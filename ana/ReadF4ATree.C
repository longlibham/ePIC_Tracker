/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		ReadF4ATree.C
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.03.28
*   	Description：	
*
=========================================================================*/

#include <iostream>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TList.h"
#include "TH1D.h"
#include "TGraph.h"
#include "Fun4AllFastTracking.h"

using namespace std;


void ReadF4ATree(){
	auto demo = new Fun4AllFastTracking();
	demo->Loop();
}
