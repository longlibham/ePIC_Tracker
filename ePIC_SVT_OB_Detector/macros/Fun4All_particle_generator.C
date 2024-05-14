/*========================================================================
*   Copyright (C) 2024 Univ. of Bham. All rights reserved.
*   
*   	FileName：		Fun4All_particle_generator.C
*   	Author：		LongLI <long.l@cern.ch>
*   	Time：			2024.05.14
*   	Description：	
*
=========================================================================*/
#ifndef MACRO_FUN4ALLPARTICLEGENERATOR_C
#define MACRO_FUN4ALLPARTICLEGENERATOR_C

#include <GlobalVariables.C>

#include <DisplayOn.C>
#include <G4_Global.C>
#include <G4_Input.C>

#include <TROOT.h>

#include <fun4all/Fun4AllServer.h>
#include <phool/recoConsts.h>
#include <RooUnblindPrecision.h>

#include <iostream>
#include <fstream>
#include <sstream>


using namespace std;

R__LOAD_LIBRARY(libfun4all.so)

namespace Enable{
	bool USE_P = false;
	bool USE_PT = true;
	//bool GENERATE_SEED = false;
}

void Fun4All_particle_generator(
		const int nEvents = 1000,
//		recoConsts* rc,
		const string& particle_name = "pi-",
		double pmin = 0.,
		double pmax = 20.,
		double etamin = -0.8,
		double etamax = 0.8,
		double phimin = -M_PI,
		double phimax = M_PI,
		const string& inputFile = "https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/sPHENIX_G4Hits_sHijing_9-11fm_00000_00010.root",
		const int skip = 0
		){

	if(Enable::USE_P && Enable::USE_PT){
		cout<<"Cant not set USE_P and USE_PT in the same time!"<<endl;
		gSystem->Exit(-1);
	}
	else if(Enable::USE_P == false && Enable::USE_PT == false){
		cout<<"None of P/PT is set, please set one of them"<<endl;
		gSystem->Exit(-1);
	}
	
//	if(Enable::GENERATE_SEED){
//		size_t findSlash = inputFile.find_last_of("/");
//		string inputFileName = inputFile.substr(findSlash+1, inputFile.size());
//
//		RooRealVar dummyVal("dummy", "", 0);
//		RooUnblindPrecision blindVal("blindVal", "blindVal", inputFileName.c_str(), nEvents, skip+1, dummyVal, kFALSE);
//		rc->set_IntFlag("RANDOMSEED", abs(ceil(blindVal.getVal()*1e2)));
//	}

	//
	//Initialize the selected Input/Event generation
	//
	InputInit();

	//Simple Input generator:
	//if you run more than one of these Input::SIMPLE_NUMBER > 1
	//add the settings for other with [1], next with [2]...
	//
	if(Input::SIMPLE){
		INPUTGENERATOR::SimpleEventGenerator[0]->add_particles(particle_name.c_str(), 1);
		if(Input::HEPMC || Input::EMBED){
			INPUTGENERATOR::SimpleEventGenerator[0]->set_reuse_existing_vertex(true);
			INPUTGENERATOR::SimpleEventGenerator[0]->set_existing_vertex_offset_vector(0., 0., 0.);
		}
		else{
			INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform, 
					PHG4SimpleEventGenerator::Uniform,
					PHG4SimpleEventGenerator::Uniform);
			INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_mean(0., 0., 0.);
			INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_width(0., 0., 0.);
		}
		INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(etamin, etamax);
		INPUTGENERATOR::SimpleEventGenerator[0]->set_phi_range(phimin, phimax);
		
		if(Enable::USE_P)
			INPUTGENERATOR::SimpleEventGenerator[0]->set_p_range(pmin, pmax);
		else if (Enable::USE_PT)
			INPUTGENERATOR::SimpleEventGenerator[0]->set_pt_range(pmin, pmax);

		//Upsilons
		//if you run more than one of these Input::UPSILON_NUMBER > 1
		//add the settings for other with [1], next with [2]...
		//
		if(Input::UPSILON){
			INPUTGENERATOR::VectorMesonGenerator[0]->add_decay_particles("mu", 0);
			INPUTGENERATOR::VectorMesonGenerator[0]->set_rapidity_range(-1, 1);
			INPUTGENERATOR::VectorMesonGenerator[0]->set_pt_range(0., 20.);
			// Y species - select only one, last one wins
			INPUTGENERATOR::VectorMesonGenerator[0]->set_upsilon_1s();

			if(Input::HEPMC || Input::EMBED){
				INPUTGENERATOR::VectorMesonGenerator[0]->set_reuse_existing_vertex(true);
				INPUTGENERATOR::VectorMesonGenerator[0]->set_existing_vertex_offset_vector(0., 0., 0.);
			}
		}
		
		//particle gun
		//if you run more than one of these Input::GUN_NUMBER > 1
		//add the settings for other with [1], next with [2]...

		if(Input::GUN){
			INPUTGENERATOR::Gun[0]->AddParticle(particle_name.c_str(), 0, 1, 0);
			INPUTGENERATOR::Gun[0]->set_vtx(0., 0., 0.);
		}

		if(Input::IONGUN){
			float theta = -25e-3;
			INPUTGENERATOR::IonGun[0]->SetA(197);
			INPUTGENERATOR::IonGun[0]->SetZ(79);
			INPUTGENERATOR::IonGun[0]->SetCharge(79);
			INPUTGENERATOR::IonGun[0]->SetMom(sin(theta)*110*197, 0, cos(theta)*110*197);

			INPUTGENERATOR::IonGun[0]->Print();
		}


		//pythia6
		if(Input::PYTHIA6){
			INPUTGENERATOR::Pythia6->set_config_file(string(getenv("CALIBRATIONROOT")) + "/Generators/phpythia6_ep.cfg");
			// apply EIC beam parameter following EIC CDR
			Input::ApplyEICBeamParameter(INPUTGENERATOR::Pythia6);
		}

		//pythia8
		if(Input::PYTHIA8){
			// apply EIC beam parameter following EIC CDR
			Input::ApplyEICBeamParameter(INPUTGENERATOR::Pythia8);
		}

		//Sartre
		if(Input::SARTRE){
			// apply EIC beam parameter following EIC CDR
			Input::ApplyEICBeamParameter(INPUTGENERATOR::Sartre);
		}

		//
		//Set Input Manager specific option
		//
		if(Input::HEPMC){
			// apply EIC beam parameter following ERIC CDR
			Input::ApplyEICBeamParameter(INPUTMANAGER::HepMCInputManager);
	
		}

		//register all input generators with Fun4All
		InputRegister();

		//Reads event generators in EIC smear files, which is registered in InputRegister
		if(Input::READEIC){
			// apply EIC beam parameter following EIC CDR
			INPUTGENERATOR::EICFileReader->SetFirstEntry(skip);
			Input::ApplyEICBeamParameter(INPUTGENERATOR::EICFileReader);
		}


	}


}


#endif
