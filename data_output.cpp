#include "data_output.h"

#include <iostream>

#include "TString.h"
#include "TVector3.h"


// constructor
data_output::data_output(const char* output_file_name, const bool include_input, const bool include_timings, const bool include_reflected, const char* original_file_name) : include_reflected{include_reflected} {

	// Get the original Input File and catch the original tree
	// create file
	output_file = new TFile(output_file_name, "RECREATE", "Output File");

	if (output_file->IsOpen()) {
		std::cout << "Output file created successfully." << std::endl << std::endl;
	}
	else {
		std::cout << "Output file could not be opened, check file name is valid." << std::endl << std::endl; exit(1);
	}

	// Include input trees in to the new output file
	if(include_input){
		std::cout << "Copying the input tree to the output file." << std::endl;
		std::cout << " May take a lot of disk space and a while - to disable this feature, set include_input to false." << std::endl;
		TFile *InputFile = new TFile(original_file_name);
		TTree *InputTree = (TTree*)InputFile->Get("event_tree");
		output_file->cd();
		InputTree->CloneTree()->Write("event_tree");
		// Get rid of the original file as copies is done
		InputFile->Close();
	}




	// create the acutall output trees we need
	test_tree = new TTree("ScintSim_tree", "Scintillation Light Simulation");
	test_tree->Branch("total_time_vuv", &total_time_vuv); //A vector of vectors ([[PMT1: t1,t2,t3,...],[PMT2: t1,t2,t3,...],[PMT3: ...],[PMT4: ...] ... ])
	test_tree->Branch("total_time_charge", &total_time_charge); //A vector of vectors ([[PMT1: t1,t2,t3,...],[PMT2: t1,t2,t3,...],[PMT3: ...],[PMT4: ...] ... ])
	test_tree->Branch("LightYield", &LightYield);
	test_tree->Branch("ChargeYield", &ChargeYield);



}

// destructor
data_output::~data_output(){
	// deleting output_file also deletes all trees properly
	delete output_file;
}

void data_output::add_maps(TH2Poly *h2charge_input, TH2Poly *h2light_input){
	output_file->cd();
	h2charge = (TH2Poly*)h2charge_input->Clone("charge_detector_map");
	h2light = (TH2Poly*)h2light_input->Clone("light_detector_map");
	h2light->Write();
	h2charge->Write();
}


void data_output::add_data_till(const std::vector<std::vector<double>> &times_vuv, const std::vector<double> &light_yield, const std::vector<double> &chargeyield){
    total_time_vuv = times_vuv;
    LightYield = light_yield;
    ChargeYield = chargeyield;
    test_tree->Fill();
}
void data_output::add_data_till(const std::vector<std::vector<double>> &times_vuv, const std::vector<std::vector<double>> &times_charge,const std::vector<double> &light_yield, const std::vector<double> &chargeyield){
    total_time_vuv = times_vuv;
    total_time_charge = times_charge;
    LightYield = light_yield;
    ChargeYield = chargeyield;
    test_tree->Fill();
}

void data_output::add_data_till(const std::vector<std::vector<double>> &times_vuv, const std::vector<std::vector<double>> &times_charge,const std::vector<double> &light_yield, const std::vector<double> &chargeyield, TH2Poly *light, TH2Poly *charge){
    TH2Poly* light_temp = (TH2Poly*)light->Clone();
    TH2Poly* charge_temp = (TH2Poly*)charge->Clone();
    total_time_vuv = times_vuv;
    total_time_charge = times_charge;
    LightYield = light_yield;
    ChargeYield = chargeyield;
    h2light = light_temp;
    h2charge = charge_temp;
    test_tree->Fill();
}
