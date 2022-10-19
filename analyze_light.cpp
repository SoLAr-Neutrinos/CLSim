// main driver for scintillation toy mc simulation code
#include<string>
#include<iostream>
#include<fstream>
#include<chrono>
#include <sstream>
#include <vector>
#include <algorithm>
#include "TH1.h"
#include "TH2Poly.h"
#include "TRandom.h"
#include "TVector3.h"
#include "data_output.h"
#include "semi_analytic_hits.h"
#include "time_parameterisation.h"
#include "utility_functions.h"
#include "radiological_parameters.h"
#include <omp.h>
#include "TROOT.h"
#include "TCanvas.h"

// File for flags (https://github.com/zhanxw/Argument)
#include "Argument.h"

// include parameter file
#include "simulation_parameters.h"
using namespace std;

// flag to enable a verbose output
bool debug = false;

int main(int argc, char* argv[]){


	// Enable Flags
	BEGIN_PARAMETER_LIST(pl)
		ADD_PARAMETER_GROUP(pl,"Input")
        	ADD_BOOL_PARAMETER(pl, isHelp, "--help", "Print this Help Message")
        	ADD_INT_PARAMETER(pl, inum, "--number", "Number of events to simulate")
        	ADD_INT_PARAMETER(pl, icores, "--core", "Number of threats to use")
        	ADD_BOOL_PARAMETER(pl, isExcl, "--exclOut", "Exclude the G4 Input from the Output files")
        	ADD_STRING_PARAMETER(pl, isCharge, "--charge", "Enable Charge Simulation")
        	ADD_DOUBLE_PARAMETER(pl, isChargeSize, "--pixSize", "Size of the charge dectector pixel")
		ADD_BOOL_PARAMETER(pl, isDiff, "--diffusion", "Enable Diffusion")
	END_PARAMETER_LIST(pl)
        ;

	pl.Read(argc, argv);


	// Get the command line arguments
	if(FLAG_isHelp || argc == 1){
		cout << "Usage: " << argv[0] << " <G4_input_file> <SiPM_placement_file> <PixelSize(cm)> <OutputFile>" << endl;
		cout << endl;
		cout << "	G4_input_file: The file containing the G4 input" << endl;
		cout << "	SiPM_placement_file: The file containing the SiPM placement" << endl;
		cout << "	PixelSize(cm): The size of the pixel in cm" << endl;
		cout << "	OutputFile: The file to write the output to" << endl;
		cout << endl;
		cout << "	Optional arguments:" << endl;
		cout << "		--number: number of events to simulate - otherwise number in root file" << endl;
		cout << "		--charge: " << " to enable charge simulation - 0 by default" << endl;
		cout << "		--diffusion: " << " If charge is enabled, diffussion is enabled. If it should be disabled set to 0" << endl;
		cout << "		--exclOut: Exclude the G4 input from the output files" << endl;
		cout << "		--number: number of threts to use in multithreded mode" << endl;
		cout << "		--help: Print this help message" << endl;
		return 1;
	}

	// Default: Include the G4 tree to the output file
	bool include_input = true;
	if(FLAG_isExcl) include_input = false;

	bool charge = false;
	bool diffusion = false;
	if(FLAG_isCharge != "" ) {
		charge = true;
		diffusion = true;
		if(FLAG_isDiff) diffusion = false;
	}



	double SiPMSize = stod(FLAG_REMAIN_ARG[2].c_str());
	cout << endl;
	cout << endl;
	cout << "Run Summary: " << endl;
	cout << "\n";
	cout << "Running with G4 file: " << FLAG_REMAIN_ARG[0].c_str() << endl;
	cout << "Number of events to simulate: " << FLAG_inum << endl;
	if(include_input){
		cout << "Including input tree event_tree in output" << endl;
	}
	else{
		cout << "Not Including input tree event_tree in output" << endl;
	}
	cout << "Output File: " << FLAG_REMAIN_ARG[3].c_str() << endl;
	cout << "Light: " << endl;
	cout << "	Running with SiPM file: " << FLAG_REMAIN_ARG[1].c_str() << endl;
	cout << "	Running with a SiPM Size of : " << stod(FLAG_REMAIN_ARG[2].c_str()) << " cm" << endl;




	std::vector<std::vector<int>> chdet_type;
	std::vector<int> chdet_direction;
	std::vector<std::vector<double>> chdet_position;
	double PixSize = 0;
	int number_chdets = 0;
	TH2Poly *chdet_map = nullptr;

	if(charge){
		if(FLAG_isChargeSize != 0){
			PixSize = FLAG_isChargeSize;
		}
		else{
			PixSize = SiPMSize;
		}

		cout << "Charge:" << "\n";
		cout << " 	Running with Pixel File : " << FLAG_isCharge << endl;
		cout << "	Running with a Pixel Size of : " << PixSize << " cm" << endl;

		if(diffusion){
			cout << "	Including diffusion" << endl;
		}
		else{
			cout << "	Not including diffusion" << endl;
		}

		// ------- Read charge detector positions and types --------
		cout << "\n";
		int max_x_ch = 0;
		int max_y_ch = 0;
		std::cout << "Loading Charge Detector positions..." << std::endl;
		std::ifstream charge_detector_positions_file;
		charge_detector_positions_file.open(FLAG_isCharge);
		if(charge_detector_positions_file.is_open()) std::cout << "File opened successfully" << std::endl;
		else {std::cout << "File not found." << std::endl; exit(1);}
		while(!charge_detector_positions_file.eof()) {
			int num_chdet, type_chdet, direction; double x_chdet, y_chdet, z_chdet;
			if(charge_detector_positions_file >> num_chdet >> x_chdet >> y_chdet >> z_chdet >> type_chdet >> direction) {
				if(x_chdet > max_x_ch) max_x_ch = x_chdet;
				if(y_chdet > max_y_ch) max_y_ch = y_chdet;
			    std::vector<int> type({num_chdet, type_chdet});
			    std::vector<double> position({x_chdet, y_chdet, z_chdet});
			    chdet_type.push_back(type);
			    chdet_position.push_back(position);
			    chdet_direction.push_back(direction);
			}
			else{ break; }
		}
		charge_detector_positions_file.close();
		number_chdets = chdet_type.size();
		cout << "Number of charge detectors: " << number_chdets << endl;

	// A little bit unconventional, but this is allows for a very quick sorting of the electrons
	// We create a variable sized historgam which corresponds to the placement of the SiPMs
	// We then sort the electrons later using the bin numbers which (correcting for a shift of 1) correspond to the SiPM number
	   chdet_map = new TH2Poly("chdet_map","chdet_map",0,max_x_ch,0,max_y_ch);
	   for(int i = 0; i < number_chdets; ++i) {
			int binN = chdet_map->AddBin(chdet_position[i][0]-PixSize/2, chdet_position[i][1]-PixSize/2, chdet_position[i][0]+PixSize/2, chdet_position[i][1]+PixSize/2);
			chdet_map->SetBinContent(binN, chdet_direction[i]+i);
			if(binN-1 != i) std::cout << "Bin number mismatch " << i << " " << binN<< std::endl;
	}




	}//end if charge
	cout << endl;
	cout << endl;

	// ------- Read photon detector positions and types --------
	std::vector<std::vector<int>> opdet_type;
	std::vector<int> opdet_direction;
	std::vector<std::vector<double>> opdet_position;
	double max_x = 0;
	double max_y = 0;

       	std::cout << "Loading Photon Detector positions..." << std::endl;
        std::ifstream detector_positions_file;
        detector_positions_file.open(argv[2]);

        if(detector_positions_file.is_open()) std::cout << "File opened successfully" << std::endl;
        else {std::cout << "File not found." << std::endl; exit(1);}
        while(!detector_positions_file.eof()) {
		int num_opdet, type_opdet, direction; double x_opdet, y_opdet, z_opdet;
		if(detector_positions_file >> num_opdet >> x_opdet >> y_opdet >> z_opdet >> type_opdet >> direction) {
		    if(x_opdet > max_x) max_x = x_opdet;
		    if(y_opdet > max_y) max_y = y_opdet;
		    std::vector<int> type({num_opdet, type_opdet});
		    std::vector<double> position({x_opdet, y_opdet, z_opdet});
		    opdet_type.push_back(type);
		    opdet_position.push_back(position);
		    opdet_direction.push_back(direction);
		}
		else{ break; }
        }
	detector_positions_file.close();
	int number_opdets = opdet_type.size();
	// Assumes that the optical detectors are also on the x-y plane. This should be generalised for field cage studies
	TH2Poly *opdet_map = nullptr;

	opdet_map = new TH2Poly("opdet_map", "opdet_map", 0, max_x, 0, max_y);
        for(int i = 0; i < number_opdets; ++i) {
			int binN = opdet_map->AddBin(opdet_position[i][0]-SiPMSize/2, opdet_position[i][1]-SiPMSize/2, opdet_position[i][0]+SiPMSize/2, opdet_position[i][1]+SiPMSize/2);
			opdet_map->SetBinContent(binN, opdet_type[i][1]);
			if(binN-1 != i) std::cout << "Bin number mismatch " << i << " " << binN<< std::endl;
	}

	std::cout << "Positions Loaded: " << number_opdets << " optical detectors." << std::endl << std::endl;



	// Set Seed for random generation - currently just take system time
	gRandom->SetSeed(0);

	// -------- Initialise semi-analytic hits class ---------
	semi_analytic_hits hits_model;
	// Sets the SiPM size for the optical detector! Do not be fooled by the naming ( I know, I know...)
	hits_model.setPixelSize(stod(FLAG_REMAIN_ARG[2].c_str()),stod(FLAG_REMAIN_ARG[2].c_str()));

	// -------- Initialise timing parametrisation class ---------
	time_parameterisation times_model(parameters::timing_discretisation_step_size);

	// -------- Initialise utility/energy spectrum class ---------
	utility_functions utility;
	utility.initalise_scintillation_functions_argon(parameters::t_singlet, parameters::t_triplet, parameters::singlet_fraction_electron, parameters::triplet_fraction_electron,
        parameters::singlet_fraction_alpha, parameters::triplet_fraction_alpha, parameters::scint_time_window);
        utility.initalise_scintillation_functions_xenon(parameters::t_singlet_Xe, parameters::t_triplet_Xe, parameters::singlet_fraction_Xe, parameters::triplet_fraction_Xe, parameters::scint_time_window);

	// ------- Read G4 simulation data --------
	// Read in the results of the qpixg4 results

	char* G4InputFileName = argv[1];
	TFile * G4InputFile = new TFile(G4InputFileName);

	TTree *G4InputTree = (TTree*)G4InputFile->Get("event_tree");

	double energy_deposit;
	vector <double> *hit_start_x = nullptr;
	vector <double> *hit_start_y = nullptr;
	vector <double> *hit_start_z = nullptr;
	vector <double> *hit_start_t = nullptr;
	vector <double> *hit_energy_deposit = nullptr;
	vector <double> *hit_length = nullptr;
  	vector<double> *particle_pdg_code = nullptr;
  	vector<double> *hit_track_id= nullptr;
  	vector<double> *hit_process_key= nullptr;
	double pixel_size;

	// Get the tree branches
	G4InputTree->SetBranchAddress("energy_deposit", &energy_deposit);
	G4InputTree->SetBranchAddress("hit_start_x", &hit_start_x);
	G4InputTree->SetBranchAddress("hit_start_y", &hit_start_y);
	G4InputTree->SetBranchAddress("hit_start_z", &hit_start_z);
	G4InputTree->SetBranchAddress("hit_length", &hit_length);
	G4InputTree->SetBranchAddress("hit_start_t", &hit_start_t);
	G4InputTree->SetBranchAddress("hit_energy_deposit", &hit_energy_deposit);
	G4InputTree->SetBranchAddress("hit_length", &hit_length);
  	G4InputTree->SetBranchAddress("hit_track_id", &hit_track_id);
  	G4InputTree->SetBranchAddress("hit_process_key", &hit_process_key);
  	G4InputTree->SetBranchAddress("particle_pdg_code", &particle_pdg_code);

	int NEventsToLoopOver = G4InputTree->GetEntries();
	if(FLAG_inum){
		NEventsToLoopOver = FLAG_inum;
	}

	if(FLAG_icores){
		omp_set_num_threads(FLAG_icores);
	}
	else{
		omp_set_num_threads(1);
	}
        cout << "Number of events to loop over: " << NEventsToLoopOver << endl;

	data_output output_file(FLAG_REMAIN_ARG[3].c_str(), include_input, parameters::include_timings, parameters::include_reflected, G4InputFileName );

	output_file.add_maps(chdet_map, opdet_map);

        std::cout << "Starting event loop" << std::endl;
	for (int EventIt=0; EventIt < NEventsToLoopOver; EventIt++)
	{
		std::cout << " --- > Event: " << EventIt << std::endl;
		G4InputTree->GetEntry(EventIt);


		//Empty chdet_map bins
		for(int i = 0; i < number_opdets; ++i) {
			opdet_map->SetBinContent(i, 0);
			chdet_map->SetBinContent(i, 0);
		}

		int max_events = hit_start_x->size();

		// Vector inlcuding all the hit positions in (x,y,z)
		std::vector<std::vector<double>> position_list(max_events, std::vector<double>(3,-999.9));

		// Vector including the current hits x,y,z positoin
		vector<TVector3> ScintPoint_array;

		// Vector of vectors.
		// [ SiPM1[t0, t1, t2,t3,...], SiPM2[t0, t1, t2,t3,...], SiPM3[t0, t1, t2,t3,...], ...]
		// For each SiPM there should be a vector. This vector will contain the time of each photon.
		//vector<std::unique_ptr<vector<double>>> total_time_vuv_array ;//(number_opdets, vector<double>);
		vector<vector<double>> total_time_vuv_array = {} ;//(number_opdets, vector<double>);
		total_time_vuv_array.clear();
		// Vector of vectors for the charge readout
		vector<vector<double>> total_time_charge;//(number_opdets, vector<double>);

		total_time_charge.clear();
		for(int i=0; i<number_opdets; i++){
			vector<double> temp;
			total_time_vuv_array.push_back(temp);
		}
		for(int i = 0; i< number_chdets; i++){
			vector<double> temp;
			total_time_charge.push_back(temp);
		}

		vector<vector<double>> op_channel_pos(number_opdets, vector<double>(3,0.0));

		// Determining the number of photons due to each hit
		std::vector<double> lightyield(hit_start_x->size());
		std::vector<double> chargeyield(hit_start_x->size());
		std::vector<double> numPhotons(hit_start_x->size());
		std::vector<double> numElectrons(hit_start_x->size());
		std::vector<TVector3> electronStartingPoints;
		std::vector<std::vector<double>> photonStartingPoints;

		std::cout << "Preparing the Energy Depositions" << std::endl;
		if(charge) std::cout << "Simulating Charge" << std::endl;

		// This loop determines the number of electrons and the light yield for each hit
		for (unsigned int i=0; i < hit_start_x->size(); i++){
				// Stores the starting points for all the electrons for hit i
				// Stores the light/charge yield for hit i
				// Check if an ionisation event
				double light_yield = hits_model.LArQL(hit_energy_deposit->at(i), hit_length->at(i), 0.5);
				lightyield[i] = light_yield;

				int number_photons = light_yield * hit_energy_deposit->at(i);
				numPhotons[i]=gRandom->Poisson(number_photons);
				double pdg = particle_pdg_code->at(hit_track_id->at(i)-1);
				photonStartingPoints.push_back({hit_start_x->at(i), hit_start_y->at(i), hit_start_z->at(i), hit_start_t->at(i), pdg, numPhotons[i]});

				chargeyield[i] = hits_model.LArQQ(hit_energy_deposit->at(i), hit_length->at(i), 0.5);
				numElectrons[i] = chargeyield[i] * hit_energy_deposit->at(i);

				if(charge){
					double z_pos = hit_start_z->at(i);
					double x_pos = hit_start_x->at(i);
					double y_pos = hit_start_y->at(i);
					double drift_time = z_pos/parameters::drift_velocity; // should be in micro seconds (drift_velcoity in us, 10e6 us = 1s)
					//convert to s
					drift_time = drift_time * 1e-6; // in seconds
					// Get the dispersion area
					double area_trans = parameters::drift_trans * drift_time; // this should be in seconds
					double area_long = parameters::drift_long * drift_time;
					// get the radius of the areas of dispersion
					double rad_trans = utility.Area2Radius(area_trans);
					double rad_long = utility.Area2Radius(area_long);

					double x_pos_final = x_pos;
					double y_pos_final = y_pos;
					double z_pos_final = z_pos;
					drift_time = drift_time * 1e6; // convert to us
					// Loop over the number of electrons and push their starting values to the hit i vector
					// Also disperse the electrons if requried
					for(int j=0; j<numElectrons[i]; j++){
						if(diffusion){
							double rand_angle = gRandom->Uniform(0,2*TMath::Pi());
							double radius = abs(gRandom->Gaus(0, rad_trans));
							x_pos_final = x_pos + radius * cos(rand_angle);
							y_pos_final = y_pos + radius * sin(rand_angle);
							z_pos_final = z_pos + gRandom->Gaus(0,rad_long);
							//Make sure that we do not create electrons on/behind the surface of the detector
							int it =0;
							while(z_pos_final  <= 0){
								z_pos_final = z_pos + gRandom->Gaus(0,rad_long/(2*it));
								it++;
							}
							drift_time = z_pos_final / parameters::drift_velocity ; // this should be in us now - think about units again
						}// End diffusion if
						// Find the bin of the TH2Poly that corresponds to the corresponding x,y position
						// The -1 is there to account for the numbering of the SiPM-IDs to start at 0, while bin numbers start at 1
						int bin = chdet_map->FindBin(x_pos_final,y_pos_final) - 1;

						if(bin < 0 || bin > number_chdets){
							continue;
						}
						total_time_charge[bin].push_back(hit_start_t->at(i)*0.001+drift_time);
					}
				} // End charge if
		}// End loop over hits

		// Go through every SiPM
		std::cout << "Simulating Light" << std::endl;


		for(int op_channel = 0; op_channel < number_opdets; op_channel++) {
			// For each photon get the time info
			std::vector<double> photon_time_vuv;
			for(auto const &photonStartingPoint : photonStartingPoints){

				TVector3 photonStartingPoint_v(photonStartingPoint[0], photonStartingPoint[1], photonStartingPoint[2]);
				double photonStartingPoint_t = photonStartingPoint[3]*0.001;
				double pdg = (int)photonStartingPoint[4];
				int number_photons = photonStartingPoint[5];

				int op_channel_type = opdet_type[op_channel][1];
				int op_direction = opdet_direction[op_channel];
				TVector3 OpDetPoint(opdet_position[op_channel][0],opdet_position[op_channel][1],opdet_position[op_channel][2]);

				// Set the ratio of fast/slow component light according to the PDG type
				// For unknown PDG types, assume electronic distribution
				if(pdg==22 || abs(pdg)==11 || abs(pdg)==12 || abs(pdg)==13 || abs(pdg)==14){
					// Gamma, e, nue, mu, numu
					double singlet_fraction = parameters::singlet_fraction_electron;
					double triplet_fraction = parameters::triplet_fraction_electron;
				}
				else if(pdg == 2212 || pdg == 2112 || pdg >= 1000000000){
					// Proton, neutron, all nuclei
					double singlet_fraction = parameters::singlet_fraction_alpha;
					double triplet_fraction = parameters::triplet_fraction_alpha;
				}
				else {
					double singlet_fraction = parameters::singlet_fraction_alpha;
					double triplet_fraction = parameters::triplet_fraction_alpha;
				}

				// gets the number of photon in the current detector
				double num_VUV;
				num_VUV = hits_model.VUVHits(number_photons, photonStartingPoint_v, OpDetPoint, op_channel_type, 0, op_direction);
				if(num_VUV < 1) continue;


				// Calculate the angle between the scinitllation point and the optical detector
				double distance_to_pmt = (OpDetPoint-photonStartingPoint_v).Mag();
				double cosine = -999.0;

				// Calculate the angle between the scinitllation point and the optical detector, this is dependent on the orientation the optical detector looks at
				if(op_direction == 1)  cosine = sqrt(pow(photonStartingPoint_v[0] - OpDetPoint[0],2)) / distance_to_pmt;
				else if(op_direction == 2)  cosine = sqrt(pow(photonStartingPoint_v[1] - OpDetPoint[1],2)) / distance_to_pmt;
				else if(op_direction == 3)  cosine = sqrt(pow(photonStartingPoint_v[2] - OpDetPoint[2],2)) / distance_to_pmt;
				else { std::cout << "Error: Optical detector direction not defined." << std::endl; }

				assert(cosine>=-1 && cosine<=1);

				double theta = acos(cosine)*180./3.14159;

				// Due to rounding errors, this can return a value like 90.000x01 which will throw an error in the angle_bin histogram.
				// Cap the angle to <90 degrees.
				if (theta >= 90){ theta=89.99999; }
				int angle_bin = theta/45;       // 45 deg bins


				// For each photon we will store its arival time in the SiPM currently being looped over.
				// Returns the transport time that the photon takes from the sicntillation point to the detector
				// Returns the time that the photon arrives at the detector in NANOSECONDS
				std::vector<double> transport_time_vuv = times_model.getVUVTime(distance_to_pmt, angle_bin, num_VUV);

				for(auto& x: transport_time_vuv) {
					// emission time similar to above
					double emission_time;
					if(pdg==22 || abs(pdg)==11 || abs(pdg)==12 || abs(pdg)==13 || abs(pdg)==14){
						emission_time = utility.get_scintillation_time_electron()*1000000.0; // in us
					}
					else if(pdg == 2212 || pdg == 2112 || pdg >= 1000000000){
						emission_time = utility.get_scintillation_time_alpha()*1000000.0; // in us
					}
					else{
					        emission_time = utility.get_scintillation_time_electron()*1000000.0; // in us
					}
					double total_time = photonStartingPoint_t+(x*0.001 + emission_time); // in microseconds
				        total_time_vuv_array[op_channel].push_back(total_time);
				}// End for transporttime
			} // End of photon starting point loop = loop over steps
	    } // End of opdet loop
	output_file.add_data_till(total_time_vuv_array,total_time_charge , lightyield, chargeyield);
	total_time_vuv_array.clear();
	lightyield.clear();
	chargeyield.clear();
	} // end event loop

     //Write to OUTPUT FILE
     output_file.write_output_file();
     sleep(2);
     std::cout << "Program finished." << std::endl;
}
