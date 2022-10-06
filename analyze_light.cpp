// main driver for scintillation toy mc simulation code
#include<string>
#include<iostream>
#include<fstream>
#include<chrono>
#include <sstream>
#include <vector>
#include <algorithm>
#include "TH1.h"
#include "TRandom.h"
#include "TVector3.h"
#include "data_output.h"
#include "semi_analytic_hits.h"
#include "time_parameterisation.h"
#include "utility_functions.h"
#include "radiological_parameters.h"

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
        	ADD_BOOL_PARAMETER(pl, isExcl, "--exclOut", "Exclude the G4 Input from the Output files")
        	ADD_BOOL_PARAMETER(pl, isCharge, "--charge", "Enable Charge Simulation")
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
		cout << "		--help: Print this help message" << endl;
		return 1;
	}

	// Default: Include the G4 tree to the output file
	bool include_input = true;
	if(FLAG_isExcl) include_input = false;

	bool charge = false;
	bool diffusion = false;
	if(FLAG_isCharge) {
		charge = true;
		diffusion = true;
		if(FLAG_isDiff) diffusion = false;
	}



	cout << endl;
	cout << endl;
	cout << "Running with G4 file: " << FLAG_REMAIN_ARG[0].c_str() << endl;
	cout << "Running with SiPM file: " << FLAG_REMAIN_ARG[1].c_str() << endl;
	cout << "Running with a Pixel Size of : " << stod(FLAG_REMAIN_ARG[2].c_str()) << " cm" << endl;
	cout << "Output File: " << FLAG_REMAIN_ARG[3].c_str() << endl;
	cout << "Number of events to simulate: " << FLAG_inum << endl;

	cout << endl;


	if(include_input){
		cout << " - Including input tree event_tree in output" << endl;
	}
	else{
		cout << " - Not Including input tree event_tree in output" << endl;
	}

	if(charge){
		cout << " - Including charge simulation" ;
		if(diffusion){
			cout << " including diffusion" << endl;
		}
		else{
			cout << " not including diffusion" << endl;
		}
	}
	cout << endl;


	// Set Seed for random generation - currently just take system time
	gRandom->SetSeed(0);

	// -------- Initialise semi-analytic hits class ---------
	semi_analytic_hits hits_model;
	hits_model.setPixelSize(stod(FLAG_REMAIN_ARG[2].c_str()),stod(FLAG_REMAIN_ARG[2].c_str()));
	double PixSize = stod(FLAG_REMAIN_ARG[2].c_str());

	// -------- Initialise timing parametrisation class ---------
	time_parameterisation times_model(parameters::timing_discretisation_step_size);

	// -------- Initialise utility/energy spectrum class ---------
	utility_functions utility;
	utility.initalise_scintillation_functions_argon(parameters::t_singlet, parameters::t_triplet, parameters::singlet_fraction_electron, parameters::triplet_fraction_electron,
        parameters::singlet_fraction_alpha, parameters::triplet_fraction_alpha, parameters::scint_time_window);
        utility.initalise_scintillation_functions_xenon(parameters::t_singlet_Xe, parameters::t_triplet_Xe, parameters::singlet_fraction_Xe, parameters::triplet_fraction_Xe, parameters::scint_time_window);

	// ------- Read photon detector positions and types --------
	std::vector<std::vector<int>> opdet_type;
	std::vector<int> opdet_direction;
	std::vector<std::vector<double>> opdet_position;

       	std::cout << "Loading Photon Detector positions..." << std::endl;
        std::ifstream detector_positions_file;
        detector_positions_file.open(argv[2]);

        if(detector_positions_file.is_open()) std::cout << "File opened successfully" << std::endl;
        else {std::cout << "File not found." << std::endl; exit(1);}
        while(!detector_positions_file.eof()) {
		int num_opdet, type_opdet, direction; double x_opdet, y_opdet, z_opdet;
		if(detector_positions_file >> num_opdet >> x_opdet >> y_opdet >> z_opdet >> type_opdet >> direction) {
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
	std::cout << "Positions Loaded: " << number_opdets << " optical detectors." << std::endl << std::endl;

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
  	G4InputTree->SetBranchAddress("particle_pdg_code", &particle_pdg_code);

	int NEventsToLoopOver = G4InputTree->GetEntries();
	if(FLAG_inum){
		NEventsToLoopOver = FLAG_inum;
	}

	cout << "Running over " << NEventsToLoopOver << " events" << endl;
	data_output output_file(FLAG_REMAIN_ARG[3].c_str(), include_input, parameters::include_timings, parameters::include_reflected, G4InputFileName );



	for (int EventIt=0; EventIt < NEventsToLoopOver; EventIt++)
	{
		//if(EventIt!=1894) continue;
		std::cout << " --- > Event: " << EventIt << std::endl;
		G4InputTree->GetEntry(EventIt);

		int max_events = hit_start_x->size();

		// Vector inlcuding all the hit positions in (x,y,z)
		std::vector<std::vector<double>> position_list(max_events, std::vector<double>(3,-999.9));

		// Vector including the number of photons per SiPM
		vector<int> num_VUV_array;

		// Vector including the current hits x,y,z positoin
		vector<TVector3> ScintPoint_array;

		// Vector of vectors.
		// [ SiPM1[t0, t1, t2,t3,...], SiPM2[t0, t1, t2,t3,...], SiPM3[t0, t1, t2,t3,...], ...]
		// For each SiPM there should be a vector. This vector will contain the time of each photon.
		vector<vector<double>> total_time_vuv_array ;//(number_opdets, vector<double>);
		total_time_vuv_array.clear();
		// Vector of vectors for the charge readout
		vector<vector<double>> total_time_charge;//(number_opdets, vector<double>);
		total_time_charge.clear();

		vector<vector<double>> op_channel_pos(number_opdets, vector<double>(3,0.0));
		std::cout << "Number of optical detectors: " << number_opdets << std::endl;

		// Determining the number of photons due to each hit
		std::vector<double> lightyield(hit_start_x->size());
		std::vector<double> chargeyield(hit_start_x->size());
		std::vector<double> numPhotons(hit_start_x->size());
		std::vector<double> numElectrons(hit_start_x->size());
		std::vector<std::vector<TVector3>> electronStartingPoints;

		std::cout << "Preparing the Energy Depositions" << std::endl;

		// This loop determines the number of electrons and the light yield for each hit

		for (int i=0; i < hit_start_x->size(); i++){
				// Stores the starting points for all the electrons for hit i
			        std::vector<TVector3> electronStartingPoints_i;
				// Stores the light/charge yield for hit i
				lightyield[i] = hits_model.LArQL(hit_energy_deposit->at(i), hit_length->at(i), 0.5);
				chargeyield[i] = hits_model.LArQQ(hit_energy_deposit->at(i), hit_length->at(i), 0.5);
				// Stores the number of photons/electrons for hit i
				numPhotons[i] = lightyield[i] * hit_energy_deposit->at(i);
				numElectrons[i] = chargeyield[i] * hit_energy_deposit->at(i);
				if(charge){
					double z_pos = hit_start_z->at(i);
					double x_pos = hit_start_x->at(i);
					double y_pos = hit_start_y->at(i);
					double drift_time = z_pos/parameters::drift_velocity;
					// Get the dispersion area
					double area_trans = parameters::drift_trans * drift_time/(1000*1000); // Think about units again
					double area_long = parameters::drift_long * drift_time/(1000*1000);
					// get the radius of the areas of dispersion
					double rad_trans = utility.Area2Radius(area_trans);
					double rad_long = utility.Area2Radius(area_long);

					double x_pos_final = x_pos;
					double y_pos_final = y_pos;
					double z_pos_final = z_pos;
					// Loop over the number of electrons and push their starting values to the hit i vector
					// Also disperse the electrons if requried
					for(int j=0; j<numElectrons[i]; j++){
						if(diffusion){
							double rand_angle = gRandom->Uniform(0,2*TMath::Pi());
							x_pos_final = x_pos + abs(gRandom->Gaus(0, rad_trans)) * cos(rand_angle);
							y_pos_final = y_pos + abs(gRandom->Gaus(0, rad_trans)) * sin(rand_angle);
							z_pos_final = z_pos + gRandom->Gaus(0,rad_long);
							drift_time = z_pos_final / parameters::drift_velocity ;
						}// End diffusion if
						// Push the starting point of the electron to the vector
						electronStartingPoints_i.push_back(TVector3(x_pos_final, y_pos_final, z_pos_final));
					}
					// Collection of all the starting points of all hits
					// [[ hit1[electron1(x,y,z), electron2(x,y,z), electron3(x,y,z), ...], hit2[electron1, electron2, electron3, ...], ...]
					electronStartingPoints.push_back(electronStartingPoints_i);
				}
		}


		std::cout << "Simulating Charge" << std::endl;
		// Here we loop over all the optical detectors and detectors.
		// We then loop over all the hits and determine if the hit is within/above the detector.
		for(int op_channel = 0; op_channel < number_opdets; op_channel++) {
			if(op_channel % (number_opdets/10) == 0)
				std::cout << "op_channel: " << op_channel
				   << " of " << number_opdets
			           << " (" << (double)op_channel*100./(double)number_opdets <<"%)"<< endl;

			std::vector<double> time_charge;
			// Position of the optical detector
			TVector3 OpDetPoint(opdet_position[op_channel][0],opdet_position[op_channel][1],opdet_position[op_channel][2]);
			for(int HitIt = 0; HitIt < hit_start_x->size(); HitIt++){
				// Get all the starting points of the electrons for hit HitIt
				std::vector<TVector3> electronStartingPoints_i = electronStartingPoints[HitIt];
				// Loop over all the electrons for hit HitIt
				for(int eIt = 0; eIt <numElectrons[HitIt] ; eIt++){
					TVector3 electronStartingPoint = electronStartingPoints_i[eIt];
					double x_pos = electronStartingPoint.X();
					double y_pos = electronStartingPoint.Y();
					double z_pos = electronStartingPoint.Z();
					// Project down on the the z=0 plane
					double drift_time = z_pos/parameters::drift_velocity;
					// Currently hard coded 10cm detectors. Need to change this to be more general using the PixSize
					if (abs(x_pos - OpDetPoint.X()) < 5 && abs(y_pos - OpDetPoint.Y()) < 5 ) {
						time_charge.push_back(drift_time);
					}
				}// end for of number electrons
			}// end for of hits
			total_time_charge.push_back(time_charge);
		}// end for of opdet

		// Get total number of electrons
		int total_num_electrons = 0;
		for(int i=0; i<numElectrons.size(); i++){
			total_num_electrons += numElectrons[i];
		}

		cout << "Total number of electrons: " << total_num_electrons << endl;
		cout << total_time_charge.size() << endl;
		for(int i=0; i<total_time_charge.size(); i++){
			cout << total_time_charge[i].size() << endl;
		}

		// Go through every SiPM
		std::cout << "Simulating Light" << std::endl;
		for(int op_channel = 0; op_channel < number_opdets; op_channel++) {
			if(op_channel % (number_opdets/10) == 0)
				std::cout << "op_channel: " << op_channel
				   << " of " << number_opdets
			           << " (" << (double)op_channel*100./(double)number_opdets <<"%)"<< endl;

			// get optical detector type - rectangular or disk aperture
			// Currently only rectangular apertures are supported
			int op_channel_type = opdet_type[op_channel][1];
			// get optical detector direction - along (x, y, z) axis facing
			int op_direction = opdet_direction[op_channel];


			// get detection channel coordinates (in cm)
			TVector3 OpDetPoint(opdet_position[op_channel][0],opdet_position[op_channel][1],opdet_position[op_channel][2]);
			vector<double> op_det_pos = {OpDetPoint.X(), OpDetPoint.Y(), OpDetPoint.Z()};

			// Final results vector
			std::vector<double> total_time_vuv;

			// Go through every hit in the event
			int num_hits = hit_start_x->size();
			for (int i=0; i < num_hits; i++)
			{
		   		position_list[i][0] = hit_start_x->at(i);
		    		position_list[i][1] = hit_start_y->at(i);
		    		position_list[i][2] = hit_start_z->at(i);

		    		double time_hit = hit_start_t->at(i); // in ns
				time_hit*=0.001; // in us

				// Light yield for the hit currently looked at
				double light_yield = lightyield[i];

				unsigned int number_photons;
		    		number_photons = gRandom->Poisson(parameters::total_QE*light_yield*hit_energy_deposit->at(i));

				assert(number_photons >= 0);

				// Set the ratio of fast/slow component light according to the PDG type
				// For unknown PDG types, assume electronic distribution
				int pdg = particle_pdg_code->at(hit_track_id->at(i)-1);
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



		    		// Setting the point of energy deposition
		    		TVector3 ScintPoint(position_list[i][0],position_list[i][1],position_list[i][2]);
				if(debug){
					cout << " ------------------------------------------ " << endl;
					cout << "ScintPoint: " << ScintPoint.X() << " " << ScintPoint.Y() << " " << ScintPoint.Z() << endl;
					cout << "OpDetPoint: " << OpDetPoint.X() << " " << OpDetPoint.Y() << " " << OpDetPoint.Z() << endl;
				}

				// gets the number of photon in the current detector
				int num_VUV = hits_model.VUVHits(number_photons, ScintPoint, OpDetPoint, op_channel_type, 0, op_direction);
				if(num_VUV== 0) { continue; } // forces the next iteration

				// Calculate the angle between the scinitllation point and the optical detector
				double distance_to_pmt = (OpDetPoint-ScintPoint).Mag();
				double cosine = -999.0;


				// Calculate the angle between the scinitllation point and the optical detector, this is dependent on the orientation the optical detector looks at
				if(op_direction == 1)  cosine = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2)) / distance_to_pmt;
				else if(op_direction == 2)  cosine = sqrt(pow(ScintPoint[1] - OpDetPoint[1],2)) / distance_to_pmt;
				else if(op_direction == 3)  cosine = sqrt(pow(ScintPoint[2] - OpDetPoint[2],2)) / distance_to_pmt;
				else { std::cout << "Error: Optical detector direction not defined." << std::endl; }

				assert(cosine>=-1 && cosine<=1);

				double theta = acos(cosine)*180./3.14159;
				if(debug){
					cout << "Distance to PMT: " << distance_to_pmt << endl;
					cout << "Theta: " << theta << endl;
				}

				// Due to rounding errors, this can return a value like 90.000x01 which will throw an error in the angle_bin histogram.
				// Cap the angle to <90 degrees.
				if (theta >= 90){ theta=89.99999; }
				int angle_bin = theta/45;       // 45 deg bins


				// For each photon we will store its arival time in the SiPM currently being looped over.
				// Returns the transport time that the photon takes from the sicntillation point to the detector
				// Returns the time that the photon arrives at the detector in NANOSECONDS
				std::vector<double> transport_time_vuv = times_model.getVUVTime(distance_to_pmt, angle_bin, num_VUV);

				// For each photon get the time info
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
					double total_time = time_hit+(x*0.001 + emission_time + 2.5*0.001); // in microseconds // WLS 2.5 ns - constant offset
					total_time_vuv.push_back(total_time);
				}// End for transporttime
		    } // end of hit
		total_time_vuv_array.push_back(total_time_vuv);
		} // End of SiPM Loop
	output_file.add_data_till(total_time_vuv_array,total_time_charge , lightyield, chargeyield);
	lightyield.clear();
	chargeyield.clear();
	} // end event loop

     //Write to OUTPUT FILE
     output_file.write_output_file();
     sleep(2);
     std::cout << "Program finished." << std::endl;
}
