#ifndef DATA_OUTPUT_H
#define DATA_OUTPUT_H

// class handling writing of event information and detected photon information to output root file

#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH2Poly.h"

class data_output {

private:
	// output file
	TFile *output_file;
	TFile *original_file;

	// branches:
	TTree *event_tree;
	TTree *data_tree;
	TTree *data_tree_vuv;
	TTree *data_tree_vis;
	TTree *test_tree;

	// tree entries:
	std::vector<std::vector<double>> total_time_vuv;
	std::vector<std::vector<double>> total_time_charge;
	std::vector<double> LightYield;
	std::vector<double> ChargeYield;

	//TH2Poly for output - currently not used
	TH2Poly *h2charge=nullptr;

	TH2Poly *h2light=nullptr;




	// event
	// most of this is not used anymore - can be deleted if this doesnt break the code (shoudlnt)
	int event_no;
	double event_E;
	double event_x_pos, event_y_pos, event_z_pos;
	// data
	int data_pmt, data_pmt_vuv, data_pmt_vis;
	int data_event, data_event_vuv, data_event_vis;
	double data_time, data_time_vuv, data_time_vis;
	double data_x_pos, data_x_pos_vuv, data_x_pos_vis;
	double data_y_pos, data_y_pos_vuv, data_y_pos_vis;
	double data_z_pos, data_z_pos_vuv, data_z_pos_vis;

	const bool include_reflected;

public:

	// constructor
	data_output(const char* output_file_name,const bool include_input, const bool include_timings, const bool include_reflected, const char* original_file_name);

	void add_maps(TH2Poly *h2charge_input, TH2Poly *h2light_input);

	void add_data_till(const int &hit_number, const std::vector<int> &detector_position,const std::vector<int> &num_VUV, const std::vector<TVector3> &ScintPoint, const std::vector<std::vector<double>> &times_vuv);
	void add_data_till(const std::vector<std::vector<double>> &times_vuv);
	void add_data_till(const std::vector<std::vector<double>> &times_vuv, const std::vector<double> &light_yield, const std::vector<double> &chargeyield);
	void add_data_till(const std::vector<std::vector<double>> &times_vuv, const std::vector<std::vector<double>> &times_charge, const std::vector<double> &light_yield, const std::vector<double> &chargeyield);
	void add_data_till(const std::vector<std::vector<double>> &times_vuv, const std::vector<std::vector<double>> &times_charge, const std::vector<double> &light_yield, const std::vector<double> &chargeyield, TH2Poly* light, TH2Poly* charge);


	// destructor
	~data_output();

	// write file
	void write_output_file() { output_file->Write(); }

};
















#endif
