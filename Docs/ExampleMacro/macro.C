// This is an overview of some operations that we perfom on the output files from the light sim
// After going through this, you should be able to perform basically any study you want using
// the light and charge outputs.
// Many of these operations you will perform very often, in different macros. So, it is a good
// idea to create a small library of functions that you can use in your macros!
// If this is your first root-macro however, this is probably not the place to start.


#include <iostream>
#include <string>
#include <fstream>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCut.h"
#include <TRandom3.h>
#include <vector>
#include "TGraph.h"
#include "TLegend.h"
#include "TMath.h"
#include "TSpectrum.h"
#include "TH2Poly.h"
#include <iterator>
#include <map>

using std::cout, std::endl, std::string, std::vector;


void macro(string CLInputFile="nuis_out", string InputDirectory="", string OutputDirectory="" )
{

//Making plots nicer
gROOT->ForceStyle(1); gStyle->SetPadTopMargin(0.07); gStyle->SetPadRightMargin(0.05); gStyle->SetPadLeftMargin(0.15); gStyle->SetPadBottomMargin(0.16); gStyle->SetLabelSize(0.06,"xyz"); gStyle->SetTitleSize(0.06,"xyz"); gStyle->SetTitleOffset(0.9,"x"); gStyle->SetTitleOffset(1.1,"y"); gStyle->SetTitleOffset(0.9,"z"); gStyle->SetStatX(0.8); gStyle->SetStatW(0.2); gStyle->SetStatY(0.85); gStyle->SetStatH(0.1); gStyle->SetOptStat(0); gStyle->SetHistLineWidth(3); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1); gStyle->SetPadGridX(kTRUE); gStyle->SetPadGridY(kTRUE);
gStyle->SetPadRightMargin(0.15);

//Defining output path for PNG and output file name for ROOT file
string PNGOutputFolderName, ROOTOutputFileName;
ROOTOutputFileName=Form("%s/ROOT/NEUTFile_%s.root", OutputDirectory.c_str(), CLInputFile.c_str());
PNGOutputFolderName=Form("%s/PNGs/%s", OutputDirectory.c_str(), CLInputFile.c_str());

char* CLFileName;
CLFileName = Form("%s/%s.root", InputDirectory.c_str(), CLInputFile.c_str());
TFile *CLFile = new TFile(CLFileName);

TTree *CLTree = (TTree*)CLFile->Get("ScintSim_tree");
TTree *G4Tree = (TTree*)CLFile->Get("event_tree");

// These are vectors of vectors, where each line corresponds to a single detection element
// Inside, the arrivial times of the phtoons/electrons gitting this element are stored as doubles.
//  [ SiPM1[t0, t1, t2,t3,...], SiPM2[t0, t1, t2,t3,...], SiPM3[t0, t1, t2,t3,...], ...]
vector<vector<double>> *total_time_vuv = nullptr;
vector<vector<double>> *total_time_charge = nullptr;

// These store the light and charge yield per hit used to generate the photons and electons
vector<double> *light_yield = 0;
vector<double> *charge_yield = 0;

// These are the placement maps for the Charge and Light detectors
TH2Poly *ChargeReadout = (TH2Poly*)CLFile->Get("charge_detector_map");
TH2Poly *LightReadout = (TH2Poly*)CLFile->Get("light_detector_map");

CLTree->SetBranchAddress("total_time_vuv", &total_time_vuv);
CLTree->SetBranchAddress("total_time_charge", &total_time_charge);
CLTree->SetBranchAddress("LightYield", &light_yield);
CLTree->SetBranchAddress("ChargeYield", &charge_yield);


// A small selection of the variables stored in the G4 part of the file
vector<double> *generator_initial_particle_energy = nullptr;
vector<double> *generator_inital_particle_x = nullptr;
vector<double> *generator_inital_particle_y = nullptr;
vector<double> *generator_inital_particle_z = nullptr;
vector<double> *hit_length = nullptr;
vector<double> *hit_energy_deposit = nullptr;
int event;

G4Tree->SetBranchAddress("generator_initial_particle_energy", &generator_initial_particle_energy);
G4Tree->SetBranchAddress("generator_initial_particle_x", &generator_inital_particle_x);
G4Tree->SetBranchAddress("generator_initial_particle_y", &generator_inital_particle_y);
G4Tree->SetBranchAddress("generator_initial_particle_z", &generator_inital_particle_z);
G4Tree->SetBranchAddress("hit_length", &hit_length);
G4Tree->SetBranchAddress("hit_energy_deposit", &hit_energy_deposit);
G4Tree->SetBranchAddress("event", &event);

int NEventsToLoopOver = 10;//CLTree->GetEntries(); // 100000
cout << "Looping over " << NEventsToLoopOver << " events..." << endl;
cout << endl;

// Get the total number of photons in the event
TH1D *totalPulseShapeLight = new TH1D("totalPulseShapeLight", "totalPulseShapeLight", 100, 0, 3); // in micro seconds
TH1D *totalPulseShapeCharge = new TH1D("totalPulseShapeCharge", "totalPulseShapeCharge", 100, 0, 600); // in micro seconds
//Looking at a specific light/charge detector
TH1D *pulseShapeLightDetector = new TH1D("pulseShapeLightDetector", "pulseShapeLightDetector", 100, 0, 3);
TH1D *pulseShapeChargeDetector = new TH1D("pulseShapeChargeDetector", "pulseShapeChargeDetector", 100, 0, 600);
//Get the total readout for the event in a TH2 - here we assume a 100x100 cm Detector with the SoLAr tiles
TH2D *totalReadoutLight = new TH2D("totalReadoutLight", "totalReadoutLight", 100, 0, 100, 100, 0, 100);
TH2D *totalReadoutCharge = new TH2D("totalReadoutCharge", "totalReadoutCharge", 300, 0, 100, 300, 0, 100);


// We need to create a lookup table for the detector positions
// For this we create two dictionaries - one for the x and one for the y positions
// CDbins?[a] returns the x/y position of the center of the bin a

CLTree->GetEntry(0);
map<int, double> CDbinsX;
map<int, double> CDbinsY;
TList * list = (TList*)ChargeReadout->GetBins();
TIter next(list);
TH2PolyBin *bin;
while((bin = (TH2PolyBin*)next())){
  CDbinsX.insert(pair<int, double>(bin->GetBinNumber(), (bin->GetXMax()+bin->GetXMin())/2));
  CDbinsY.insert(pair<int, double>(bin->GetBinNumber(), (bin->GetYMax()+bin->GetYMin())/2));
}

map<int, double> LDbinsX;
map<int, double> LDbinsY;
list = (TList*)LightReadout->GetBins();
next = TIter(list);
while((bin = (TH2PolyBin*)next())){
  LDbinsX.insert(pair<int, double>(bin->GetBinNumber(), (bin->GetXMax()+bin->GetXMin())/2));
  LDbinsY.insert(pair<int, double>(bin->GetBinNumber(), (bin->GetYMax()+bin->GetYMin())/2));
}



for (unsigned int EventIt=9; EventIt < NEventsToLoopOver; EventIt++){
	CLTree->GetEntry(EventIt);
	G4Tree->GetEntry(EventIt);


	// We want to store the event display in these TH2Poly's
	TH2Poly *Charge = (TH2Poly*)ChargeReadout->Clone();
	TH2Poly *Light = (TH2Poly*)LightReadout->Clone();

  	cout << "Event " << EventIt << '\n';

	//Get the number of photo detectors
  	int numberPD = total_time_vuv->size();
	//Get the number of charge detectors
  	int numberCD = total_time_charge->size();

	//Get total number of photons
	int detPhotons = 0;
	// We also get the CD with the highest photon count.
	// This is done to later plot the time distirbution of it as a function of time
	int maxPhotons = 0;
	for(int i = 0; i < total_time_vuv->size(); i++){
		detPhotons += total_time_vuv->at(i).size();
		if(total_time_vuv->at(i).size() > total_time_vuv->at(maxPhotons).size()){
			maxPhotons = i;
		}
	}

	cout << "Total number of photons: " << detPhotons << '\n';

	//Get total number of electrons
	int detElectrons = 0;
	// We also get the CD with the highest electron count
	int maxElectrons = 0;
	for(int i = 0; i < total_time_charge->size(); i++){
		detElectrons += total_time_charge->at(i).size();
		if(total_time_charge->at(i).size() > total_time_charge->at(maxElectrons).size()){
			maxElectrons = i;
		}
	}

	cout << "Total number of electrons: " << detElectrons << '\n';

	// Filling the total pulse shape
	// This loops over all the photo detectors
	cout << "Filling the total pulse shape" << '\n';
	for(auto PD: *total_time_vuv){
		// In each photo detector we loop over all entries, which are hit-times
		for(auto time: PD){
			totalPulseShapeLight->Fill(time);
		}
	}

	// Filling the total pulse shape - old school version
	// Does the same as above, but slower
	/* for( int PDit = 0; PDit < numberPD; PDit++){ */
	/* 	for(int PhotonIt = 0; PhotonIt < total_time_vuv->at(PDit).size(); PhotonIt++){ */
	/*		double time = total_time_vuv->at(PDit).at(PhotonIt);*/
	/* 		totalPulseShapeLight->Fill(time); */
	/* 	} */
	/* } */

	// Filling the total pulse shape - charge
	cout << "Filling the total pulse shape - charge" << '\n';
	for(auto CD: *total_time_charge){
		// In each charge pixel detector, we loop over all entries, which are arrival times of electrons
		for(auto time: CD){
			totalPulseShapeCharge->Fill(time);
		}
	}

	// Filling the pulse shape for a specific detector
	cout << "Filling the pulse shape for a specific detector" << '\n';
	for(auto time: total_time_vuv->at(maxPhotons)){
		pulseShapeLightDetector->Fill(time);
	}

	// Filling the pulse shape for a specific detector - charge
	cout << "Filling the pulse shape for a specific detector - charge" << '\n';
	for(auto time: total_time_charge->at(maxElectrons)){
		pulseShapeChargeDetector->Fill(time);
	}



	// Here we fill the total readout , once as a simple TH2, and once as a TH2Poly
	cout << "Filling Event Displays" << '\n';
	for( int vecIt = 0; vecIt < numberPD; vecIt++){
		//Get the position of the detector
		double x = LDbinsX.at(vecIt+1);
		double y = LDbinsY.at(vecIt+1);
		//Get the number of photons
		int nPhotons = total_time_vuv->at(vecIt).size();
		//Fill the histogram
		totalReadoutLight->Fill(x, y, nPhotons);

		//Fill the TH2Poly
		int  bin = Light->FindBin(x, y);
		Light->SetBinContent(bin, nPhotons);
	}

	for(int vecIt = 0; vecIt < numberCD; vecIt++){
		//Get the position of the detector
		double x = CDbinsX.at(vecIt+1);
		double y = CDbinsY.at(vecIt+1);
		//Get the number of electrons
		int nElectrons = total_time_charge->at(vecIt).size();
		//Fill the histogram
		totalReadoutCharge->Fill(x, y, nElectrons);

		//Fill the TH2Poly
		int  bin = Charge->FindBin(x, y);
		Charge->SetBinContent(bin, nElectrons);
	}


	//Plot the total pulse shape
	TCanvas *cTotalPulseShapeLight = new TCanvas("cTotalPulseShapeLight", "cTotalPulseShapeLight", 800,800);
	totalPulseShapeLight->Draw("HIST");
	cTotalPulseShapeLight->SaveAs("cTotalPulseShapeLight.png");

	//Plot the total pulse shape - charge
	TCanvas *cTotalPulseShapeCharge = new TCanvas("cTotalPulseShapeCharge", "cTotalPulseShapeCharge", 800,800);
	totalPulseShapeCharge->Draw("HIST");
	cTotalPulseShapeCharge->SaveAs("cTotalPulseShapeCharge.png");

	//Plot the pulse shape for a specific detector
	TCanvas *cPulseShapeLightDetector17 = new TCanvas("cPulseShapeLightDetector17", "cPulseShapeLightDetector17", 800,800);
	pulseShapeLightDetector->Draw("HIST");
	cPulseShapeLightDetector17->SaveAs("cPulseShapeLightDetector17.png");

	//Plot the pulse shape for a specific detector - charge
	TCanvas *cPulseShapeChargeDetector17 = new TCanvas("cPulseShapeChargeDetector17", "cPulseShapeChargeDetector17", 800,800);
	pulseShapeChargeDetector->Draw("HIST");
	cPulseShapeChargeDetector17->SaveAs("cPulseShapeChargeDetector17.png");

	//Plot the total readout
	TCanvas *cTotalReadoutCharge = new TCanvas("cTotalReadoutCharge", "cTotalReadoutCharge", 10000,10000);
	totalReadoutCharge->Draw("COLZ");
	cTotalReadoutCharge->SaveAs("cTotalReadoutCharge.png");

	//Plot the total readout - light
	TCanvas *cTotalReadoutLight = new TCanvas("cTotalReadoutLight", "cTotalReadoutLight", 10000,10000);
	totalReadoutLight->Draw("COLZ");
	cTotalReadoutLight->SaveAs("cTotalReadoutLight.png");

	TCanvas *cChargeReadout = new TCanvas("cChargeReadout", "cChargeReadout", 10000,10000);
	Charge->Draw("COLZ");
	cChargeReadout->SaveAs("cChargeReadout.png");

	TCanvas *cLightReadout = new TCanvas("cLightReadout", "cLightReadout", 10000,10000);
	Light->Draw("COLZ");
	cLightReadout->SaveAs("cLightReadout.png");

	cout << "Done with event " << EventIt << endl;
}//End for EventIt - loop over events in the tree

return 0;

}//End main
