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

using namespace std;

void macro(string File="nuis_out", string InputDirectory="", string OutputDirectory="" )
{

  //Making plots nicer
  gROOT->ForceStyle(1); gStyle->SetPadTopMargin(0.07); gStyle->SetPadRightMargin(0.05); gStyle->SetPadLeftMargin(0.15); gStyle->SetPadBottomMargin(0.16); gStyle->SetLabelSize(0.06,"xyz"); gStyle->SetTitleSize(0.06,"xyz"); gStyle->SetTitleOffset(0.9,"x"); gStyle->SetTitleOffset(1.1,"y"); gStyle->SetTitleOffset(0.9,"z"); gStyle->SetStatX(0.8); gStyle->SetStatW(0.2); gStyle->SetStatY(0.85); gStyle->SetStatH(0.1); gStyle->SetOptStat(0); gStyle->SetHistLineWidth(3); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1); gStyle->SetPadGridX(kTRUE); gStyle->SetPadGridY(kTRUE);


  //Defining output path for PNG and output file name for ROOT file
  string PNGOutputFolderName, ROOTOutputFileName;
  ROOTOutputFileName=Form("%s/ROOT/File_%s.root", OutputDirectory.c_str(), File.c_str());
  PNGOutputFolderName=Form("%s/PNGs/%s", OutputDirectory.c_str(), File.c_str());

  char* InputFileName;
  InputFileName = Form("%s/%s.root", InputDirectory.c_str(), File.c_str());
  TFile *InputFile = new TFile(InputFileName);

  TTree *SimInputTree = (TTree*)InputFile->Get("ScintSim_tree");
  TTree *G4InputTree = (TTree*)InputFile->Get("event_tree");

  vector<vector<double>> *total_time_vuv = 0;
  vector<vector<double>> *total_time_charge = 0;
  vector<double> *light_yield = 0;
  vector<double> *charge_yield = 0;
  SimInputTree->SetBranchAddress("total_time_vuv", &total_time_vuv);
  SimInputTree->SetBranchAddress("total_time_charge", &total_time_charge);
  SimInputTree->SetBranchAddress("LightYield", &light_yield);
  SimInputTree->SetBranchAddress("ChargeYield", &charge_yield);

  vector<double> *generator_initial_particle_energy = nullptr;
  vector<double> *generator_inital_particle_x = nullptr;
  vector<double> *generator_inital_particle_y = nullptr;
  vector<double> *generator_inital_particle_z = nullptr;
  vector<double> *hit_length = nullptr;
  vector<double> *hit_energy_deposit = nullptr;
  int event = 0;

  G4InputTree->SetBranchAddress("generator_initial_particle_energy", &generator_initial_particle_energy);
  G4InputTree->SetBranchAddress("generator_initial_particle_x", &generator_inital_particle_x);
  G4InputTree->SetBranchAddress("generator_initial_particle_y", &generator_inital_particle_y);
  G4InputTree->SetBranchAddress("generator_initial_particle_z", &generator_inital_particle_z);
  G4InputTree->SetBranchAddress("hit_length", &hit_length);
  G4InputTree->SetBranchAddress("hit_energy_deposit", &hit_energy_deposit);
  G4InputTree->SetBranchAddress("event", &event);

  vector<vector<double>> *pixel_reset = nullptr;


  int NEventsToLoopOver = SimInputTree->GetEntries();
  cout << "Looping over " << NEventsToLoopOver << " events..." << endl;
  cout << endl;

std::vector<std::vector<int>> opdet_type;
std::vector<int> opdet_direction;
std::vector<int> opdet_id;
std::vector<std::vector<double>> opdet_position;


std::cout << "Loading Photon Detector positions..." << std::endl;
std::ifstream detector_positions_file;
detector_positions_file.open("full_dune.txt_");
if(detector_positions_file.is_open()) std::cout << "File opened successfully" << std::endl;
else {std::cout << "File not found." << std::endl; exit(1);}
while(!detector_positions_file.eof()) {
	int num_opdet, type_opdet, direction; double x_opdet, y_opdet, z_opdet;
	if(detector_positions_file >> num_opdet >> x_opdet >> y_opdet >> z_opdet >> type_opdet >> direction) {
	    std::vector<int> type({num_opdet, type_opdet});
	    std::vector<double> position({x_opdet, y_opdet, z_opdet});

	    opdet_type.push_back(type);
	    opdet_id.push_back(num_opdet);
	    opdet_position.push_back(position);
	    opdet_direction.push_back(direction);
	}
	else{ break; }
}
detector_positions_file.close();
cout << "Dine reading file" << endl;


TH1F *h_PulseShape = new TH1F("h_PulseShape", "h_PulseShape", 1000, 0, 4);
TH2F *h_Detector = new TH2F("h_Detector", "h_Detector", 140, 0, 1400, 60, 0, 600);
TH1F *h_NDetPhot = new TH1F("h_NDetPhot", "h_NDetPhot", 10, 15e3, 20e3);

for (unsigned int EventIt=0; EventIt < NEventsToLoopOver; EventIt++){
	  SimInputTree->GetEntry(EventIt);
	  G4InputTree->GetEntry(EventIt);

	  // Loop over all PD
	  int pdIT = 0;
	  int NDetPhot = 0;
	  for(auto pmt:*total_time_vuv){
		  for(auto time:pmt){
			  h_PulseShape->Fill(time);
		  }
		  NDetPhot += pmt.size();

		  // Get position of current PD
 		  int opdet = opdet_id[pdIT];
		  double x = opdet_position[opdet][0];
		  double y = opdet_position[opdet][1];
		  double z = opdet_position[opdet][2];
		  int bin = h_Detector->FindBin(x,y);
		  h_Detector->SetBinContent(bin, h_Detector->GetBinContent(bin)+pmt.size());
		  pdIT+=1;
	  }
	  h_NDetPhot->Fill(NDetPhot);

}



h_Detector->Scale(1./NEventsToLoopOver);

gStyle->SetPadRightMargin(0.15);
TCanvas *c1 = new TCanvas("c1", "c1", 1920, 1080);
h_Detector->SetContour(500);
h_Detector->Draw("colz");
c1->SaveAs(Form("%s/Detector_Average.png", PNGOutputFolderName.c_str()));

TCanvas *c2 = new TCanvas("c2", "c2", 1920, 1080);
h_PulseShape->Draw("HIST");
c2->SetLogy();
c2->SaveAs(Form("%s/PulseShape.png", PNGOutputFolderName.c_str()));

TCanvas *c3 = new TCanvas("c3", "c3", 1920, 1080);
h_NDetPhot->GetXaxis()->SetNdivisions(505);
h_NDetPhot->Draw("HIST");
c3->SaveAs(Form("%s/NDetPhot.png", PNGOutputFolderName.c_str()));



}
