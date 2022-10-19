// This script creates events from a signal and a background file
// It implements certain detector effects like the darkcount and efficiency
// Further a very primitve trigger is implemented

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



void macro(string CLSignal="nuis_out", string CLBackground = "",  string InputDirectory="", string OutputDirectory="" )
{

//Making plots nicer
gROOT->ForceStyle(1); gStyle->SetPadTopMargin(0.07); gStyle->SetPadRightMargin(0.05); gStyle->SetPadLeftMargin(0.15); gStyle->SetPadBottomMargin(0.16); gStyle->SetLabelSize(0.06,"xyz"); gStyle->SetTitleSize(0.06,"xyz"); gStyle->SetTitleOffset(0.9,"x"); gStyle->SetTitleOffset(1.1,"y"); gStyle->SetTitleOffset(0.9,"z"); gStyle->SetStatX(0.8); gStyle->SetStatW(0.2); gStyle->SetStatY(0.85); gStyle->SetStatH(0.1); gStyle->SetOptStat(0); gStyle->SetHistLineWidth(3); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1); gStyle->SetPadGridX(kTRUE); gStyle->SetPadGridY(kTRUE);
gStyle->SetPadRightMargin(0.15);

//Defining output path for PNG and output file name for ROOT file
string PNGOutputFolderName, ROOTOutputFileName;
ROOTOutputFileName=Form("%s/ROOT/NEUTFile_%s.root", OutputDirectory.c_str(), CLSignal.c_str());
PNGOutputFolderName=Form("%s/PNGs/%s", OutputDirectory.c_str(), CLSignal.c_str());





TRandom *r = new TRandom3(0);


char* CLFileName;
CLFileName = Form("%s/%s.root", InputDirectory.c_str(), CLSignal.c_str());
TFile *CLFile = new TFile(CLFileName);
TTree *CLTree = (TTree*)CLFile->Get("ScintSim_tree");

// These are vectors of vectors, where each line corresponds to a single detection element
// Inside, the arrivial times of the phtoons/electrons gitting this element are stored as doubles.
//  [ SiPM1[t0, t1, t2,t3,...], SiPM2[t0, t1, t2,t3,...], SiPM3[t0, t1, t2,t3,...], ...]
vector<vector<double>> *signal_light = nullptr;
vector<vector<double>> *signal_charge = nullptr;

// These are the placement maps for the Charge and Light detectors
TH2Poly *ChargeReadout = (TH2Poly*)CLFile->Get("charge_detector_map");
TH2Poly *LightReadout = (TH2Poly*)CLFile->Get("light_detector_map");

CLTree->SetBranchAddress("total_time_vuv", &signal_light);
CLTree->SetBranchAddress("total_time_charge", &signal_charge);


char* CLFileName2;
CLFileName2 = Form("%s/%s.root", InputDirectory.c_str(), CLBackground.c_str());
TFile *CLBkg = new TFile(CLFileName2);
TTree *CLBkgTree = (TTree*)CLBkg->Get("ScintSim_tree");

vector<vector<double>> *background_light = nullptr;
vector<vector<double>> *background_charge = nullptr;

CLBkgTree->SetBranchAddress("total_time_vuv", &background_light);
CLBkgTree->SetBranchAddress("total_time_charge", &background_charge);


int NEventsToLoopOver = CLTree->GetEntries(); // 100000
cout << "Looping over " << NEventsToLoopOver << " events..." << endl;
cout << endl;


// Building the look up tables
// See other example code
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

int nhist = 100*100;


TH1F *hTotalPhot = new TH1F("hTotalPhot","test",20,0,20000);



// Build a single event

//Find number of background evnets expected

// We have 1500 kg of lAr, with 1.01 Bq/kg of Ar39 backgorund, meaning in 20us, we ahve around 0.03 events.
// Lets be very pessimistic and say we have 3 decays per recording window
double BkgEventsPerReadoutWindow = 3;
int NBkg = r->Poisson(BkgEventsPerReadoutWindow);
//We want to take random events from the background file
int startingBkgEvent = r->Uniform(0, CLBkgTree->GetEntries()-(2*NBkg));

//For each PD we build a histogram of arrival times (time resoultion of 0.2 microseconds = 200ns)
TH1F** arrayBkg = new TH1F*[nhist];
for (int i=0;i<nhist;i++){
    arrayBkg[i] = new TH1F(Form("hBkg%d",i),"test",100,0,20); // 20 micro seconds
}
//Build the total accumulated background histogram over all PDs
TH1F *hTotalPhotBkg = new TH1F("hTotalPhotBkg","test",100,0,20);

//Loop over the background events
for(int i=0; i<NBkg; i++){
  //Get the corresponding background event
  CLBkgTree->GetEntry(startingBkgEvent+i);
  // The background events are generated in a time between -10s and +10s.
  // This is due ot how G4 does things, so we want to shift them first to zero and then distribute them normally through the readout window
  double time_offset = r->Uniform(0, 20);
  // Loop over all PDs
  for(int j=0; j<background_light->size(); j++){
    // Sort the arrival times to be able to zero them
    sort(background_light->at(j).begin(), background_light->at(j).end());
    // Loop over all photons in this PD
    for(int k=0; k<background_light->at(j).size(); k++){
      // Assume a 40 % efficiency, also for Ar30 events
      if(r->Uniform(0,1)<0.4){
	// Fill the histograms with the arrival times, shifted to 0 and then shifted by the random offset
      	arrayBkg[j]->Fill((background_light->at(j).at(k)-background_light->at(j).at(0))+time_offset);
      	hTotalPhotBkg->Fill(abs(background_light->at(j).at(k)-background_light->at(j).at(0))+time_offset);
      }//End 40 eff
    }//End loop over photons
  }//End loop over SiPMs
}//End loop over background events

// Now we have histograms with Ar39 background events


// We now want to get the dark counts
CLTree->GetEntry(0);
//We assume 3 coutns/ micro second per PD
double darkCountRatePerUS = 3; // 3 per micro second
double darkCountRatePerReadoutWindow = darkCountRatePerUS*20; // 20 microseconds as readout winodw
// Build our favourite vector
vector<vector<double>> darkcounttimes;
// Loop over all PDs
for (int i=0;i<signal_light->size();i++){
  // Uniformly distribute the dark counts in the readout window
  // Sample the number of dark counts from a Poisson distribution
  vector<double> darkcounttimes_i;
  int Ndarkcounts = r->Poisson(darkCountRatePerReadoutWindow);
  for (int j=0;j<Ndarkcounts;j++){
    darkcounttimes_i.push_back(r->Uniform(0,20));
  }//End loop over dark counts
  darkcounttimes.push_back(darkcounttimes_i);
}//End loop over PDs

// We now sort the dark counts in to the same histograms as the background events
// Loop over all PDs
for(int j = 0; j<darkcounttimes.size(); j++){
  // Loop over all dark counts in this PD
  for(int k = 0; k<darkcounttimes.at(j).size(); k++){
    // Fill the histograms with the arrival times.
    // This time no shifts are required, as we dsitrbuted them inside the readout window
    arrayBkg[j]->Fill(darkcounttimes.at(j).at(k));
    hTotalPhotBkg->Fill(darkcounttimes.at(j).at(k));
  }
}

// Build the signal histograms
TH1F** arraySig = new TH1F*[nhist];
for (int i=0;i<nhist;i++){
    arraySig[i] = new TH1F(Form("hSig%d",i),"test",100,0,20); // 2000 micro seconds
}
TH1F *hTotalPhotSig = new TH1F("hTotalPhotSig","test",100,0,20);

//Our favourite event number
int eventNumber = r->Uniform(0, CLTree->GetEntries());
CLTree->GetEntry(eventNumber);

// We want to record the 5 microseconds before the event
double time_offset = 5;
// Loop over all PDs
for(int j=0; j<signal_light->size(); j++){
   // Loop over all photons in this PD
  for(int k=0; k<signal_light->at(j).size(); k++){
      // Assume a 40 % efficiency
      if(r->Uniform(0,1)<0.4){
	//Fill histogram
    	arraySig[j]->Fill(signal_light->at(j).at(k)+time_offset);
    	hTotalPhotSig->Fill(signal_light->at(j).at(k)+time_offset);
      }
  }
}

// Combine the signal and background histograms
TH1F* hTotal = (TH1F*)hTotalPhotSig->Clone("hTotal");
hTotal->Add(hTotalPhotBkg);
TH1F** array = new TH1F*[nhist];
for (int i=0;i<nhist;i++){
    array[i] = (TH1F*)arraySig[i]->Clone(Form("h%d",i));
    array[i]->Add(arrayBkg[i]);
}


// We now have the background and signal histograms.
// We now can build event displays and try to trigger

// Lets get things ready for event displays
TH2Poly* Sig = (TH2Poly*)LightReadout->Clone("Sig");
TH2Poly* Bkg = (TH2Poly*)LightReadout->Clone("Bkg");
TH2Poly* Total = (TH2Poly*)LightReadout->Clone("Total");

// We fill the TH2Polys with the events in the histograms
for(int i=0; i<nhist; i++){
	int Nsig = arraySig[i]->GetEntries();
	int Nbkg = arrayBkg[i]->GetEntries();
	Sig->SetBinContent(i+1, Nsig);
	Bkg->SetBinContent(i+1, Nbkg);
	Total->SetBinContent(i+1, Nsig+Nbkg);
}

// We disect the detector in to 100 10*10 squares and check if we can trigger on those
int Nsquares = 100;
int Nbins = 10;
// Each element correpsonds to a sqaure of 10x10 PDs
// It holds the corresponding bin numbers
vector<vector<int>> bins;
// This is the overall bin we want to look at: (binItX, binItY) is a bin of 100 PD
for(int binItX = 0; binItX < 10; binItX ++){
	for(int binItY = 0; binItY < 10; binItY++){
		// In this one bin, we are grouping together 100 PD.
		// We store those bins in a vector - because vectors are cool
		vector<int> bin_;
		for(int squareItY = 0; squareItY < Nbins; squareItY++){
			for(int squareItX = 0; squareItX < Nbins; squareItX++){
				// Numbering:
				// 100 PD in Y axis, before a new bin starts.
				// Going down one line requries 1000 PD
				int bin = squareItY * 100 + squareItX + 1000*binItY + 10*binItX;
				// We add the bin to the vector
				bin_.push_back(bin);
			}
		}
		// We add the vector to the vector of vectors,
		bins.push_back(bin_);
	}
}

// We now check if we have more than 6 squares giving more than 500 photons in 1.2 micro seconds
int TriggerTimeWindow = 6; // In units of bins
int TotalTriggerThreshold = 500; // In units of photons
int NumberTriggeredSquares = 6;
//This sores how often a specifc time was triggered
vector<int> triggered(100,0);

for(int i = 0; i < bins.size(); i++){
	// All of these bins sit in the same square
	// Now we loop over different time windows and check if we can trigger on this square
	for(int timeIt = 0; timeIt + TriggerTimeWindow < array[0]->GetNbinsX(); timeIt++){
		int TriggerCounter = 0;
		int totalPhotons = 0;
		for(int j = 0; j < bins[i].size(); j++){
			totalPhotons+=array[bins[i][j]]->Integral(timeIt, timeIt+TriggerTimeWindow);
		}
		int allPDPhtons = hTotal->Integral(timeIt, timeIt+TriggerTimeWindow);
		if(totalPhotons > TotalTriggerThreshold ){
			triggered[timeIt]++;
		}
	}//end of time loop
}//end of square loop


// Draw the trigger window
TH1F* hTriggerTime = new TH1F("hTriggerTime", "hTriggerTime", 100, 0, 20);
for(int i = 0; i < triggered.size(); i++){
	if(triggered[i] > NumberTriggeredSquares){
		cout << "Time " << array[0]->GetBinCenter(i) << endl;
		hTriggerTime->SetBinContent(i, 7000);
	}
}

TCanvas *c1 = new TCanvas("c1","c1",800,600);
c1->Divide(2,2);
c1->cd(1);
hTotalPhotBkg->SetTitle("Background + Darkcount");
hTotalPhotBkg->Draw("HIST");
c1->cd(2);
hTotalPhotSig->SetTitle("Signal");
hTotalPhotSig->Draw("HIST");
c1->cd(3);
hTotal->SetTitle("Total");
hTotal->Draw("HIST");
c1->cd(4);
hTotal->Draw("HIST");
hTriggerTime->SetLineColor(kRed);
hTriggerTime->Draw("SAME");
c1->SaveAs(Form("%s/Total.png",PNGOutputFolderName.c_str()));

TCanvas *c2 = new TCanvas("c2","c2",1000,1000);
c2->Divide(2,2);
c2->cd(1);
Sig->SetTitle("Signal");
Sig->Draw("colz");
c2->cd(2);
Bkg->SetTitle("Background");
Bkg->Draw("colz");
c2->cd(3);
Total->SetTitle("Total");
Total->Draw("colz");
c2->SaveAs(Form("%s/Total2D.png",PNGOutputFolderName.c_str()));

return 0;

}//End main
