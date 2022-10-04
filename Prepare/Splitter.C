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

using namespace std;


void Splitter(string inputfile = "")
{
  TFile *f = new TFile(inputfile.c_str());
  TTree *T = (TTree*)f->Get("event_tree");

  for(int i = 1; i<100; i++){
    TFile *f1 = new TFile(Form("%s_%d",inputfile.c_str(),i),"RECREATE");
    TTree *T1 = T->CopyTree(Form("(event<=%d*100)*(event>=%d*100)",i, i-1));
    f1->cd();
    T1->Write("event_tree");
    f1->Close();
  }
}
