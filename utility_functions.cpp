// implementation of utility functions class

#include "utility_functions.h"

#include <cmath>
#include "TMath.h" //chrisflynn
#include "TVector3.h" //chriflynn
#include "TF1.h" //chriflynn
#include "TRandom.h" //chriflynn
#include <chrono> //chriflynn
#include <iostream> //chriflynn
#include "TSpline.h" //chriflynn
#include "TFile.h" //chriflynn

// destructor
utility_functions::~utility_functions(){
	//if fScintillation_function != nullptr {
		delete fScintillation_function_electron;
		delete fScintillation_function_alpha;
		delete fScintillation_function_prompt;
		delete fScintillation_function_xenon;
	//}
}

// poisson distribution function
int utility_functions::poisson(const double mean, const double draw, const double eng) const {
	int number = 0;
	const int border = 16;
	double limit = 2e9;

	if(mean <= border) {
		double position = draw;
		double poissonValue = std::exp(-mean);
		double poissonSum = poissonValue;

		while(poissonSum <= position) {
			number++;
			poissonValue *= mean/number;
			poissonSum += poissonValue;
		}
		return number;
	} // the case of mean <= 16

	double value, t, y;
	t = std::sqrt(-2*std::log(draw));
	y = 2*3.141592654*eng;
	t *= std::cos(y);
	value = mean + t*std::sqrt(mean) + 0.5;
	if(value <= 0) {return 0; }
	if(value >= limit) { return limit; }
	return value;
}

// scintillation time spectrum function
double utility_functions::scintillation_function(const double *t, const double *par) {
	double time = *t;
	double t_singlet = par[0];
	double t_triplet = par[1];
	double singlet_part = par[2];
	double triplet_part = par[3];

	double Scintillation = exp(-(time/t_singlet))*singlet_part/t_singlet + exp(-(time/t_triplet))*triplet_part/t_triplet;

	return Scintillation;
}

// function to create scintillation function TF1 with required parameters
void utility_functions::initalise_scintillation_functions_argon(const double t_singlet, const double t_triplet, const double singlet_fraction_electron,
					const double triplet_fraction_electron, const double singlet_fraction_alpha, const double triplet_fraction_alpha,
		                	const double scint_time_window) {

	// create scintillation spectrum for electrons
	fScintillation_function_electron = new TF1("Scintillation Timing", scintillation_function, 0, scint_time_window, 4);
	fScintillation_function_electron->SetParameter(0, t_singlet);
    fScintillation_function_electron->SetParameter(1, t_triplet);
    fScintillation_function_electron->SetParameter(2, singlet_fraction_electron);
    fScintillation_function_electron->SetParameter(3, triplet_fraction_electron);
    fScintillation_function_electron->SetNpx(50000);

    // create scintillation spectrum for alphas
	fScintillation_function_alpha = new TF1("Scintillation Timing", scintillation_function, 0, scint_time_window, 4);
	fScintillation_function_alpha->SetParameter(0, t_singlet);
    fScintillation_function_alpha->SetParameter(1, t_triplet);
    fScintillation_function_alpha->SetParameter(2, singlet_fraction_alpha);
    fScintillation_function_alpha->SetParameter(3, triplet_fraction_alpha);
    fScintillation_function_alpha->SetNpx(50000);

    // create scintillation spectrum for prompt only case
	fScintillation_function_prompt = new TF1("Scintillation Timing", scintillation_function, 0, scint_time_window, 4);
	fScintillation_function_prompt->SetParameter(0, t_singlet);
    fScintillation_function_prompt->SetParameter(1, t_triplet);
    fScintillation_function_prompt->SetParameter(2, 1.00);
    fScintillation_function_prompt->SetParameter(3, 0.00);
    fScintillation_function_prompt->SetNpx(50000);
}

void utility_functions::initalise_scintillation_functions_xenon(const double t_singlet_Xe, const double t_triplet_Xe, const double singlet_fraction_Xe,
																const double triplet_fraction_Xe, const double scint_time_window) {

    // create scintillation spectrum for xenon
    fScintillation_function_xenon = new TF1("Scintillation Timing", scintillation_function, 0, scint_time_window, 4);
    fScintillation_function_xenon->SetParameter(0, t_singlet_Xe);
    fScintillation_function_xenon->SetParameter(1, t_triplet_Xe);
    fScintillation_function_xenon->SetParameter(2, singlet_fraction_Xe);
    fScintillation_function_xenon->SetParameter(3, triplet_fraction_Xe);
    fScintillation_function_xenon->SetNpx(50000);
}

//Beta Decay Function ///chrisflynn
double utility_functions::SpectrumFunction(double *x, double *par)
{
       double KE = *x;
       double Q = par[0];
       double MassE = 0.510998910; //mass electron - Mev/c^2

       double N = std::sqrt(pow(KE,2) + 2*KE*MassE) * std::pow((Q-KE),2) * (KE+MassE);

       return N;
}

//Supernova neutrino energy spectrum ///chrisflynn
//Spectrum energy comes froma 1kpc galactic supernova
double utility_functions::fsn(double *x, double *par)
{
       double E = *x;
       double Eav = par[0];

       double f_nu = std::pow(E,3)*std::exp(-4*E/Eav);

       return f_nu;
}

//Solar neutrino energy spectrum ////chrisflynn
double utility_functions::fso(double *x, double *par)
{
       double E = *x;
       double Eav = par[0];

       TFile *f = new TFile("../spectra/solar_neu_energy_spline.root");
       TSpline3 *spline = (TSpline3*)f->Get("../spectra/Spline3");
       double f_s_neu = spline->Eval(E);
       f->Close();

       return f_s_neu;
}


//hep spectrum ///chrisflynn
double utility_functions::fhep(double *x, double *par)
{
      double E = *x;
      double Eav = par[0];

      TFile *f1 = new TFile("../spectra/hep_neu_energy_spline.root");
      TSpline3 *spline = (TSpline3*)f1->Get("../spectra/Spline3");
      double f_hep_neu = spline->Eval(E);
      f1->Close();

      return f_hep_neu;
}

//Radon-22 decay energy spectrum (Gaussian about the alpha decay energy) //chrisflynn
double utility_functions::Rn_function(double *x, double *par)
{
      double E = *x;
      double Q_Rn = par[0];

      double sigma = 0.01;
      //double sigma = (E-Q_Rn) / std::sgrt(1.3863); //1.3863 = ln(4)
      double sigma_sq = sigma * sigma;
      double gauss = 1/(sigma*std::sqrt(2*3.1416)) * std::exp((-1*std::pow((E-Q_Rn),2))/(2*sigma_sq));

      return gauss;
      }
