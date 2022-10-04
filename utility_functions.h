#ifndef UTILITY_FUNCTION_H
#define UTILITY_FUNCTION_H

// class containing energy spectrums and other utility functions required by main code

#include "TF1.h"

class utility_functions {

private:
	TF1 *fScintillation_function_electron = nullptr;
        TF1 *fScintillation_function_alpha = nullptr;
        TF1 *fScintillation_function_prompt = nullptr;
        TF1 *fScintillation_function_xenon = nullptr;

public:

	// constructor
	utility_functions(){};

	// destructor
	~utility_functions();

	// poisson distribution function
	int poisson(const double mean, const double draw, const double eng) const;

	// scintillation function
	static double scintillation_function( const double *t, const double *par);
        void initalise_scintillation_functions_argon(const double t_singlet, const double t_triplet, const double singlet_fraction_electron, const double triplet_fraction_electron,
                                                        const double singlet_fraction_alpha, const double triplet_fraction_alpha, const double scint_time_window);
        void initalise_scintillation_functions_xenon(const double t_singlet_Xe, const double t_triplet_Xe, const double singlet_fraction_Xe, const double triplet_fraction_Xe,
                                                        const double scint_time_window);

	// RETURNS THE SCINTILLATION TIME IN SECONDS
        double get_scintillation_time_electron() { return fScintillation_function_electron->GetRandom(); }

        double get_scintillation_time_alpha() { return fScintillation_function_alpha->GetRandom(); }
        double get_scintillation_time_prompt() { return fScintillation_function_prompt->GetRandom(); }
        double get_scintillation_time_xenon() { return fScintillation_function_xenon->GetRandom(); }

        // Spectrum Function (Beta Decay) /////chrisflynn
        static double SpectrumFunction(double *x, double *par);

        //Supernova Spectrum
        static double fsn(double *x, double *par);

        // HEP spectrum ////chrisflynn
        static double fhep(double *x, double *par);

        //Solar spectrum ////chrisflynn
        static double fso(double *x, double *par);

        //Rn Spectrum ////chrisflynn
        static double Rn_function(double *x, double *par);

};

#endif
