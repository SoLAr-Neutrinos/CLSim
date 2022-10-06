// this file contains parameters used in simulation

#include "TRandom3.h"
namespace parameters {

	// output file
	const char *output_file_name = "testing_xenon.root";

    ///-------------------------------------
    //--------Simulation Settings-----------
    ///-------------------------------------
    // Xenon doping
    const bool simulate_xenon = false;   // enables xenon doping, false = pure argon

    // WLS Reflective foils
    const bool include_reflected = false; // enables WLS relfective foils on the cathode (visible light)

    // Timings
    const bool include_timings = true;   // enables timings (emission and transport)

    // Electric Field in kV/cm on the detector
    const double electric_field = 0.5;

    // Drift values
    double drift_velocity = 0.16; // cm/us
    double drift_long = 6.627; //cm^2/s
    double drift_trans = 13.237; //cm^2/s


    // Detection Efficency
    double SiPM_QE = 0.4;
    double total_QE = SiPM_QE;

    // timing parametersiation properties
    const double timing_discretisation_step_size = 2.12; // cm

    /* // scintillation timing properties */
    const double scint_time_window = 0.00001;   // 10 us
    // timing
    const double t_singlet = 6e-9;//0.000000006;       // 6ns
    const double t_triplet = 1.5e-6;//0.0000015;         // 1.5 us
    // prompt/late ratio
    // electron-like
    const double singlet_fraction_electron = 0.30;
    const double triplet_fraction_electron = 0.70;
    // alpha
    const double singlet_fraction_alpha = 0.75;
    const double triplet_fraction_alpha = 0.25;

    // xenon -- Aprile, Doke https://arxiv.org/abs/0910.4956
    // timing
    const double t_singlet_Xe = 0.0000000042;      // 4.2 ns
    const double t_triplet_Xe = 0.000000022;       // 22 ns
    // prompt/late ratio
    // minimal dependence on particle type, approximate as single value
    const double singlet_fraction_Xe = 0.30;
    const double triplet_fraction_Xe = 0.70;

}
