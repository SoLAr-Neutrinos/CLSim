// this file contains parameters used in simulation MERGED FILE
#include "TRandom3.h"
namespace radiological {
	
        ///-------------------------------------
        ////----LAr & Ar39 properties--------------------
        /////-------------------------------------
        const double MassE = 0.510998910;       // mass electron - MeV/c^2
        const double Q_Ar = 0.565;                      //Q value of decay - Ar39
        const double Q_Co60B = 0.318;           //End point of beta decay of Co60[MeV]
        const double Q_Ar42 = 0.599;            //End point of beta decay of Ar42[MeV]
        const double Q_K42 = 3.525;
        const double Q_40KB = 1.35;             //End point of beta decay of 40K[MeV]
        const double Q_Kr85B1 = 0.84;           //End point of beta decay of Kr85[MeV]
        const double Q_Kr85B2 = 0.687;
        const double Q_Kr85G1 = 0.151;
        const double Q_Kr85G2 = 0.305;
        const double Q_Pb214 = 1.03;
        const double Q_Bi214 = 3.2;

        ///-------------------------------------
        ////----SN properties---------------------
        /////-------------------------------------
        const double Eav = 20.;                         // Average energy for SN spectrum
        const double expected_sn = 2.8;         // For poisson weighting

        ///-------------------------------------
        ////----Radon properties---------------------
        /////-------------------------------------
        const double Q_Rn222 = 5.490;                              // deposited energy from a radon decay - Rn-222 --> Po-218
        const double Q_Po218 = 6.0;
        const double Q_Po214 = 7.7;
        const double Q_Po210 = 5.3;

	//----------------------------------------------
	//-------Positions of radioactive material------
	//----------------------------------------------
	//------[these come from Jingyuan/LArSoft]------
        //40KB,40KG 
        const double K_x_position_range[2] {3.495e2,3.505e2};      // cm
        const double K_y_position_range[2] {-600,600};    // cm
        const double K_z_position_range[2] {0,1395};    // cm
        //Co60B
        const double Co_x_position_range[2] {-5e-1,5e-1};      // cm
        const double Co_y_position_range[2] {-600,600};    // cm
        const double Co_z_position_range[2] {0,1395};    // cm
        //Po210
        const double Po_x_position_range[2] {4.77e-1,1.477};      // cm
        const double Po_y_position_range[2] {-600,600};    // cm
        const double Po_z_position_range[2] {0,1395};    // cm
        
}
