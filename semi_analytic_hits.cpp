#include "semi_analytic_hits.h"

// implementation of semi-analytic model for number of incident photons

#include <iostream>
#include<fstream>
#include <cmath>

#include "TRandom.h"
#include "TSystem.h"
#include "TMath.h"
#include "TFormula.h"
#include "Math/SpecFuncMathMore.h"

using namespace std;
bool debug_2 = false;
// constructor
semi_analytic_hits::semi_analytic_hits() {

  // load mathmore library
  gSystem->Load("libMathMore.so");
  if(gSystem->Load("libMathMore.so") < 0) {
      throw(std::runtime_error("Unable to load MathMore library"));
    }
  _mathmore_loaded_ = true;

  std::cout << "Light simulation for DUNE Single Phase detector." << std::endl;
  std::cout << std::endl;

}


// Function to set the SiPM size
void semi_analytic_hits::setPixelSize(const double y,const double z){
	y_dimension_detector = y;	// cm
	z_dimension_detector = z;	// cm

}


// Funciton for the LArQL modeling
double QBirks(double dE, double dx, double E){

	double k_Birks=0.0486;
	double A_Birks=0.800;
	double Wion = 23.6 ;
	double Nex_Ni = 0.29;
	double Ni = 1.0e6/Wion;
	double pLAr = 1.30;

    double edep = dE/dx;
    if(edep<1) edep = 1.;
    double nom = A_Birks/Wion;
    double denom = 1+ k_Birks*edep /(pLAr * E);
    return nom/denom*1.0e6;
}

double Corr(double dE, double dx, double E){
    double alpha = 0.032;
    double beta = 0.008258;
    double edep = dE/dx;
    if(edep<1) edep = 1.;
    return exp(-E/(alpha *log(edep) + beta));
}

double QChi(double dE, double dx, double E){
    vector<double> chi_param = {2.151572666e-5, -3.988504, 1.38343421, 1.9919521e-6};
    double edep = dE/dx;
    if(edep<1) edep = 1.;
    return chi_param[0]/(chi_param[1]+exp(chi_param[2] + chi_param[3]*edep));
}

double Qinf(double dE, double dx,double  E){

double k_Birks=0.0486;
double A_Birks=0.800;
double Wion = 23.6 ;
double Nex_Ni = 0.29;
double Ni = 1.0e6/Wion;
double pLAr = 1.30;

    return 1.0e6/Wion;
}

double semi_analytic_hits::LArQL(const double energy_deposit, const double hit_distance, const double electric_field){
    return Nex_Ni * Ni + Ni - semi_analytic_hits::LArQQ(energy_deposit,hit_distance, electric_field);
}

double semi_analytic_hits::LArQQ(const double energy_deposit, const double hit_distance, const double electric_field){
    return QBirks(energy_deposit, hit_distance, electric_field)+Corr(energy_deposit, hit_distance, electric_field)*QChi(energy_deposit, hit_distance, electric_field)*Qinf(energy_deposit, hit_distance, electric_field);
  }


// VUV hits calculation
int semi_analytic_hits::VUVHits(const int &Nphotons_created, const TVector3 &ScintPoint, const TVector3 &OpDetPoint, const int &optical_detector_type, const int &scintillation_type, const int &optical_direction) {

  // distance and angle between ScintPoint and OpDetPoint
  double distance = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2) + pow(ScintPoint[1] - OpDetPoint[1],2) + pow(ScintPoint[2] - OpDetPoint[2],2));
  double cosine;
  double theta;


  // distance from center for border corrections
  double r_distance = -1;
  // calculate solid angle:
  double solid_angle = -1;
  // rectangular aperture
  if (optical_detector_type == 1) {
    // set Arapuca geometry struct for solid angle function
    acc detPoint;
    TVector3 ScintPoint_Temp = ScintPoint;
    TVector3 OpDetPoint_Temp = OpDetPoint;
    if(optical_direction == 1){ // x == const plane - this is how the parametrisation was done
	    detPoint.ax = OpDetPoint[0]; detPoint.ay = OpDetPoint[1]; detPoint.az = OpDetPoint[2];
  	    cosine = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2)) / distance;
	    r_distance = sqrt( pow(ScintPoint[1] - center_y_1, 2) + pow(ScintPoint[2] - center_z_1, 2));
    }  // centre coordinates of optical detector
    else if(optical_direction == 2){ // y == const plane - need to excange x and y
	    detPoint.ax = OpDetPoint[1]; detPoint.ay = OpDetPoint[0]; detPoint.az = OpDetPoint[2];
	    OpDetPoint_Temp[0] = detPoint.ax; OpDetPoint_Temp[1] = detPoint.ay; OpDetPoint_Temp[2] = detPoint.az;
	    ScintPoint_Temp[0] = ScintPoint[1]; ScintPoint_Temp[1] = ScintPoint[0];
  	    cosine = sqrt(pow(ScintPoint[1] - OpDetPoint[1],2)) / distance;
	    r_distance = sqrt( pow(ScintPoint[0] - center_x_2, 2) + pow(ScintPoint[2] - center_z_2, 2));
    }  // centre coordinates of optical detector
    else if(optical_direction == 3){// z == const plane - need to excange x and z
	    detPoint.ax = OpDetPoint[2]; detPoint.ay = OpDetPoint[1]; detPoint.az = OpDetPoint[0];
	    OpDetPoint_Temp[0] = detPoint.ax; OpDetPoint_Temp[1] = detPoint.ay; OpDetPoint_Temp[2] = detPoint.az;
	    ScintPoint_Temp[0] = ScintPoint[2]; ScintPoint_Temp[2] = ScintPoint[0];
  	    cosine = sqrt(pow(ScintPoint[2] - OpDetPoint[2],2)) / distance;
	    //cout << "ScintPoint[0]" << ScintPoint[0] << "ScintPoint[1]" << ScintPoint[1] << "ScintPoint[2]" << ScintPoint[2] << endl;
	    //cout << center_x_3 << " " << center_y_3 << " " << center_z_3 << endl;
	    r_distance = sqrt( pow(ScintPoint[0] - center_x_3, 2) + pow(ScintPoint[1] - center_y_3, 2));
	    //cout << "r_distance" << r_distance << endl;
    }  // centre coordinates of optical detector
    else {std::cout << "Error: optical_direction not set correctly" << std::endl; exit(1);}




    detPoint.w = y_dimension_detector; detPoint.h = z_dimension_detector; // width and height in cm of arapuca active window
    TVector3 ScintPoint_rel = ScintPoint_Temp - OpDetPoint_Temp;

    if (cosine < 0.001)
	    solid_angle = 0;
    else
    	solid_angle = solid(detPoint, ScintPoint_rel);

    // calculate solid angle
   if(solid_angle < 0){
	    std::cout << "Error: solid angle is negative" << std::endl;
	    exit(1);
    }
    if(r_distance < 0){
	    std::cout << "Error: r_distance is negative" << std::endl;
	    exit(1);
    }
  }
   else {
    std::cout << "Error: Invalid optical detector type." << endl;
    exit(1);
  }

  theta = acos(cosine)*180./pi;

  // calculate number of photons hits by geometric acceptance: accounting for solid angle and LAr absorbtion length
  double hits_geo = exp(-1.*distance/L_abs) * (solid_angle / (4*pi)) * Nphotons_created;
  if(debug_2){
	  cout << "solid_angle: " << solid_angle << endl;
	  cout << "Nphotons_created: " << Nphotons_created << endl;
	  cout << "Geometric accepted gammas: " << hits_geo << endl;
  }


  // determine Gaisser-Hillas correction for Rayleigh scattering distance and angular dependence, accounting for border effects
  // offset angle bin
  if(theta>89.0){
	  return hits_geo;
  }
  int j = (theta/delta_angle);

  /* std::cout << " theta = " << theta << " delta_angle = " << delta_angle << " j = " << j << std::endl; */
    // identify GH parameters and border corrections by optical detector type and scintillation type
  double pars_ini[4] = {0,0,0,0};
  double s1, s2, s3;
  // determine initial parameters and border corrections by optical detector type and scintillation type
  // flat PDs
  if (optical_detector_type == 0 || optical_detector_type == 1){
    if (scintillation_type == 0) { // argon
      if (j >= 8){
        pars_ini[0] = fGHVUVPars_flat_argon[0][j];
      }
      else{
	double temp1 = fGHVUVPars_flat_argon[0][j];
      	double temp2 = fGHVUVPars_flat_argon[0][j+1];
      	pars_ini[0] = temp1 + (temp2-temp1)*(theta-j*delta_angle)/delta_angle;
      }

      if(j >= 8){
	pars_ini[1] = fGHVUVPars_flat_argon[1][j];
      }
      else{
	double temp1 = fGHVUVPars_flat_argon[1][j];
      	double temp2 = fGHVUVPars_flat_argon[1][j+1];
      	pars_ini[1] = temp1 + (temp2-temp1)*(theta-j*delta_angle)/delta_angle;
      }

      if(j >= 8){
	pars_ini[2] = fGHVUVPars_flat_argon[2][j];
      }
      else{
	double temp1 = fGHVUVPars_flat_argon[2][j];
      	double temp2 = fGHVUVPars_flat_argon[2][j+1];
      	pars_ini[2] = temp1 + (temp2-temp1)*(theta-j*delta_angle)/delta_angle;
      }

      if(j >= 8){
	pars_ini[3] = fGHVUVPars_flat_argon[3][j];
      }
      else{
	double temp1 = fGHVUVPars_flat_argon[3][j];
      	double temp2 = fGHVUVPars_flat_argon[3][j+1];
      	pars_ini[3] = temp1 + (temp2-temp1)*(theta-j*delta_angle)/delta_angle;
      }
      s1 = interpolate( angulo, slopes1_flat_argon, theta, true);
      s2 = interpolate( angulo, slopes2_flat_argon, theta, true);
      s3 = interpolate( angulo, slopes3_flat_argon, theta, true);
    }
    else if (scintillation_type == 1) { // xenon
      pars_ini[0] = fGHVUVPars_flat_xenon[0][j];
      pars_ini[1] = fGHVUVPars_flat_xenon[1][j];
      pars_ini[2] = fGHVUVPars_flat_xenon[2][j];
      pars_ini[3] = fGHVUVPars_flat_xenon[3][j];
      s1 = interpolate( angulo, slopes1_flat_xenon, theta, true);
      s2 = interpolate( angulo, slopes2_flat_xenon, theta, true);
      s3 = interpolate( angulo, slopes3_flat_xenon, theta, true);
    }
    else {
      std::cout << "Error: Invalid scintillation type configuration." << endl;
      exit(1);
    }
  }
  // add border correction
  pars_ini[0] = pars_ini[0] + s1 * r_distance;
  pars_ini[1] = pars_ini[1] + s2 * r_distance;
  pars_ini[2] = pars_ini[2] + s3 * r_distance;
  pars_ini[3] = pars_ini[3];

  // calculate correction factor
  double GH_correction = GaisserHillas(distance, pars_ini);





  // apply correction
  double hits_vuv = 0 ;
  hits_vuv = gRandom->Poisson(GH_correction*hits_geo/cosine);
  return hits_vuv;
  if(GH_correction*hits_geo/cosine < 1) {
  	//hits_vuv = gRandom->Poisson(GH_correction*hits_geo/cosine);
	double rnd = gRandom->Rndm();
	if(rnd < GH_correction*hits_geo/cosine){
	 hits_vuv = 1;
	}
	else {
		hits_vuv = 0;
	}
  }
  else {
  	hits_vuv = (int) (GH_correction*hits_geo/cosine);
  }

  if(debug_2){
	  cout << "VUV hits: " << hits_vuv << endl;
	  std::cout<< " Nphotons_created = " << Nphotons_created << std::endl;
	  std::cout<< " hits_geo = " << hits_geo << std::endl;
	  std::cout << " distance = " << distance << std::endl;
	  std::cout<< " GH_correction/cosine = " << GH_correction/cosine << std::endl;
	  std::cout << " hits_vuv = " << hits_vuv << std::endl;
	  std::cout<< " GH_correction = " << GH_correction << std::endl;
	  std::cout<< " cosine = " << cosine << std::endl;
	  std::cout<< " hits_geo/Nphotons_created = " << hits_geo/Nphotons_created << std::endl;
	  std::cout<<" Solid Anlge = " << solid_angle << std::endl;
	  if(hits_vuv < 0) {return hits_geo;}
	  if(hits_vuv >= Nphotons_created || hits_vuv>= 1.2*hits_geo){return hits_geo;}
  }
  return hits_vuv;
}

// gaisser-hillas function definition
Double_t semi_analytic_hits::GaisserHillas(double x,double *par) {
  //This is the Gaisser-Hillas function
  Double_t X_mu_0=par[3];
  Double_t Normalization=par[0];
  Double_t Diff=par[1]-X_mu_0;
  Double_t Term=pow((x-X_mu_0)/Diff,Diff/par[2]);
  Double_t Exponential=TMath::Exp((par[1]-x)/par[2]);

  return ( Normalization*Term*Exponential);
}


// solid angle of rectanglular aperture calculation functions

double semi_analytic_hits::omega(const double &a, const double &b, const double &d) const{

  double aa = a/(2.0*d);
  double bb = b/(2.0*d);
  double aux = (1+aa*aa+bb*bb)/((1.+aa*aa)*(1.+bb*bb));
  return 4*std::acos(std::sqrt(aux));

}

double semi_analytic_hits::solid(const SiPM& out, const TVector3 &v) const{

  //v is the position of the track segment with respect to
  //the center position of the SiPM face
 // https://vixra.org/abs/2001.0603

   // We need to change this

   TVector3 v_copy = v;
   cout << "DET SIZE: " << out.h << " " << out.w << endl;

   // The hit is directly above the SiPM
  if( v.Y()==0.0 && v.Z()==0.0){
    return omega(out.w,out.h,v.X());
  }

  if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.Y())-out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(A+a),2*(B+b),d)-omega(2*A,2*(B+b),d)-omega(2*(A+a),2*B,d)+omega(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.Y())+out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(a-A),2*(b-B),d)+omega(2*A,2*(b-B),d)+omega(2*(a-A),2*B,d)+omega(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.Y())-out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(A+a),2*(b-B),d)-omega(2*A,2*(b-B),d)+omega(2*(A+a),2*B,d)-omega(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.Y())+out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(a-A),2*(B+b),d)-omega(2*(a-A),2*B,d)+omega(2*A,2*(B+b),d)-omega(2*A,2*B,d))/4.0;
    return to_return;
  }
  // error message if none of these cases, i.e. something has gone wrong!
  std::cout << "Warning: invalid solid angle call." << std::endl;
  return 0.0;
}


double semi_analytic_hits::solid(const acc& out, const TVector3 &v) const{

  //v is the position of the track segment with respect to
  //the center position of the arapuca window

  // arapuca plane fixed in x direction

  if( v.Y()==0.0 && v.Z()==0.0){
    return omega(out.w,out.h,v.X());
  }

  if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.Y())-out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    //std::cout << A <<" " << B << " " <<a <<" " << b <<" " << d << std::endl;
    double to_return = (omega(2*(A+a),2*(B+b),d)-omega(2*A,2*(B+b),d)-omega(2*(A+a),2*B,d)+omega(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.Y())+out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(a-A),2*(b-B),d)+omega(2*A,2*(b-B),d)+omega(2*(a-A),2*B,d)+omega(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.Y())-out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(A+a),2*(b-B),d)-omega(2*A,2*(b-B),d)+omega(2*(A+a),2*B,d)-omega(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.Y())+out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(a-A),2*(B+b),d)-omega(2*(a-A),2*B,d)+omega(2*A,2*(B+b),d)-omega(2*A,2*B,d))/4.0;
    return to_return;
  }
  // error message if none of these cases, i.e. something has gone wrong!
  std::cout << "Warning: invalid solid angle call." << std::endl;
  return 0.0;
}


// solid angle of circular aperture
double semi_analytic_hits::Disk_SolidAngle(double *x, double *p) {
  const double d = x[0];
  const double h = x[1];
  const double b = p[0];
  if(b <= 0. || d < 0. || h <= 0.) return 0.;
  const double aa = TMath::Sqrt(h*h/(h*h+(b+d)*(b+d)));
  if(d == 0) {
    return 2.*TMath::Pi()*(1.-aa);
  }
  const double bb = TMath::Sqrt(4*b*d/(h*h+(b+d)*(b+d)));
  const double cc = 4*b*d/((b+d)*(b+d));

  if(!_mathmore_loaded_) {
    if(gSystem->Load("libMathMore.so") < 0) {
      throw(std::runtime_error("Unable to load MathMore library"));
    }
    _mathmore_loaded_ = true;
  }
  if(TMath::Abs(ROOT::Math::comp_ellint_1(bb) - bb) < 1e-10 && TMath::Abs(ROOT::Math::comp_ellint_3(cc,bb) - cc) <1e-10) {
    throw(std::runtime_error("please do gSystem->Load(\"libMathMore.so\") before running Disk_SolidAngle for the first time!"));
  }
  if(d < b) {
    return 2.*TMath::Pi() - 2.*aa*(ROOT::Math::comp_ellint_1(bb) + TMath::Sqrt(1.-cc)*ROOT::Math::comp_ellint_3(cc,bb));
  }
  if(d == b) {
    return TMath::Pi() - 2.*aa*ROOT::Math::comp_ellint_1(bb);
  }
  if(d > b) {
    return 2.*aa*(TMath::Sqrt(1.-cc)*ROOT::Math::comp_ellint_3(cc,bb) - ROOT::Math::comp_ellint_1(bb));
  }

  return 0.;
}

double semi_analytic_hits::Disk_SolidAngle(double d, double h, double b) {
  double x[2] = { d, h };
  double p[1] = { b };
  if(!_mathmore_loaded_) {
    if(gSystem->Load("libMathMore.so") < 0) {
      throw(std::runtime_error("Unable to load MathMore library"));
    }
    _mathmore_loaded_ = true;
  }
  return Disk_SolidAngle(x,p);
}

double semi_analytic_hits::Omega_Dome_Model(const double distance, const double theta) const {
  // this function calculates the solid angle of a semi-sphere of radius b,
  // as a correction to the analytic formula of the on-axix solid angle,
  // as we move off-axis an angle theta. We have used 9-angular bins
  // with delta_theta width.

  // par0 = Radius correction close
  // par1 = Radius correction far
  // par2 = breaking distance betwween "close" and "far"

  double par0[9] = {0., 0., 0., 0., 0., 0.597542, 1.00872, 1.46993, 2.04221};
  double par1[9] = {0, 0, 0.19569, 0.300449, 0.555598, 0.854939, 1.39166, 2.19141, 2.57732};
  const double delta_theta = 10.;
  int j = int(theta/delta_theta);
  // 8" PMT radius
  const double b = 8*2.54/2.;
  // distance form which the model parameters break (empirical value)
  const double d_break = 5*b;//par2

  if(distance >= d_break) {
    double R_apparent_far = b - par1[j];
    return  (2*3.1416 * (1 - sqrt(1 - pow(R_apparent_far/distance,2))));

  }
  else {
    double R_apparent_close = b - par0[j];
    return (2*3.1416 * (1 - sqrt(1 - pow(R_apparent_close/distance,2))));
  }
}

double semi_analytic_hits::interpolate( const std::vector<double> &xData, const std::vector<double> &yData, double x, bool extrapolate ) {
  int size = xData.size();
  int i = 0;                                          // find left end of interval for interpolation
  if ( x >= xData[size - 2] )                         // special case: beyond right end
    {
      i = size - 2;
    }
  else
    {
      while ( x > xData[i+1] ) i++;
    }
  double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1]; // points on either side (unless beyond ends)
  if ( !extrapolate )                                                    // if beyond ends of array and not extrapolating
    {
      if ( x < xL ) yR = yL;
      if ( x > xR ) yL = yR;
    }
  double dydx = ( yR - yL ) / ( xR - xL );            // gradient
  return yL + dydx * ( x - xL );                      // linear interpolation
}
