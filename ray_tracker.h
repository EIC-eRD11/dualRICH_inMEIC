//#######2D optical ray tracker C++/ROOT code#########
//Code version 1.1
//Author Alessio Del Dotto 

//To compile and run in root framework: 
//.L ray_tracker.cpp++
//the file ray_tracker.h has to be in linked with ray_tracker.cpp

#ifndef __RAY_TRACKER_H_
#define __RAY_TRACKER_H_

#include <Riostream.h>
#include <TCut.h>
#include <TTree.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <TGraph.h>
#include <TFile.h>
#include <TMatrixD.h>
#include <TVectorD.h>

class ray_tracing {

 private:

 public:

  Double_t sp[5],sm[5],ang[5];
  Double_t center;
  Double_t yr, xr;

  Double_t xm_d, xp_d, ym_d, yp_d;
  Double_t m1,q1,mf,qf;
  Double_t gamma1;

  Double_t m1_d1, q1_d1, m2_d1, q2_d1;

  Double_t xm, xp, ym, yp;
  Double_t phi, alpha, alpha1;

  Double_t xm_1, xp_1, ym_1, yp_1;
  Double_t phi_1, alpha_1, alpha1_1;

  Double_t x1_det, y1_det, x2_det, y2_det;

  void set_detector();
  void spherical_mirror(Double_t xc, Double_t yc, Double_t r, Double_t m, Double_t q);
  void spherical_mirror1(Double_t xc, Double_t yc, Double_t r, Double_t m, Double_t q);
  void fresnel(Double_t xc, Double_t yc, Double_t r, Double_t m, Double_t q, Double_t n, Double_t n1, Double_t d1);
  void fresnel1(Double_t xc, Double_t yc, Double_t r, Double_t m, Double_t q, Double_t n, Double_t n1, Double_t d1);
  void trac(Double_t xc = 155., Double_t yc = 110., Double_t r = 245., Double_t xc1 = 155., Double_t yc1 = 95., Double_t r1 = 245., Double_t y_sep = 100., Double_t ang = 10.);
  void trac_sep(Double_t xc = 155., Double_t yc = 110., Double_t r = 245., Double_t xc1 = 155., Double_t yc1 = 95., Double_t r1 = 245., Double_t y_sep = 100., Double_t ang = 10., Int_t radiator = 0);
  void trac_fres(Double_t xc = 155., Double_t yc = 110., Double_t r = 245., Double_t xc1 = 155., Double_t yc1 = 95., Double_t r1 = 245., Double_t y_sep = 100., Double_t ang = 10.);
  void trac_fres_sep(Double_t xc = 155., Double_t yc = 110., Double_t r = 245., Double_t xc1 = 155., Double_t yc1 = 95., Double_t r1 = 245., Double_t y_sep = 100., Double_t ang = 10., Int_t radiator = 0);

};

#endif
