#include <Riostream.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>

#include <TMath.h>
#include <TEnv.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TTree.h>
#include <TDatabasePDG.h>
#include <TStyle.h>
#include <TLine.h>
#include <TEllipse.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TH1D.h>
#include <TH3F.h>
#include <TF1.h>
#include <TText.h>
#include <TLegend.h>
#include <TChain.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <TRandom.h>
#include <TMatrixTBase.h>
#include <TMatrixT.h>
#include "TPaveText.h"
#include <TGraphErrors.h>
#include <TFile.h>
#include <TTree.h>
#include <TSpline.h>
#include <TFrame.h>

using namespace std;

//*******eic_dual_rich class*******//

class eic_dual_rich {

 private:

 public:

  LongDouble_t sx,sy,sz;

  void acq(string input_filename, int pixel_eff, int emi_eff, int qe_eff, int tr_eff, int mag_lable, int shield_lable);
  void ind_ray(double Ex, double Ey, double Ez, double Dx, double Dy, double Dz, double cx, double cy, double cz);
  void acq_QE(string input_filename, int select);
  void err_plots(Int_t points, Int_t label);
  void nph_plots(Int_t points);

  Int_t phd_count = 0;
  Double_t ev_a[500], nm_a[500], qe_a[500];

  Double_t ch_er_a[100], ch_er_g[100], emi_er_a[100], emi_er_g[100], pix_er_a[100], pix_er_g[100], mg_er_a[100], mg_er_g[100], tr_er_a[100], tr_er_g[100];   
  Double_t angl[100];

  Double_t Nph_a[100], err_Nph_a[100], Nph_g[100], err_Nph_g[100];
  Double_t Nph_a_shield[100], err_Nph_a_shield[100];

  //double newp(double x, double th, double a, double d, double R);
  //double newpp(double x, double th, double a, double d, double R);

};


void eic_dual_rich::ind_ray(double Ex, double Ey, double Ez, double Dx, double Dy, double Dz, double cx, double cy, double cz){

  //Double_t cx = 100.;
  //Double_t cy = 0.;
  //Double_t cz = 134.5;

  //Double_t cx,cy,cz;
  //Double_t sx,sy,sz;
  LongDouble_t cex,cey,cez;
  LongDouble_t cdx,cdy,cdz;

  Int_t i,iwhere;

  LongDouble_t th,a,d;
  LongDouble_t x,dx;

  LongDouble_t y,y1;

  LongDouble_t eps = 0.00000000001;
  LongDouble_t R = 289.9;  //rich
  //LongDouble_t R = 299.9;    //rich3

  //eic_dual_rich *f = new eic_dual_rich();

  cex = -cx+Ex;
  cey = -cy+Ey;
  cez = -cz+Ez;

  cdx = -cx+Dx;
  cdy = -cy+Dy;
  cdz = -cz+Dz;

  //cout<<"ce is: "<<cex<<"  "<<cey<<"  "<<cez<<endl;
  //cout<<"cd is: "<<cdx<<"  "<<cdy<<"  "<<cdz<<endl;

  a = TMath::Sqrt(cex*cex+cey*cey+cez*cez);
  d = TMath::Sqrt(cdx*cdx+cdy*cdy+cdz*cdz);
  th = TMath::ACos((cdx*cex+cdy*cey+cdz*cez)/a/d);
  
  //cout<<"a,d,th is: "<<a<<"  "<<d<<"  "<<th<<endl;

  i = 0;
  x = th/2.;  
  //cout<<"x, sinx, sin(th-x) is:  "<<x<<"  "<<sin(x)<<"  "<<sin(th-x)<<endl;
  y = R*(a*sin(x)-d*sin(th-x))+a*d*sin(th-2*x);
  y1 = R*(a*cos(x)+d*cos(th-x))-2*a*d*cos(th-2*x);
  //cout<<"y, y1 is:  "<<y<<"  "<<y1<<endl;
  //dx = -(f->newp(x, th, a, d, R)/f->newpp(x, th, a, d, R));
  dx = -y/y1;

  while(TMath::Abs(dx)>eps && i<100){

    x+=dx;
    y = R*(a*sin(x)-d*sin(th-x))+a*d*sin(th-2*x);
    y1 = R*(a*cos(x)+d*cos(th-x))-2*a*d*cos(th-2*x);
    //dx = -(f->newp(x, th, a, d, R)/f->newpp(x, th, a, d, R));
    dx = -y/y1;
    i++;

  }

  //if(i>=100) cout<<"Not convergent"<<endl;

  if(i<100){
    sx = cx + (R*cos(x)/a-R*sin(x)/tan(th)/a)*cex + (R*sin(x)/d/sin(th))*cdx;
    sy = cy + (R*cos(x)/a-R*sin(x)/tan(th)/a)*cey + (R*sin(x)/d/sin(th))*cdy;
    sz = cz + (R*cos(x)/a-R*sin(x)/tan(th)/a)*cez + (R*sin(x)/d/sin(th))*cdz;

  }
}
/*
double newp(double x, double th, double a, double d, double R){

  Double_t y;

  y = R*(a*sin(x)-d*sin(th-x))+a*d*sin(th-2*x);

  return y;

}

double newpp(double x, double th, double a, double d, double R){

  Double_t y;

  y = R*(a*cos(x)+d*cos(th-x))-2*a*d*cos(th-2*x);

  return y;

}
*/

//*******acq data function*******//


void eic_dual_rich::acq(string input_filename, int pixel_eff, int emi_eff, int qe_eff, int tr_eff, int mag_lable, int shield_lable){

  gStyle->SetCanvasColor(0);
  //gStyle->SetTitleW(1.2);
  //gStyle->SetTitleH(0.2);
  gStyle->SetTitleSize(0.065,"xyzt");
  gStyle->SetTitleXOffset(1.);
  gStyle->SetTitleYOffset(0.7);
  gStyle->SetLabelSize(0.065,"XYZ");
  gStyle->SetPadTopMargin(.10);
  gStyle->SetPadLeftMargin(.14);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadBottomMargin(.15);
  gStyle->SetPalette(1);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetTitleFillColor(0); 
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.3); 
  //gStyle->SetOptStat(0);
  gStyle->SetStatBorderSize(0);
  gStyle->SetStatColor(0);
  gStyle->SetStatW(0.35);
  gStyle->SetStatH(0.35);


  LongDouble_t emi_x=0.;
  LongDouble_t emi_y=0.;
  LongDouble_t emi_z=0.;

  LongDouble_t dir_x=0.;
  LongDouble_t dir_y=0.;
  LongDouble_t dir_z=0.;

  LongDouble_t in_px=0.;
  LongDouble_t in_py=0.;
  LongDouble_t in_pz=0.;

  LongDouble_t theta_t=0.;
  LongDouble_t phi_t=0.;

  LongDouble_t emig_x=0.;
  LongDouble_t emig_y=0.;
  LongDouble_t emig_z=0.;

  //LongDouble_t mir_x = 100.; //rich
  LongDouble_t mir_x = 145.; //rich2
  //LongDouble_t mir_x = 130.; //rich3
  LongDouble_t mir_y = 0.;
  //LongDouble_t mir_z = 134.5; //rich
  LongDouble_t mir_z = 122.5; //rich2 in gemc coordinates 250+82.5-210.
  //LongDouble_t mir_z = 114.5;  //rich3
  LongDouble_t mir_xx = mir_x;

  Double_t mir_R = 289.9;

  LongDouble_t det_x=0.;
  LongDouble_t det_y=0.;
  LongDouble_t det_z=0.;

  Double_t theta1 = 0.;
  Double_t theta2 = 0.;

  Double_t ref_frac = 0.;
  Double_t ph_energy = 0.;

  Double_t det_xl=0.;
  Double_t det_yl=0.;
  Double_t det_zl=0.;

  Double_t det_xr=0.;
  Double_t det_yr=0.;
  Double_t det_zr=0.;

  Double_t det_xr1 = 0.; 
  Double_t det_yr1 = 0.;

  eic_dual_rich *f = new eic_dual_rich();

  Double_t ph_E_v[5],ph_E_v_g[11],ref_v[5],ref_v_g[11];
  ref_v[0]=1.01963;   
  ref_v[1]=1.01992; 
  ref_v[2]=1.02029; 
  ref_v[3]=1.02074; 
  ref_v[4]=1.02128;

  ref_v_g[0]=1.00048; 
  ref_v_g[1]=1.00048; 
  ref_v_g[2]=1.00049; 
  ref_v_g[3]=1.00049; 
  ref_v_g[4]=1.00050;
  ref_v_g[5]=1.00050;
  ref_v_g[6]=1.00051;
  ref_v_g[7]=1.00052;
  ref_v_g[8]=1.00052;
  ref_v_g[9]=1.00053;
  ref_v_g[10]=1.00054;

	for(Int_t p=0;p<5;p++){
	  ph_E_v[p] = 2.+p*0.5;
	}
	for(Int_t p=0;p<11;p++){
	  ph_E_v_g[p] = 2.+p*0.5;
	}
	

  TFile *file=new TFile(input_filename.c_str());
  if (file->IsZombie()) {
    cout << "Error opening file" << input_filename << endl;
    exit(-1);
  }
  else cout << "open file " << input_filename << endl;

  Double_t *Ne_ph_a = new Double_t[1];
  Double_t *Npi_ph_a = new Double_t[1];
  Double_t *Nk_ph_a = new Double_t[1];
  Double_t *Np_ph_a = new Double_t[1];

  Double_t *Ne_ph_g = new Double_t[1];
  Double_t *Npi_ph_g = new Double_t[1];
  Double_t *Nk_ph_g = new Double_t[1];
  Double_t *Np_ph_g = new Double_t[1];

  for(Int_t n=0;n<1;n++){
    Ne_ph_a[n] = 0;
    Npi_ph_a[n] = 0;
    Nk_ph_a[n] = 0;
    Np_ph_a[n] = 0;

    Ne_ph_g[n] = 0;
    Npi_ph_g[n] = 0;
    Nk_ph_g[n] = 0;
    Np_ph_g[n] = 0;
  }

  TCanvas **cch_ph;
  cch_ph = new TCanvas*[64];

    TH1F **nph_e_aerogel_h;
  nph_e_aerogel_h=new TH1F*[64];
  for(Int_t i=0;i<64;i++){
    nph_e_aerogel_h[i] = new TH1F(Form("nph_e_aerogel_h%d",i),"",2000,0,200.);
  }

    TH1F **nph_pi_aerogel_h;
  nph_pi_aerogel_h=new TH1F*[64];
  for(Int_t i=0;i<64;i++){
    nph_pi_aerogel_h[i] = new TH1F(Form("nph_pi_aerogel_h%d",i),"",2000,0,200.);
  }

    TH1F **nph_k_aerogel_h;
  nph_k_aerogel_h=new TH1F*[64];
  for(Int_t i=0;i<64;i++){
    nph_k_aerogel_h[i] = new TH1F(Form("nph_k_aerogel_h%d",i),"",2000,0,200.);
  }

    TH1F **nph_p_aerogel_h;
  nph_p_aerogel_h=new TH1F*[64];
  for(Int_t i=0;i<64;i++){
    nph_p_aerogel_h[i] = new TH1F(Form("nph_p_aerogel_h%d",i),"",2000,0,200.);
  }

    TH1F **nph_e_gas_h;
  nph_e_gas_h=new TH1F*[64];
  for(Int_t i=0;i<64;i++){
    nph_e_gas_h[i] = new TH1F(Form("nph_e_gas_h%d",i),"",1000,0,400.);
  }

  TH1F **nph_pi_gas_h;
  nph_pi_gas_h=new TH1F*[64];
  for(Int_t i=0;i<64;i++){
    nph_pi_gas_h[i] = new TH1F(Form("nph_pi_gas_h%d",i),"",1000,0,1000.);
  }

    TH1F **nph_k_gas_h;
  nph_k_gas_h=new TH1F*[64];
  for(Int_t i=0;i<64;i++){
    nph_k_gas_h[i] = new TH1F(Form("nph_k_gas_h%d",i),"",1000,0,400.);
  }

    TH1F **nph_p_gas_h;
  nph_p_gas_h=new TH1F*[64];
  for(Int_t i=0;i<64;i++){
    nph_p_gas_h[i] = new TH1F(Form("nph_p_gas_h%d",i),"",1000,0,400.);
  }

  TH1F *id_h = new TH1F("id_h","",1000,0.,100.);
  TH1F *hit_h = new TH1F("hit_h","",1000,0.,100.);

  TH1F *x_h = new TH1F("x_h","",8000,-4000.,4000.);
  TH1F *y_h = new TH1F("y_h","",8000,-4000.,4000.);
  TH1F *z_h = new TH1F("z_h","",8000,-4000.,4000.);

  //TH1F *tc_pi_aerogel_h = new TH1F("tc_pi_aerogel_h","",1000,0.,0.3);
  //TH1F *tc_pi_gas_h = new TH1F("tc_pi_gas_h","",1000,0.,0.3);
  //TH1F *tc_k_aerogel_h = new TH1F("tc_k_aerogel_h","",1000,0.,0.3);
  //TH1F *tc_k_gas_h = new TH1F("tc_k_gas_h","",1000,0.,0.3);

  TCanvas **cch;
  cch = new TCanvas*[64];
  //cch[0] = new TCanvas("cch","",1000,500); 
  TH1F **tc_e_aerogel_h;
  tc_e_aerogel_h=new TH1F*[64];
  for(Int_t i=0;i<64;i++){
    tc_e_aerogel_h[i] = new TH1F(Form("tc_e_aerogel_h%d",i),"",10000,0,1.);
  }
  TH1F **tc_e_gas_h;
  tc_e_gas_h=new TH1F*[64];
  for(Int_t i=0;i<64;i++){
    tc_e_gas_h[i] = new TH1F(Form("tc_e_gas_h%d",i),"",10000,0,0.3);
  }
  TH1F **tc_pi_aerogel_h;
  tc_pi_aerogel_h=new TH1F*[64];
  for(Int_t i=0;i<64;i++){
    tc_pi_aerogel_h[i] = new TH1F(Form("tc_pi_aerogel_h%d",i),"",10000,0,1.);
  }
  TH1F **tc_pi_gas_h;
  tc_pi_gas_h=new TH1F*[64];
  for(Int_t i=0;i<64;i++){
    tc_pi_gas_h[i] = new TH1F(Form("tc_pi_gas_h%d",i),"",10000,0.,0.3);
  }
  TH1F **tc_k_aerogel_h;
  tc_k_aerogel_h=new TH1F*[64];
  for(Int_t i=0;i<64;i++){
    tc_k_aerogel_h[i] = new TH1F(Form("tc_k_aerogel_h%d",i),"",10000,0,1.);
  }
  TH1F **tc_k_gas_h;
  tc_k_gas_h=new TH1F*[64];
  for(Int_t i=0;i<64;i++){
    tc_k_gas_h[i] = new TH1F(Form("tc_k_gas_h%d",i),"",10000,0,0.3);
  }
  TH1F **tc_p_aerogel_h;
  tc_p_aerogel_h=new TH1F*[64];
  for(Int_t i=0;i<64;i++){
    tc_p_aerogel_h[i] = new TH1F(Form("tc_p_aerogel_h%d",i),"",10000,0,1.);
  }
  TH1F **tc_p_gas_h;
  tc_p_gas_h=new TH1F*[64];
  for(Int_t i=0;i<64;i++){
    tc_p_gas_h[i] = new TH1F(Form("tc_p_gas_h%d",i),"",10000,0,0.3);
  }

  TH1F **tc_pi_aerogel_ene_h;
  tc_pi_aerogel_ene_h=new TH1F*[64];
  for(Int_t i=0;i<64;i++){
    tc_pi_aerogel_ene_h[i] = new TH1F(Form("tc_pi_aerogel_ene_h%d",i),"",1000,0,1000.);
  }
  TH1F **tc_pi_aerogel_ene_bk_h;
  tc_pi_aerogel_ene_bk_h=new TH1F*[64];
  for(Int_t i=0;i<64;i++){
    tc_pi_aerogel_ene_bk_h[i] = new TH1F(Form("tc_pi_aerogel_ene_bk_h%d",i),"",1000,0,1000.);
  }
  

  TCanvas **cch_d;
  cch_d = new TCanvas*[64];

  TCanvas **dist_ph;
  dist_ph = new TCanvas*[64];

  TH1F **dist_pi_aerogel_h;
  dist_pi_aerogel_h=new TH1F*[64];
  for(Int_t i=0;i<64;i++){
    dist_pi_aerogel_h[i] = new TH1F(Form("dist_pi_aerogel_h%d",i),"",1000000,0.,350.);
  }
  TH1F **dist_pi_gas_h;
  dist_pi_gas_h=new TH1F*[64];
  for(Int_t i=0;i<64;i++){
    dist_pi_gas_h[i] = new TH1F(Form("dist_pi_gas_h%d",i),"",1000000,0.,280.);
  }

  //TH2F *ph_det = new TH2F("ph_det","", 1666, -250., 250., 1666, -250., 250.);
  //TH2F *ph_det_sph = new TH2F("ph_det_sph","", 1666, -3.14, 3.14, 1666, -3.14, 3.14);
  //TH2F *ph_det1 = new TH2F("ph_det1","", 1200, -200., 200., 1200, -200., 200.);
  //TH1F *ph_det_x = new TH1F("ph_det_x","", 1000, -1000., 1000.);
  //TH1F *ph_det_y = new TH1F("ph_det_y","", 1000, -1000., 1000.);
  /*
  TH1F *emi_x_h = new TH1F("emi_x_h","", 10000, -1000., 1000.);
  TH1F *emi_y_h = new TH1F("emi_y_h","", 10000, -1000., 1000.);
  TH1F *emi_z_h = new TH1F("emi_z_h","", 10000, 0., 1000.);
  TH1F *emig_x_h = new TH1F("emig_x_h","", 10000, -1000., 1000.);
  TH1F *emig_y_h = new TH1F("emig_y_h","", 10000, -1000., 1000.);
  TH1F *emig_z_h = new TH1F("emig_z_h","", 10000, 0., 1000.);
  TH1F *emig1_x_h = new TH1F("emig1_x_h","", 10000, -1000., 1000.);
  TH1F *emig1_y_h = new TH1F("emig1_y_h","", 10000, -1000., 1000.);
  TH1F *emig1_z_h = new TH1F("emig1_z_h","", 10000, 0., 1000.);
  */
  Double_t xbin = 0.;
  Double_t ybin = 0.;

  TTree *generated = (TTree*) file->Get("generated");
  vector <double> *gen_pid=0;
  vector <double> *gen_px=0,*gen_py=0,*gen_pz=0,*gen_vx=0,*gen_vy=0,*gen_vz=0;
  generated->SetBranchAddress("pid",&gen_pid);
  generated->SetBranchAddress("px",&gen_px);
  generated->SetBranchAddress("py",&gen_py);
  generated->SetBranchAddress("pz",&gen_pz);
  generated->SetBranchAddress("vx",&gen_vx);
  generated->SetBranchAddress("vy",&gen_vy);
  generated->SetBranchAddress("vz",&gen_vz);

  //generated->Draw("pz");

  TTree *eic_rich = (TTree*) file->Get("eic_rich");

  vector<double> *eic_rich_id=0,*eic_rich_hitn=0;
  vector<double> *eic_rich_pid=0,*eic_rich_mpid=0,*eic_rich_tid=0,*eic_rich_mtid=0,*eic_rich_otid=0;
  vector<double> *eic_rich_trackE=0,*eic_rich_totEdep=0,*eic_rich_avg_x=0,*eic_rich_avg_y=0,*eic_rich_avg_z=0,*eic_rich_avg_lx=0,*eic_rich_avg_ly=0,*eic_rich_avg_lz=0,*eic_rich_px=0,*eic_rich_py=0,*eic_rich_pz=0,*eic_rich_vx=0,*eic_rich_vy=0,*eic_rich_vz=0,*eic_rich_mvx=0,*eic_rich_mvy=0,*eic_rich_mvz=0,*eic_rich_avg_t=0;
  vector<double> *eic_rich_in_px=0,*eic_rich_in_py=0,*eic_rich_in_pz=0,*eic_rich_in_x=0,*eic_rich_in_y=0,*eic_rich_in_z=0,*eic_rich_in_t=0,*eic_rich_out_px=0,*eic_rich_out_py=0,*eic_rich_out_pz=0,*eic_rich_out_x=0,*eic_rich_out_y=0,*eic_rich_out_z=0,*eic_rich_out_t=0,*eic_rich_nsteps=0;
  eic_rich->SetBranchAddress("hitn",&eic_rich_hitn);  
  eic_rich->SetBranchAddress("id",&eic_rich_id);
  eic_rich->SetBranchAddress("pid",&eic_rich_pid);
  eic_rich->SetBranchAddress("mpid",&eic_rich_mpid);
  eic_rich->SetBranchAddress("tid",&eic_rich_tid);
  eic_rich->SetBranchAddress("mtid",&eic_rich_mtid);
  eic_rich->SetBranchAddress("otid",&eic_rich_otid);
  eic_rich->SetBranchAddress("trackE",&eic_rich_trackE);
  eic_rich->SetBranchAddress("totEdep",&eic_rich_totEdep);
  eic_rich->SetBranchAddress("avg_x",&eic_rich_avg_x);
  eic_rich->SetBranchAddress("avg_y",&eic_rich_avg_y);
  eic_rich->SetBranchAddress("avg_z",&eic_rich_avg_z);
  eic_rich->SetBranchAddress("avg_lx",&eic_rich_avg_lx);
  eic_rich->SetBranchAddress("avg_ly",&eic_rich_avg_ly);
  eic_rich->SetBranchAddress("avg_lz",&eic_rich_avg_lz);
  eic_rich->SetBranchAddress("px",&eic_rich_px);
  eic_rich->SetBranchAddress("py",&eic_rich_py);
  eic_rich->SetBranchAddress("pz",&eic_rich_pz);
  eic_rich->SetBranchAddress("vx",&eic_rich_vx);
  eic_rich->SetBranchAddress("vy",&eic_rich_vy);
  eic_rich->SetBranchAddress("vz",&eic_rich_vz);
  eic_rich->SetBranchAddress("mvx",&eic_rich_mvx);
  eic_rich->SetBranchAddress("mvy",&eic_rich_mvy);
  eic_rich->SetBranchAddress("mvz",&eic_rich_mvz);
  eic_rich->SetBranchAddress("avg_t",&eic_rich_avg_t);
  eic_rich->SetBranchAddress("in_px",&eic_rich_in_px);
  eic_rich->SetBranchAddress("in_py",&eic_rich_in_py);
  eic_rich->SetBranchAddress("in_pz",&eic_rich_in_pz);
  eic_rich->SetBranchAddress("in_x",&eic_rich_in_x);
  eic_rich->SetBranchAddress("in_y",&eic_rich_in_y);
  eic_rich->SetBranchAddress("in_z",&eic_rich_in_z);
  eic_rich->SetBranchAddress("in_t",&eic_rich_in_t);
  eic_rich->SetBranchAddress("out_px",&eic_rich_out_px);
  eic_rich->SetBranchAddress("out_py",&eic_rich_out_py);
  eic_rich->SetBranchAddress("out_pz",&eic_rich_out_pz);
  eic_rich->SetBranchAddress("out_x",&eic_rich_out_x);
  eic_rich->SetBranchAddress("out_y",&eic_rich_out_y);
  eic_rich->SetBranchAddress("out_z",&eic_rich_out_z);
  eic_rich->SetBranchAddress("out_t",&eic_rich_out_t);
  eic_rich->SetBranchAddress("nsteps",&eic_rich_nsteps);

  cout<<"eic_rich entries: "<<eic_rich->GetEntries()<<endl;
  //eic_rich->Draw("id");

  TRandom rran;
  rran.SetSeed(0);

  Double_t qe_p = 0.;

  Double_t esx,esy,esz;
  Double_t es;
  Double_t dis_aerogel, dis_gas;
  Double_t p;
  Double_t theta_c;
  Int_t index = 0;
  Double_t temp_index;
  Int_t temp = 0;
  Double_t m[4]={0.000511, 0.13957018, 0.493677, 0.938272};
  Double_t the_a[4], the_g[4];
  Double_t mom_p;
		
  TH1F *mmm    = new TH1F("mmm","",1000,12,14);		

  for(Int_t i=0;(i<eic_rich->GetEntries());i++){

    TH2F *ph_det = new TH2F("ph_det","", 1666, -250., 250., 1666, -250., 250.);

    generated->GetEntry(i);
    eic_rich->GetEntry(i);

    if(sqrt(gen_px->at(0)*gen_px->at(0)+gen_py->at(0)*gen_py->at(0)+gen_pz->at(0)*gen_pz->at(0))<8000.){
      //index = (sqrt(gen_px->at(0)*gen_px->at(0)+gen_py->at(0)*gen_py->at(0)+gen_pz->at(0)*gen_pz->at(0))-999.)/2000.;  //it depends on the imput file momentum spacing!!!
    }
    else{
      //index = 4. + (sqrt(gen_px->at(0)*gen_px->at(0)+gen_py->at(0)*gen_py->at(0)+gen_pz->at(0)*gen_pz->at(0))-9999.)/3000.; 
    }
    //if(i<2600) index = i/100;
    //if(i>=2600 && i<5200) index = (i-2600)/100;
    //if(i>=5200 && i<7800) index = (i-5200)/100;
    //if(i>=7800 && i<10400) index = (i-7800)/100;
    temp_index = (TMath::ACos((gen_pz->at(0))/sqrt(gen_px->at(0)*gen_px->at(0)+gen_py->at(0)*gen_py->at(0)+gen_pz->at(0)*gen_pz->at(0)))*180./TMath::Pi());

    if(index<4) mom_p = 1. + index*2.;
    if(index>=4) mom_p = 10. + (index-4)*3.;
    if(index==5 && i>=5200 && i<7800) {
      mmm->Fill(sqrt(gen_px->at(0)*gen_px->at(0)+gen_py->at(0)*gen_py->at(0)+gen_pz->at(0)*gen_pz->at(0))/1000.);
    }
    mom_p = 30.;
    
    for(Int_t idp=0;idp<4;idp++){
      the_a[idp]=acos(sqrt(mom_p*mom_p+m[idp]*m[idp])/mom_p/1.02);
      the_g[idp]=acos(sqrt(mom_p*mom_p+m[idp]*m[idp])/mom_p/1.00086);
      //if(index==0)cout<<the_a[idp]<<"  "<<the_g[idp]<<endl;
    }
    
    //index=0;

    /*
    if(temp_index<7) index = 0;
    if(temp_index<11 && temp_index>6) index = 1;
    if(temp_index<15 && temp_index>10) index = 2;
    if(temp_index<19 && temp_index>14) index = 3;
    if(temp_index<23 && temp_index>18) index = 4;
    */
    /*if(temp_index<5) index = 0;
    if(temp_index<7 && temp_index>4) index = 1;
    if(temp_index<9 && temp_index>6) index = 2;
    if(temp_index<11 && temp_index>8) index = 3;
    if(temp_index<13 && temp_index>10) index = 4;
    if(temp_index<15 && temp_index>12) index = 5;
    if(temp_index<17 && temp_index>14) index = 6;
    if(temp_index<19 && temp_index>16) index = 7;
    if(temp_index<21 && temp_index>18) index = 8;	
    if(temp_index<23 && temp_index>20) index = 9;
    if(temp_index<25 && temp_index>22) index = 10;
    if(temp_index<27 && temp_index>24) index = 11;*/
    
    if(temp_index<7) index = 0;
    if(temp_index<11 && temp_index>6) index = 1;
    if(temp_index<16 && temp_index>11) index = 2;
    if(temp_index<21 && temp_index>17) index = 3;
    if(temp_index<26 && temp_index>22) index = 4;
    
    if(temp_index <= 3.1 + 2*temp) temp++;
    
    //index = temp -1; 
    //if(index<1) index = 0;

    //cout<<index<<endl;

    //cout<<"px,py,pz,i, index: "<<gen_px->at(0)<<"  "<<gen_py->at(0)<<"  "<<gen_pz->at(0)<<"  "<<i<<"  "<<index<<"  "<<gen_pid->at(0)<<endl;
    //cout<<eic_rich_id->size()<<endl;

    //cout<<eic_rich_id->size()<<"  "<<eic_rich_pid->size()<<endl;
    if((gen_pid->at(0)) == 211 || (gen_pid->at(0)) == 321  || (gen_pid->at(0)) == 2212 || (gen_pid->at(0)) == 11){
      for(Int_t j=0;j<eic_rich_id->size();j++){

	if(((eic_rich_pid->at(j))==211 || (eic_rich_pid->at(j))==321 || (eic_rich_pid->at(j))==2212 || (eic_rich_pid->at(j))==11) && (eic_rich_mtid->at(j))==0  && (eic_rich_id->at(j))>=11 && (eic_rich_id->at(j))<20){
	  //if(i==0 || i==30 || i==60 || i==90 || i==120 || i==150 || i==180 || i==210 || i==240 || i==270 || i==300 || i==330){
	  //if(i==0 || i==5 || i==10 || i==15 || i==20){
	  in_px = eic_rich_in_px->at(j)/sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j));
	  in_py = eic_rich_in_py->at(j)/sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j));
	  in_pz = eic_rich_in_pz->at(j)/sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j));

	  if(tr_eff!=0){
	  theta_t = TMath::ACos(in_pz);
	  phi_t = TMath::ATan2(in_py,in_px);
	  theta_t = rran.Gaus(theta_t,0.0005);
	  phi_t = rran.Gaus(phi_t,0.0005);
	  in_pz = TMath::Cos(theta_t);
	  in_px = TMath::Sin(theta_t)*TMath::Cos(phi_t);
	  in_py = TMath::Sin(theta_t)*TMath::Sin(phi_t);
	  }
	  //}
	  //cout<<in_px<<"  "<<in_py<<"  "<<in_pz<<endl;
	  //cout<<eic_rich_tid->at(j)<<"  "<<eic_rich_mtid->at(j)<<endl;

	}

      //id_h->Fill(eic_rich_id->at(j));
      if((eic_rich_pid->at(j))==0 && (eic_rich_id->at(j))>=11 && (eic_rich_id->at(j))<20 && emi_eff==0){
	//cout<<eic_rich_vx->at(j)*0.1<<"  "<<eic_rich_vy->at(j)*0.1<<"  "<<eic_rich_vz->at(j)*0.1<<endl;
	emi_x = eic_rich_vx->at(j)*0.1;
	emi_y = eic_rich_vy->at(j)*0.1;
	emi_z = eic_rich_vz->at(j)*0.1;

	/*if(emi_z<254.){
	emi_x_h->Fill(emi_x);
	emi_y_h->Fill(emi_y);
	emi_z_h->Fill(emi_z);
	//tc_pi_aerogel_ene_bk_h[index]->Fill(eic_rich_trackE->at(j)*1000000.);
	}*/
	//emi_x = eic_rich_out_x->at(j)*0.1;
	//emi_y = eic_rich_out_y->at(j)*0.1;
	//emi_z = eic_rich_out_z->at(j)*0.1;

	//cout<<"emission aerogel pre: "<<emi_x<<"  "<<emi_y<<"  "<<emi_z<<"  "<<j<<endl;

      }
      if((eic_rich_pid->at(j))==0 && (eic_rich_id->at(j))>=21 && (eic_rich_id->at(j))<30 && emi_eff==0){
	//cout<<eic_rich_vx->at(j)*0.1<<"  "<<eic_rich_vy->at(j)*0.1<<"  "<<eic_rich_vz->at(j)*0.1<<endl;
	emig_x = eic_rich_vx->at(j)*0.1;
	emig_y = eic_rich_vy->at(j)*0.1;
	emig_z = eic_rich_vz->at(j)*0.1;
	/*if(emig_z > 254.){
	emig_x_h->Fill(emig_x);
	emig_y_h->Fill(emig_y);
	emig_z_h->Fill(emig_z);
	}*/
      }
      if((eic_rich_pid->at(j))==(gen_pid->at(0)) && (eic_rich_id->at(j))>=11 && (eic_rich_id->at(j))<20 && emi_eff!=0 && (eic_rich_mtid->at(j))==0){
	//cout<<"pid is: "<<eic_rich_pid->at(j)<<"  "<<j<<endl;
	//cout<<"id,hitn: "<<eic_rich_id->at(j)<<"  "<<eic_rich_hitn->at(j)<<endl;
	//cout<<"mpid,tid,otid: "<<eic_rich_mpid->at(j)<<"  "<<eic_rich_tid->at(j)<<"  "<<eic_rich_otid->at(j)<<endl;
	//cout<<"px,py,pz: "<<eic_rich_px->at(j)<<"  "<<eic_rich_py->at(j)<<"  "<<eic_rich_pz->at(j)<<endl;
	//cout<<"av_x,av_y,av_z: "<<eic_rich_avg_x->at(j)<<"  "<<eic_rich_avg_y->at(j)<<"  "<<eic_rich_avg_z->at(j)<<endl;
	//if(i==0 || i==30 || i==60 || i==90 || i==120 || i==150 || i==180 || i==210 || i==240 || i==270 || i==300 || i==330){
	emi_x = rran.Gaus(eic_rich_avg_x->at(j)*0.1,0.01);
	emi_y = rran.Gaus(eic_rich_avg_y->at(j)*0.1,0.01);
	emi_z = eic_rich_avg_z->at(j)*0.1;
	//}
	//cout<<"emission aerogel: "<<emi_x<<"  "<<emi_y<<"  "<<emi_z<<endl;

      }
      if((eic_rich_pid->at(j))==(gen_pid->at(0)) && (eic_rich_id->at(j))>=21 && (eic_rich_id->at(j))<30 && emi_eff!=0 && (eic_rich_mtid->at(j))==0){
	//cout<<"pid is: "<<eic_rich_pid->at(j)<<"  "<<j<<endl;
	//cout<<"id,hitn: "<<eic_rich_id->at(j)<<"  "<<eic_rich_hitn->at(j)<<endl;
	//cout<<"mpid,tid,otid: "<<eic_rich_mpid->at(j)<<"  "<<eic_rich_tid->at(j)<<"  "<<eic_rich_otid->at(j)<<endl;
	//cout<<"px,py,pz: "<<eic_rich_px->at(j)<<"  "<<eic_rich_py->at(j)<<"  "<<eic_rich_pz->at(j)<<endl;
	//cout<<"av_x,av_y,av_z: "<<eic_rich_avg_x->at(j)<<"  "<<eic_rich_avg_y->at(j)<<"  "<<eic_rich_avg_z->at(j)<<endl;
	//if(i==0 || i==30 || i==60 || i==90 || i==120 || i==150 || i==180 || i==210 || i==240 || i==270 || i==300 || i==330){
	emig_x = rran.Gaus(eic_rich_avg_x->at(j)*0.1,0.01);
	emig_y = rran.Gaus(eic_rich_avg_y->at(j)*0.1,0.01);
	emig_z = eic_rich_avg_z->at(j)*0.1;
	//}
	//emig_x = 19.3;
	//emig_y = 0.;
	//emig_z = 340.1;	

	//cout<<"emission gas: "<<emig_x<<"  "<<emig_y<<"  "<<emig_z<<endl;
      }
      //cout<<eic_rich_id->at(j)<<endl;
      //if((eic_rich_id->at(j))==41 && (eic_rich_vz->at(j))<2540){   

      //&& eic_rich_vz->at(j)<2540.
      if((eic_rich_id->at(j)>40) && (eic_rich_id->at(j)<50) && (eic_rich_mpid->at(j)==211 || eic_rich_mpid->at(j)==321 || eic_rich_mpid->at(j)==2212 || eic_rich_mpid->at(j)==11)){
	Int_t phi_count = eic_rich_id->at(j);
	phi_count = phi_count -41;
        mir_x = mir_xx*TMath::Cos(TMath::Pi()*-60.*phi_count/180.);
	mir_y = mir_xx*TMath::Sin(TMath::Pi()*-60.*phi_count/180.);
 
	//cout<<eic_rich_id->at(j)<<endl;
	//cout<<"phi_count is: "<<phi_count<<" mir x "<<mir_x<<" mir y "<<mir_y<<" mir_z  "<<mir_z<<endl;
	//cout<<"vx,vy,vz: "<<eic_rich_vx->at(j)*0.1<<"  "<<eic_rich_vy->at(j)*0.1<<"  "<<eic_rich_vz->at(j)*0.1<<endl;
	//cout<<eic_rich_hitn->at(j)<<endl;
	//}
	//if((eic_rich_hitn->at(j))==255){
	//if((eic_rich_mpid->at(j)) != 211)cout<<eic_rich_mpid->at(j)<<endl;

	x_h->Fill(eic_rich_avg_x->at(j));
	y_h->Fill(eic_rich_avg_y->at(j));
	z_h->Fill(eic_rich_avg_z->at(j));

	ph_det->Fill(eic_rich_avg_lx->at(j)*0.1,eic_rich_avg_ly->at(j)*0.1);
	//ph_det_sph->Fill(TMath::ACos(eic_rich_avg_lz->at(j)*0.1/279.9),TMath::ATan2(eic_rich_avg_ly->at(j)*0.1,eic_rich_avg_lx->at(j)*0.1));

	xbin = ph_det->GetXaxis()->FindBin(eic_rich_avg_lx->at(j)*0.1);
	ybin = ph_det->GetYaxis()->FindBin(eic_rich_avg_ly->at(j)*0.1);

	//ph_det1->Fill(sqrt(eic_rich_avg_lx->at(j)*0.1*eic_rich_avg_lx->at(j)*0.1+eic_rich_avg_ly->at(j)*0.1*eic_rich_avg_ly->at(j)*0.1),TMath::ATan2(eic_rich_avg_ly->at(j),eic_rich_avg_lx->at(j))*180./TMath::Pi());
	//ph_det_x->Fill(eic_rich_avg_lx->at(j)*0.1);
	//ph_det_y->Fill(eic_rich_avg_ly->at(j)*0.1);

	/*emi_x = eic_rich_vx->at(j)*0.1;  //emission corrispondence with the previous above
	emi_y = eic_rich_vy->at(j)*0.1;
	emi_z = eic_rich_vz->at(j)*0.1;*/

	det_x = eic_rich_in_x->at(j)*0.1;
	det_y = eic_rich_in_y->at(j)*0.1;
	det_z = eic_rich_in_z->at(j)*0.1;

	//det_xl = eic_rich_avg_lx->at(j)*0.1;
	//det_yl = eic_rich_avg_ly->at(j)*0.1;

	det_xl = ph_det->GetXaxis()->GetBinCenter(xbin);
	det_yl = ph_det->GetYaxis()->GetBinCenter(ybin);
	det_zl = eic_rich_avg_lz->at(j)*0.1;

	//det_xr = det_x*TMath::Cos(TMath::Pi()*10./180.) + (det_z-254.5)*TMath::Sin(TMath::Pi()*10./180.);
	//det_zr = -det_x*TMath::Sin(TMath::Pi()*10./180.) + (det_z-254.5)*TMath::Cos(TMath::Pi()*10./180.);
	//**********Plane*********************//
	/*det_xr = det_xl*TMath::Cos(-TMath::Pi()*0./180.) + (det_zl)*TMath::Sin(-TMath::Pi()*0./180.);
	det_zr = -det_xl*TMath::Sin(-TMath::Pi()*0./180.) + (det_zl)*TMath::Cos(-TMath::Pi()*0./180.);

	det_xr1 = det_xr*TMath::Cos(TMath::Pi()*60.*phi_count/180.) + det_yl*TMath::Sin(TMath::Pi()*60.*phi_count/180.);
	det_yr1 = -det_xr*TMath::Sin(TMath::Pi()*60.*phi_count/180.) + det_yl*TMath::Cos(TMath::Pi()*60.*phi_count/180.);*/
	//**********Mirror*********************//
	det_xr = det_xl*TMath::Cos(-TMath::ASin(mir_xx/mir_R)) + (det_zl)*TMath::Sin(-TMath::ASin(mir_xx/mir_R));
	det_zr = -det_xl*TMath::Sin(-TMath::ASin(mir_xx/mir_R)) + (det_zl)*TMath::Cos(-TMath::ASin(mir_xx/mir_R))+mir_z;

	det_xr1 = det_xr*TMath::Cos(TMath::Pi()*60.*phi_count/180.) + det_yl*TMath::Sin(TMath::Pi()*60.*phi_count/180.)+mir_x;
	det_yr1 = -det_xr*TMath::Sin(TMath::Pi()*60.*phi_count/180.) + det_yl*TMath::Cos(TMath::Pi()*60.*phi_count/180.)+mir_y;

	//cout<<det_x<<"  "<<det_y<<"  "<<det_z<<endl;
	//cout<<det_xl<<"  "<<det_yl<<"  "<<det_zl<<endl;
	//cout<<det_xr1<<"  "<<det_yr1<<endl;

	//cout<<eic_rich_avg_z->at(j)<<endl;
	if(pixel_eff == 0) f->ind_ray(emi_x,emi_y,emi_z,det_x,det_y,det_z,mir_x,mir_y,mir_z);
	if(pixel_eff != 0) f->ind_ray(emi_x,emi_y,emi_z,det_xr1,det_yr1,det_z,mir_x,mir_y,mir_z);

	//cout<<"emission aerogel: "<<emi_x<<"  "<<emi_y<<"  "<<emi_z<<"  "<<j<<endl;
	//cout<<"detector point: "<<det_x<<"  "<<det_y<<"  "<<det_z<<endl;
	//cout<<"mirror center: "<<mir_x<<"  "<<mir_y<<"  "<<mir_z<<endl;
	//cout<<emi_x<<"  "<<emi_y<<"  "<<emi_z<<endl;
	//if((eic_rich_id->at(j))==46)cout<<det_x<<"  "<<det_y<<"  "<<det_z<<endl;
	//cout<<f->sx<<"  "<<f->sy<<"  "<<f->sz<<endl;

	esx = f->sx - emi_x;
	esy = f->sy - emi_y;
	esz = f->sz - emi_z;

	es = sqrt(esx*esx+esy*esy+esz*esz);
	dis_aerogel = es+sqrt((f->sx-det_x)*(f->sx-det_x)+(f->sy-det_y)*(f->sy-det_y)+(f->sz-det_z)*(f->sz-det_z));

	if((gen_pid->at(0))==211 && eic_rich_vz->at(j)*0.1<254.){
	dist_pi_aerogel_h[index]->Fill(dis_aerogel);
	}

	esx = esx/es;
	esy = esy/es;
	esz = esz/es;

	//cout<<eic_rich_trackE->at(j)<<endl;

	ph_energy = eic_rich_trackE->at(j)*1000000.;

	ref_frac = 1.000823/1.02;

	//TSpline3 *ref_ind;
	//ref_ind = new TSpline3("",ph_E_v,ref_v,5);

	//TSpline3 *ref_ind_g;
	//ref_ind_g = new TSpline3("",ph_E_v_g,ref_v_g,11);

	//ref_frac = 1.00048/ref_ind->Eval(ph_energy);

	//cout<<ref_ind_g->Eval(ph_energy)<<"  "<<ph_energy<<endl;
	
	//if(ph_energy>=2. && ph_energy<2.5 && (eic_rich_vz->at(j)*0.1)<254.) ref_frac = 1.00048/1.01963;
	//if(ph_energy>=2.5 && ph_energy<3 && (eic_rich_vz->at(j)*0.1)<254.) ref_frac = 1.00048/1.01992;
	//if(ph_energy>=3 && ph_energy<3.5 && (eic_rich_vz->at(j)*0.1)<254.) ref_frac = 1.00048/1.02029;
	//if(ph_energy>=3.5 && ph_energy<4 && (eic_rich_vz->at(j)*0.1)<254.) ref_frac = 1.00048/1.02074;


	theta2 = TMath::ACos(esz);
	theta1 = TMath::ASin(TMath::Sin(theta2)*ref_frac);

	//cout<<theta2<<"  "<<theta1<<endl;

	esx = esx*TMath::Sin(theta1)/TMath::Sin(theta2);
	esy = esy*TMath::Sin(theta1)/TMath::Sin(theta2);
	esz = TMath::Cos(theta1);

	p = sqrt(in_px*in_px+in_py*in_py+in_pz*in_pz);
	//p = sqrt(gen_px->at(0)*gen_px->at(0)+gen_py->at(0)*gen_py->at(0)+gen_pz->at(0)*gen_pz->at(0));
	//p = sqrt(emi_x*emi_x+emi_y*emi_y+emi_z*emi_z);

	//cout<<esx<<"  "<<esy<<"  "<<esz<<"  "<<es<<"  "<<p<<endl;

	theta_c = TMath::ACos((esx*in_px+esy*in_py+esz*in_pz)/p);
	//theta_c = TMath::ACos((esx*gen_px->at(0)+esy*gen_py->at(0)+esz*gen_pz->at(0))/p);
	//theta_c = TMath::ACos((esx*emi_x+esy*emi_y+esz*emi_z)/es/p);

	//cout<<"primo: "<<in_px<<"  "<<in_py<<"  "<<in_pz<<endl;
	//cout<<"secondo: "<<gen_px->at(0)<<"  "<<gen_py->at(0)<<"  "<<gen_pz->at(0)<<endl;

	//cout<<eic_rich_trackE->at(j)<<"  "<<rran.Uniform(0,1)<<endl;

	for(Int_t k=0;k<phd_count-1;k++){

	  if((eic_rich_trackE->at(j)*1000000.)>=ev_a[k] && (eic_rich_trackE->at(j)*1000000.)<ev_a[k+1]){
	    qe_p = qe_a[k];
	    //cout<<eic_rich_trackE->at(j)*1000000.<<"  "<<qe_a[k]<<endl;
	  }
	}

	if((eic_rich_trackE->at(j)*1000000.)<ev_a[0] || (eic_rich_trackE->at(j)*1000000.)>ev_a[phd_count-1]) qe_p = 0.;

	//if(rran.Uniform(0,1)<qe_p && qe_eff != 0 && ph_energy<=4.2 && ph_energy>=2.04358 && (ph_det->GetBinContent(xbin,ybin) == 1)){
	if(rran.Uniform(0,1)<qe_p && qe_eff != 0 && ph_energy<=10. && ph_energy>=2.04358 && (ph_det->GetBinContent(xbin,ybin) == 1) && (shield_lable != 0)){

	  if(theta_c>(the_a[0]-0.02) && theta_c<(the_a[0]+0.02) && (gen_pid->at(0))==11){
	    tc_e_aerogel_h[index]->Fill(theta_c);
	    Ne_ph_a[0]++;
	  }
	  if(theta_c>(the_a[1]-0.02) && theta_c<(the_a[1]+0.02) && (gen_pid->at(0))==211){
	    tc_pi_aerogel_h[index]->Fill(theta_c);
	    /*if(pixel_eff != 0 && theta_c>(the_a[1]-0.003) && theta_c<(the_a[1]+0.003)) tc_pi_aerogel_h[index]->Fill(theta_c);
	    if(emi_eff != 0 && theta_c>(the_a[1]-0.01) && theta_c<(the_a[1]+0.01)) tc_pi_aerogel_h[index]->Fill(theta_c);
	    if(tr_eff != 0 && theta_c>(the_a[1]-0.01) && theta_c<(the_a[1]+0.01)) tc_pi_aerogel_h[index]->Fill(theta_c);
	    if(mag_lable != 0 && theta_c>(the_a[1]-0.02) && theta_c<(the_a[1]+0.02)) tc_pi_aerogel_h[index]->Fill(theta_c);
	    else tc_pi_aerogel_h[index]->Fill(theta_c);*/
	    Npi_ph_a[0]++;

	    tc_pi_aerogel_ene_h[index]->Fill(1239./ph_energy);
	    //if(i==1)cout<<theta_c<<endl;
	  }
	  else if((theta_c<(the_a[1]-0.02) || theta_c>(the_a[1]+0.02)) && (gen_pid->at(0))==211){
	    tc_pi_aerogel_ene_bk_h[index]->Fill(1239./ph_energy);
	    //tc_pi_aerogel_h[index]->Fill(theta_c);
	    //if(i==1)cout<<eic_rich_in_z->at(j)<<"  "<<eic_rich_out_z->at(j)<<endl;
	  }
	  if(theta_c>(the_a[2]-0.02) && theta_c<(the_a[2]+0.02) && (gen_pid->at(0))==321){
	    tc_k_aerogel_h[index]->Fill(theta_c);
	    Nk_ph_a[0]++;
	  }
	  if(theta_c>((the_a[3]-0.02)) && theta_c<((the_a[3]+0.02))  && (gen_pid->at(0))==2212){
	    tc_p_aerogel_h[index]->Fill(theta_c);
	    Np_ph_a[0]++;
	  }

	}
	else if(rran.Uniform(0,1)<qe_p && qe_eff != 0 && ph_energy<=4.2 && ph_energy>=2.04358 && (ph_det->GetBinContent(xbin,ybin) == 1) && (shield_lable == 0)){

	  if(theta_c>(the_a[0]-0.02) && theta_c<(the_a[0]+0.02) && (gen_pid->at(0))==11){
	    tc_e_aerogel_h[index]->Fill(theta_c);
	    Ne_ph_a[0]++;
	  }
	  if(theta_c>(the_a[1]-0.02) && theta_c<(the_a[1]+0.02) && (gen_pid->at(0))==211){
	    tc_pi_aerogel_h[index]->Fill(theta_c);
	    Npi_ph_a[0]++;

	    tc_pi_aerogel_ene_h[index]->Fill(1239./ph_energy);
	    //if(i==1)cout<<theta_c<<endl;
	  }
	  else if((theta_c<(the_a[1]-0.02) || theta_c>(the_a[1]+0.02)) && (gen_pid->at(0))==211){
	    tc_pi_aerogel_ene_bk_h[index]->Fill(1239./ph_energy);
	    //tc_pi_aerogel_h[index]->Fill(theta_c);
	    //if(i==1)cout<<eic_rich_in_z->at(j)<<"  "<<eic_rich_out_z->at(j)<<endl;
	  }
	  if(theta_c>(the_a[2]-0.02) && theta_c<(the_a[2]+0.02) && (gen_pid->at(0))==321){
	    tc_k_aerogel_h[index]->Fill(theta_c);
	    Nk_ph_a[0]++;
	  }
	  if(theta_c>((the_a[3]-0.02)) && theta_c<((the_a[3]+0.02))  && (gen_pid->at(0))==2212){
	    tc_p_aerogel_h[index]->Fill(theta_c);
	    Np_ph_a[0]++;
	  }

	}
	else if(qe_eff == 0 && ph_energy<=10. && ph_energy>=2.04358){

	  if(theta_c>(the_a[0]-0.02) && theta_c<(the_a[0]+0.02) && (gen_pid->at(0))==11){
	    tc_e_aerogel_h[index]->Fill(theta_c);
	    Ne_ph_a[0]++;
	  }
	  if(theta_c>(the_a[1]-0.02) && theta_c<(the_a[1]+0.02) && (gen_pid->at(0))==211){
	    tc_pi_aerogel_h[index]->Fill(theta_c);
	    Npi_ph_a[0]++;

	    tc_pi_aerogel_ene_h[index]->Fill(1239./ph_energy);
	  }
	  else if((theta_c<(the_a[1]-0.02) || theta_c>(the_a[1]+0.02)) && (gen_pid->at(0))==211){
	    tc_pi_aerogel_ene_bk_h[index]->Fill(1239./ph_energy);
	    //if(i==1) cout<<eic_rich_id->at(j)<<endl;
	  }
	  if(theta_c>(the_a[2]-0.02) && theta_c<(the_a[2]+0.02) && (gen_pid->at(0))==321){
	    tc_k_aerogel_h[index]->Fill(theta_c);
	    Nk_ph_a[0]++;
	  }
	  if(theta_c>(the_a[3]-0.02) && theta_c<(the_a[3]+0.02) && (gen_pid->at(0))==2212){
	    tc_p_aerogel_h[index]->Fill(theta_c);
	    Np_ph_a[0]++;
	  }

	}
	//cout<<"theta_c: "<<theta_c<<endl;

      }
      //if((eic_rich_id->at(j))==41 && (eic_rich_vz->at(j))>2540){
      if((eic_rich_id->at(j))>40  && (eic_rich_id->at(j))<50  && (eic_rich_mpid->at(j)==211 || eic_rich_mpid->at(j)==321 || eic_rich_mpid->at(j)==2212 || eic_rich_mpid->at(j)==11)){
	Int_t phi_count = eic_rich_id->at(j);
	phi_count = phi_count - 41;
        mir_x = mir_xx*TMath::Cos(TMath::Pi()*-60.*phi_count/180.);
	mir_y = mir_xx*TMath::Sin(TMath::Pi()*-60.*phi_count/180.);

	hit_h->Fill(eic_rich_hitn->at(j));
	//cout<<"vx,vy,vz: "<<eic_rich_vx->at(j)*0.1<<"  "<<eic_rich_vy->at(j)*0.1<<"  "<<eic_rich_vz->at(j)*0.1<<endl;
	//cout<<eic_rich_hitn->at(j)<<endl;
	//}
	//if((eic_rich_hitn->at(j))==255){
	x_h->Fill(eic_rich_avg_x->at(j));
	y_h->Fill(eic_rich_avg_y->at(j));
	z_h->Fill(eic_rich_avg_z->at(j));

	//ph_det1->Fill(eic_rich_avg_lx->at(j)*0.1,eic_rich_avg_ly->at(j)*0.1);
	
	xbin = ph_det->GetXaxis()->FindBin(eic_rich_avg_lx->at(j)*0.1);
	ybin = ph_det->GetYaxis()->FindBin(eic_rich_avg_ly->at(j)*0.1);

	det_x = eic_rich_in_x->at(j)*0.1;
	det_y = eic_rich_in_y->at(j)*0.1;
	det_z = eic_rich_in_z->at(j)*0.1;

	det_xl = ph_det->GetXaxis()->GetBinCenter(xbin);
	det_yl = ph_det->GetYaxis()->GetBinCenter(ybin);
	det_zl = eic_rich_avg_lz->at(j)*0.1;
	//**********Plane*********************//
	/*det_xr = det_xl*TMath::Cos(-TMath::Pi()*0./180.) + (det_zl)*TMath::Sin(-TMath::Pi()*0./180.);
	det_zr = -det_xl*TMath::Sin(-TMath::Pi()*0./180.) + (det_zl)*TMath::Cos(-TMath::Pi()*0./180.);

	det_xr1 = det_xr*TMath::Cos(TMath::Pi()*60.*phi_count/180.) + det_yl*TMath::Sin(TMath::Pi()*60.*phi_count/180.);
	det_yr1 = -det_xr*TMath::Sin(TMath::Pi()*60.*phi_count/180.) + det_yl*TMath::Cos(TMath::Pi()*60.*phi_count/180.);*/

	//**********Mirror*********************//
	det_xr = det_xl*TMath::Cos(-TMath::ASin(mir_xx/mir_R)) + (det_zl)*TMath::Sin(-TMath::ASin(mir_xx/mir_R));
	det_zr = -det_xl*TMath::Sin(-TMath::ASin(mir_xx/mir_R)) + (det_zl)*TMath::Cos(-TMath::ASin(mir_xx/mir_R))+mir_z;

	det_xr1 = det_xr*TMath::Cos(TMath::Pi()*60.*phi_count/180.) + det_yl*TMath::Sin(TMath::Pi()*60.*phi_count/180.)+mir_x;
	det_yr1 = -det_xr*TMath::Sin(TMath::Pi()*60.*phi_count/180.) + det_yl*TMath::Cos(TMath::Pi()*60.*phi_count/180.)+mir_y;

	//cout<<det_x<<"  "<<det_y<<"  "<<det_z<<endl;
	//cout<<det_xl<<"  "<<det_yl<<"  "<<det_zl<<endl;
	//cout<<det_xr1<<"  "<<det_yr1<<endl;

	/*emig_x = eic_rich_vx->at(j)*0.1;  //emission corrispondence with the previous above
	emig_y = eic_rich_vy->at(j)*0.1;
	emig_z = eic_rich_vz->at(j)*0.1;

	if(emig_z>254){
	emig1_x_h->Fill(emig_x);
	emig1_y_h->Fill(emig_y);
	emig1_z_h->Fill(emig_z);
	}*/

	//cout<<eic_rich_avg_z->at(j)<<endl;
	if(pixel_eff == 0) f->ind_ray(emig_x,emig_y,emig_z,det_x,det_y,det_z,mir_x,mir_y,mir_z);
	if(pixel_eff != 0) f->ind_ray(emig_x,emig_y,emig_z,det_xr1,det_yr1,det_z,mir_x,mir_y,mir_z);
	//cout<<"emission gas: "<<emig_x<<"  "<<emig_y<<"  "<<emig_z<<endl;
	//cout<<"detector point: "<<det_x<<"  "<<det_y<<"  "<<det_z<<endl;
	//cout<<emi_x<<"  "<<emi_y<<"  "<<emi_z<<endl;
	//cout<<det_x<<"  "<<det_y<<"  "<<det_z<<endl;
	//cout<<f->sx<<"  "<<f->sy<<"  "<<f->sz<<endl;

	esx = f->sx - emig_x;
	esy = f->sy - emig_y;
	esz = f->sz - emig_z;

	es = sqrt(esx*esx+esy*esy+esz*esz);
	dis_gas = es+sqrt((f->sx-det_x)*(f->sx-det_x)+(f->sy-det_y)*(f->sy-det_y)+(f->sz-det_z)*(f->sz-det_z));

	if((gen_pid->at(0))==211 && eic_rich_vz->at(j)*0.1>254.){
	  dist_pi_gas_h[index]->Fill(dis_gas);
	}

	p = sqrt(in_px*in_px+in_py*in_py+in_pz*in_pz);

	theta_c = TMath::ACos((esx*in_px+esy*in_py+esz*in_pz)/es/p);



	//p = sqrt(gen_px->at(0)*gen_px->at(0)+gen_py->at(0)*gen_py->at(0)+gen_pz->at(0)*gen_pz->at(0));

	//cout<<in_px<<"  "<<in_py<<"  "<<in_pz<<endl;

	//theta_c = TMath::ACos((esx*gen_px->at(0)+esy*gen_py->at(0)+esz*gen_pz->at(0))/es/p);

	for(Int_t k=0;k<phd_count-1;k++){

	  if((eic_rich_trackE->at(j)*1000000.)>=ev_a[k] && (eic_rich_trackE->at(j)*1000000.)<ev_a[k+1]){
	    qe_p = qe_a[k];
	    //cout<<eic_rich_trackE->at(j)*1000000.<<"  "<<qe_a[k]<<endl;
	  }
	}
	
	if((eic_rich_trackE->at(j)*1000000.)<ev_a[0] || (eic_rich_trackE->at(j)*1000000.)>ev_a[phd_count-1]) qe_p = 0.;
 

	if(rran.Uniform(0,1)<qe_p && qe_eff != 0 && ph_energy>=2.04358 && ph_energy<=10. && (ph_det->GetBinContent(xbin,ybin) == 1)){

	  if(theta_c>(the_g[0]-0.01) && theta_c<(the_g[0]+0.01) && (gen_pid->at(0))==11){
	    tc_e_gas_h[index]->Fill(theta_c);
	    //cout<< eic_rich_mpid->at(j)<<"  "<<eic_rich_mtid->at(j)<<"  "<<eic_rich_vz->at(j)*0.1<<endl;
	    Ne_ph_g[0]++;
	  }
	  if(theta_c>(the_g[1]-0.01) && theta_c<(the_g[1]+0.01) && (gen_pid->at(0))==211){
	    tc_pi_gas_h[index]->Fill(theta_c);
	    /*if(pixel_eff != 0 && theta_c>(the_g[1]-0.003) && theta_c<(the_g[1]+0.003)) tc_pi_gas_h[index]->Fill(theta_c); 
	    if(emi_eff != 0 && theta_c>(the_g[1]-0.01) && theta_c<(the_g[1]+0.01)) tc_pi_gas_h[index]->Fill(theta_c); 
	    if(tr_eff != 0 && theta_c>(the_g[1]-0.01) && theta_c<(the_g[1]+0.01)) tc_pi_gas_h[index]->Fill(theta_c); 
	    if(mag_lable != 0 && theta_c>(the_g[1]-0.02) && theta_c<(the_g[1]+0.02)) tc_pi_gas_h[index]->Fill(theta_c); 
	    else tc_pi_gas_h[index]->Fill(theta_c);*/
	    Npi_ph_g[0]++;
	  }
	  if(theta_c>(the_g[2]-0.01) && theta_c<(the_g[2]+0.01) && (gen_pid->at(0))==321){
	    tc_k_gas_h[index]->Fill(theta_c);
	    Nk_ph_g[0]++;
	  }
	  if(theta_c>(the_g[3]-0.01) && theta_c<(the_g[3]+0.01) && (gen_pid->at(0))==2212){
	    tc_p_gas_h[index]->Fill(theta_c);
	    Np_ph_g[0]++;
	  }

	}
	else if(qe_eff == 0){

	  if(theta_c>(the_g[0]-0.009) && theta_c<(the_g[0]+0.009) && (gen_pid->at(0))==11){
	    tc_e_gas_h[index]->Fill(theta_c);	    //cout<< eic_rich_mpid->at(j)<<"  "<<eic_rich_mtid->at(j)<<"  "<<eic_rich_vz->at(j)*0.1<<endl;
	    Ne_ph_g[0]++;
	  }
	  if(theta_c>(the_g[1]-0.009) && theta_c<(the_g[1]+0.009) && (gen_pid->at(0))==211){
	    tc_pi_gas_h[index]->Fill(theta_c);
	    Npi_ph_g[0]++;
	  }
	  if(theta_c>(the_g[2]-0.009) && theta_c<(the_g[2]+0.009) && (gen_pid->at(0))==321){
	    tc_k_gas_h[index]->Fill(theta_c);
	    Nk_ph_g[0]++;
	  }
	  if(theta_c>(the_g[3]-0.009) && theta_c<(the_g[3]+0.009) && (gen_pid->at(0))==2212){
	    tc_p_gas_h[index]->Fill(theta_c);
	    Np_ph_g[0]++;
	  }

	}
	//cout<<"theta_c: "<<theta_c<<endl;

      }

    }
    }

    //if(Ne_ph_g[0] == 0) cout<<i<<"----------------------"<<endl;
    //ph_det->Draw("col z");
    delete ph_det;

    if(Ne_ph_a[0]>0)nph_e_aerogel_h[index]->Fill(Ne_ph_a[0]);
    if(Npi_ph_a[0]>0)nph_pi_aerogel_h[index]->Fill(Npi_ph_a[0]);
    if(Nk_ph_a[0]>0)nph_k_aerogel_h[index]->Fill(Nk_ph_a[0]);
    if(Np_ph_a[0]>0)nph_p_aerogel_h[index]->Fill(Np_ph_a[0]);

    if(Ne_ph_g[0]>0)nph_e_gas_h[index]->Fill(Ne_ph_g[0]*0.7);
    if(Npi_ph_g[0]>0)nph_pi_gas_h[index]->Fill(Npi_ph_g[0]*0.7);
    if(Nk_ph_g[0]>0)nph_k_gas_h[index]->Fill(Nk_ph_g[0]*0.7);
    if(Np_ph_g[0]>0)nph_p_gas_h[index]->Fill(Np_ph_g[0]*0.7);

    for(Int_t n=0;n<1;n++){
      Ne_ph_a[n] = 0;
      Npi_ph_a[n] = 0;
      Nk_ph_a[n] = 0;
      Np_ph_a[n] = 0;
      
      Ne_ph_g[n] = 0;
      Npi_ph_g[n] = 0;
      Nk_ph_g[n] = 0;
      Np_ph_g[n] = 0;
    }

  }
  
  eic_rich->ResetBranchAddresses();
  generated->ResetBranchAddresses();

  delete eic_rich_hitn;  
  delete eic_rich_id;
  delete eic_rich_pid;
  delete eic_rich_mpid;
  delete eic_rich_tid;
  delete eic_rich_mtid;
  delete eic_rich_otid;
  delete eic_rich_trackE;
  delete eic_rich_totEdep;
  delete eic_rich_avg_x;
  delete eic_rich_avg_y;
  delete eic_rich_avg_z;
  delete eic_rich_avg_lx;
  delete eic_rich_avg_ly;
  delete eic_rich_avg_lz;
  delete eic_rich_px;
  delete eic_rich_py;
  delete eic_rich_pz;
  delete eic_rich_vx;
  delete eic_rich_vy;
  delete eic_rich_vz;
  delete eic_rich_mvx;
  delete eic_rich_mvy;
  delete eic_rich_mvz;
  delete eic_rich_avg_t;
  delete eic_rich_in_px;
  delete eic_rich_in_py;
  delete eic_rich_in_pz;
  delete eic_rich_in_x;
  delete eic_rich_in_y;
  delete eic_rich_in_z;
  delete eic_rich_in_t;
  delete eic_rich_out_px;
  delete eic_rich_out_py;
  delete eic_rich_out_pz;
  delete eic_rich_out_x;
  delete eic_rich_out_y;
  delete eic_rich_out_z;
  delete eic_rich_out_t;
  delete eic_rich_nsteps;

  delete gen_pid;
  delete gen_px;
  delete gen_py;
  delete gen_pz;
  delete gen_vx;
  delete gen_vy;
  delete gen_vz;

  delete generated;
  delete eic_rich;

  
  /*cout<<" pi_mean "<<tc_pi_gas_h[0]->GetMean()<<" k_mean "<<tc_k_gas_h[0]->GetMean()<<endl;
  cout<<" pi_sigma "<<tc_pi_gas_h[0]->GetRMS()<<" k_sigma "<<tc_k_gas_h[0]->GetRMS()<<endl;

  Double_t mean_diff_g = TMath::Abs(tc_pi_gas_h[0]->GetMean() -  tc_k_gas_h[0]->GetMean());
  Double_t mean_sigma_g = (tc_pi_gas_h[0]->GetRMS() +  tc_k_gas_h[0]->GetRMS())/2.;

  Double_t n_sigmas = mean_diff_g/mean_sigma_g;

  cout<<"n sigmas is: "<<n_sigmas<<endl;*/
  /*
  TCanvas *c1 = new TCanvas("c1","c1",1000,500);
  c1->Divide(2,1);
  c1->cd(1);
  id_h->Draw();
  c1->cd(2);
  hit_h->Draw();

  TCanvas *c2 = new TCanvas("c2","c2",1000,500);
  c2->Divide(3,1);
  c2->cd(1);
  x_h->Draw();
  c2->cd(2);
  y_h->Draw();
  c2->cd(3);
  z_h->Draw();
  */

  

  Double_t n_sigmas_pi_k_gas[100],n_sigmas_k_p_gas[100],n_sigmas_e_pi_gas[100],n_sigmas_pi_k_aerogel[100],n_sigmas_k_p_aerogel[100],n_sigmas_e_pi_aerogel[100];
  Double_t momentum1[100],momentum2[100],momentum3[100],momentum4[100],momentum5[100], momentum6[100];
  Double_t momentum[100];
  Int_t points=0;
  Int_t i1=0;
  Int_t i2=0;
  Int_t i3=0;
  Int_t i4=0;
  Int_t i5=0;
  Int_t i6=0;

  Double_t n_th_pi_k_gas[100],n_th_k_p_gas[100],n_th_e_pi_gas[100],n_th_pi_k_aerogel[100],n_th_k_p_aerogel[100],n_th_e_pi_aerogel[100];
  Double_t momentum1_th[100],momentum2_th[100],momentum3_th[100],momentum4_th[100],momentum5_th[100], momentum6_th[100];
  Int_t i1_th=0;
  Int_t i2_th=0;
  Int_t i3_th=0;
  Int_t i4_th=0;
  Int_t i5_th=0;
  Int_t i6_th=0;

  Double_t emi_err_a[100], emi_err_g[100], dist_a[100], dist_g[100], pol_angle[100];

  Double_t bgs_cont[100]; 

  //TCanvas *cch = new TCanvas("cch","cch",1000,500);
  for(Int_t i=0;i<=index;i++){


    //cout<<" pi_mean "<<tc_pi_gas_h[i]->GetMean()<<" k_mean "<<tc_k_gas_h[i]->GetMean()<<endl;
    //cout<<" pi_sigma "<<tc_pi_gas_h[i]->GetRMS()<<" k_sigma "<<tc_k_gas_h[i]->GetRMS()<<endl;


    Double_t mean_diff_pi_k_a = TMath::Abs(tc_pi_aerogel_h[i]->GetMean() -  tc_k_aerogel_h[i]->GetMean());
    Double_t mean_sigma_pi_k_a = (tc_pi_aerogel_h[i]->GetRMS() +  tc_k_aerogel_h[i]->GetRMS())/2.;
    Double_t mean_diff_k_p_a = TMath::Abs(tc_k_aerogel_h[i]->GetMean() -  tc_p_aerogel_h[i]->GetMean());
    Double_t mean_sigma_k_p_a = (tc_k_aerogel_h[i]->GetRMS() +  tc_p_aerogel_h[i]->GetRMS())/2.;
    Double_t mean_diff_e_pi_a = TMath::Abs(tc_e_aerogel_h[i]->GetMean() -  tc_pi_aerogel_h[i]->GetMean());
    Double_t mean_sigma_e_pi_a = (tc_e_aerogel_h[i]->GetRMS() +  tc_pi_aerogel_h[i]->GetRMS())/2.;

    Double_t mean_diff_pi_k_g = TMath::Abs(tc_pi_gas_h[i]->GetMean() -  tc_k_gas_h[i]->GetMean());
    Double_t mean_sigma_pi_k_g = (tc_pi_gas_h[i]->GetRMS() +  tc_k_gas_h[i]->GetRMS())/2.;
    Double_t mean_diff_k_p_g = TMath::Abs(tc_k_gas_h[i]->GetMean() -  tc_p_gas_h[i]->GetMean());
    Double_t mean_sigma_k_p_g = (tc_k_gas_h[i]->GetRMS() +  tc_p_gas_h[i]->GetRMS())/2.;
    Double_t mean_diff_e_pi_g = TMath::Abs(tc_e_gas_h[i]->GetMean() -  tc_pi_gas_h[i]->GetMean());
    Double_t mean_sigma_e_pi_g = (tc_e_gas_h[i]->GetRMS() +  tc_pi_gas_h[i]->GetRMS())/2.;

    Double_t Nph_a_pi_k = (tc_pi_aerogel_h[i]->GetEntries() + tc_k_aerogel_h[i]->GetEntries())/2./30.;  //remember the 30--->number of generated events!
    Double_t Nph_a_k_p = (tc_k_aerogel_h[i]->GetEntries() + tc_p_aerogel_h[i]->GetEntries())/2./30.;
    Double_t Nph_a_e_pi = (tc_e_aerogel_h[i]->GetEntries() + tc_pi_aerogel_h[i]->GetEntries())/2./30.; 

    Double_t Nph_g_pi_k = (tc_pi_gas_h[i]->GetEntries() + tc_k_gas_h[i]->GetEntries())/2./30.;  //remember the 30--->number of generated events!
    Double_t Nph_g_k_p = (tc_k_gas_h[i]->GetEntries() + tc_p_gas_h[i]->GetEntries())/2./30.;
    Double_t Nph_g_e_pi = (tc_e_gas_h[i]->GetEntries() + tc_pi_gas_h[i]->GetEntries())/2./30.; 

    if(i<4) momentum[i] = 1.+2.*i;
    if(i>=4) momentum[i] = 10.+3.*(i-4);
    pol_angle[i] = 5.+5.*i;

    angl[i] = pol_angle[i];

   if(momentum[i]>=0.5 && momentum[i]<=15.){
      momentum5[i5]=momentum[i];
      n_sigmas_e_pi_aerogel[i5] = (sqrt((nph_e_aerogel_h[i]->GetMean()+nph_pi_aerogel_h[i]->GetMean())/2.))*mean_diff_e_pi_a/mean_sigma_e_pi_a;
      cout<<"1 "<<" 1 "<<n_sigmas_e_pi_aerogel[i5]<<" "<<momentum[i]<<endl;
      //cout<<"e/pi(aerogel):          "<<n_sigmas_e_pi_aerogel[i5]<< "          "<<momentum[i]<<endl;
      i5++;
    }
    if(momentum[i]>=2.5 && momentum[i]<=21.){
      momentum1[i1]=momentum[i];
      n_sigmas_pi_k_aerogel[i1] = (sqrt((nph_pi_aerogel_h[i]->GetMean()+nph_k_aerogel_h[i]->GetMean())/2.))*mean_diff_pi_k_a/mean_sigma_pi_k_a;
      cout<<"1 "<<" 2 "<<n_sigmas_pi_k_aerogel[i1]<<" "<<momentum[i]<<endl;
      //cout<<"pi/k(aerogel):           "<<n_sigmas_pi_k_aerogel[i1]<< "          "<<momentum[i]<<endl;
      i1++;
    }
    if(momentum[i]>=5. && momentum[i]<=27.){
      momentum2[i2]=momentum[i];
      n_sigmas_k_p_aerogel[i2] = (sqrt((nph_k_aerogel_h[i]->GetMean()+nph_p_aerogel_h[i]->GetMean())/2.))*mean_diff_k_p_a/mean_sigma_k_p_a;
      cout<<"1 "<<" 3 "<<n_sigmas_k_p_aerogel[i2]<<" "<<momentum[i]<<endl;
      //cout<<"k/p(aerogel):           "<<n_sigmas_k_p_aerogel[i2]<< "          "<<momentum[i]<<endl;
      i2++;
    }
    if(momentum[i]>=3.5 && momentum[i]<=30.){
      momentum6[i6]=momentum[i];
      n_sigmas_e_pi_gas[i6] = (sqrt((nph_e_gas_h[i]->GetMean()+nph_pi_gas_h[i]->GetMean())/2.))*mean_diff_e_pi_g/mean_sigma_e_pi_g;
      cout<<"2 "<<" 1 "<<n_sigmas_e_pi_gas[i6]<<" "<<momentum[i]<<endl;
      //cout<<"e/pi(gas):          "<<n_sigmas_e_pi_gas[i6]<< "          "<<momentum[i]<<endl;
      i6++;
    }
    if(momentum[i]>=12.5){
      momentum3[i3]=momentum[i];
      n_sigmas_pi_k_gas[i3] = (sqrt((nph_pi_gas_h[i]->GetMean()+nph_k_gas_h[i]->GetMean())/2.))*mean_diff_pi_k_g/mean_sigma_pi_k_g;
      cout<<"2 "<<" 2 "<<n_sigmas_pi_k_gas[i3]<<" "<<momentum[i]<<endl;
      //cout<<"pi/k(gas):           "<<n_sigmas_pi_k_gas[i3]<< "          "<<momentum[i]<<endl;
      i3++;
    }
    if(momentum[i]>=23.5){
      momentum4[i4]=momentum[i];
      n_sigmas_k_p_gas[i4] = (sqrt((nph_k_gas_h[i]->GetMean()+nph_p_gas_h[i]->GetMean())/2.))*mean_diff_k_p_g/mean_sigma_k_p_g;
      cout<<"2 "<<" 3 "<<n_sigmas_k_p_gas[i4]<<" "<<momentum[i]<<endl;
      //cout<<"k/p(gas):           "<<n_sigmas_k_p_gas[i4]<< "          "<<momentum[i]<<endl;
      i4++;
    }



   if(momentum[i]>=1. && momentum[i]<=20.){
      momentum5_th[i5_th]=momentum[i];
      n_th_e_pi_aerogel[i5_th] = (nph_e_aerogel_h[i]->GetMean()-nph_pi_aerogel_h[i]->GetMean())/((nph_e_aerogel_h[i]->GetRMS()+nph_pi_aerogel_h[i]->GetRMS())/2.);
      i5_th++;
    }
    if(momentum[i]>=1. && momentum[i]<=20.){
      momentum1_th[i1_th]=momentum[i];
      n_th_pi_k_aerogel[i1_th] = (nph_pi_aerogel_h[i]->GetMean()-nph_k_aerogel_h[i]->GetMean())/((nph_pi_aerogel_h[i]->GetRMS()+nph_k_aerogel_h[i]->GetRMS())/2.);
      i1_th++;
    }
    if(momentum[i]>=2.5 && momentum[i]<=20.){
      momentum2_th[i2_th]=momentum[i];
      n_th_k_p_aerogel[i2_th] = (nph_k_aerogel_h[i]->GetMean()-nph_p_aerogel_h[i]->GetMean())/((nph_k_aerogel_h[i]->GetRMS()+nph_p_aerogel_h[i]->GetRMS())/2.);
      i2_th++;
    }
    if(momentum[i]>=1.){
      momentum6_th[i6_th]=momentum[i];
      n_th_e_pi_gas[i6_th] = (nph_e_gas_h[i]->GetMean()-nph_pi_gas_h[i]->GetMean())/((nph_e_gas_h[i]->GetRMS()+nph_pi_gas_h[i]->GetRMS())/2.);
      i6_th++;
    }
    if(momentum[i]>=3.5){
      momentum3_th[i3_th]=momentum[i];
      n_th_pi_k_gas[i3_th] = (nph_pi_gas_h[i]->GetMean()-nph_k_gas_h[i]->GetMean())/((nph_pi_gas_h[i]->GetRMS()+nph_k_gas_h[i]->GetRMS())/2.);
      i3_th++;
    }
    if(momentum[i]>=12.5){
      momentum4_th[i4_th]=momentum[i];
      n_th_k_p_gas[i4_th] = (nph_k_gas_h[i]->GetMean()-nph_p_gas_h[i]->GetMean())/((nph_k_gas_h[i]->GetRMS()+nph_p_gas_h[i]->GetRMS())/2.);
      i4_th++;
    }


    emi_err_a[i] = tc_pi_aerogel_h[i]->GetRMS();
    emi_err_g[i] = tc_pi_gas_h[i]->GetRMS();
    dist_a[i] = dist_pi_aerogel_h[i]->GetRMS();
    dist_g[i] = dist_pi_gas_h[i]->GetRMS();

    //...errors for plots.... 
    
    if(mag_lable!=0 && emi_eff == 0 && pixel_eff == 0) mg_er_a[i] = tc_pi_aerogel_h[i]->GetRMS();
    if(mag_lable!=0 && emi_eff == 0 && pixel_eff == 0) mg_er_g[i] = tc_pi_gas_h[i]->GetRMS();
    if(emi_eff == 0 && pixel_eff == 0 && mag_lable==0 && tr_eff==0) ch_er_a[i] = tc_pi_aerogel_h[i]->GetRMS();
    if(emi_eff == 0 && pixel_eff == 0 && mag_lable==0 && tr_eff==0) ch_er_g[i] = tc_pi_gas_h[i]->GetRMS();
    if(emi_eff != 0) emi_er_a[i] = tc_pi_aerogel_h[i]->GetRMS();
    if(emi_eff != 0) emi_er_g[i] = tc_pi_gas_h[i]->GetRMS();
    if(pixel_eff != 0) pix_er_a[i] = tc_pi_aerogel_h[i]->GetRMS();
    if(pixel_eff != 0) pix_er_g[i] = tc_pi_gas_h[i]->GetRMS();
    if(tr_eff != 0) tr_er_a[i] = tc_pi_aerogel_h[i]->GetRMS();
    if(tr_eff != 0) tr_er_g[i] = tc_pi_gas_h[i]->GetRMS();

    //emi_err_a[i] = tc_k_aerogel_h[i]->GetRMS();
    //emi_err_g[i] = tc_k_gas_h[i]->GetRMS();

    points = i+1;
    
    //cout<<"n sigmas is: "<<n_sigmas_pi_k_gas[i]<<endl;
    
    
    cch[i] = new TCanvas(Form("cch%d",i),"",1000,800); 
    cch[i]->Divide(2,4);
    cch[i]->cd(1);
    tc_pi_gas_h[i]->SetFillColor(2);
    tc_pi_gas_h[i]->Draw();
    //tc_aerogel_h->Draw("same");
    cch[i]->cd(2);
    tc_pi_aerogel_h[i]->Draw();
    cch[i]->cd(3);
    tc_k_gas_h[i]->SetFillColor(2);
    tc_k_gas_h[i]->Draw();
    cch[i]->cd(4);
    tc_k_aerogel_h[i]->Draw();
    cch[i]->cd(5);
    tc_p_gas_h[i]->SetFillColor(2);
    tc_p_gas_h[i]->Draw();
    cch[i]->cd(6);
    tc_p_aerogel_h[i]->Draw();
    cch[i]->cd(7);
    tc_e_gas_h[i]->SetFillColor(2);
    tc_e_gas_h[i]->Draw();
    cch[i]->cd(8);
    tc_e_aerogel_h[i]->Draw();
    
    /*
    cch_ph[i] = new TCanvas(Form("cch_ph%d",i),"",1000,800); 
    cch_ph[i]->Divide(2,4);
    cch_ph[i]->cd(1);
    nph_pi_gas_h[i]->SetFillColor(2);
    nph_pi_gas_h[i]->Draw();
    cch_ph[i]->cd(2);
    nph_pi_aerogel_h[i]->Draw();
    cch_ph[i]->cd(3);
    nph_k_gas_h[i]->SetFillColor(2);
    nph_k_gas_h[i]->Draw();
    cch_ph[i]->cd(4);
    nph_k_aerogel_h[i]->Draw();
    cch_ph[i]->cd(5);
    nph_p_gas_h[i]->SetFillColor(2);
    nph_p_gas_h[i]->Draw();
    cch_ph[i]->cd(6);
    nph_p_aerogel_h[i]->Draw();
    cch_ph[i]->cd(7);
    nph_e_gas_h[i]->SetFillColor(2);
    nph_e_gas_h[i]->Draw();
    cch_ph[i]->cd(8);
    nph_e_aerogel_h[i]->Draw();
    */
    if(shield_lable==0){
      Nph_a[i] = nph_pi_aerogel_h[i]->GetMean();
      err_Nph_a[i] = nph_pi_aerogel_h[i]->GetRMS();
      Nph_g[i] = nph_pi_gas_h[i]->GetMean();
      err_Nph_g[i] = nph_pi_gas_h[i]->GetRMS();
    }
    else{
      Nph_a_shield[i] = nph_pi_aerogel_h[i]->GetMean();
      err_Nph_a_shield[i] = nph_pi_aerogel_h[i]->GetRMS();
    } 

    //cout<<pol_angle[i]<<"  "<<nph_pi_aerogel_h[i]->GetMean()<<endl;

    /*cch_d[i] = new TCanvas(Form("cch_d%d",i),"",1000,800); 
    cch_d[i]->Divide(2,2);
    cch_d[i]->cd(1);
    dist_pi_gas_h[i]->SetFillColor(2);
    dist_pi_gas_h[i]->Draw();
    cch_d[i]->cd(2);
    dist_pi_aerogel_h[i]->Draw();*/    
    /*
    cch[i] = new TCanvas(Form("cch%d",i),"",500,500); 
    cch[i]->Divide(1,3);
    cch[i]->cd(1);
    //tc_pi_gas_h[i]->SetFillColor(2);
    //tc_pi_gas_h[i]->Draw();
    //tc_aerogel_h->Draw("same");
    //cch[i]->cd(2);
    tc_pi_aerogel_h[i]->Draw();
    //cch[i]->cd(3);
    //tc_k_gas_h[i]->SetFillColor(2);
    //tc_k_gas_h[i]->Draw();
    cch[i]->cd(2);
    tc_k_aerogel_h[i]->Draw();
    //cch[i]->cd(5);
    //tc_p_gas_h[i]->SetFillColor(2);
    //tc_p_gas_h[i]->Draw();
    cch[i]->cd(3);
    tc_p_aerogel_h[i]->Draw();
    */

    dist_ph[i] = new TCanvas(Form("dist_ph%d",i),"",500,500);
    dist_ph[i]->Divide(1,2);
    dist_ph[i]->cd(1);
    tc_pi_aerogel_ene_h[i]->Draw();
    tc_pi_aerogel_ene_h[i]->GetXaxis()->SetTitle("#lambda [nm]");
    dist_ph[i]->cd(2);
    tc_pi_aerogel_ene_bk_h[i]->Draw();
    tc_pi_aerogel_ene_bk_h[i]->SetFillColor(2);
    tc_pi_aerogel_ene_bk_h[i]->GetXaxis()->SetTitle("#lambda [nm]");

    bgs_cont[i]=tc_pi_aerogel_ene_bk_h[i]->GetEntries()*100./tc_pi_aerogel_ene_h[i]->GetEntries();

    }

  TGraph *a_sig_e_pi = new TGraph(i5,momentum5,n_sigmas_e_pi_aerogel);
  TGraph *a_sig_pi_k = new TGraph(i1,momentum1,n_sigmas_pi_k_aerogel);
  TGraph *a_sig_k_p = new TGraph(i2,momentum2,n_sigmas_k_p_aerogel);
  TGraph *g_sig_e_pi = new TGraph(i6,momentum6,n_sigmas_e_pi_gas);
  TGraph *g_sig_pi_k = new TGraph(i3,momentum3,n_sigmas_pi_k_gas);
  TGraph *g_sig_k_p = new TGraph(i4,momentum4,n_sigmas_k_p_gas);

  TGraph *a_th_e_pi = new TGraph(i5_th,momentum5_th,n_th_e_pi_aerogel);
  TGraph *a_th_pi_k = new TGraph(i1_th,momentum1_th,n_th_pi_k_aerogel);
  TGraph *a_th_k_p = new TGraph(i2_th,momentum2_th,n_th_k_p_aerogel);
  TGraph *g_th_e_pi = new TGraph(i6_th,momentum6_th,n_th_e_pi_gas);
  TGraph *g_th_pi_k = new TGraph(i3_th,momentum3_th,n_th_pi_k_gas);
  TGraph *g_th_k_p = new TGraph(i4_th,momentum4_th,n_th_k_p_gas);


  Double_t fx[11] , fy[11], f3sig[11];
  for(Int_t e=0;e<11;e++){
    fx[e] = e*7;
    fy[e] = 200;   
    f3sig[e]=3;
  }
  TGraph *fake = new TGraph(11,fx,fy);
  TGraph *fake1 = new TGraph(11,fx,f3sig);

  /* 
  TCanvas *c_sig_aerogel = new TCanvas("c_sig_aerogel","",1000,500);
  c_sig_aerogel->Divide(1,3);
  c_sig_aerogel->cd(1);
  a_sig_pi_k->Draw("AC*");
  a_sig_pi_k->SetMarkerStyle(20);
  a_sig_pi_k->SetMarkerColor(1);
  a_sig_pi_k->SetLineStyle(1);
  a_sig_pi_k->SetLineWidth(1);
  a_sig_pi_k->SetTitle("Aerogel");
  a_sig_pi_k->GetXaxis()->SetTitle("momentum [GeV/c]");
  a_sig_pi_k->GetYaxis()->SetTitle("N_{#sigma}^{#pi,K}");
  c_sig_aerogel->cd(2);
  a_sig_k_p->Draw("AC*");
  a_sig_k_p->SetMarkerStyle(20);
  a_sig_k_p->SetMarkerColor(1);
  a_sig_k_p->SetLineStyle(1);
  a_sig_k_p->SetLineWidth(1);
  a_sig_k_p->SetTitle("Aerogel");
  a_sig_k_p->GetXaxis()->SetTitle("momentum [GeV/c]");
  a_sig_k_p->GetYaxis()->SetTitle("N_{#sigma}^{K,p}");
  c_sig_aerogel->cd(3);
  a_sig_e_pi->Draw("AC*");
  a_sig_e_pi->SetMarkerStyle(20);
  a_sig_e_pi->SetMarkerColor(1);
  a_sig_e_pi->SetLineStyle(1);
  a_sig_e_pi->SetLineWidth(1);
  a_sig_e_pi->SetTitle("Aerogel");
  a_sig_e_pi->GetXaxis()->SetTitle("momentum [GeV/c]");
  a_sig_e_pi->GetYaxis()->SetTitle("N_{#sigma}^{e,#pi}");

  TCanvas *c_sig_gas = new TCanvas("c_sig_gas","",1000,500);
  c_sig_gas->Divide(1,3);
  c_sig_gas->cd(1);
  g_sig_pi_k->Draw("AC*");
  g_sig_pi_k->SetMarkerStyle(20);
  g_sig_pi_k->SetMarkerColor(1);
  g_sig_pi_k->SetLineStyle(1);
  g_sig_pi_k->SetLineWidth(1);
  g_sig_pi_k->SetTitle("CF_{4} gas");
  g_sig_pi_k->GetXaxis()->SetTitle("momentum [GeV/c]");
  g_sig_pi_k->GetYaxis()->SetTitle("N_{#sigma}^{#pi,K}");
  c_sig_gas->cd(2);
  g_sig_k_p->Draw("AC*");
  g_sig_k_p->SetMarkerStyle(20);
  g_sig_k_p->SetMarkerColor(1);
  g_sig_k_p->SetLineStyle(1);
  g_sig_k_p->SetLineWidth(1);
  g_sig_k_p->SetTitle("CF_{4} gas");
  g_sig_k_p->GetXaxis()->SetTitle("momentum [GeV/c]");
  g_sig_k_p->GetYaxis()->SetTitle("N_{#sigma}^{K,p}");
  c_sig_gas->cd(3);
  g_sig_e_pi->Draw("AC*");
  g_sig_e_pi->SetMarkerStyle(20);
  g_sig_e_pi->SetMarkerColor(1);
  g_sig_e_pi->SetLineStyle(1);
  g_sig_e_pi->SetLineWidth(1);
  g_sig_e_pi->SetTitle("CF_{4} gas");
  g_sig_e_pi->GetXaxis()->SetTitle("momentum [GeV/c]");
  g_sig_e_pi->GetYaxis()->SetTitle("N_{#sigma}^{e,#pi}");
  */

  
  TCanvas *c_sig_same = new TCanvas("c_sig_same","",1000,500);
  c_sig_same->Divide(1,1);
  c_sig_same->cd(1);
  fake->Draw("AC*");
  fake->SetTitle("");
  fake->GetXaxis()->SetTitle("momentum [GeV/c]");
  fake->GetYaxis()->SetTitle("N_{#sigma}^{ring}");
  fake->SetMarkerColor(1);
  fake->SetMarkerStyle(1);
  fake->SetLineColor(1);
  fake1->Draw("C*");
  fake1->SetMarkerColor(1);
  fake1->SetMarkerStyle(1);
  fake1->SetLineColor(1);
  fake1->SetLineStyle(1);
  fake1->SetLineWidth(2);
  g_sig_pi_k->Draw("C*");
  g_sig_pi_k->SetMarkerStyle(20);
  g_sig_pi_k->SetMarkerColor(4);
  g_sig_pi_k->SetLineColor(4);
  g_sig_pi_k->SetLineStyle(1);
  g_sig_pi_k->SetLineWidth(2);
  //g_sig_pi_k->SetTitle("");
  //g_sig_pi_k->GetXaxis()->SetTitle("momentum [GeV/c]");
  //g_sig_pi_k->GetYaxis()->SetTitle("N_{#sigma}");
  g_sig_k_p->Draw("C*");
  g_sig_k_p->SetMarkerStyle(20);
  g_sig_k_p->SetMarkerColor(4);
  g_sig_k_p->SetLineColor(4);
  g_sig_k_p->SetLineStyle(2);
  g_sig_k_p->SetLineWidth(2);
  //g_sig_k_p->SetTitle("CF_{4} gas");
  //g_sig_k_p->GetXaxis()->SetTitle("momentum [GeV/c]");
  //g_sig_k_p->GetYaxis()->SetTitle("N_{#sigma}^{K,p}");
  g_sig_e_pi->Draw("C*");
  g_sig_e_pi->SetMarkerStyle(20);
  g_sig_e_pi->SetMarkerColor(4);
  g_sig_e_pi->SetLineStyle(8);
  g_sig_e_pi->SetLineColor(4);
  g_sig_e_pi->SetLineWidth(2);

  a_sig_pi_k->Draw("C*");
  a_sig_pi_k->SetMarkerStyle(20);
  a_sig_pi_k->SetMarkerColor(2);
  a_sig_pi_k->SetLineColor(2);
  a_sig_pi_k->SetLineStyle(1);
  a_sig_pi_k->SetLineWidth(2);
  //a_sig_pi_k->SetTitle("Aerogel");
  //a_sig_pi_k->GetXaxis()->SetTitle("momentum [GeV/c]");
  //a_sig_pi_k->GetYaxis()->SetTitle("N_{#sigma}^{#pi,K}");
  a_sig_k_p->Draw("C*");
  a_sig_k_p->SetMarkerStyle(20);
  a_sig_k_p->SetMarkerColor(2);
  a_sig_k_p->SetLineColor(2);
  a_sig_k_p->SetLineStyle(2);
  a_sig_k_p->SetLineWidth(2);
  //a_sig_k_p->SetTitle("Aerogel");
  //a_sig_k_p->GetXaxis()->SetTitle("momentum [GeV/c]");
  //a_sig_k_p->GetYaxis()->SetTitle("N_{#sigma}^{K,p}");
  a_sig_e_pi->Draw("C*");
  a_sig_e_pi->SetMarkerStyle(20);
  a_sig_e_pi->SetMarkerColor(2);
  a_sig_e_pi->SetLineStyle(8);
  a_sig_e_pi->SetLineColor(2);
  a_sig_e_pi->SetLineWidth(2);


  TCanvas *c_th_same = new TCanvas("c_th_same","",1000,500);
  c_th_same->Divide(1,1);
  c_th_same->cd(1);
  fake->Draw("AC*");
  fake->SetTitle("");
  fake->GetXaxis()->SetTitle("momentum [GeV/c]");
  fake->GetYaxis()->SetTitle("N_{#sigma}");
  fake->SetMarkerColor(1);
  fake->SetMarkerStyle(1);
  fake->SetLineColor(1);
  fake1->Draw("C*");
  fake1->SetMarkerColor(1);
  fake1->SetMarkerStyle(1);
  fake1->SetLineColor(1);
  fake1->SetLineStyle(1);
  fake1->SetLineWidth(2);
  g_th_pi_k->Draw("C*");
  g_th_pi_k->SetMarkerStyle(20);
  g_th_pi_k->SetMarkerColor(4);
  g_th_pi_k->SetLineColor(4);
  g_th_pi_k->SetLineStyle(1);
  g_th_pi_k->SetLineWidth(2);
  g_th_k_p->Draw("C*");
  g_th_k_p->SetMarkerStyle(20);
  g_th_k_p->SetMarkerColor(4);
  g_th_k_p->SetLineColor(4);
  g_th_k_p->SetLineStyle(2);
  g_th_k_p->SetLineWidth(2);
  g_th_e_pi->Draw("C*");
  g_th_e_pi->SetMarkerStyle(20);
  g_th_e_pi->SetMarkerColor(4);
  g_th_e_pi->SetLineStyle(8);
  g_th_e_pi->SetLineColor(4);
  g_th_e_pi->SetLineWidth(2);
  
  a_th_pi_k->Draw("C*");
  a_th_pi_k->SetMarkerStyle(20);
  a_th_pi_k->SetMarkerColor(2);
  a_th_pi_k->SetLineColor(2);
  a_th_pi_k->SetLineStyle(1);
  a_th_pi_k->SetLineWidth(2);
  a_th_k_p->Draw("C*");
  a_th_k_p->SetMarkerStyle(20);
  a_th_k_p->SetMarkerColor(2);
  a_th_k_p->SetLineColor(2);
  a_th_k_p->SetLineStyle(2);
  a_th_k_p->SetLineWidth(2);
  a_th_e_pi->Draw("C*");
  a_th_e_pi->SetMarkerStyle(20);
  a_th_e_pi->SetMarkerColor(2);
  a_th_e_pi->SetLineStyle(8);
  a_th_e_pi->SetLineColor(2);
  a_th_e_pi->SetLineWidth(2);  

  /*
  TCanvas *ph_space = new TCanvas("ph_space","",600,600);
  ph_space->Divide(1,2);
  ph_space->cd(1);
  ph_det->Draw("COLZ");
  ph_det->GetXaxis()->SetTitle("cm");
  ph_det->GetYaxis()->SetTitle("cm");
  ph_space->cd(2);
  ph_det_sph->Draw("COLZ");
  ph_det_sph->GetXaxis()->SetTitle("rad");
  ph_det_sph->GetYaxis()->SetTitle("rad");
  */
  /*ph_space->cd(2);
  ph_det_x->Draw();
  ph_det_x->GetXaxis()->SetTitle("cm");
  ph_space->cd(3);
  ph_det_y->Draw();
  ph_det_y->GetXaxis()->SetTitle("cm");
  ph_space->cd(4);
  ph_det1->Draw("COLZ");
  ph_det1->GetXaxis()->SetTitle("cm");
  ph_det1->GetYaxis()->SetTitle("cm");*/
  /*
  TCanvas *c_emi = new TCanvas("c_emi","",1000,1500);
  c_emi->Divide(2,3);
  c_emi->cd(1);
  emi_x_h->Draw();
  c_emi->cd(3);
  emi_y_h->Draw();
  c_emi->cd(5);
  emi_z_h->Draw();
  c_emi->cd(2);
  emig_x_h->Draw();
  c_emi->cd(4);
  emig_y_h->Draw();
  c_emi->cd(6);
  emig_z_h->Draw();
  */
  TGraph *a_emi = new TGraph(points,pol_angle,emi_err_a);
  TGraph *g_emi = new TGraph(points,pol_angle,emi_err_g);
  TGraph *a_dist = new TGraph(points,pol_angle,dist_a);
  TGraph *g_dist = new TGraph(points,pol_angle,dist_g);

  TGraph *contamination = new TGraph(points,pol_angle,bgs_cont);

  TGraphErrors *a_nph = new TGraphErrors(points,pol_angle,Nph_a,0,err_Nph_a);

  //TGraph *g_emi = new TGraph(points,momentum,emi_err_g);  

  TCanvas *emi_e = new TCanvas("emi_e","",800,800);
  emi_e->Divide(1,2);
  emi_e->cd(1);
  a_emi->Draw("AC*");
  a_emi->SetMarkerStyle(20);
  a_emi->SetMarkerColor(1);
  a_emi->SetLineStyle(1);
  a_emi->SetLineWidth(1);
  a_emi->SetTitle("Aerogel");
  a_emi->GetXaxis()->SetTitle("polar angle [deg]");
  if(emi_eff!=0) a_emi->GetYaxis()->SetTitle("#sigma_{emission (1 p.e.)}^{Aerogel} [rad]");
  if(pixel_eff!=0) a_emi->GetYaxis()->SetTitle("#sigma_{pixel (1 p.e.)}^{Aerogel} [rad]"); 
  if(pixel_eff==0 && emi_eff==0) a_emi->GetYaxis()->SetTitle("#sigma_{chromatic (1 p.e.)}^{Aerogel} [rad]");
  emi_e->cd(2);
  g_emi->Draw("AC*");
  g_emi->SetMarkerStyle(20);
  g_emi->SetMarkerColor(1);
  g_emi->SetLineStyle(1);
  g_emi->SetLineWidth(1);
  g_emi->SetTitle("C_{2}F_{6} gas");
  g_emi->GetXaxis()->SetTitle("polar angle [deg]");
  if(emi_eff!=0) g_emi->GetYaxis()->SetTitle("#sigma_{emission (1 p.e.)}^{C_{2}F_{6}} [rad]");
  if(pixel_eff!=0) g_emi->GetYaxis()->SetTitle("#sigma_{pixel (1 p.e.)}^{CF_{4}} [rad]");
  if(pixel_eff==0 && emi_eff==0) g_emi->GetYaxis()->SetTitle("#sigma_{chromatic (1 p.e.)}^{CF_{4}} [rad]");

  TCanvas *mmm_c = new TCanvas("mmm_c","",800,800);
  mmm_c->Divide(1,1);
  //mmm->Draw();
  contamination->Draw("AC*");
  contamination->SetMarkerStyle(20);
  contamination->SetMarkerColor(1);
  contamination->SetLineStyle(1);
  contamination->SetLineWidth(1);
  contamination->SetTitle("Contamination %");
  contamination->GetXaxis()->SetTitle("polar angle [deg]");

  /*TCanvas *mean_nph = new TCanvas("mean_nph","",800,800);
  mean_nph->Divide(1,1);
  a_nph->Draw("AB1");
  a_nph->SetMarkerStyle(20);
  a_nph->SetMarkerColor(1);
  a_nph->SetLineStyle(1);
  a_nph->SetLineWidth(1);
  a_nph->SetTitle("");
  a_nph->GetXaxis()->SetTitle("polar angle [deg]");
  a_nph->GetYaxis()->SetTitle("Mean number of p.e.");*/

  /*TCanvas *dist_e = new TCanvas("dist_e","",1000,500);
  dist_e->Divide(1,2);
  dist_e->cd(1);
  a_dist->Draw("AC*");
  a_dist->SetMarkerStyle(20);
  a_dist->SetMarkerColor(1);
  a_dist->SetLineStyle(1);
  a_dist->SetLineWidth(1);
  a_dist->SetTitle("Aerogel");
  a_dist->GetXaxis()->SetTitle("polar angle [deg]");
  a_dist->GetYaxis()->SetTitle("#sigma_{dist}^{a}");
  dist_e->cd(2);
  g_dist->Draw("AC*");
  g_dist->SetMarkerStyle(20);
  g_dist->SetMarkerColor(1);
  g_dist->SetLineStyle(1);
  g_dist->SetLineWidth(1);
  g_dist->SetTitle("Gas");
  g_dist->GetXaxis()->SetTitle("polar angle [deg]");
  g_dist->GetYaxis()->SetTitle("#sigma_{dist}^{g}");*/

  /*
  for(Int_t i=0;i<64;i++){
    delete nph_e_aerogel_h[i];
    delete nph_pi_aerogel_h[i];
    delete nph_k_aerogel_h[i];
    delete nph_p_aerogel_h[i];

    delete nph_e_gas_h[i];
    delete nph_pi_gas_h[i];
    delete nph_k_gas_h[i];
    delete nph_p_gas_h[i];

    delete tc_e_aerogel_h[i];
    delete tc_pi_aerogel_h[i];
    delete tc_k_aerogel_h[i];
    delete tc_p_aerogel_h[i];

    delete tc_e_gas_h[i];
    delete tc_pi_gas_h[i];
    delete tc_k_gas_h[i];
    delete tc_p_gas_h[i];

    delete dist_pi_aerogel_h[i];
    delete dist_pi_gas_h[i];

    delete tc_pi_aerogel_ene_h[i];
    delete tc_pi_aerogel_ene_bk_h[i];

    //if(i<=index)delete cch[i];
    //if(i<=index)delete cch_ph[i];
    //if(i<=index)delete cch_d[i];

  }

  delete[] nph_e_aerogel_h;
  delete[] nph_pi_aerogel_h;
  delete[] nph_k_aerogel_h;
  delete[] nph_p_aerogel_h;

  delete[] nph_e_gas_h;
  delete[] nph_pi_gas_h;
  delete[] nph_k_gas_h;
  delete[] nph_p_gas_h;

  delete[] tc_e_aerogel_h;
  delete[] tc_pi_aerogel_h;
  delete[] tc_k_aerogel_h;
  delete[] tc_p_aerogel_h;

  delete[] tc_e_gas_h;
  delete[] tc_pi_gas_h;
  delete[] tc_k_gas_h;
  delete[] tc_p_gas_h;
  
  delete[] dist_pi_aerogel_h;
  delete[] dist_pi_gas_h;
  
  delete[] tc_pi_aerogel_ene_h;
  delete[] tc_pi_aerogel_ene_bk_h;


  delete[] cch;
  delete[] cch_ph;
  delete[] cch_d;
  delete dist_ph;
  delete file;
  */  

}


void eic_dual_rich::acq_QE(string input_filename, int select){

  ifstream inputFile;
  string line;

  Double_t ev,nm,qe;

  Double_t QE_H12700[41];
  Double_t QE_H12700_03_WLS_meas[41] = {
    0.016,0.02,0.0243455,0.0349796,0.0400769,0.0495496,0.054666,0.0612895,0.0758019,0.0853365,0.100662,0.121331,0.144678,0.162644,0.180719,0.194414,0.202599,0.224051,0.235051,0.253334,0.268143,0.285398,0.30002,0.309013,0.319247,0.328839,0.333333,0.335,0.33337,0.327161,0.321697,0.328776,0.333637,0.318123,0.313051,0.326953,0.331335,0.331335,0.331335,0.331335,0.331335
  };

  inputFile.open(input_filename.c_str());
  if (!inputFile) {
    cerr << "Error opening imput file";
    return;
  }
  else{
    cout<<input_filename.c_str()<<endl;
  }

  while(inputFile.good()){
    getline(inputFile, line);
    istringstream iss(line);
    while(iss>>ev>>nm>>qe){
      
      ev_a[phd_count]=ev;
      nm_a[phd_count]=nm;
      if(select == 0) qe_a[phd_count]=qe;
      if(select == 1) qe_a[phd_count]=QE_H12700_03_WLS_meas[phd_count];
      QE_H12700[phd_count]=qe;

      //cout<<ev<<"  "<<nm<<"  "<<qe<<"  "<<qe_a[phd_count]<<"  "<<phd_count<<endl;

      phd_count++;

    }

  }

  TGraph *gqe = new TGraph(phd_count,nm_a,QE_H12700);
  TGraph *gqe1 = new TGraph(phd_count,nm_a,QE_H12700_03_WLS_meas);

  TCanvas *qeff = new TCanvas("qeff","",800,800);
  qeff->Divide(1,1);
  qeff->cd(1);
  gqe->Draw("AC*");
  gqe->SetMarkerStyle(20);
  gqe->SetMarkerColor(1);
  gqe->SetLineStyle(1);
  gqe->SetLineWidth(2);
  gqe->SetLineColor(1);
  gqe->SetTitle("");
  gqe->GetXaxis()->SetTitle("#lambda [nm]");
  gqe->GetYaxis()->SetTitle("QE");
  gqe1->Draw("C*");
  gqe1->SetMarkerStyle(20);
  gqe1->SetMarkerColor(2);
  gqe1->SetLineStyle(1);
  gqe1->SetLineWidth(2);
  gqe1->SetLineColor(2);

}

void eic_dual_rich::err_plots(Int_t points, Int_t label){

  TGraph *a_ch = new TGraph(points,angl,ch_er_a);
  TGraph *g_ch = new TGraph(points,angl,ch_er_g);
  TGraph *a_emi = new TGraph(points,angl,emi_er_a);
  TGraph *g_emi = new TGraph(points,angl,emi_er_g);
  TGraph *a_emi1 = new TGraph(points,angl,emi_er_a);
  TGraph *g_emi1 = new TGraph(points,angl,emi_er_g);
  TGraph *a_pix = new TGraph(points,angl,pix_er_a);
  TGraph *g_pix = new TGraph(points,angl,pix_er_g);
  TGraph *a_mag = new TGraph(points,angl,mg_er_a);
  TGraph *g_mag = new TGraph(points,angl,mg_er_g);
  TGraph *a_tr = new TGraph(points,angl,tr_er_a);
  TGraph *g_tr = new TGraph(points,angl,tr_er_g);

  TCanvas *errors = new TCanvas("errors","",800,800);
  errors->Divide(1,2);
  errors->cd(1);
  a_ch->Draw("AC*");
  a_ch->SetMarkerStyle(20);
  a_ch->SetMarkerColor(1);
  a_ch->SetLineStyle(1);
  a_ch->SetLineWidth(1);
  a_ch->SetTitle("Aerogel");
  a_ch->GetXaxis()->SetTitle("polar angle [deg]");
  a_ch->GetYaxis()->SetTitle("#sigma_{(1 p.e.)}[rad]");
  a_emi->Draw("C*");
  a_emi->SetMarkerStyle(20);
  a_emi->SetMarkerColor(2);
  a_emi->SetLineStyle(1);
  a_emi->SetLineColor(2);
  a_emi->SetLineWidth(1);
  if(label =! 0) a_emi1->Draw("C*");
  a_emi1->SetMarkerStyle(20);
  a_emi1->SetMarkerColor(2);
  a_emi1->SetLineStyle(2);
  a_emi1->SetLineColor(2);
  a_emi1->SetLineWidth(1);
  a_pix->Draw("C*");
  a_pix->SetMarkerStyle(20);
  a_pix->SetMarkerColor(3);
  a_pix->SetLineStyle(1);
  a_pix->SetLineColor(3);
  a_pix->SetLineWidth(1);
  a_mag->Draw("C*");
  a_mag->SetMarkerStyle(20);
  a_mag->SetMarkerColor(4);
  a_mag->SetLineStyle(1);
  a_mag->SetLineColor(4);
  a_mag->SetLineWidth(1);
  a_tr->Draw("C*");
  a_tr->SetMarkerStyle(20);
  a_tr->SetMarkerColor(6);
  a_tr->SetLineStyle(1);
  a_tr->SetLineColor(6);
  a_tr->SetLineWidth(1);
  errors->cd(2);
  g_ch->Draw("AC*");
  g_ch->SetMarkerStyle(20);
  g_ch->SetMarkerColor(1);
  g_ch->SetLineStyle(1);
  g_ch->SetLineWidth(1);
  g_ch->SetTitle("C_{2}F_{6} gas");
  g_ch->GetXaxis()->SetTitle("polar angle [deg]");
  g_ch->GetYaxis()->SetTitle("#sigma_{(1 p.e.)}[rad]");
  g_emi->Draw("C*");
  g_emi->SetMarkerStyle(20);
  g_emi->SetMarkerColor(2);
  g_emi->SetLineStyle(1);
  g_emi->SetLineColor(2);
  g_emi->SetLineWidth(1);
  if(label =! 0) g_emi1->Draw("C*");
  g_emi1->SetMarkerStyle(20);
  g_emi1->SetMarkerColor(2);
  g_emi1->SetLineStyle(2);
  g_emi1->SetLineColor(2);
  g_emi1->SetLineWidth(1);
  g_pix->Draw("C*");
  g_pix->SetMarkerStyle(20);
  g_pix->SetMarkerColor(3);
  g_pix->SetLineStyle(1);
  g_pix->SetLineColor(3);
  g_pix->SetLineWidth(1);
  g_mag->Draw("C*");
  g_mag->SetMarkerStyle(20);
  g_mag->SetMarkerColor(4);
  g_mag->SetLineStyle(1);
  g_mag->SetLineColor(4);
  g_mag->SetLineWidth(1);
  g_tr->Draw("C*");
  g_tr->SetMarkerStyle(20);
  g_tr->SetMarkerColor(6);
  g_tr->SetLineStyle(1);
  g_tr->SetLineColor(6);
  g_tr->SetLineWidth(1);

}

void eic_dual_rich::nph_plots(Int_t points){

  Double_t ratio[points];
  Double_t eff[points], eff_shield[points];

  for(Int_t i=0;i<points;i++){

    ratio[i] = (1.-Nph_a_shield[i]/Nph_a[i])*100.; 
    eff[i] = TMath::Exp(-Nph_a[i])*(1.+Nph_a[i]+Nph_a[i]*Nph_a[i]/2.)*100.;
    eff_shield[i] = TMath::Exp(-Nph_a_shield[i])*(1.+Nph_a_shield[i]+Nph_a_shield[i]*Nph_a_shield[i]/2.)*100.;

  }

  TGraphErrors *a_nph = new TGraphErrors(points,angl,Nph_a,0,err_Nph_a);
  TGraphErrors *a_nph_shield = new TGraphErrors(points,angl,Nph_a_shield,0,err_Nph_a_shield);

  TGraphErrors *a_nph_ratio = new TGraphErrors(points,angl,ratio,0,0);

  TGraphErrors *a_nph_eff = new TGraphErrors(points,angl,eff,0,0);
  TGraphErrors *a_nph_eff_s = new TGraphErrors(points,angl,eff_shield,0,0);

  
  TCanvas *c_nph = new TCanvas("c_nph","",800,800);
  c_nph->Divide(1,2);
  TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.2,1.0,1.0);
  TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.2);
  pad1->Draw();
  //pad2->Draw();
  pad1->cd();
  //a_nph->Draw("AB1");
  a_nph->SetFillColor(kOrange+9);
  a_nph->SetFillStyle(3001);
  a_nph->SetLineColor(kOrange+7);
  //a_nph->SetMarkerStyle(20);
  //a_nph->SetMarkerColor(4);
  //a_nph->SetLineStyle(1);
  //a_nph->SetLineWidth(1);
  a_nph_shield->SetTitle("");
  a_nph_shield->GetXaxis()->SetTitle("polar angle [deg]");
  a_nph_shield->GetYaxis()->SetTitle("Mean number of p.e.");
  a_nph_shield->Draw("AB1");
  a_nph_shield->SetFillColor(kBlue+3);
  a_nph_shield->SetFillStyle(3001);
  a_nph_shield->SetLineColor(kBlue+2);
  /*pad2->cd();
  a_nph_ratio->Draw("AB1");
  a_nph_ratio->GetXaxis()->SetLabelSize(0);
  a_nph_ratio->GetYaxis()->SetTitle("Diff. %");
  a_nph_ratio->GetYaxis()->SetRangeUser(0,15);
  a_nph_ratio->GetYaxis()->SetNdivisions(5, 5, 0, kTRUE);
  a_nph_ratio->GetYaxis()->SetTitleSize(0.18);
  a_nph_ratio->GetYaxis()->SetTitleOffset(0.25);
  a_nph_ratio->GetYaxis()->SetLabelSize(0.15);
  a_nph_ratio->SetFillColor(38);
  a_nph_ratio->SetFillStyle(3001);
  a_nph_ratio->SetTitle("");*/
  
  TCanvas *c_eff = new TCanvas("c_eff","",800,800);
  c_eff->Divide(1,1);
  c_eff->cd(1);
  a_nph_eff_s->Draw("AB1");
  a_nph_eff_s->SetFillColor(kBlue+3);
  a_nph_eff_s->SetFillStyle(3001);
  a_nph_eff_s->SetLineColor(kBlue+2);
  //a_nph->SetMarkerStyle(20);
  //a_nph->SetMarkerColor(4);
  //a_nph->SetLineStyle(1);
  //a_nph->SetLineWidth(1);
  a_nph_eff_s->SetTitle("");
  a_nph_eff_s->GetXaxis()->SetTitle("polar angle [deg]");
  a_nph_eff_s->GetYaxis()->SetTitle("Ineff. %");
  //a_nph_eff->Draw("B1");
  a_nph_eff->SetFillColor(kOrange+9);
  a_nph_eff->SetFillStyle(3001);
  a_nph_eff->SetLineColor(kOrange+7);
  
}
