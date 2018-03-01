/*******************
Usage in ROOT:
.L rich_analyzer_pythia.cpp++
eic_dual_rich *a = new eic_dual_rich()
a->acq_QE("h12700_pmt.txt",0)
a->acq_phytia("out ... .root",numberofevents,pid_table(from 0 to 3, e,pi,K,p))
********************/

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

  void acq_phytia(string input_filename, int event_n, int pid_type);

  Int_t phd_count = 0;
  Double_t ev_a[500], nm_a[500], qe_a[500];
  
  Double_t ind_ray1(double Ex, double Ey, double Ez, double Dx, double Dy, double Dz, double cx, double cy, double cz, double vx, double vy, double vz, int select_radiator);
  void acq_QE(string input_filename, int select);
  double resolution(double pol, int radiator);
  
  Double_t result;
  

};

double eic_dual_rich::resolution(double pol, int radiator){

  Double_t *aerogel = new Double_t[5];
  aerogel[0]=0.00269856; aerogel[1]=0.00247247; aerogel[2]=0.00243567; aerogel[3]=0.00253132; aerogel[4]=0.00204384;
  Double_t *gas = new Double_t[5];
  gas[0]=0.00190304; gas[1]=0.00137526; gas[2]=0.0015026; gas[3]=0.00182912; gas[4]=0.00220638;
  Double_t *mom = new Double_t[5];
  mom[0]=5; mom[1]=10; mom[2]=15; mom[3]=20; mom[4]=25;  
  
  TSpline3 *res_a;
  res_a = new TSpline3("",mom,aerogel,5);
  TSpline3 *res_g;
  res_g = new TSpline3("",mom,gas,5);
  
  if(radiator==1) result=res_a->Eval(pol);
  if(radiator==2) result=res_g->Eval(pol);

  delete[] aerogel;
  delete[] gas;
  delete[] mom;
  delete res_a;
  delete res_g;

  return result;
  
  
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

Double_t eic_dual_rich::ind_ray1(double Ex, double Ey, double Ez, double Dx, double Dy, double Dz, double cx, double cy, double cz, double vx, double vy, double vz, int select_radiator){

  //Int_t sel_radiator = 2;

  //if(select_radiator == 1) sel_radiator = 1;

  LongDouble_t cex,cey,cez;
  LongDouble_t cdx,cdy,cdz;

  Int_t i,iwhere;

  LongDouble_t th,a,d;
  LongDouble_t x,dx;

  LongDouble_t y,y1;

  LongDouble_t eps = 0.00000000001;
  LongDouble_t R = 0;  

  R = 289.9;
  Double_t refidx1 = 1.02;
  Double_t refidx2 = 1.00086;

  Double_t esx, esy, esz, es;
  Double_t ref_frac, theta1, theta2;

  Double_t theta_c = 0.;

  cex = -cx+Ex;
  cey = -cy+Ey;
  cez = -cz+Ez;

  cdx = -cx+Dx;
  cdy = -cy+Dy;
  cdz = -cz+Dz;

  a = TMath::Sqrt(cex*cex+cey*cey+cez*cez);
  d = TMath::Sqrt(cdx*cdx+cdy*cdy+cdz*cdz);
  th = TMath::ACos((cdx*cex+cdy*cey+cdz*cez)/a/d);

  i = 0;
  x = th/2.;  
  y = R*(a*sin(x)-d*sin(th-x))+a*d*sin(th-2*x);
  y1 = R*(a*cos(x)+d*cos(th-x))-2*a*d*cos(th-2*x);
  dx = -y/y1;

  while(TMath::Abs(dx)>eps && i<100){

    x+=dx;
    y = R*(a*sin(x)-d*sin(th-x))+a*d*sin(th-2*x);
    y1 = R*(a*cos(x)+d*cos(th-x))-2*a*d*cos(th-2*x);
    dx = -y/y1;
    i++;

  }

  if(i>=100) cout<<"Not convergent"<<endl;

  if(i<100){
    sx = cx + (R*cos(x)/a-R*sin(x)/tan(th)/a)*cex + (R*sin(x)/d/sin(th))*cdx;
    sy = cy + (R*cos(x)/a-R*sin(x)/tan(th)/a)*cey + (R*sin(x)/d/sin(th))*cdy;
    sz = cz + (R*cos(x)/a-R*sin(x)/tan(th)/a)*cez + (R*sin(x)/d/sin(th))*cdz;
  }

  esx = sx - Ex;
  esy = sy - Ey;
  esz = sz - Ez;

  es = sqrt(esx*esx+esy*esy+esz*esz);

  esx = esx/es;
  esy = esy/es;
  esz = esz/es;

  if(select_radiator == 1){
    ref_frac = refidx2/refidx1;
  
    theta2 = TMath::ACos(esz);
    theta1 = TMath::ASin(TMath::Sin(theta2)*ref_frac);
    
    esx = esx*TMath::Sin(theta1)/TMath::Sin(theta2);
    esy = esy*TMath::Sin(theta1)/TMath::Sin(theta2);
    esz = TMath::Cos(theta1);
  }
  theta_c = TMath::ACos((esx*vx+esy*vy+esz*vz));

  return theta_c;

}


//*******acq data function*******//

void eic_dual_rich::acq_phytia(string input_filename, int event_n, int pid_type){

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
  gStyle->SetPalette(53);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  //gStyle->SetTextColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetTitleFillColor(0); 
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.3); 
  //gStyle->SetOptStat(0);
  gStyle->SetStatBorderSize(0);
  gStyle->SetStatColor(0);
  gStyle->SetStatW(0.35);
  gStyle->SetStatH(0.35);

  LongDouble_t dir_x=0.;
  LongDouble_t dir_y=0.;
  LongDouble_t dir_z=0.;

  LongDouble_t in_px=0.;
  LongDouble_t in_py=0.;
  LongDouble_t in_pz=0.;

  LongDouble_t theta_t=0.;
  LongDouble_t phi_t=0.;


  LongDouble_t mir_x = 145.; //rich2
  LongDouble_t mir_y = 0.;
  LongDouble_t mir_z = 122.5; //rich2 in gemc coordinates 250+82.5-210.
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

  TFile *file=new TFile(input_filename.c_str());
  if (file->IsZombie()) {
    cout << "Error opening file" << input_filename << endl;
    exit(-1);
  }
  else cout << "open file " << input_filename << endl;

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
  
  Int_t NPART = 30;
  Int_t *molteplicity = new Int_t[eic_rich->GetEntries()];
  for(Int_t i=0;i<eic_rich->GetEntries();i++)molteplicity[i]=0;
  TH1F *molt_h = new TH1F("","",20,0,20);
  
  TH2F *ph_det = new TH2F("ph_det","", 1666, -250., 250., 1666, -250., 250.);
  TH1F **ch_h;
  ch_h=new TH1F*[NPART];
  for(Int_t i=0;i<NPART;i++){
    ch_h[i] = new TH1F(Form("track%d",i),"",1000,0.,1.);
  }
  TH1F **chb_h;
  chb_h=new TH1F*[NPART];
  for(Int_t i=0;i<NPART;i++){
    chb_h[i] = new TH1F(Form("track%d",i),"",1000,0.,1.);
  }
  TH1F **ch_part;
  ch_part=new TH1F*[4];
  for(Int_t i=0;i<4;i++){
    ch_part[i] = new TH1F(Form("part%d",i),"",100,0.,100.);
  }
  TH1F **pi_id_h;
  pi_id_h=new TH1F*[4];
  for(Int_t i=0;i<4;i++){
    pi_id_h[i] = new TH1F(Form("1vs%d",i),"",50,0.,50.);
  }
  TCanvas **ch;
  ch = new TCanvas*[NPART];
  
  Int_t track_number = 0;
  Int_t p_id[NPART];
  Double_t likel_a = 0.;
  Double_t likel_g = 0.;
  LongDouble_t emi_x[NPART];
  LongDouble_t emi_y[NPART];
  LongDouble_t emi_z[NPART];
  LongDouble_t emig_x[NPART];
  LongDouble_t emig_y[NPART];
  LongDouble_t emig_z[NPART];
  Double_t mom_p[NPART];
  Int_t sector_id[NPART];
  Double_t mom_v[NPART][3];
  Double_t pol_angle[NPART];
  Double_t phi_angle[NPART];
  Int_t p_type[NPART];
  Double_t the_a[NPART][4];
  Double_t the_g[NPART][4];
  Double_t tch_dist[NPART][4], tch_dist1[NPART][4];
  Double_t a_like[NPART][4], g_like[NPART][4];
  Int_t a_counts[NPART][4], g_counts[NPART][4], sub_cycle[NPART];
  Double_t refer[4];
  Double_t n_a[NPART][4], n_g[NPART][4];

  for(Int_t ii=0;ii<NPART;ii++){
    mom_v[ii][0] = 0.;
    mom_v[ii][1] = 0.;
    mom_v[ii][2] = 0.;
    pol_angle[ii] = 0.;
    phi_angle[ii] = 0.; 
    mom_p[ii] = 0.;
    sector_id[ii] = 0.;

    the_a[ii][0] = 0.;
    the_a[ii][1] = 0.;
    the_a[ii][2] = 0.;
    the_a[ii][3] = 0.;
      
    the_g[ii][0] = 0.;
    the_g[ii][1] = 0.;
    the_g[ii][2] = 0.;
    the_g[ii][3] = 0.;

    n_a[ii][0] = 0.;
    n_a[ii][1] = 0.;
    n_a[ii][2] = 0.;
    n_a[ii][3] = 0.;
      
    n_g[ii][0] = 0.;
    n_g[ii][1] = 0.;
    n_g[ii][2] = 0.;
    n_g[ii][3] = 0.;      
      
    tch_dist[ii][0] = 0.;
    tch_dist[ii][1] = 0.;
    tch_dist[ii][2] = 0.;
    tch_dist[ii][3] = 0.;

    tch_dist1[ii][0] = 0.;
    tch_dist1[ii][1] = 0.;
    tch_dist1[ii][2] = 0.;
    tch_dist1[ii][3] = 0.;
      
    a_like[ii][0] = 0.;
    a_like[ii][1] = 0.;
    a_like[ii][2] = 0.;
    a_like[ii][3] = 0.;
      
    a_counts[ii][0] = 0.;
    a_counts[ii][1] = 0.;
    a_counts[ii][2] = 0.;
    a_counts[ii][3] = 0.;

    g_like[ii][0] = 0.;
    g_like[ii][1] = 0.;
    g_like[ii][2] = 0.;
    g_like[ii][3] = 0.;
      
    g_counts[ii][0] = 0.;
    g_counts[ii][1] = 0.;
    g_counts[ii][2] = 0.;
    g_counts[ii][3] = 0.;
      
    emi_x[ii] = 0.;
    emi_y[ii] = 0.;
    emi_z[ii] = 0.;

    emig_x[ii] = 0.;
    emig_y[ii] = 0.;
    emig_z[ii] = 0.;

    sub_cycle[ii] = 0;
    p_id[ii] = 999.;
  }

  track_number = 0.; 
  
  //Int_t event_n = 2;
  Double_t qe_p = 0.;
  Double_t c_angle = 0.;
  Double_t c_angle_a = 0.;

  Double_t m[4]={0.000511, 0.13957018, 0.493677, 0.938272};

  Double_t xbin = 0.;
  Double_t ybin = 0.;

  Int_t idx0 = 0;
  Int_t idx1 = 1;
  Int_t evnt = 0;
  
  for(Int_t i=0;(i<eic_rich->GetEntries());i++){

    generated->GetEntry(i);
    eic_rich->GetEntry(i);
    evnt = i;
    
    TH2F *ph_det_s = new TH2F("ph_det","", 1666, -250., 250., 1666, -250., 250.);

    for(Int_t j=0;j<eic_rich_id->size();j++){

      if(i<=event_n && eic_rich_pid->at(j)!=0 && eic_rich_pid->at(j)!=22 && eic_rich_pid->at(j)!=2112 && eic_rich_vz->at(j)<2500 && eic_rich_out_z->at(j)==2540) {
	
	if(eic_rich_id->at(j)<20){
	  //cout<<"PID NUMBER: "<<eic_rich_pid->at(j)<<endl;
	  if(eic_rich_pid->at(j)==11){
	    //ch_part[0]->Fill(sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j))/1000.);
	    p_id[track_number] = 0;
	  }
	  if(eic_rich_pid->at(j)==211){
	    //ch_part[1]->Fill(sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j))/1000.);
	    p_id[track_number] = 1;
	  }
	  if(eic_rich_pid->at(j)==321){
	    //ch_part[2]->Fill(sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j))/1000.);
	    p_id[track_number] = 2;
	  }
	  if(eic_rich_pid->at(j)==2212){
	    //ch_part[3]->Fill(sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j))/1000.);
	    p_id[track_number] = 3;
	  }
	  
	  //molteplicity[i]++;
	  //cout<<"MOMENTUM: "<<sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j))/1000.<<endl;
	}

	mom_v[track_number][0] = eic_rich_in_px->at(j)/sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j));
	mom_v[track_number][1] = eic_rich_in_py->at(j)/sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j));
	mom_v[track_number][2] = eic_rich_in_pz->at(j)/sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j));

	pol_angle[track_number] = TMath::ACos(mom_v[track_number][2])*57.2958;
	phi_angle[track_number] = TMath::ATan2(mom_v[track_number][1],mom_v[track_number][0])*57.2958;
	
	mom_p[track_number] = sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j))/1000.;

	if(pol_angle[track_number]>3. && pol_angle[track_number]<25.){
	  ch_part[p_id[track_number]]->Fill(mom_p[track_number]);
	}
	
	for(Int_t idp=0;idp<4;idp++){
	  the_a[track_number][idp]=acos(sqrt(mom_p[track_number]*mom_p[track_number]+m[idp]*m[idp])/mom_p[track_number]/1.02);
	  refer[idp] = 8./TMath::Sin(acos(sqrt(50.*50.+m[idp]*m[idp])/50./1.02))/TMath::Sin(acos(sqrt(50.*50.+m[idp]*m[idp])/50./1.02));
	  n_a[track_number][idp]=207.*sin(the_a[track_number][idp])*sin(the_a[track_number][idp]);

	  the_g[track_number][idp]=acos(sqrt(mom_p[track_number]*mom_p[track_number]+m[idp]*m[idp])/mom_p[track_number]/1.00082);
	  refer[idp] = 20./TMath::Sin(acos(sqrt(50.*50.+m[idp]*m[idp])/50./1.00086))/TMath::Sin(acos(sqrt(50.*50.+m[idp]*m[idp])/50./1.00086));
	  n_g[track_number][idp]=11600.*sin(the_g[track_number][idp])*sin(the_g[track_number][idp]);
	  //cout<<n_a[track_number][idp]<<"  "<<n_g[track_number][idp]<<"  "<<idp<<endl;
	}

	//cout<<"Momentum_versor: "<<mom_v[track_number][0]<<"  "<<mom_v[track_number][1]<<"  "<<mom_v[track_number][2]<<endl;
	
	emi_x[track_number] = eic_rich_avg_x->at(j)*0.1;
	emi_y[track_number] = eic_rich_avg_y->at(j)*0.1;
	emi_z[track_number] = eic_rich_avg_z->at(j)*0.1;

	//cout<<"Emission: "<<emi_x[track_number]<<"  "<<emi_y[track_number]<<"  "<<emi_z[track_number]<<endl;

	p_type[track_number] = eic_rich_tid->at(j);
	sector_id[track_number] = eic_rich_id->at(j)-10;
	
	//cout<<"Track number:  "<<p_type[track_number]<<endl;
	//cout<<"Sector id:  "<<sector_id[track_number]<<endl;
	
	if(pol_angle[track_number]>3. && pol_angle[track_number]<25.){
	  track_number++;
	  molteplicity[i]++;
	}
	//track_number++;
	//molteplicity[i]++;
      }
      if(i<=event_n && eic_rich_pid->at(j)!=0 && eic_rich_pid->at(j)!=22 && eic_rich_vz->at(j)<2500 && eic_rich_id->at(j)<30 && eic_rich_id->at(j)>20 &&  eic_rich_tid->at(j)==p_type[track_number-1]) {
       	emig_x[track_number-1] = eic_rich_avg_x->at(j)*0.1;
	emig_y[track_number-1] = eic_rich_avg_y->at(j)*0.1;
	emig_z[track_number-1] = eic_rich_avg_z->at(j)*0.1;
      }

      
      //if(i<=event_n && eic_rich_pid->at(j)==0 && eic_rich_id->at(j)<50 && eic_rich_id->at(j)>40 && (eic_rich_mtid->at(j)==p_type[0] || eic_rich_mtid->at(j)==p_type[1]) && pol_angle[0]>3. && pol_angle[0]<25. && p_id[0]==1){ //********************uncomment if multiplicity 1
      if(i<=event_n && eic_rich_pid->at(j)==0 && eic_rich_id->at(j)<50 && eic_rich_id->at(j)>40){ //to be c ommented if multiplicity 1
	ph_energy = eic_rich_trackE->at(j)*1000000.;
	
	for(Int_t k=0;k<phd_count-1;k++){

	  if((eic_rich_trackE->at(j)*1000000.)>=ev_a[k] && (eic_rich_trackE->at(j)*1000000.)<ev_a[k+1]){
	    qe_p = qe_a[k];
	  }
	}

	if((eic_rich_trackE->at(j)*1000000.)<ev_a[0] || (eic_rich_trackE->at(j)*1000000.)>ev_a[phd_count-1]) qe_p = 0.;

	Int_t phi_count = eic_rich_id->at(j);
	phi_count = phi_count -41;
        mir_x = mir_xx*TMath::Cos(TMath::Pi()*-60.*phi_count/180.);
	mir_y = mir_xx*TMath::Sin(TMath::Pi()*-60.*phi_count/180.);
	
	if(rran.Uniform(0,1)<qe_p && ph_energy<=10. && ph_energy>=2.04358){

	  det_x = eic_rich_in_x->at(j)*0.1;
	  det_y = eic_rich_in_y->at(j)*0.1;
	  det_z = eic_rich_in_z->at(j)*0.1;
	  
	  ph_det_s->Fill(eic_rich_in_x->at(j)*0.1,eic_rich_in_y->at(j)*0.1);
	  
	  if(evnt==0)ph_det->Fill(eic_rich_in_x->at(j)*0.1,eic_rich_in_y->at(j)*0.1);

	  xbin = ph_det->GetXaxis()->FindBin(eic_rich_in_x->at(j)*0.1);
	  ybin = ph_det->GetYaxis()->FindBin(eic_rich_in_y->at(j)*0.1);
	  
	  if(ph_det_s->GetBinContent(xbin,ybin) == 1){
	    for(Int_t i=0;i<track_number;i++){

	      c_angle_a = ind_ray1(emi_x[i],emi_y[i],emi_z[i],ph_det_s->GetXaxis()->GetBinCenter(xbin),ph_det_s->GetYaxis()->GetBinCenter(ybin),det_z,mir_x,mir_y,mir_z,mom_v[i][0],mom_v[i][1],mom_v[i][2],1);
	      c_angle = ind_ray1(emig_x[i],emig_y[i],emig_z[i],ph_det_s->GetXaxis()->GetBinCenter(xbin),ph_det_s->GetYaxis()->GetBinCenter(ybin),det_z,mir_x,mir_y,mir_z,mom_v[i][0],mom_v[i][1],mom_v[i][2],2);
	      if(eic_rich_mtid->at(j)==p_type[i] && evnt==366)ch_h[i]->Fill(c_angle);
	      else if(evnt==366) chb_h[i]->Fill(c_angle);
	      
	      for(Int_t idp=0;idp<4;idp++){
		if(the_a[i][idp]>0. && c_angle_a<(the_a[i][idp]+resolution(pol_angle[i],1)*5) && c_angle_a>(the_a[i][idp]-resolution(pol_angle[i],1)*5)){
		  //if(the_a[i][idp]>0. && c_angle_a<(the_a[i][idp]+0.012) && c_angle_a>(the_a[i][idp]-0.012)){
		  tch_dist[i][idp] = tch_dist[i][idp] + c_angle_a;
		  a_counts[i][idp]++;
		  //if(i==1)cout<<c_angle<<"  "<<the_g[i][idp]+0.006<<endl;
		}
		if(the_g[i][idp]>0. && c_angle<(the_g[i][idp]+resolution(pol_angle[i],2)*5) && c_angle>(the_g[i][idp]-resolution(pol_angle[i],2)*5)){
		  //if(the_g[i][idp]>0. && c_angle<(the_g[i][idp]+0.0025*5) && c_angle>(the_g[i][idp]-0.0025*5)){
		  tch_dist1[i][idp] = tch_dist1[i][idp] + c_angle;
		  g_counts[i][idp]++;
		}
	      }
	    }
	  } //if ph_det_s ****************************
	}
      }
      
    } //end of j cycle 

    for(Int_t j=0;j<track_number;j++){
      //cout<<"TRACK NUMBER: "<<j<<endl;
      for(Int_t idp=0;idp<4;idp++){
	a_like[j][idp] = n_a[j][idp]*(the_a[j][idp]-(tch_dist[j][idp]/a_counts[j][idp]))*(the_a[j][idp]-(tch_dist[j][idp]/a_counts[j][idp]))/resolution(pol_angle[j],1)/resolution(pol_angle[j],1);
	g_like[j][idp] = n_g[j][idp]*(the_g[j][idp]-(tch_dist1[j][idp]/g_counts[j][idp]))*(the_g[j][idp]-(tch_dist1[j][idp]/g_counts[j][idp]))/resolution(pol_angle[j],2)/resolution(pol_angle[j],2);
	if(a_like[j][idp]>0 && a_counts[j][idp]>3){likel_a = a_like[j][idp];}
	else {a_like[j][idp] = 999.;}
	if(g_like[j][idp]>0 && g_counts[j][idp]>3){likel_g = g_like[j][idp];}
	else {g_like[j][idp] = 999.;}
	//if(a_like[j][idp] != 999 && g_like[j][idp] != 999) cout<<"LogLikelihood: "<<likel_a+likel_g<<"  Particle type: "<<idp<<endl;
	//if(a_like[j][idp] != 999 && g_like[j][idp] == 999) cout<<"LogLikelihood (no gas): "<<likel_a<<"  Particle type: "<<idp<<endl;
	//if(a_like[j][idp] == 999 && g_like[j][idp] == 999) cout<<"LogLikelihood (no aerogel and no gas)"<<endl;
	//g_like[j][idp] = exp(-g_like[j][idp]);
	//if(j==0)cout<<"AEROGEL"<<" Sector: "<<sector_id[j]<<"  "<<pol_angle[j]<<endl;
	//if(j==0)cout<<tch_dist[j][idp]/a_counts[j][idp]<<"  "<<a_counts[j][idp]<<"  "<<the_a[j][idp]<<"  "<<idp<<endl;
	//if(j==0)cout<<a_like[j][idp]<<" p_type  "<<p_id[j]<<" Nph:  "<<n_a[j][idp]<<" Mom: "<<mom_p[j]<<endl;
	//if(j==0 && evnt==366)cout<<"GAS"<<endl;
	//if(j==0 && evnt==366)cout<<tch_dist1[j][idp]/g_counts[j][idp]<<"  "<<g_counts[j][idp]<<"  "<<the_g[j][idp]<<"  "<<mom_p[j]<<endl;
	//if(j==0 && evnt==366)cout<<likel_g<<"  "<<idp<<endl;
	//if(j==2)cout<<"Sum: "<<a_like[j][idp]+g_like[j][idp]<<endl;
      }
      //cout<<"PID+++++++++++++++++++++++++++++++++++++++++"<<endl;
      idx0=0;
      for(Int_t idp=1;idp<4;idp++){
	if((a_like[j][idp]!=999 && g_like[j][idp]==999) || (a_like[j][idp]!=999 && g_counts[j][idp]<9)){
	  if((a_like[j][idp])<(a_like[j][idx0])) idx0 = idp;
	  //if(p_id[j]==2)cout<<"Test0: "<<a_like[j][idx0]<<"  "<<a_like[j][idp]<<endl;
	}
	else if(g_like[j][idp]!=999 && g_counts[j][idp]>=9 && mom_p[j]<30.){
	  if((a_like[j][idp]+g_like[j][idp])<(a_like[j][idx0]+g_like[j][idx0])) idx0 = idp;
	  //if(p_id[j]==2)cout<<"Test1: "<<a_like[j][idx0]+g_like[j][idx0]<<"  "<<a_like[j][idp]+g_like[j][idp]<<endl;
	}
	else if(g_like[j][idp]!=999 && mom_p[j]>=30.){
	  if((g_like[j][idp])<(g_like[j][idx0])) idx0 = idp;
	  //if(p_id[j]==2)cout<<"Test1: "<<a_like[j][idx0]+g_like[j][idx0]<<"  "<<a_like[j][idp]+g_like[j][idp]<<endl;
	}
	else if(a_like[j][idp]==999 && a_counts[j][idp-1]>4 && idp!=1){
	  //if(a_counts[j][idp-1]>a_counts[j][idp-2]+2) idx0=idp-1;
	  //if(p_id[j]==2)cout<<"Particle above th a: "<<idp<<endl;
	  //cout<<"Selected "<<idx0<<endl;
	  break;
	}
	else if(g_like[j][idp]==999 && g_counts[j][idp-1]>9 && idp!=1){
	  //if(g_counts[j][idp-1]>g_counts[j][idp-2]+4) idx0=idp-1;
	  //if(p_id[j]==2)cout<<"Particle above th g: "<<idp<<endl;
	  //if(evnt==3)cout<<"Selected "<<idx0<<endl;
	  break;
	}
	/*else if(g_like[j][idp]!=999 && (the_a[j][idx0]-the_a[j][idp])<0.009){
	//cout<<"Thetadiff: "<<(the_a[j][idx0]-the_a[j][idp])<<endl;
	if(p_id[j]==2)cout<<"Test1: "<<g_like[j][idx0]<<"  "<<g_like[j][idp]<<" "<<idx0<<"  "<<idp<<endl;
	if((g_like[j][idp])<(g_like[j][idx0])) idx0 = idp;
	}
	else if(g_like[j][idp]==999 && g_counts[j][idp-1]>9 && idp!=1){
	if(p_id[j]==2)cout<<"Particle above th: "<<idp<<endl;
	//if(evnt==3)cout<<"Selected "<<idx0<<endl;
	break;
	}
	else if(g_like[j][idp]==999 && a_like[j][idp]!=999){
	if(p_id[j]==2)cout<<"Test2: "<<a_like[j][idx0]<<"  "<<a_like[j][idp]<<" "<<idx0<<"  "<<idp<<endl;
	if((a_like[j][idp])<(a_like[j][idx0])) idx0 = idp;
	}*/
	//if(p_id[j]==2)cout<<"Selected "<<idx0<<endl;	
      }
      //if(p_id[j]==2)cout<<"Selected "<<idx0<<endl;
      if(mom_p[j]<=50 && mom_p[j]>=3. && pol_angle[j]>3. && pol_angle[j]<25.){ // Here PID evaluation .........................................
	if(p_id[j]==pid_type && idx0==0)pi_id_h[0]->Fill(mom_p[j]);
	if(p_id[j]==pid_type && idx0==1)pi_id_h[1]->Fill(mom_p[j]);
	if(p_id[j]==pid_type && idx0==2)pi_id_h[2]->Fill(mom_p[j]);
	if(p_id[j]==pid_type && idx0==3)pi_id_h[3]->Fill(mom_p[j]);
	//cout<<mom_p[j]<<endl;
      }
    }

    delete ph_det_s;

    //for(Int_t i=0;i<NPART;i++){
    //ch[i] = new TCanvas(Form("ch%d",i),"",700,500);
    //}
    if(evnt==0){
      
      for(Int_t ii=0;ii<track_number;ii++){
	ch[ii] = new TCanvas(Form("ch%d",ii),"",700,500);
	ch_h[ii]->SetFillColor(3);
	ch_h[ii]->Draw("");
	chb_h[ii]->Draw("same");
      }
    }
    
    for(Int_t ii=0;ii<NPART;ii++){
      mom_v[ii][0] = 0.;
      mom_v[ii][1] = 0.;
      mom_v[ii][2] = 0.;
      pol_angle[ii] = 0.;
      phi_angle[ii] = 0.; 
      mom_p[ii] = 0.;
      sector_id[ii] = 0.;

      the_a[ii][0] = 0.;
      the_a[ii][1] = 0.;
      the_a[ii][2] = 0.;
      the_a[ii][3] = 0.;
      
      the_g[ii][0] = 0.;
      the_g[ii][1] = 0.;
      the_g[ii][2] = 0.;
      the_g[ii][3] = 0.;

      n_a[ii][0] = 0.;
      n_a[ii][1] = 0.;
      n_a[ii][2] = 0.;
      n_a[ii][3] = 0.;
      
      n_g[ii][0] = 0.;
      n_g[ii][1] = 0.;
      n_g[ii][2] = 0.;
      n_g[ii][3] = 0.;  
      
      tch_dist[ii][0] = 0.;
      tch_dist[ii][1] = 0.;
      tch_dist[ii][2] = 0.;
      tch_dist[ii][3] = 0.;

      tch_dist1[ii][0] = 0.;
      tch_dist1[ii][1] = 0.;
      tch_dist1[ii][2] = 0.;
      tch_dist1[ii][3] = 0.;
      
      a_like[ii][0] = 0.;
      a_like[ii][1] = 0.;
      a_like[ii][2] = 0.;
      a_like[ii][3] = 0.;
      
      a_counts[ii][0] = 0.;
      a_counts[ii][1] = 0.;
      a_counts[ii][2] = 0.;
      a_counts[ii][3] = 0.;

      g_like[ii][0] = 0.;
      g_like[ii][1] = 0.;
      g_like[ii][2] = 0.;
      g_like[ii][3] = 0.;
      
      g_counts[ii][0] = 0.;
      g_counts[ii][1] = 0.;
      g_counts[ii][2] = 0.;
      g_counts[ii][3] = 0.;
      
      emi_x[ii] = 0.;
      emi_y[ii] = 0.;
      emi_z[ii] = 0.;

      emig_x[ii] = 0.;
      emig_y[ii] = 0.;
      emig_z[ii] = 0.;

      sub_cycle[ii] = 0;
      p_id[ii] = 0;
    }
    track_number = 0;

    if(i<=event_n)molt_h->Fill(molteplicity[i]);
    //cout<<"Molteplicity: "<<molteplicity[i]<<" Event n.: "<<i<<endl;
    
  } //end loop in i


  Double_t cont[4][100];
  Double_t p0[100],p1[100],p2[100],p3[100],momentum[100];

  for(Int_t k=0;k<100;k++){
    momentum[k] = k + 1.5;
    p0[k]=-99;
    p1[k]=-99;
    p2[k]=-99;
    p3[k]=-99;
    for(Int_t y=0;y<4;y++){
      cont[y][k] = pi_id_h[y]->GetBinContent(pi_id_h[0]->FindBin(k+1));
      //cout<<"Counts "<<y<<"  "<<k+1<<"  "<<cont[y][k]<<endl;
    }
    if(cont[0][k]!=0)p0[k] = cont[0][k]/(cont[0][k]+cont[1][k]+cont[2][k]+cont[3][k]);
    if(cont[1][k]!=0)p1[k] = cont[1][k]/(cont[0][k]+cont[1][k]+cont[2][k]+cont[3][k]);
    if(cont[2][k]!=0)p2[k] = cont[2][k]/(cont[0][k]+cont[1][k]+cont[2][k]+cont[3][k]);
    if(cont[3][k]!=0)p3[k] = cont[3][k]/(cont[0][k]+cont[1][k]+cont[2][k]+cont[3][k]);
  }

  TGraph *g0 = new TGraph(100,momentum,p0);
  TGraph *g1 = new TGraph(100,momentum,p1);
  TGraph *g2 = new TGraph(100,momentum,p2);
  TGraph *g3 = new TGraph(100,momentum,p3);
    
  //ph_det->Draw("");
  /*
    TCanvas *ch = new TCanvas("ch","",700,500);
    ch_h[0]->SetFillColor(3);
    chb_h[0]->Draw("");
    ch_h[0]->Draw("same");
    TCanvas *ch1 = new TCanvas("ch1","",700,500);
    ch_h[1]->SetFillColor(3);
    ch_h[1]->Draw("");
    chb_h[1]->Draw("same");*/

  TCanvas *ph = new TCanvas("ph","",500,500);
  ph_det->Draw("");
  //molt_h->Draw("");

  TCanvas *pa = new TCanvas("pa","",700,700);
  pa->Divide(2,2);
  pa->cd(1);
  ch_part[0]->Draw();
  for(Int_t k=0;k<4;k++)ch_part[k]->GetXaxis()->SetTitle("Momentum [GeV/c]");
  pa->cd(2);
  ch_part[1]->Draw();
  pa->cd(3);
  ch_part[2]->Draw();
  pa->cd(4);
  ch_part[3]->Draw();

  TCanvas *pa_pi = new TCanvas("pa_pi","",700,700);
  pa_pi->Divide(2,2);
  pa_pi->cd(1);
  pi_id_h[0]->Draw();
  //pi_id_h[0]->DrawNormalized("",pi_id_h[0]->GetEntries()/(pi_id_h[0]->GetEntries()+pi_id_h[1]->GetEntries()+pi_id_h[2]->GetEntries()+pi_id_h[3]->GetEntries()));
  for(Int_t k=0;k<4;k++)pi_id_h[k]->GetXaxis()->SetTitle("Momentum [GeV/c]");
  pa_pi->cd(2);
  pi_id_h[1]->Draw();
  //pi_id_h[1]->DrawNormalized("",pi_id_h[1]->GetEntries()/(pi_id_h[0]->GetEntries()+pi_id_h[1]->GetEntries()+pi_id_h[2]->GetEntries()+pi_id_h[3]->GetEntries()));
  pa_pi->cd(3);
  pi_id_h[2]->Draw();
  //pi_id_h[2]->DrawNormalized("",pi_id_h[2]->GetEntries()/(pi_id_h[0]->GetEntries()+pi_id_h[1]->GetEntries()+pi_id_h[2]->GetEntries()+pi_id_h[3]->GetEntries()));
  pa_pi->cd(4);
  pi_id_h[3]->Draw();
  //pi_id_h[3]->DrawNormalized("",pi_id_h[2]->GetEntries()/(pi_id_h[0]->GetEntries()+pi_id_h[1]->GetEntries()+pi_id_h[2]->GetEntries()+pi_id_h[3]->GetEntries()));

  TCanvas *cgr = new TCanvas("cgr","",700,700);
  cgr->Divide(2,2);
  cgr->cd(1);
  g0->GetXaxis()->SetRangeUser(2,50);
  g0->GetYaxis()->SetRangeUser(0,1);
  //g0->Draw("AP*");
  g0->SetFillColor(38);
  g0->Draw("AB");
  cgr->cd(2);
  g1->GetXaxis()->SetRangeUser(2,50);
  g1->GetYaxis()->SetRangeUser(0,1);
  //g1->Draw("AP*");
  g1->SetFillColor(38);
  g1->Draw("AB");
  cgr->cd(3);
  g2->GetXaxis()->SetRangeUser(2,50);
  g2->GetYaxis()->SetRangeUser(0,1);
  //g2->Draw("AP*");
  g2->SetFillColor(38);
  g2->Draw("AB");
  cgr->cd(4);
  g3->GetXaxis()->SetRangeUser(2,50);
  g3->GetYaxis()->SetRangeUser(0,1);
  g3->SetFillColor(38);
  g3->Draw("AB");
  //g3->Draw("AP*");
    
}

 
