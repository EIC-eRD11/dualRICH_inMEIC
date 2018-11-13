/*******************
Usage in ROOT:
.L rich_analyzer_pythia.cpp++
eic_dual_rich *a = new eic_dual_rich()
a->acq_QE("h12700_pmt.txt",0)
a->acq_phytia("out ... .root",numberofevents,pid_table(from 0 to 3, e,pi,K,p))
a->acq_pnew("../outputs/out.1.root",200)
********************/

#include <Riostream.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include <Math/ProbFunc.h>
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
#include <TGraph2D.h>

using namespace std;

//*******eic_dual_rich class*******//

class eic_dual_rich {

private:

public:

  LongDouble_t sx,sy,sz;

  void acq_pnew(string input_filename, int event_n);
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
  /*
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
  */

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

void eic_dual_rich::acq_pnew(string input_filename, int event_n){

  gStyle->SetCanvasColor(0);
  //gStyle->SetTitleW(1.2);
  //gStyle->SetTitleH(0.2);
  gStyle->SetTitleSize(0.055,"xyzt");
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleYOffset(1.4);
  gStyle->SetLabelSize(0.055,"XYZ");
  gStyle->SetPadTopMargin(.05);
  gStyle->SetPadLeftMargin(0.25);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(.25);
  gStyle->SetPalette(55);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  //gStyle->SetTextColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetTitleFillColor(0); 
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.3); 
  gStyle->SetOptStat(0);
  gStyle->SetStatBorderSize(0);
  gStyle->SetStatColor(0);
  gStyle->SetStatW(0.35);
  gStyle->SetStatH(0.35);

  Double_t dir_x=0.;
  Double_t dir_y=0.;
  Double_t dir_z=0.;

  Double_t in_px=0.;
  Double_t in_py=0.;
  Double_t in_pz=0.;

  Double_t theta_t=0.;
  Double_t phi_t=0.;


  Double_t mir_x = 145.; //rich2
  Double_t mir_y = 0.;
  Double_t mir_z = 122.5; //rich2 in gemc coordinates 250+82.5-210.
  Double_t mir_xx = mir_x;

  Double_t mir_R = 289.9;

  Double_t det_x=0.;
  Double_t det_y=0.;
  Double_t det_z=0.;

  Double_t ref_frac = 0.;

  Double_t det_xl=0.;
  Double_t det_yl=0.;
  Double_t det_zl=0.;

  Double_t det_xr=0.;
  Double_t det_yr=0.;
  Double_t det_zr=0.;

  Double_t det_xr1 = 0.; 
  Double_t det_yr1 = 0.;

  Int_t evnt = 0;

  TRandom rran;
  rran.SetSeed(0);
  
  Int_t NPART = 30;
  Int_t NHITS = 400;
  Int_t multiplicity[event_n];
  for(Int_t j=0;j<event_n;j++)multiplicity[j]=0;
  TH1F *molt_h = new TH1F("","",20,0,20);
  
  Int_t track_number = 0;
  Int_t p_id[NPART];
  Int_t t_id[NPART];
  Double_t mom_p[NPART];
  Int_t sector_id[NPART];
  Double_t mom_v[NPART][3];
  Double_t pol_angle[NPART];
  Double_t phi_angle[NPART];
  Int_t p_type[NPART];
  Double_t the_a[NPART][4];
  Double_t the_g[NPART][4];
  Double_t n_a[NPART][4], n_g[NPART][4];
  Double_t emia_x[NPART];
  Double_t emia_y[NPART];
  Double_t emia_z[NPART];
  Double_t emig_x[NPART];
  Double_t emig_y[NPART];
  Double_t emig_z[NPART];
  Double_t nph_a[NPART];  Double_t nph_g[NPART], chsum_a[NPART];  Double_t chsum_g[NPART];
  Double_t nph_a_bar[NPART];  Double_t nph_g_bar[NPART];

  Double_t m[4]={0.000511, 0.13957018, 0.493677, 0.938272};
  Int_t nsigma = 5;
  Double_t aidx = 1.019;
  Double_t gidx = 1.00086;

  Double_t ph_energy = 0.;
  Double_t qe_p = 0.;
  Int_t phi_count = 0;
  Int_t ph_count = 0;
  Double_t c_angle_a = 0.;
  Double_t c_angle_g = 0.;
  Double_t c_a[NPART][NHITS], c_g[NPART][NHITS];
  Double_t m_tr[NHITS];
  Double_t nph_cor[NHITS], nph_wro[NHITS];
  Int_t h_idx[NHITS];
  Float_t rdummy[NHITS];
  Double_t xbin = 0.;
  Double_t ybin = 0.;

  Double_t min_angle = 5.;
  Double_t max_angle = 25.;
  Double_t min_mom = 2.5;
  Double_t max_mom = 50.;

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
  
  TH1F **p_id_h;
  p_id_h=new TH1F*[16];
  for(Int_t i=0;i<16;i++){
    p_id_h[i] = new TH1F(Form("1vs%d",i),"",10,0.,50.);
  }


  TH1F **p_id_tot;
  TH1F **ch_part;
  ch_part=new TH1F*[4];
  p_id_tot=new TH1F*[4];
  for(Int_t i=0;i<4;i++){
    ch_part[i] = new TH1F(Form("part%d",i),"",100,0.,100.);
    p_id_tot[i] = new TH1F("","",10,0.,50.);
  }
  

  TH2F *ph_det = new TH2F("ph_det","", 1666, -250., 250., 1666, -250., 250.);
  TCanvas *c_ph = new TCanvas("c_ph","c_ph",600,600);

  TH2F *h2_bad = new TH2F("h2_bad","", 10, -0.5, 9.5, 10, -0.5, 9.5);
  TH1F *h_tot = new TH1F("h_tot","", 10, -0.5, 9.5);

  TH1F *tot_ph = new TH1F("tot_ph","",10,-0.5,9.5);
  TH1F *good_ph = new TH1F("good_ph","",10,-0.5,9.5);
  TH1F *bad_ph = new TH1F("bad_ph","",10,-0.5,9.5);
  TH1F *bg_ph  = new TH1F("bg_ph","",10,-0.5,9.5);
  TH1F *good_ph_m = new TH1F("good_ph_m","",50,-0.5,49.5);
  TH1F *bad_ph_m  = new TH1F("bad_ph_m","",50,-0.5,49.5);
  TH1F *tot_ph_m  = new TH1F("tot_ph_m","",50,-0.5,49.5);
  
  TChain *generated = new TChain("generated");
  generated->Add(input_filename.c_str());
  vector <double> *gen_pid=0;
  vector <double> *gen_px=0,*gen_py=0,*gen_pz=0,*gen_vx=0,*gen_vy=0,*gen_vz=0;
  generated->SetBranchAddress("pid",&gen_pid);
  generated->SetBranchAddress("px",&gen_px);
  generated->SetBranchAddress("py",&gen_py);
  generated->SetBranchAddress("pz",&gen_pz);
  generated->SetBranchAddress("vx",&gen_vx);
  generated->SetBranchAddress("vy",&gen_vy);
  generated->SetBranchAddress("vz",&gen_vz);

  TChain *flux = new TChain("flux");
  flux->Add(input_filename.c_str());

  vector<double> *flux_id=0,*flux_hitn=0;
  vector<double> *flux_pid=0,*flux_mpid=0,*flux_tid=0,*flux_mtid=0,*flux_otid=0;
  vector<double> *flux_trackE=0,*flux_totEdep=0,*flux_avg_x=0,*flux_avg_y=0,*flux_avg_z=0,*flux_avg_lx=0,*flux_avg_ly=0,*flux_avg_lz=0,*flux_px=0,*flux_py=0,*flux_pz=0,*flux_vx=0,*flux_vy=0,*flux_vz=0,*flux_mvx=0,*flux_mvy=0,*flux_mvz=0,*flux_avg_t=0;
  vector<double> *flux_procID=0,*flux_nsteps=0;
  flux->SetBranchAddress("hitn",&flux_hitn);  
  flux->SetBranchAddress("id",&flux_id);
  flux->SetBranchAddress("pid",&flux_pid);
  flux->SetBranchAddress("mpid",&flux_mpid);
  flux->SetBranchAddress("tid",&flux_tid);
  flux->SetBranchAddress("mtid",&flux_mtid);
  flux->SetBranchAddress("otid",&flux_otid);
  flux->SetBranchAddress("trackE",&flux_trackE);
  flux->SetBranchAddress("totEdep",&flux_totEdep);
  flux->SetBranchAddress("avg_x",&flux_avg_x);
  flux->SetBranchAddress("avg_y",&flux_avg_y);
  flux->SetBranchAddress("avg_z",&flux_avg_z);
  flux->SetBranchAddress("avg_lx",&flux_avg_lx);
  flux->SetBranchAddress("avg_ly",&flux_avg_ly);
  flux->SetBranchAddress("avg_lz",&flux_avg_lz);
  flux->SetBranchAddress("px",&flux_px);
  flux->SetBranchAddress("py",&flux_py);
  flux->SetBranchAddress("pz",&flux_pz);
  flux->SetBranchAddress("vx",&flux_vx);
  flux->SetBranchAddress("vy",&flux_vy);
  flux->SetBranchAddress("vz",&flux_vz);
  flux->SetBranchAddress("mvx",&flux_mvx);
  flux->SetBranchAddress("mvy",&flux_mvy);
  flux->SetBranchAddress("mvz",&flux_mvz);
  flux->SetBranchAddress("avg_t",&flux_avg_t);
  flux->SetBranchAddress("nsteps",&flux_nsteps);
  flux->SetBranchAddress("procID",&flux_procID);

  TChain *eic_rich = new TChain("eic_rich");
  eic_rich->Add(input_filename.c_str());

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

  if (event_n>eic_rich->GetEntries()) {
    cout <<" Error! Number of events > events in the Tree" << endl;
    exit(-1);
  }

  TRandom *rnd = new TRandom();
  rnd->SetSeed(0);
  
  for(Int_t i=0;i<event_n;i++){ //loop on the events
    
    generated->GetEntry(i);
    flux->GetEntry(i);
    eic_rich->GetEntry(i);
    
    evnt = i;
    track_number = 0;

    for(Int_t k=0;k<NPART;k++){
      for(Int_t kk=0;kk<400;kk++){
	c_a[k][kk]=0.;
	c_g[k][kk]=0.;

	m_tr[kk]=0.;
      }
    }
    
    //if(evnt%1000==0)cout<<evnt<<" events processed "<<endl;
    for(Int_t ii=0;ii<NPART;ii++){
    
      mom_v[ii][0] = 0.;
      mom_v[ii][1] = 0.;
      mom_v[ii][2] = 0.;
      pol_angle[ii] = 0.;
      phi_angle[ii] = 0.; 
      mom_p[ii] = 0.;
      sector_id[ii] = 0;

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

      emia_x[ii] = 0.;
      emia_y[ii] = 0.;
      emia_z[ii] = 0.;
    
      emig_x[ii] = 0.;
      emig_y[ii] = 0.;
      emig_z[ii] = 0.;
    
      p_id[ii] = 999;
      t_id[ii] = 999;
      nph_cor[ii] = 0;
      nph_wro[ii] = 0;
      p_type[ii] = 99;
    }

    ph_count = 0;

    for(Int_t j=0;j<eic_rich_id->size();j++){ //loop on the sub-events
      
      if((eic_rich_pid->at(j)==11 || abs(eic_rich_pid->at(j))==211  || eic_rich_pid->at(j)==2212 || abs(eic_rich_pid->at(j))==321) && eic_rich_vz->at(j)<2500 && eic_rich_out_z->at(j)==2540 && eic_rich_id->at(j)<20) {

	p_id[track_number] = eic_rich_pid->at(j);
	t_id[track_number] = eic_rich_tid->at(j);
	
	//if(eic_rich_id->at(j)<20){
	  if(eic_rich_pid->at(j)==11){
	    p_id[track_number] = 0;
	  }
	  if(abs(eic_rich_pid->at(j))==211){
	    p_id[track_number] = 1;
	  }
	  if(abs(eic_rich_pid->at(j))==321){
	    p_id[track_number] = 2;
	  }
	  if(eic_rich_pid->at(j)==2212){
	    p_id[track_number] = 3;
	  }
	  //if(evnt==0)cout<<"MOMENTUM: "<<sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j))/1000.<<endl;
	  //}

	mom_v[track_number][0] = eic_rich_in_px->at(j)/sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j));
	mom_v[track_number][1] = eic_rich_in_py->at(j)/sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j));
	mom_v[track_number][2] = eic_rich_in_pz->at(j)/sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j));

	//cout<<mom_v[track_number][0]<<" "<<mom_v[track_number][1]<<" "<<mom_v[track_number][2]<<endl;
	
	pol_angle[track_number] = TMath::ACos(mom_v[track_number][2])*57.2958;
	phi_angle[track_number] = TMath::ATan2(mom_v[track_number][1],mom_v[track_number][0])*57.2958;
	
	mom_p[track_number] = sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j))/1000.;

	if(pol_angle[track_number]>min_angle && pol_angle[track_number]<max_angle){
	  ch_part[p_id[track_number]]->Fill(mom_p[track_number]);
	}
	
	for(Int_t idp=0;idp<4;idp++){
	  the_a[track_number][idp]=acos(sqrt(mom_p[track_number]*mom_p[track_number]+m[idp]*m[idp])/mom_p[track_number]/aidx);
	  n_a[track_number][idp]=207.*sin(the_a[track_number][idp])*sin(the_a[track_number][idp]);
	  if(std::isnan(the_a[track_number][idp])){
	    n_a[track_number][idp]=0;
	    the_a[track_number][idp]=999.;
	  }

	  the_g[track_number][idp]=acos(sqrt(mom_p[track_number]*mom_p[track_number]+m[idp]*m[idp])/mom_p[track_number]/gidx);
	  n_g[track_number][idp]=11600.*2.7*sin(the_g[track_number][idp])*sin(the_g[track_number][idp]);
	  if(std::isnan(the_g[track_number][idp])){
	    n_g[track_number][idp]=0;
	    the_g[track_number][idp]=999.;
	  }
	  //cout<<"Nph: "<<n_a[track_number][idp]<<"  "<<n_g[track_number][idp]<<"  "<<idp<<"  "<<track_number<<endl;
	}

	//cout<<"Momentum_versor: "<<mom_v[track_number][0]<<"  "<<mom_v[track_number][1]<<"  "<<mom_v[track_number][2]<<endl;
	
	emia_x[track_number] = eic_rich_avg_x->at(j)*0.1;
	emia_y[track_number] = eic_rich_avg_y->at(j)*0.1;
	emia_z[track_number] = eic_rich_avg_z->at(j)*0.1;

	//cout<<"Emission: "<<emi_x[track_number]<<"  "<<emi_y[track_number]<<"  "<<emi_z[track_number]<<endl;

	p_type[track_number] = eic_rich_tid->at(j);
	sector_id[track_number] = eic_rich_id->at(j)-10;
	
	//cout<<"Track number:  "<<p_type[track_number]<<endl;
	//cout<<"Sector id:  "<<sector_id[track_number]<<endl;
	
	if(pol_angle[track_number]>min_angle && pol_angle[track_number]<max_angle && mom_p[track_number]>min_mom && mom_p[track_number]<max_mom){
	  track_number++;
	  multiplicity[i]++;
	  //cout<<track_number<<" "<<evnt<<endl;
	}
      }
      if(track_number>0){
	if((eic_rich_pid->at(j)==11 || abs(eic_rich_pid->at(j))==211 || eic_rich_pid->at(j)==2212 || abs(eic_rich_pid->at(j))==321) && eic_rich_vz->at(j)<2500 && eic_rich_id->at(j)<30 && eic_rich_id->at(j)>20 &&  eic_rich_tid->at(j)==p_type[track_number-1]) {
	  emig_x[track_number-1] = eic_rich_avg_x->at(j)*0.1;
	  emig_y[track_number-1] = eic_rich_avg_y->at(j)*0.1;
	  emig_z[track_number-1] = eic_rich_avg_z->at(j)*0.1;
	  //cout<<"G: "<<track_number<<" "<<evnt<<endl;
	}

	if(i<=event_n && eic_rich_pid->at(j)==0 && eic_rich_id->at(j)<50 && eic_rich_id->at(j)>40){ //********************reduced multiplicity
	  for(Int_t k=0;k<track_number;k++){
	    if(eic_rich_mtid->at(j)==p_type[k] && pol_angle[k]>min_angle && pol_angle[k]<max_angle && mom_p[k]>min_mom && mom_p[k]<max_mom){
	      //if(eic_rich_pid->at(j)==0 && eic_rich_id->at(j)<50 && eic_rich_id->at(j)>40){ //here photons are collected & cherenkov angles calculated
	
	      ph_energy = eic_rich_trackE->at(j)*1000000.;
	
	      for(Int_t k=0;k<phd_count-1;k++){

		if((eic_rich_trackE->at(j)*1000000.)>=ev_a[k] && (eic_rich_trackE->at(j)*1000000.)<ev_a[k+1]){
		  qe_p = qe_a[k];
		}
	      }

	      if((eic_rich_trackE->at(j)*1000000.)<ev_a[0] || (eic_rich_trackE->at(j)*1000000.)>ev_a[phd_count-1]) qe_p = 0.;

	      phi_count = eic_rich_id->at(j);
	      phi_count = phi_count -41;
	      mir_x = mir_xx*TMath::Cos(TMath::Pi()*-60.*phi_count/180.);
	      mir_y = mir_xx*TMath::Sin(TMath::Pi()*-60.*phi_count/180.);
	
	      if(rran.Uniform(0,1)<qe_p && ph_energy<=10. && ph_energy>=2.04358){
	  
		det_x = eic_rich_in_x->at(j)*0.1;
		det_y = eic_rich_in_y->at(j)*0.1;
		det_z = eic_rich_in_z->at(j)*0.1;

		m_tr[ph_count] = eic_rich_mtid->at(j);
	  
		ph_det->Fill(eic_rich_in_x->at(j)*0.1,eic_rich_in_y->at(j)*0.1);

		xbin = ph_det->GetXaxis()->FindBin(eic_rich_in_x->at(j)*0.1);
		ybin = ph_det->GetYaxis()->FindBin(eic_rich_in_y->at(j)*0.1);

		if(ph_det->GetBinContent(xbin,ybin) == 1){
		  for(Int_t i=0;i<track_number;i++){
		    c_angle_a = ind_ray1(emia_x[i],emia_y[i],emia_z[i],ph_det->GetXaxis()->GetBinCenter(xbin),ph_det->GetYaxis()->GetBinCenter(ybin),det_z,mir_x,mir_y,mir_z,mom_v[i][0],mom_v[i][1],mom_v[i][2],1);
		    c_a[i][ph_count] = c_angle_a;
		    //cout<<mom_p[0]<<"  "<<det_x<<"  "<<det_y<<"  "<<c_angle_a<<"  "<<evnt<<endl;
		    c_angle_g = ind_ray1(emig_x[i],emig_y[i],emig_z[i],ph_det->GetXaxis()->GetBinCenter(xbin),ph_det->GetYaxis()->GetBinCenter(ybin),det_z,mir_x,mir_y,mir_z,mom_v[i][0],mom_v[i][1],mom_v[i][2],2);
		    c_g[i][ph_count] = c_angle_g;
		    if(eic_rich_mtid->at(j)==p_type[i] && evnt==3)ch_h[i]->Fill(c_angle_a);
		    else if(evnt==3) chb_h[i]->Fill(c_angle_a);
		    //cout<<c_a[i][ph_count]<<" track number: "<<i<<" ph_number: "<<ph_count<<" event: "<<evnt<<endl;
	      
		  }

		  ph_count++;
	    
		} // if BinContent end
	      }

	    }
	  }
	
	} // end if photons in the detector
      }
    } //end loop sub-events  
    
    c_ph->cd(1);
    if(evnt==3){
      cout<<"Multiplicity is:  "<<track_number<<endl;
      ph_det->Draw("");
      ph_det->SetMarkerStyle(7);
      ph_det->GetXaxis()->SetRangeUser(-250,200);
      ph_det->GetYaxis()->SetRangeUser(-99.9,250);
      ph_det->GetXaxis()->SetTitle("cm");
      ph_det->GetYaxis()->SetTitle("cm");
      c_ph->SaveAs("figs/ph_det.png");

      //c_ph->Update();
      ch_h[0]->SetFillColor(12);
      chb_h[0]->SetFillColor(11);
      ch_h[0]->SetLineColor(12);
      chb_h[0]->SetLineColor(11);
      ch_h[0]->Draw("");
      chb_h[0]->Draw("same");
      ch_h[0]->GetXaxis()->SetTitle("#theta_{Ch} [rad]");
      ch_h[0]->GetYaxis()->SetTitle("");
      c_ph->SaveAs("figs/ph_histo.png");
      
    }
    ph_det->Reset();

    //cout<<"Nph: "<<ph_count<<" Nevent: "<<evnt<<" Ntracks: "<<track_number<<endl;
    
    if(track_number==0) continue;
    molt_h->Fill(track_number);
    
    Double_t lh_bar = -1E99;
    Int_t c_bar = -1;

    for(int c=0;c<TMath::Power(4,track_number);c++){
      //for(int c=(TMath::Power(4,track_number)-1);c>=0;c--){
      //cout<<c<<endl;
      Double_t lh = 0.;
      for(Int_t t=0;t<track_number;t++)nph_a[t]=nph_g[t]=chsum_a[t]=chsum_g[t]=0;
      Double_t nph_b = 0;
      rnd->RndmArray(ph_count,rdummy);
      TMath::Sort(ph_count, rdummy, h_idx);
      for(Int_t p=0;p<ph_count;p++){
	Double_t dmin = 0; Int_t rad_bar = -1; Int_t t_bar = -1; Double_t dmin_bar = 0;
	//Double_t pnot_a[4][track_number], pnot_g[4][track_number];
	dmin = 0.003*ROOT::Math::poisson_cdf_c(nph_b+1,1);
	for(Int_t t=0;t<track_number;t++){
	  Int_t part=3&(c>>(2*t));
	  //Double_t dnorm = TMath::Exp(-(the_a[t][part]-c_a[t][h_idx[p])*(the_a[t][part]-c_a[t][h_idx[p])/resolution(pol_angle[t],1)/resolution(pol_angle[t],1)/2.)/sqrt(2*TMath::Pi()*resolution(pol_angle[t],1)*resolution(pol_angle[t],1));
	  Double_t dnorm = TMath::Gaus(c_a[t][h_idx[p]],the_a[t][part],resolution(pol_angle[t],1),0);
	  //pnot_a[part][t]=dnorm;
	  Double_t pois = ROOT::Math::poisson_cdf_c(nph_a[t]+1,n_a[t][part]);
	  //cout<<p<<" "<<dnorm<<"  "<<pois<<" a "<<part<<"  "<<t<<" v  "<<c_a[t][h_idx[p]<<"  "<<the_a[t][part]<<"  "<<resolution(pol_angle[t],1)<<endl;
	  if((dnorm*pois)>dmin){
	    dmin=dnorm*pois;
	    rad_bar = 1;
	    t_bar = t;
	    dmin_bar = dnorm;
	  }
	  //dnorm = (the_g[t][part]-c_g[t][h_idx[p])*(the_g[t][part]-c_g[t][h_idx[p])/resolution(pol_angle[t],2)/resolution(pol_angle[t],2);
	  //dnorm = TMath::Exp(-(the_g[t][part]-c_g[t][h_idx[p])*(the_g[t][part]-c_g[t][h_idx[p])/resolution(pol_angle[t],2)/resolution(pol_angle[t],2)/2.)/sqrt(2*TMath::Pi()*resolution(pol_angle[t],2)*resolution(pol_angle[t],2));
	  dnorm = TMath::Gaus(c_g[t][h_idx[p]],the_g[t][part],resolution(pol_angle[t],2),0);
	  //pnot_g[part][t] = dnorm;
	  //if(dnorm>1) cout<<dnorm<<" g "<<t<<endl;
	  pois = ROOT::Math::poisson_cdf_c(nph_g[t]+1,n_g[t][part]);
	  if((dnorm*pois)>dmin){
	    dmin=dnorm*pois;
	    rad_bar = 2;
	    t_bar = t;
	    dmin_bar = dnorm;
	  }
	  //if(p==0)cout<<(the_a[t][part]-c_a[t][h_idx[p])*(the_a[t][part]-c_a[t][h_idx[p])<<" "<<(the_g[t][part]-c_g[t][h_idx[p])*(the_g[t][part]-c_g[t][h_idx[p])<<"  "<<part<<"  "<<t<<endl;
	  //cout<<dmin<<"  "<<p<<"  "<<part<<"  "<<c_a[t][h_idx[p]<<"  "<<c_g[t][h_idx[p]<<endl;
	  
	} // end track
	
	if(rad_bar==1){ nph_a[t_bar]++;  chsum_a[t_bar]+=c_a[t_bar][h_idx[p]]; }
	if(rad_bar==2){ nph_g[t_bar]++;  chsum_g[t_bar]+=c_g[t_bar][h_idx[p]]; }
	if(rad_bar==-1){ nph_b++; bg_ph->Fill(track_number);}
	if((rad_bar==1 || rad_bar==2) && t_id[t_bar] == m_tr[h_idx[p]]) {
	  nph_cor[t_bar]++;
	  good_ph->Fill(track_number);
	  tot_ph->Fill(track_number);
	  good_ph_m->Fill(pol_angle[t_bar]);
	  tot_ph_m->Fill(pol_angle[t_bar]);
	}
	if((rad_bar==1 || rad_bar==2) && t_id[t_bar] != m_tr[h_idx[p]]) {
	  nph_wro[t_bar]++;
	  tot_ph->Fill(track_number);
	  bad_ph->Fill(track_number);
	  bad_ph_m->Fill(pol_angle[t_bar]);
	  tot_ph_m->Fill(pol_angle[t_bar]);
	}
	Int_t p_bar = -1; 
	/*if(rad_bar==1 || rad_bar==2){
	  lh = lh + TMath::Log(dmin_bar);
	  p_bar =3&(c>>(2*t_bar));
	  }*/
	/*for(Int_t pa=0;pa<4;pa++){
	  if(pa!=p_bar) {
	    lh = lh +
	      TMath::Log(1.-TMath::Gaus(c_a[t_bar][h_idx[p],the_a[t_bar][pa],resolution(pol_angle[t_bar],1),0)) +
	      TMath::Log(1.-TMath::Gaus(c_g[t_bar][h_idx[p],the_g[t_bar][pa],resolution(pol_angle[t_bar],2),0));
	  }
	  else if(rad_bar==1)lh = lh + TMath::Log(1.-TMath::Gaus(c_g[t_bar][h_idx[p],the_g[t_bar][pa],resolution(pol_angle[t_bar],2),0));
	  else if(rad_bar==2)lh = lh + TMath::Log(1.-TMath::Gaus(c_a[t_bar][h_idx[p],the_a[t_bar][pa],resolution(pol_angle[t_bar],1),0));
	    
	}*/
      }  // end photons
      //cout<<" C: "<<(hex)<<c<<" L: "<<(dec)<<lh<<endl;

      Double_t np_ass = 0;
      lh = 0;
      
      for(Int_t t=0;t<track_number;t++){
	Int_t part=3&(c>>(2*t));
	if(nph_a[t]>0){
	  lh += TMath::Log(TMath::Gaus(chsum_a[t]/nph_a[t],the_a[t][part],resolution(pol_angle[t],1)/sqrt(nph_a[t]),0));
	}
	if(nph_g[t]>0){
	  lh += TMath::Log(TMath::Gaus(chsum_g[t]/nph_g[t],the_g[t][part],resolution(pol_angle[t],2)/sqrt(nph_g[t]),0));
	}
	
	lh += TMath::Log(TMath::Poisson(nph_a[t],n_a[t][part])) + TMath::Log(TMath::Poisson(nph_g[t],n_g[t][part]));

	lh += TMath::Log(TMath::Poisson(nph_b,1))+nph_b*TMath::Log(0.003);
	
	np_ass+=nph_a[t]+nph_g[t];
	//cout<<nph_a[t]<<"  "<<nph_g[t]<<"  "<<t<<endl;
	//cout<<part<<" "<<p_id[t]<<"  "<<mom_p[t]<<"  "<< nph_a_bar[t]<<"  "<< nph_g_bar[t]<<"  "<<n_a[t][part]<<"  "<<n_g[t][part]<<"  "<<n_a[t][p_id[t]]<<"  "<<n_g[t][p_id[t]]<<"  "<<t<<endl;
	}

      lh+=TMath::Log(TMath::Gaus(ph_count-np_ass,1.,1.,0)); //probabilita' dei fotoni non asseganti (fondo)

      //cout<<"lh:  "<<lh<<" lh_bar: "<<lh_bar<<" ph_count:  "<<ph_count<<" np_ass:  "<<np_ass<<" ++  "<<TMath::Log(TMath::Gaus(ph_count-np_ass,1.,1.,0))<<endl;
      
      if(lh>lh_bar){
	lh_bar = lh;
	c_bar = c;
	//cout<<" C_bar: "<<(hex)<<c_bar<<" L_bar: "<<(dec)<<lh_bar<<"  "<<np_ass<<endl;
	for(Int_t t=0;t<track_number;t++){
	  nph_a_bar[t]=nph_a[t];
	  nph_g_bar[t]=nph_g[t];
	}
	//cout<<c<<"  "<<evnt<<endl;
      }
 
    } //end of c combination loop 

    

    
    //cout<<" C_bar: "<<(hex)<<c_bar<<" L_bar: "<<(dec)<<lh_bar<<endl;

    Int_t pid_type = 1; Int_t flag_p = 0;

    for(Int_t t=0;t<track_number;t++){
      Int_t part=3&(c_bar>>(2*t));

      //if(p_id[t]==pid_type && part!=0) {flag_p = 1;
      //cout<<"Nph: "<<ph_count<<" Nevent: "<<evnt<<" Ntracks: "<<track_number<<"################"<<endl;
      //break;}
      
    }

    Int_t bad=0;
    for(Int_t t=0;t<track_number;t++){
      Int_t part=3&(c_bar>>(2*t));

      if(p_id[t]==0)p_id_tot[0]->Fill(mom_p[t]);
      if(p_id[t]==1)p_id_tot[1]->Fill(mom_p[t]);
      if(p_id[t]==2)p_id_tot[2]->Fill(mom_p[t]);
      if(p_id[t]==3 && mom_p[t]>5.)p_id_tot[3]->Fill(mom_p[t]);
      
      if(p_id[t]==0 && part==0)p_id_h[0]->Fill(mom_p[t]);
      if(p_id[t]==0 && part==1)p_id_h[1]->Fill(mom_p[t]);
      if(p_id[t]==0 && part==2)p_id_h[2]->Fill(mom_p[t]);
      if(p_id[t]==0 && part==3)p_id_h[3]->Fill(mom_p[t]);
      if(p_id[t]==1 && part==0)p_id_h[4]->Fill(mom_p[t]);
      if(p_id[t]==1 && part==1)p_id_h[5]->Fill(mom_p[t]);
      if(p_id[t]==1 && part==2)p_id_h[6]->Fill(mom_p[t]);
      if(p_id[t]==1 && part==3)p_id_h[7]->Fill(mom_p[t]);
      if(p_id[t]==2 && part==0)p_id_h[8]->Fill(mom_p[t]);
      if(p_id[t]==2 && part==1)p_id_h[9]->Fill(mom_p[t]);
      if(p_id[t]==2 && part==2)p_id_h[10]->Fill(mom_p[t]);
      if(p_id[t]==2 && part==3)p_id_h[11]->Fill(mom_p[t]);
      if(p_id[t]==3 && part==0 && mom_p[t]>5.)p_id_h[12]->Fill(mom_p[t]);
      if(p_id[t]==3 && part==1 && mom_p[t]>5.)p_id_h[13]->Fill(mom_p[t]);
      if(p_id[t]==3 && part==2 && mom_p[t]>5.)p_id_h[14]->Fill(mom_p[t]);
      if(p_id[t]==3 && part==3 && mom_p[t]>5.)p_id_h[15]->Fill(mom_p[t]);

      if(p_id[t]==0 && part!=0) bad++;
      if(p_id[t]==1 && part!=1) bad++;
      if(p_id[t]==2 && part!=2) bad++;
      if(p_id[t]==3 && part!=3) bad++;
  
      if(flag_p == 1){cout<<part<<" "<<p_id[t]<<"  "<<mom_p[t]<<"  "<< nph_a_bar[t]<<"  "<< nph_g_bar[t]<<"  "<<n_a[t][part]<<"  "<<n_g[t][part]<<"  "<<n_a[t][p_id[t]]<<"  "<<n_g[t][p_id[t]]<<"  "<<t<<endl;
	cout<<" C_bar: "<<(hex)<<c_bar<<" L_bar: "<<(dec)<<lh_bar<<endl;
      }
      //if (part != p_id[t]) cout <<" ****************** ";
      //cout <<endl;
    }

    //cout<<bad<<"  "<<track_number<<endl;
    h2_bad->Fill(track_number,bad);
    h_tot->Fill(track_number);

    bad = 0;
    
  } //end loop on the events

  Int_t NC = 16;
  Double_t cont[NC][100];
  Double_t p[NC][100], perr[NC][100];
  Double_t momentum[100];

  for(Int_t k=0;k<10;k++){
    momentum[k] = k*5 + 5.5;
    for(Int_t y=0;y<NC;y++)p[y][k]=0;
    for(Int_t y=0;y<NC;y++){
      cont[y][k] = p_id_h[y]->GetBinContent(p_id_h[0]->FindBin(k+1));
      //cout<<"Counts "<<y<<"  "<<k+1<<"  "<<cont[y][k]<<endl;
    }
    for(Int_t y=0;y<4;y++){
    if(cont[0+4*y][k]!=0)p[0+4*y][k] = cont[0+4*y][k]/(cont[0+4*y][k]+cont[1+4*y][k]+cont[2+4*y][k]+cont[3+4*y][k]);
    if(cont[1+4*y][k]!=0)p[1+4*y][k] = cont[1+4*y][k]/(cont[0+4*y][k]+cont[1+4*y][k]+cont[2+4*y][k]+cont[3+4*y][k]);
    if(cont[2+4*y][k]!=0)p[2+4*y][k] = cont[2+4*y][k]/(cont[0+4*y][k]+cont[1+4*y][k]+cont[2+4*y][k]+cont[3+4*y][k]);
    if(cont[3+4*y][k]!=0)p[3+4*y][k] = cont[3+4*y][k]/(cont[0+4*y][k]+cont[1+4*y][k]+cont[2+4*y][k]+cont[3+4*y][k]);

    if(cont[0+4*y][k]!=0)perr[0+4*y][k] = 1./cont[0+4*y][k];
    if(cont[1+4*y][k]!=0)perr[1+4*y][k] = 1./cont[1+4*y][k];
    if(cont[2+4*y][k]!=0)perr[2+4*y][k] = 1./cont[2+4*y][k];
    if(cont[3+4*y][k]!=0)perr[3+4*y][k] = 1./cont[3+4*y][k];
    }
  }

  TGraphErrors **gre;
  gre=new TGraphErrors*[NC];
  for(Int_t y=0;y<NC;y++){
    gre[y] = new TGraphErrors(10,momentum,p[y],0,0);
  }


  TCanvas *cgr = new TCanvas("cgr","",700,700);
  cgr->Divide(4,4);
  for(Int_t y=0;y<NC;y++){
    cgr->cd(y+1);
    /*gre[y]->SetTitle("");
    gre[y]->GetXaxis()->SetRangeUser(3,55);
    gre[y]->GetXaxis()->SetTitle("Momentum [GeV/c]");
    gre[y]->GetYaxis()->SetRangeUser(0,1);
    gre[y]->GetYaxis()->SetTitle("");
    gre[y]->SetFillColor(38);
    gre[y]->Draw("AB");*/

    p_id_h[y]->SetTitle("");
    p_id_h[y]->GetXaxis()->SetRangeUser(3,55);
    p_id_h[y]->GetXaxis()->SetTitle("Momentum [GeV/c]");
    if(y>=1 && y<=3)p_id_h[y]->GetYaxis()->SetRangeUser(0,0.1);
    if(y==6 || y==7)p_id_h[y]->GetYaxis()->SetRangeUser(0,0.1);
    if(y==4)p_id_h[y]->GetYaxis()->SetRangeUser(0,0.4);
    if(y==11)p_id_h[y]->GetYaxis()->SetRangeUser(0,0.2);
    if(y==12 || y==13 || y==14 || y==8 || y==9)p_id_h[y]->GetYaxis()->SetRangeUser(0,0.1);
    p_id_h[y]->GetYaxis()->SetTitle("");
    p_id_h[y]->SetFillColor(38);
    if(y>=0 && y<=3)p_id_h[y]->Divide(p_id_tot[0]);
    if(y>=4 && y<=7)p_id_h[y]->Divide(p_id_tot[1]);
    if(y>=8 && y<=11)p_id_h[y]->Divide(p_id_tot[2]);
    if(y>=12 && y<=15)p_id_h[y]->Divide(p_id_tot[3]);
    p_id_h[y]->Draw();
  }
  cgr->SaveAs("figs/pid_table.png");
  
  /*for(Int_t y=0;y<4;y++){
    cgr->cd(1+y);
    gre[0+4*y]->GetXaxis()->SetRangeUser(2,50);
    gre[0+4*y]->GetYaxis()->SetRangeUser(0,1);
    gre[0+4*y]->SetMarkerColor(1);
    gre[1+4*y]->SetMarkerColor(2);
    gre[2+4*y]->SetMarkerColor(3);
    gre[3+4*y]->SetMarkerColor(4);
    gre[0+4*y]->SetMarkerStyle(21);
    gre[1+4*y]->SetMarkerStyle(21);
    gre[2+4*y]->SetMarkerStyle(21);
    gre[3+4*y]->SetMarkerStyle(21);
    gre[0+4*y]->Draw("AP");
    gre[1+4*y]->Draw("P");
    gre[2+4*y]->Draw("P");
    gre[3+4*y]->Draw("P");
    }*/
  TCanvas *cbad = new TCanvas("cbad","",700,700);
  //h2_bad->Scale(1./h2_bad->GetEntries());
  for(Int_t i=0;i<10;i++){
    //cout<<h_tot->GetBinContent(i+1)<<"  "<<i+1<<endl;
    for(Int_t j=0;j<10;j++){
      if(h2_bad->GetBinContent(i+1,j+1)!=0){
	Double_t bc = h2_bad->GetBinContent(i+1,j+1);
	h2_bad->SetBinContent(i+1,j+1,bc/h_tot->GetBinContent(i+1));
      }
    }
      //cout<<h2_bad->GetBinContent(i+1,j+1)<<"  "<<i+1<<"  "<<j+1<<endl;
  }
  h2_bad->GetZaxis()->SetRangeUser(0., 1.);
  h2_bad->Draw("colz");
  cbad->SaveAs("figs/eff_table.png");
  
  TCanvas *cph = new TCanvas("cph","",700,700);
  //bad_ph->Scale(1./(bad_ph->GetEntries()+good_ph->GetEntries()));
  bad_ph->Divide(tot_ph);
  bad_ph->Draw("bar");
  cph->SaveAs("figs/bad_ph.png");

  TCanvas *cph_bg = new TCanvas("cph_bg","",700,700);
  bg_ph->Scale(1./(bad_ph->GetEntries()+good_ph->GetEntries()+bg_ph->GetEntries()));
  bg_ph->Draw("bar");
  cph_bg->SaveAs("figs/bg_ph.png");

  TCanvas *cph_m = new TCanvas("cph_m","",700,700);
  cph_m->Divide(1,3);
  cph_m->cd(1);
  tot_ph_m->Draw("");
  cph_m->cd(2);
  bad_ph_m->Draw("");
  cph_m->cd(3);
  bad_ph_m->Divide(tot_ph_m);
  bad_ph_m->Draw("bar");

  TCanvas *pa = new TCanvas("pa","",900,900);
  pa->Divide(2,2);
  pa->cd(1);
  ch_part[0]->Draw();
  pa->cd(2);
  ch_part[1]->Draw();
  pa->cd(3);
  ch_part[2]->Draw();
  pa->cd(4);
  ch_part[3]->Draw();
  for(Int_t k=0;k<4;k++){
    ch_part[k]->GetXaxis()->SetTitle("Momentum [GeV/c]");
    pa->cd(k+1)->SetLogy();
  }
  pa->SaveAs("figs/phase_space.png");

  TCanvas *molt = new TCanvas("molt","",500,500);
  molt_h->Draw("");
  molt->SaveAs("figs/mult.png");
  
} //end of class member acq_pnew


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

  /*TFile *file=new TFile(input_filename.c_str());
  if (file->IsZombie()) {
    cout << "Error opening file" << input_filename << endl;
    exit(-1);
  }
  else cout << "open file " << input_filename << endl;*/

  //TTree *generated = (TTree*) file->Get("generated");
  TChain *generated = new TChain("generated");
  generated->Add(input_filename.c_str());
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

  //TTree *flux = (TTree*) file->Get("flux");
  TChain *flux = new TChain("flux");
  flux->Add(input_filename.c_str());

  vector<double> *flux_id=0,*flux_hitn=0;
  vector<double> *flux_pid=0,*flux_mpid=0,*flux_tid=0,*flux_mtid=0,*flux_otid=0;
  vector<double> *flux_trackE=0,*flux_totEdep=0,*flux_avg_x=0,*flux_avg_y=0,*flux_avg_z=0,*flux_avg_lx=0,*flux_avg_ly=0,*flux_avg_lz=0,*flux_px=0,*flux_py=0,*flux_pz=0,*flux_vx=0,*flux_vy=0,*flux_vz=0,*flux_mvx=0,*flux_mvy=0,*flux_mvz=0,*flux_avg_t=0;
  vector<double> *flux_procID=0,*flux_nsteps=0;
  flux->SetBranchAddress("hitn",&flux_hitn);  
  flux->SetBranchAddress("id",&flux_id);
  flux->SetBranchAddress("pid",&flux_pid);
  flux->SetBranchAddress("mpid",&flux_mpid);
  flux->SetBranchAddress("tid",&flux_tid);
  flux->SetBranchAddress("mtid",&flux_mtid);
  flux->SetBranchAddress("otid",&flux_otid);
  flux->SetBranchAddress("trackE",&flux_trackE);
  flux->SetBranchAddress("totEdep",&flux_totEdep);
  flux->SetBranchAddress("avg_x",&flux_avg_x);
  flux->SetBranchAddress("avg_y",&flux_avg_y);
  flux->SetBranchAddress("avg_z",&flux_avg_z);
  flux->SetBranchAddress("avg_lx",&flux_avg_lx);
  flux->SetBranchAddress("avg_ly",&flux_avg_ly);
  flux->SetBranchAddress("avg_lz",&flux_avg_lz);
  flux->SetBranchAddress("px",&flux_px);
  flux->SetBranchAddress("py",&flux_py);
  flux->SetBranchAddress("pz",&flux_pz);
  flux->SetBranchAddress("vx",&flux_vx);
  flux->SetBranchAddress("vy",&flux_vy);
  flux->SetBranchAddress("vz",&flux_vz);
  flux->SetBranchAddress("mvx",&flux_mvx);
  flux->SetBranchAddress("mvy",&flux_mvy);
  flux->SetBranchAddress("mvz",&flux_mvz);
  flux->SetBranchAddress("avg_t",&flux_avg_t);
  flux->SetBranchAddress("nsteps",&flux_nsteps);
  flux->SetBranchAddress("procID",&flux_procID);

  //flux->Draw("avg_z");

  //TTree *eic_rich = (TTree*) file->Get("eic_rich");
  TChain *eic_rich = new TChain("eic_rich");
  eic_rich->Add(input_filename.c_str());

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

  /***flux def***/

  Int_t fcount = 0;
  Double_t fmomv[NPART][3];
  
  for(Int_t ii=0;ii<NPART;ii++){

    fmomv[ii][0] = 0.;
    fmomv[ii][1] = 0.;
    fmomv[ii][2] = 0.;
    
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
  Int_t nsigma = 5;
  Double_t aidx = 1.019;
  Double_t gidx = 1.00086;
    

  
  for(Int_t i=0;(i<eic_rich->GetEntries());i++){
    
    generated->GetEntry(i);
    flux->GetEntry(i);
    eic_rich->GetEntry(i);
    
    evnt = i;
    
    TH2F *ph_det_s = new TH2F("ph_det","", 1666, -250., 250., 1666, -250., 250.);

    for(Int_t k=0;k<flux_id->size();k++){
      
      //cout<<flux_id->at(k)<<endl;

      if(i<=event_n && flux_pid->at(k)!=0 && flux_pid->at(k)!=22 && flux_pid->at(k)!=2112 && flux_avg_z->at(k)>0 && flux_vz->at(k)==0){

	fmomv[fcount][0] = flux_px->at(k)/sqrt(flux_px->at(k)*flux_px->at(k)+flux_py->at(k)*flux_py->at(k)+flux_pz->at(k)*flux_pz->at(k));
	fmomv[fcount][1] = flux_py->at(k)/sqrt(flux_px->at(k)*flux_px->at(k)+flux_py->at(k)*flux_py->at(k)+flux_pz->at(k)*flux_pz->at(k));
	fmomv[fcount][2] = flux_pz->at(k)/sqrt(flux_px->at(k)*flux_px->at(k)+flux_py->at(k)*flux_py->at(k)+flux_pz->at(k)*flux_pz->at(k));

	//cout<<" * "<<fmomv[fcount][0]<<" "<<fmomv[fcount][1]<<" "<<fmomv[fcount][2]<<endl;
	//cout<<TMath::ACos(fmomv[fcount][2])*57.2958<<endl;
	//cout<<"PID NUMBER*: "<<flux_pid->at(k)<<endl;

	if((TMath::ACos(fmomv[fcount][2])*57.2958)>3 && (TMath::ACos(fmomv[fcount][2])*57.2958)<25)fcount++;
	
      }
      
    }

    //if(i==1)cout<<fcount<<endl;
    
    for(Int_t j=0;j<eic_rich_id->size();j++){

      if(i<=event_n && eic_rich_pid->at(j)!=0 && eic_rich_pid->at(j)!=22 && eic_rich_pid->at(j)!=2112 && eic_rich_vz->at(j)<2500 && eic_rich_out_z->at(j)==2540) {
	
	if(eic_rich_id->at(j)<20){
	  //if(evnt==69)cout<<"PID NUMBER: "<<eic_rich_pid->at(j)<<endl;
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
	  if(evnt<5)cout<<"MOMENTUM: "<<sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j))/1000.<<endl;
	}

	mom_v[track_number][0] = eic_rich_in_px->at(j)/sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j));
	mom_v[track_number][1] = eic_rich_in_py->at(j)/sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j));
	mom_v[track_number][2] = eic_rich_in_pz->at(j)/sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j));

	//cout<<mom_v[track_number][0]<<" "<<mom_v[track_number][1]<<" "<<mom_v[track_number][2]<<endl;
	
	pol_angle[track_number] = TMath::ACos(mom_v[track_number][2])*57.2958;
	phi_angle[track_number] = TMath::ATan2(mom_v[track_number][1],mom_v[track_number][0])*57.2958;
	
	mom_p[track_number] = sqrt(eic_rich_in_px->at(j)*eic_rich_in_px->at(j)+eic_rich_in_py->at(j)*eic_rich_in_py->at(j)+eic_rich_in_pz->at(j)*eic_rich_in_pz->at(j))/1000.;

	if(pol_angle[track_number]>0. && pol_angle[track_number]<25.){
	  ch_part[p_id[track_number]]->Fill(mom_p[track_number]);
	}
	
	for(Int_t idp=0;idp<4;idp++){
	  the_a[track_number][idp]=acos(sqrt(mom_p[track_number]*mom_p[track_number]+m[idp]*m[idp])/mom_p[track_number]/aidx);
	  refer[idp] = 8./TMath::Sin(acos(sqrt(50.*50.+m[idp]*m[idp])/50./aidx))/TMath::Sin(acos(sqrt(50.*50.+m[idp]*m[idp])/50./aidx));
	  n_a[track_number][idp]=207.*sin(the_a[track_number][idp])*sin(the_a[track_number][idp]);

	  the_g[track_number][idp]=acos(sqrt(mom_p[track_number]*mom_p[track_number]+m[idp]*m[idp])/mom_p[track_number]/gidx);
	  refer[idp] = 20./TMath::Sin(acos(sqrt(50.*50.+m[idp]*m[idp])/50./gidx))/TMath::Sin(acos(sqrt(50.*50.+m[idp]*m[idp])/50./gidx));
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
	
	if(pol_angle[track_number]>0. && pol_angle[track_number]<25.){
	  track_number++;
	  molteplicity[i]++;
	  //track_number=2;
	  //molteplicity[i]=2;
	}
	//track_number++;
	//molteplicity[i]++;
      }
      if(i<=event_n && eic_rich_pid->at(j)!=0 && eic_rich_pid->at(j)!=22 && eic_rich_vz->at(j)<2500 && eic_rich_id->at(j)<30 && eic_rich_id->at(j)>20 &&  eic_rich_tid->at(j)==p_type[track_number-1]) {
       	emig_x[track_number-1] = eic_rich_avg_x->at(j)*0.1;
	emig_y[track_number-1] = eic_rich_avg_y->at(j)*0.1;
	emig_z[track_number-1] = eic_rich_avg_z->at(j)*0.1;
      }
      
      //if(i<=event_n && eic_rich_pid->at(j)==0 && eic_rich_id->at(j)<50 && eic_rich_id->at(j)>40 && (eic_rich_mtid->at(j)==p_type[0] || eic_rich_mtid->at(j)==p_type[1]) && pol_angle[0]>5. && pol_angle[0]<25. && p_id[0]==pid_type && p_id[1]==pid_type){ //********************uncomment if multiplicity 1
	if(i<=event_n && eic_rich_pid->at(j)==0 && eic_rich_id->at(j)<50 && eic_rich_id->at(j)>40){ //to be c ommented if multiplicity 1
	ph_energy = eic_rich_trackE->at(j)*1000000.;
	
	for(Int_t k=0;k<phd_count-1;k++){

	  if((eic_rich_trackE->at(j)*1000000.)>=ev_a[k] && (eic_rich_trackE->at(j)*1000000.)<ev_a[k+1]){
	    qe_p = qe_a[k];
	    //qe_p = 1;
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
	  
	  if(evnt==69)ph_det->Fill(eic_rich_in_x->at(j)*0.1,eic_rich_in_y->at(j)*0.1);

	  xbin = ph_det->GetXaxis()->FindBin(eic_rich_in_x->at(j)*0.1);
	  ybin = ph_det->GetYaxis()->FindBin(eic_rich_in_y->at(j)*0.1);
	  
	  if(ph_det_s->GetBinContent(xbin,ybin) == 1){
	    for(Int_t i=0;i<track_number;i++){

	      c_angle_a = ind_ray1(emi_x[i],emi_y[i],emi_z[i],ph_det_s->GetXaxis()->GetBinCenter(xbin),ph_det_s->GetYaxis()->GetBinCenter(ybin),det_z,mir_x,mir_y,mir_z,mom_v[i][0],mom_v[i][1],mom_v[i][2],1);
	      //cout<<mom_p[0]<<"  "<<det_x<<"  "<<det_y<<"  "<<c_angle_a<<"  "<<evnt<<endl;
	      c_angle = ind_ray1(emig_x[i],emig_y[i],emig_z[i],ph_det_s->GetXaxis()->GetBinCenter(xbin),ph_det_s->GetYaxis()->GetBinCenter(ybin),det_z,mir_x,mir_y,mir_z,mom_v[i][0],mom_v[i][1],mom_v[i][2],2);
	      if(eic_rich_mtid->at(j)==p_type[i] && evnt==69)ch_h[i]->Fill(c_angle_a);
	      else if(evnt==69) chb_h[i]->Fill(c_angle_a);
	      
	      for(Int_t idp=0;idp<4;idp++){
		if(the_a[i][idp]>0. && c_angle_a<(the_a[i][idp]+resolution(pol_angle[i],1)*nsigma) && c_angle_a>(the_a[i][idp]-resolution(pol_angle[i],1)*nsigma)){
		  //if(the_a[i][idp]>0. && c_angle_a<(the_a[i][idp]+0.012) && c_angle_a>(the_a[i][idp]-0.012)){
		  tch_dist[i][idp] = tch_dist[i][idp] + c_angle_a;
		  a_counts[i][idp]++;
		  //if(i==1)cout<<c_angle<<"  "<<the_g[i][idp]+0.006<<endl;
		}
		if(the_g[i][idp]>0. && c_angle<(the_g[i][idp]+resolution(pol_angle[i],2)*nsigma) && c_angle>(the_g[i][idp]-resolution(pol_angle[i],2)*nsigma)){
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
	//if(p_id[0]==pid_type && mom_p[0]<5.)cout<<"AEROGEL"<<" Sector: "<<sector_id[j]<<"  "<<pol_angle[j]<<endl;
	//if(p_id[0]==pid_type && mom_p[0]<5.)cout<<tch_dist[j][idp]/a_counts[j][idp]<<"  "<<a_counts[j][idp]<<"  "<<the_a[j][idp]<<"  "<<evnt<<endl;
	//if(p_id[0]==pid_type && mom_p[0]<5.)cout<<a_like[j][idp]<<" p_type  "<<p_id[j]<<" Nph:  "<<n_a[j][idp]<<" Mom: "<<mom_p[j]<<endl;
	//if(j==0 && evnt==366)cout<<"GAS"<<endl;
	//if(j==0 && evnt==366)cout<<tch_dist1[j][idp]/g_counts[j][idp]<<"  "<<g_counts[j][idp]<<"  "<<the_g[j][idp]<<"  "<<mom_p[j]<<endl;
	//if(j==0 && evnt==366)cout<<likel_g<<"  "<<idp<<endl;
	//if(j==2)cout<<"Sum: "<<a_like[j][idp]+g_like[j][idp]<<endl;
      }
      //cout<<"PID+++++++++++++++++++++++++++++++++++++++++"<<endl;
      idx0=0;
      for(Int_t idp=1;idp<4;idp++){
	if(a_like[j][0]==999 && a_like[j][1]==999 && a_like[j][2]==999 && a_like[j][3]==999){
	  idx0 = 99;
	  break;
	}
	else if((a_like[j][idp]!=999 && g_like[j][idp]==999) || (a_like[j][idp]!=999 && g_counts[j][idp]<9)){
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
      if(mom_p[j]<=50 && mom_p[j]>=3. && pol_angle[j]>5. && pol_angle[j]<25.){ // Here PID evaluation .........................................
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
    fcount = 0;

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
  //ph_det->Draw("");
  molt_h->Draw("");

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

 
