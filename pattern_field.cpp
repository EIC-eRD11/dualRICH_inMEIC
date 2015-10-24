//#######Calculate the magnetic field bending effect for the dual RICH in MEIC concept#########
//Code version 1.1
//Author Alessio Del Dotto 

//To compile and run in root framework: 
//.L pattern_field.cpp++
//magnetic_field *a = new magnetic_field()
//a->acq() 


#include <Riostream.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>

#include <TMath.h>
#include <TEnv.h>
#include <TH2F.h>
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
#include <vector>
#include <string>
#include <TGraphErrors.h>
//#include <pattern_field.h>

class magnetic_field {

 private:

 public:
  void acq();

};

//*******magnetic_field class acq*********//

void magnetic_field::acq(){

  gStyle->SetCanvasColor(0);
  //gStyle->SetTitleW(0.4);
  //gStyle->SetTitleH(0.095);
  gStyle->SetTitleSize(0.05,"xyzt");
  gStyle->SetTitleXOffset(1);
  gStyle->SetTitleYOffset(1.4);
  gStyle->SetLabelSize(0.03,"XYZ");
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

  Double_t theta, omega, p_z, p_perp;
  Double_t theta_min, theta_max, eta_min, eta_max;

  Double_t s, sv, xx, x[500], dphi[500], dphi1[500], dphis[500], dphis1[500], dphi_mean, delta_ring;

  Double_t dphi_rms[100], mag_mean[100], mag_int[100], eta_v[100], p_v[100], pk_v[100], ppr_v[100];
  Double_t thetac[100], thetac_p[100], thetac_m[100], thetack[100], thetacpr[100];
  Int_t eta_count = 0;
  Int_t p_count = 0;
  Int_t pk_count = 0;
  Int_t ppr_count = 0;

  Double_t h[3], th[3];
  Double_t xs[3][500];
  xs[0][0]=0.;
  xs[1][0]=0.;
  xs[2][0]=0.;
  Double_t t[3][500];
  Double_t turn_ang, rho, thetai, lambdai;
  Double_t ang_diff, phi, ang_r, dphii;

  Double_t r,l,B_r,B_z, theta_B;
  Double_t B_perp[500], zz[500], B_par[500];
  Double_t p_perpp[500], p_par[500];
  Int_t y,z;
  ifstream inputFile;
  Double_t p = 15.;
  //Double_t eta = 1.;
  //Double_t p = 31.2429;
  Int_t start_v = 220;

  
  for(Double_t eta=1.;eta<=4.;eta=(eta+0.5)){
    //for(Double_t p=1.;p<=50.;p=(p+0.5)){

      Int_t count = 0;

      theta_min = TMath::ATan(25./416.);
      theta_max = TMath::ATan(225./416.);
      
      eta_min = -TMath::Log(TMath::Tan(theta_max/2.));
      eta_max = -TMath::Log(TMath::Tan(theta_min/2.));

      cout<<eta<<endl;

      theta = 2*TMath::ATan(exp(-eta));
      //theta = 0.41627;
      omega = TMath::Pi()/2.-theta;
      p_perp = p*TMath::Cos(omega);
      p_z = p*TMath::Sin(omega);
      s = 0.01/TMath::Sin(omega);
      xx = 0.;
      x[0]=0.;
      dphi[0]=0.;
      dphi1[0]=0.;
      dphis[0]=0.;
      dphis1[0]=0.;

      cout<<"theta:  "<<theta*180./TMath::Pi()<<endl;
      cout<<"eta_min:  "<<eta_min<<endl;
      cout<<"eta_max:  "<<eta_max<<endl;

      TH1F *t1_dis = new TH1F("#theta","",1000,-50.,50.);
      TH1F *t2_dis = new TH1F("#Delta#phi","",1000,-50.,50.);
      TH1F *ang_r_dis = new TH1F("#lambda","",1000,1.47354-0.01,1.47354+0.01);
      TH1F *phi_dis = new TH1F("#phi","",1000,-0.1,0.1);

      TH1F *mag_per = new TH1F("","",10000,0.,10000.);

      t1_dis->SetFillColor(4);
      t2_dis->SetFillColor(4);
      ang_r_dis->SetFillColor(4);
      phi_dis->SetFillColor(4);

      mag_per->SetFillColor(4);

      inputFile.open("meic_det1_solenoid_dual_v3_inter.dat");
      if (!inputFile) {
	cerr << "Error in opening imput file";
	exit(1);
      }
    
      while (inputFile>>r>>l>>B_r>>B_z) {
      
	z=l;
	y=l*TMath::Tan(theta);

	//if(z>=220. && z<=385.) B_r = B_r*0.1;
	//if(z>=220. && z<=385.) B_z = B_z*0.1;
	
	if(l>=0. && l<500. && r==y){ //prima era 261.  416.

	  //versors systems 

	  h[0] = B_r/sqrt(B_r*B_r+B_z*B_z);
	  h[1] = 0.;
	  h[2] = B_z/sqrt(B_r*B_r+B_z*B_z);

	  t[0][0] = TMath::Cos(omega);
	  t[1][0] = 0;
	  t[2][0] = TMath::Sin(omega);

	  //cout<<t[2][z]<<"  "<<h[2]<<"  "<<z<<endl;

	  B_par[z] = B_z*t[2][z]+B_r*t[0][z];
	  //B_perp[z] = sqrt(B_r*B_r+B_z*B_z-B_par[z]*B_par[z]);
	  theta_B =  TMath::ACos(B_par[z]/sqrt(B_r*B_r+B_z*B_z));
	  B_perp[z] = TMath::Abs(sqrt(B_r*B_r+B_z*B_z)*TMath::Sin(theta_B));
	  zz[z]=z;
	  p_par[z] = p*t[2][z]*h[2]+p*t[0][z]*h[0];
	  p_perpp[z] = sqrt(p*p-p_par[z]*p_par[z]);

	  if(z>220. && z<300.) mag_per->Fill(TMath::Abs(B_perp[z]));

	  //cout<<p_par[z]<<"  "<<p_perpp[z]<<"  "<<B_par[z]<<"  "<<B_perp[z]<<z<<endl;
	  //cout<<B_z<<"  "<<B_r<<"  "<<z<<"  "<<l*TMath::Tan(theta)<<endl;

	  thetai = TMath::ACos(h[0]*t[0][z]+h[1]*t[1][z]+h[2]*t[2][z]);
	  lambdai = TMath::Pi()/2.-thetai;;

	  th[0] = t[1][z]*h[2]-t[2][z]*h[1];
	  th[1] = t[2][z]*h[0]-t[0][z]*h[2];
	  th[2] = t[0][z]*h[1]-t[1][z]*h[0];

	  sv = 0.01/t[2][z];
	  rho = p*TMath::Cos(lambdai)*10000./0.3/sqrt(B_r*B_r+B_z*B_z); 
	  turn_ang = sqrt(B_r*B_r+B_z*B_z)*0.3*s/p/10000.;


	  xs[0][z+1]=xs[0][z]+(rho/TMath::Cos(lambdai))*((1.-TMath::Cos(turn_ang))*th[0]+TMath::Sin(turn_ang)*t[0][z]+(turn_ang-TMath::Sin(turn_ang))*TMath::Sin(lambdai)*h[0]);
	  xs[1][z+1]=xs[1][z]+(rho/TMath::Cos(lambdai))*((1.-TMath::Cos(turn_ang))*th[1]+TMath::Sin(turn_ang)*t[1][z]+(turn_ang-TMath::Sin(turn_ang))*TMath::Sin(lambdai)*h[1]);
	  xs[2][z+1]=xs[2][z]+(rho/TMath::Cos(lambdai))*((1.-TMath::Cos(turn_ang))*th[2]+TMath::Sin(turn_ang)*t[2][z]+(turn_ang-TMath::Sin(turn_ang))*TMath::Sin(lambdai)*h[2]);

	  t[0][z+1] = TMath::Sin(turn_ang)*th[0]+TMath::Cos(turn_ang)*t[0][z]+(1.-TMath::Cos(turn_ang))*TMath::Sin(lambdai)*h[0];
	  t[1][z+1] = TMath::Sin(turn_ang)*th[1]+TMath::Cos(turn_ang)*t[1][z]+(1.-TMath::Cos(turn_ang))*TMath::Sin(lambdai)*h[1];
	  t[2][z+1] = TMath::Sin(turn_ang)*th[2]+TMath::Cos(turn_ang)*t[2][z]+(1.-TMath::Cos(turn_ang))*TMath::Sin(lambdai)*h[2];


	  if(l>=220. && l<=300.){

	    if(z>=0.){ ang_diff = TMath::ACos(t[0][start_v]*t[0][z]+t[1][start_v]*t[1][z]+t[2][start_v]*t[2][z]);
	    phi = TMath::ATan(t[1][z]/t[0][z]);
	    ang_r = TMath::ASin(t[2][z]);
	    dphii = (TMath::ATan(t[1][z]/t[0][z])-TMath::ATan(t[1][start_v]/t[0][start_v]));
	    
	    }
	    //cout<<z<<"  "<<phi<<endl;
	    
	    if(z>=0.) t1_dis->Fill(ang_diff*1000.);
	    if(z>=0.) t2_dis->Fill(dphii*1000.);
	    ang_r_dis->Fill(ang_r);
	    phi_dis->Fill(phi);

      
	  }
	  
	  count++;
	}
      }

      TGraph *graph_mag_abs = new TGraph(count,zz,B_perp);

      TF1 f1("f",[&](double *zz, double *){ return graph_mag_abs->Eval(zz[0]); },0,500,0);
      Double_t integral = f1.Integral(220,385);

      //eta_v[eta_count] = eta;
      eta_v[eta_count] = theta*180./TMath::Pi();
      //eta_v[eta_count] = p;      //case of for on p
      dphi_rms[eta_count] = t1_dis->GetRMS();
      mag_mean[eta_count] = mag_per->GetMean();
      
      mag_int[eta_count] = integral;

      Double_t m[4]={0.000511, 0.13957018, 0.493677, 0.938272};
      //Double_t ng = 1.00137; //C4F10
      Double_t ng = 1.000482;  //CF4
      Double_t p_th[4]; 
      Double_t dR;
      Double_t Nph = 15.;  //C4F10 about 24  

      dR = dphi_rms[eta_count]/sqrt(2*Nph)*10./p;

      p_th[1] = m[1]/sqrt(ng*ng-1);
      p_th[2] = m[2]/sqrt(ng*ng-1);
      p_th[3] = m[3]/sqrt(ng*ng-1);

      if(p>p_th[1]) {
	thetac[p_count]=acos(sqrt(p*p+m[1]*m[1])/p/ng)*1000.;
	thetac_m[p_count]=acos(sqrt(p*p+m[1]*m[1])/p/ng)*1000.-dR;
	thetac_p[p_count]=acos(sqrt(p*p+m[1]*m[1])/p/ng)*1000.+dR;
	p_v[p_count] = p;

	//cout<<thetac[p_count]<<"  "<<thetac_m[p_count]<<"  "<<thetac_p[p_count]<<endl; 

	p_count++;
      }

      if(p>p_th[2]) {
	thetack[pk_count]=acos(sqrt(p*p+m[2]*m[2])/p/ng)*1000.;
	pk_v[pk_count] = p; 

	pk_count++;
      }

      if(p>p_th[3]) {
	thetacpr[ppr_count]=acos(sqrt(p*p+m[3]*m[3])/p/ng)*1000.;
	ppr_v[ppr_count] = p; 

	ppr_count++;
      }
      
      cout<<eta_v[eta_count]<<"  "<<dphi_rms[eta_count]<<"  "<<eta_count<<endl; 


      if(eta<4.) delete t1_dis;
      if(eta<4.) delete t2_dis;
      if(eta<4.) delete ang_r_dis;
      if(eta<4.) delete phi_dis;
      if(eta<4.) delete mag_per;

      /*if(p<50.) delete t1_dis;             //case of for on p
      if(p<50.) delete t2_dis;
      if(p<50.) delete ang_r_dis;
      if(p<50.) delete phi_dis;*/
    
      eta_count++;
      
      inputFile.close();

      if(eta == 4.){
	//if(p == 50.){

	TCanvas *cang = new TCanvas("cang","cang",1000,500);
	cang->Divide(2,1);
	cang->cd(1);
	ang_r_dis->Draw("BAR1");
	ang_r_dis->GetXaxis()->SetTitle("#lambda [rad]");
	cang->cd(2);
	phi_dis->Draw("BAR1");
	phi_dis->GetXaxis()->SetTitle("#phi [rad]");

	cang->Update();

	TCanvas *cang1 = new TCanvas("cang1","cang1",1000,1000);
	cang1->Divide(1,1);
	cang1->cd(1);
	t1_dis->Draw("BAR1");
	t1_dis->SetTitle("p = 30 GeV/c, #eta = 3");
	t1_dis->GetXaxis()->SetTitle("#theta [mrad]");
	//t2_dis->Draw("BAR1");
	//t2_dis->GetXaxis()->SetTitle("#Delta#phi [mrad]");

	cang1->Update();

	TGraph *graph_mag = new TGraph(count,zz,B_perp);
	graph_mag->SetMarkerStyle(20);
	graph_mag->SetMarkerColor(1);
	graph_mag->SetLineStyle(1);
	graph_mag->SetLineWidth(1);
	graph_mag->SetTitle("B#perp to #vec{p}");
	graph_mag->GetXaxis()->SetTitle("z [cm]");
	graph_mag->GetYaxis()->SetTitle("B_{#perp} [G]");
	TGraph *graph_p = new TGraph(count,zz,p_perpp);
	graph_p->SetMarkerStyle(20);
	graph_p->SetMarkerColor(1);
	graph_p->SetLineStyle(1);
	graph_p->SetLineWidth(1);
	graph_p->SetTitle("p#perp to #vec{B}");
	graph_p->GetXaxis()->SetTitle("z [cm]");
	graph_p->GetYaxis()->SetTitle("|p_{#perp}| [GeV]");
	TGraph *graph_mag1 = new TGraph(count,zz,B_par);
	graph_mag1->SetMarkerStyle(20);
	graph_mag1->SetMarkerColor(1);
	graph_mag1->SetLineStyle(1);
	graph_mag1->SetLineWidth(1);
	graph_mag1->SetTitle("B#parallel to #vec{p}");
	graph_mag1->GetXaxis()->SetTitle("z [cm]");
	graph_mag1->GetYaxis()->SetTitle("B_{#parallel} [G]");
	TGraph *graph_p1 = new TGraph(count,zz,p_par);
	graph_p1->SetMarkerStyle(20);
	graph_p1->SetMarkerColor(1);
	graph_p1->SetLineStyle(1);
	graph_p1->SetLineWidth(1);
	graph_p1->SetTitle("p#parallel to #vec{B}");
	graph_p1->GetXaxis()->SetTitle("z [cm]");
	graph_p1->GetYaxis()->SetTitle("|p_{#parallel}| [GeV]");;

	TCanvas *cgra = new TCanvas("cgra","cgra",500,1000);
	cgra->Divide(1,2);
	cgra->cd(1);
	graph_mag->Draw("AC");
	//cgra->cd(3);
	//graph_p->Draw("AC");
	cgra->cd(2);
	graph_mag1->Draw("AC");
	//cgra->cd(4);
	//graph_p1->Draw("AC");

	cgra->Update();

	TCanvas *ch = new TCanvas("ch","ch",500,500);
	ch->Divide(1,1);
	ch->cd(1);
	mag_per->Draw();
	ch->Update();

      }      
  } //end for eta

  TGraphErrors *graph = new TGraphErrors(eta_count,eta_v,dphi_rms);
  TGraphErrors *graph1 = new TGraphErrors(eta_count,eta_v,mag_mean);
  TGraphErrors *graph2 = new TGraphErrors(eta_count,eta_v,mag_int);

  /*TGraphErrors *graph_pi = new TGraphErrors(p_count,p_v,thetac);       //for plot vs p
  TGraphErrors *graph_pim = new TGraphErrors(p_count,p_v,thetac_m);
  TGraphErrors *graph_pip = new TGraphErrors(p_count,p_v,thetac_p);

  TGraphErrors *graph_k = new TGraphErrors(pk_count,pk_v,thetack);
  TGraphErrors *graph_pr = new TGraphErrors(ppr_count,ppr_v,thetacpr);*/

  //TLegend *leg = new TLegend(0.2,0.7,0.4,0.89);
  //leg->SetFillStyle(0);
  //leg->SetFillColor(0);
  //leg->SetShadowColor(0);
  //leg->AddEntry(graph,"","l");

  TCanvas *gr = new TCanvas("gr","gr",1500,500);

  gr->Divide(3,1);

  gr->cd(1);
  graph->Draw("AC");
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(1);
  graph->SetLineStyle(1);
  graph->SetLineWidth(3);
  //graph->SetTitle("#Delta#theta for p = 30 GeV/c");
  graph->SetTitle("");
  //graph->GetXaxis()->SetTitle("#eta");
  graph->GetXaxis()->SetTitle("#theta [deg]");
  graph->GetYaxis()->SetTitle("#Delta#theta [mrad]");
  //graph->GetYaxis()->SetTitle("|B|_{#perp} [G]");
  gr->Update();
  gr->cd(2);
  graph1->Draw("AC");
  graph1->SetMarkerStyle(20);
  graph1->SetMarkerColor(1);
  graph1->SetLineStyle(1);
  graph1->SetLineWidth(3);
  graph1->SetTitle("Mean |B|_{#perp} to #vec{p}");
  //graph1->GetXaxis()->SetTitle("#eta");
  graph1->GetXaxis()->SetTitle("#theta [deg]");
  graph1->GetYaxis()->SetTitle("|B|_{#perp} [G]");
  gr->Update();
  gr->cd(3);
  graph2->Draw("AC");
  graph2->SetMarkerStyle(20);
  graph2->SetMarkerColor(1);
  graph2->SetLineStyle(1);
  graph2->SetLineWidth(3);
  graph2->SetTitle("Integral |B|_{#perp} to #vec{p}");
  //graph2->GetXaxis()->SetTitle("#eta");
  graph2->GetXaxis()->SetTitle("#theta [deg]");
  graph2->GetYaxis()->SetTitle("|B|_{#perp} [G]");
  gr->Update();


  /*TCanvas *gr_pi = new TCanvas("gr_pi","gr_pi",500,500);

  gr_pi->Divide(1,1);

  gr_pi->cd(1);
  graph_pi->Draw("AC");
  graph_pi->SetMarkerStyle(20);
  graph_pi->SetMarkerColor(1);
  graph_pi->SetLineStyle(1);
  graph_pi->SetLineWidth(2);
  graph_pi->SetTitle("");
  graph_pi->GetXaxis()->SetTitle("p [GeV/c]");
  graph_pi->GetYaxis()->SetTitle("#theta_{C} [mrad]");
  graph_pim->Draw("C");
  graph_pim->SetMarkerStyle(20);
  graph_pim->SetMarkerColor(1);
  graph_pim->SetLineStyle(1);
  graph_pim->SetLineColor(1);
  graph_pim->SetLineWidth(1);
  graph_pip->Draw("C");
  graph_pip->SetMarkerStyle(20);
  graph_pip->SetMarkerColor(1);
  graph_pip->SetLineStyle(1);
  graph_pip->SetLineColor(1);
  graph_pip->SetLineWidth(1);
  graph_k->Draw("C");
  graph_k->SetMarkerStyle(20);
  graph_k->SetMarkerColor(1);
  graph_k->SetLineStyle(1);
  graph_k->SetLineWidth(2);
  graph_pr->Draw("C");
  graph_pr->SetMarkerStyle(20);
  graph_pr->SetMarkerColor(1);
  graph_pr->SetLineStyle(1);
  graph_pr->SetLineWidth(2);
  gr_pi->Update();*/

}

