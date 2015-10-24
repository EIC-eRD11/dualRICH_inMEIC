//#######2D optical ray tracker C++/ROOT code#########
//Code version 1.1
//Author Alessio Del Dotto 

//To compile and run in root framework: 
//.L ray_tracker.cpp++
//ray_tracing *a = new ray_tracing()
//a->trac() 
//a->trac_fres() 
//a->trac_sep() 
//a->trac_fres_sep() 

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
#include <TArc.h>
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
#include <ray_tracker.h>

void ray_tracing::spherical_mirror(Double_t xc, Double_t yc, Double_t r, Double_t m, Double_t q){

  //Double_t xc, yc, r, m, q;
  Double_t a,b,c;

  //xc = 0.;
  //yc = 0.;
  //r = 4.;
  //m = TMath::Sqrt(3.)/3.;
  //q = -2.;

  a = m*m +1.;
  b = 2*(m*q-m*yc-xc);
  c = yc*yc-r*r+xc*xc-2*q*yc+q*q;

  xm = (-b-TMath::Sqrt(b*b-4*a*c))/(2*a);
  xp = (-b+TMath::Sqrt(b*b-4*a*c))/(2*a);

  ym = m*xm+q;
  yp = m*xp+q;

  phi = TMath::ATan2(yp-yc,xp-xc);
  alpha = TMath::ATan2(yp-q,xp);
  alpha1 = TMath::Pi()-(alpha-2*phi);

  //Double_t q1 = yp-TMath::Tan(alpha1)*xp;
  //Double_t m1 = TMath::Tan(alpha1);
  
  //cout<<xm<<"  "<<ym<<"  "<<xp<<"  "<<yp<<endl;
  //cout<<phi*180./3.14<<"  "<<alpha*180/3.14<<"  "<<alpha1*180./3.14<<endl;

  /*TH2F *plane = new TH2F("plane","",500.,0.,500.,500,0.,500.);
  plane->GetXaxis()->SetTitle("cm");
  plane->GetYaxis()->SetTitle("cm");

  TCanvas *cc1 = new TCanvas("cc1","cc1",700,700);
  cc1->Divide(1,1);
  cc1->cd(1);
  plane->Draw();

  TArc *mirror = new TArc(xc,yc,r,-60,60);
  mirror->SetFillStyle(0);
  mirror->Draw("only");

  TLine *ray1 = new TLine(0,q,xp,yp);
  TLine *ray2 = new TLine(0,q1,xp,yp);
  //TLine *ray2 = new TLine(xp,yp,500.,400.);
  ray1->SetLineColor(kRed);
  ray2->SetLineColor(kBlue);
  TLine *radius = new TLine(xc,yc,xp,yp);

  ray1->Draw();
  ray2->Draw();
  radius->Draw();

  Double_t mr = (yp-yc)/(xp-xc);

  Double_t ang = TMath::ATan(TMath::Abs((m-mr)/(1.+m*mr)))*180/3.14;
  Double_t ang1 = TMath::ATan(TMath::Abs((mr-m1)/(1.+mr*m1)))*180/3.14;
  
  cout<<ang<<"  "<<ang1<<endl;*/

}

void ray_tracing::spherical_mirror1(Double_t xc, Double_t yc, Double_t r, Double_t m, Double_t q){

  //Double_t xc, yc, r, m, q;
  Double_t a,b,c;

  //xc = 0.;
  //yc = 0.;
  //r = 4.;
  //m = TMath::Sqrt(3.)/3.;
  //q = -2.;

  a = m*m +1.;
  b = 2*(m*q-m*yc-xc);
  c = yc*yc-r*r+xc*xc-2*q*yc+q*q;

  xm_1 = (-b-TMath::Sqrt(b*b-4*a*c))/(2*a);
  xp_1 = (-b+TMath::Sqrt(b*b-4*a*c))/(2*a);

  ym_1 = m*xm_1+q;
  yp_1 = m*xp_1+q;

  phi_1 = TMath::ATan2(yp_1-yc,xp_1-xc);
  alpha_1 = TMath::ATan2(yp_1-q,xp_1);
  alpha1_1 = TMath::Pi()-(alpha_1-2*phi_1);
  
  //cout<<xm<<"  "<<ym<<"  "<<xp<<"  "<<yp<<endl;
  //cout<<phi*180./3.14<<"  "<<alpha*180/3.14<<"  "<<alpha1*180./3.14<<endl;

}


void ray_tracing::fresnel(Double_t xc, Double_t yc, Double_t r, Double_t m, Double_t q, Double_t n, Double_t n1, Double_t d1){

  Double_t aa,bb,cc;
  Double_t xm_f[4],xp_f[4],ym_f[4],yp_f[4],xci[4];
  Int_t count = 0;

  Double_t gamma, omega, theta_i, theta_r, theta_r1;
  Double_t d;

  d = xc-r;

  gamma = TMath::ATan(m);

  Double_t a,b,c;


  for(Int_t i=0; i<5; i++){

    aa = 1.;
    bb = -2.*yc;
    cc = yc*yc+(d-(xc-d1*i))*(d-(xc-d1*i))-r*r;

    sm[i] = (-bb-TMath::Sqrt(bb*bb-4*aa*cc))/(2*aa);
    sp[i] = (-bb+TMath::Sqrt(bb*bb-4*aa*cc))/(2*aa);

    ang[i]=180.-TMath::ASin((sp[i]-yc)/r)*180./TMath::Pi();

    //cout<<"Intersezioni perpendicolare: "<<sm[i]<<"  "<<sp[i]<<"  "<<i<<endl;
  }
  for(Int_t i=0; i<4; i++){
    a = m*m +1.;
    b = 2*(m*q-m*yc-(xc-d1*i));
    c = yc*yc-r*r+(xc-d1*i)*(xc-d1*i)-2*q*yc+q*q;

    xci[i] = xc-d1*i;
    
    xm_f[i] = (-b-TMath::Sqrt(b*b-4*a*c))/(2*a);
    xp_f[i] = (-b+TMath::Sqrt(b*b-4*a*c))/(2*a);
    
    ym_f[i] = m*xm_f[i]+q;
    yp_f[i] = m*xp_f[i]+q;

    //cout<<a<<"  "<<b<<"  "<<c<<"  "<<i<<endl;
    //cout<<xm_f[i]<<"  "<<ym_f[i]<<"  "<<xp_f[i]<<"  "<<yp_f[i]<<"  "<<i<<endl;
    //cout<<ym_f[i]<<"  "<<sm[i]<<"  "<<sm[i+1]<<"  "<<sp[i]<<"  "<<sp[i+1]<<endl;

    if(((ym_f[i]<sm[i]) && (ym_f[i]>sm[i+1])) || ((ym_f[i]>sp[i]) && (ym_f[i]<sp[i+1])) &&  (xm_f[i]>=d)  && (count == 0)){
      //if((TMath::Abs(ym_f[i])<sp[i+1]) &&  (xm_f[i]>=d)  && (count == 0)){ 
      ym_d = ym_f[i];  
      xm_d = xm_f[i];
      center = xci[i];
      //cout<<"Punti presi:  "<<ym_f[i]<<"  "<<sp[i+1]<<"  "<<xm_d<<"  "<<ym_d<<"  "<<count<<endl;
      count++;
    }
    else {
      //cout<<"Lens "<<i+1<<" not intercepted: "<<count<<endl;
    }
    
  }

  if(count > 0){
  if((m>=0) && ((ym_d -yc)>=0)){
  omega = TMath::Pi()-TMath::ATan2((ym_d-yc),(xm_d-center));
  theta_i = gamma+omega;
  theta_r = TMath::ASin(n*TMath::Sin(theta_i)/n1);
  gamma1 = omega - theta_r;
  m1_d1 = TMath::ATan(-gamma1);
  q1_d1 = ym_d-m1_d1*xm_d;
  theta_r1 = TMath::ASin(n1*TMath::Sin(gamma1)/n);
  m2_d1 = TMath::ATan(-theta_r1);

  q2_d1 = m1_d1*(xc-r+d1)+q1_d1 - m2_d1*(xc-r+d1);
  }
  if(m<0 && ((ym_d -yc)>=0)){
  omega = TMath::Pi()-TMath::ATan2((ym_d-yc),(xm_d-center));
  theta_i = gamma+omega;
  theta_r = TMath::ASin(n*TMath::Sin(theta_i)/n1);
  gamma1 = omega - theta_r;
  m1_d1 = TMath::ATan(-gamma1);
  q1_d1 = ym_d-m1_d1*xm_d;
  theta_r1 = TMath::ASin(n1*TMath::Sin(gamma1)/n);
  m2_d1 = TMath::ATan(-theta_r1);

  q2_d1 = m1_d1*(xc-r+d1)+q1_d1 - m2_d1*(xc-r+d1);
  }

  if((m<=0) && ((ym_d -yc)<0)){
  gamma = -gamma;  
  omega = TMath::Pi()+TMath::ATan2((ym_d-yc),(xm_d-center));
  theta_i = gamma+omega;
  theta_r = TMath::ASin(n*TMath::Sin(theta_i)/n1);
  gamma1 = omega - theta_r;
  m1_d1 = TMath::ATan(gamma1);
  q1_d1 = ym_d-m1_d1*xm_d;
  theta_r1 = TMath::ASin(n1*TMath::Sin(gamma1)/n);
  m2_d1 = TMath::ATan(theta_r1);

  q2_d1 = m1_d1*(xc-r+d1)+q1_d1 - m2_d1*(xc-r+d1);
  }
  if(m>0 && ((ym_d -yc)<=0)){
  gamma = -gamma;
  omega = TMath::Pi()+TMath::ATan2((ym_d-yc),(xm_d-center));
  theta_i = gamma+omega;
  theta_r = TMath::ASin(n*TMath::Sin(theta_i)/n1);
  gamma1 = omega - theta_r;
  m1_d1 = TMath::ATan(gamma1);
  q1_d1 = ym_d-m1_d1*xm_d;
  theta_r1 = TMath::ASin(n1*TMath::Sin(gamma1)/n);
  m2_d1 = TMath::ATan(theta_r1);

  q2_d1 = m1_d1*(xc-r+d1)+q1_d1 - m2_d1*(xc-r+d1);
  }
  

  //cout<<xm_d<<"  "<<ym_d<<"  "<<m2_d1<<"  "<<q2_d1<<endl;

  //cout<<gamma*180./3.14<<"  "<<omega*180./3.14<<"  "<<theta_i*180./3.14<<"  "<<theta_r*180./3.14<<"  "<<gamma1*180./3.14<<endl;
  //cout<<xm_d<<"  "<<ym_d<<"  "<<xp_d<<"  "<<yp_d<<endl;

  /*TH2F *plane = new TH2F("plane","",500.,0.,1000.,500,-1000.,1000.);
  plane->GetXaxis()->SetTitle("cm");
  plane->GetYaxis()->SetTitle("cm");

  TCanvas *cc1 = new TCanvas("cc1","cc1",700,700);
  cc1->Divide(1,1);
  cc1->cd(1);
  plane->Draw();
  
  TLine *diot1 = new TLine(xc-r+d1,sm[4],xc-r+d1,sp[4]);
  TLine *diotr = new TLine(xm_d,ym_d,xc,yc);
  TArc *diot2 = new TArc(xc,yc,r,ang[1],ang[0]);
  TArc *diot3 = new TArc(xc-d1,yc,r,ang[2],ang[1]);
  TArc *diot4 = new TArc(xc-d1*2,yc,r,ang[3],ang[2]);
  TArc *diot22 = new TArc(xc,yc,r,-ang[0],-ang[1]);
  TArc *diot33 = new TArc(xc-d1,yc,r,-ang[1],-ang[2]);
  TArc *diot44 = new TArc(xc-d1*2,yc,r,-ang[2],-ang[3]);
  //TArc *diot5 = new TArc(xc-d1*3,yc,r,ang[4],ang[3]);
  TArc *diot55 = new TArc(xc-d1*3,yc,r,-ang[4],-ang[3]);
  TArc *diot5 = new TArc(xc-d1*3,yc,r,0,360);
  TLine *diot6 = new TLine(xc-r,-800,xc-r,800);
  TLine *diot7 = new TLine(xc-r+20,-800,xc-r+20,800);
  TLine *diot8 = new TLine(xc-r-20,-800,xc-r-20,800);
  diot1->Draw();
  diotr->Draw();
  diot2->SetFillStyle(0);
  diot3->SetFillStyle(0);
  diot4->SetFillStyle(0);
  diot22->SetFillStyle(0);
  diot33->SetFillStyle(0);
  diot44->SetFillStyle(0);
  diot5->SetFillStyle(0);
  diot55->SetFillStyle(0);
  diot2->Draw("only");
  diot3->Draw("only");
  diot4->Draw("only");
  diot22->Draw("only");
  diot33->Draw("only");
  diot44->Draw("only");
  diot5->Draw("only");
  diot55->Draw("only");
  //diot6->Draw();
  //diot7->Draw();
  //diot8->Draw();

  TLine *ray1 = new TLine(-1000,m*(-1000)+q,xm_d,ym_d);
  TLine *ray2 = new TLine(xm_d,ym_d,(xc-r+d1),m1_d1*(xc-r+d1)+q1_d1);
  TLine *ray3 = new TLine((xc-r+d1),m1_d1*(xc-r+d1)+q1_d1,900,m2_d1*900+q2_d1);
  ray1->SetLineColor(kRed);
  ray2->SetLineColor(kBlue);
  ray3->SetLineColor(26);
  TLine *radius = new TLine(xm_d,ym_d,xc,yc);

  ray1->Draw();
  ray2->Draw();
  ray3->Draw();
  //radius->Draw();

  //cout<<xm_d<<"  "<<ym_d<<"  "<<xp_d<<"  "<<yp_d<<endl;
  //cout<<gamma<<"  "<<omega<<"  "<<theta_r<<"  "<<gamma1<<endl;*/
  }
  else{cout<<"no lens intercepted"<<endl;}

  count = 0;

}


void ray_tracing::fresnel1(Double_t xc, Double_t yc, Double_t r, Double_t m, Double_t q, Double_t n, Double_t n1, Double_t d1){

  Double_t aa,bb,cc;
  Double_t xm_f[4],xp_f[4],ym_f[4],yp_f[4],xci[4];
  Int_t count = 0;

  Double_t gamma, omega, theta_i, theta_r, theta_r1;
  Double_t d;

  d = xc+r;

  theta_r = TMath::ASin(n*TMath::Sin(TMath::ATan(m))/n1);
  m1_d1 = TMath::Tan(theta_r);
  q1_d1 = m*(d-d1)+q-m1_d1*(d-d1);

  gamma = TMath::ATan(m1_d1);

  yr=m*(d-d1)+q;
  xr=d-d1;

  Double_t a,b,c;


  for(Int_t i=0; i<5; i++){

    aa = 1.;
    bb = -2.*yc;
    cc = yc*yc+(d-(xc+d1*i))*(d-(xc+d1*i))-r*r;

    sm[i] = (-bb-TMath::Sqrt(bb*bb-4*aa*cc))/(2*aa);
    sp[i] = (-bb+TMath::Sqrt(bb*bb-4*aa*cc))/(2*aa);

    ang[i]= TMath::ASin((sp[i]-yc)/r)*180./TMath::Pi();

    //cout<<"Intersezioni perpendicolare: "<<sm[i]<<"  "<<sp[i]<<"  "<<i<<endl;
    //cout<<"Angoli:  "<<ang[i]<<endl;
  }
  for(Int_t i=0; i<4; i++){
    a = m1_d1*m1_d1 +1.;
    b = 2*(m1_d1*q1_d1-m1_d1*yc-(xc+d1*i));
    c = yc*yc-r*r+(xc+d1*i)*(xc+d1*i)-2*q1_d1*yc+q1_d1*q1_d1;

    xci[i] = xc-d1*i;
    
    xm_f[i] = (-b-TMath::Sqrt(b*b-4*a*c))/(2*a);
    xp_f[i] = (-b+TMath::Sqrt(b*b-4*a*c))/(2*a);
    
    ym_f[i] = m1_d1*xm_f[i]+q1_d1;
    yp_f[i] = m1_d1*xp_f[i]+q1_d1;

    //cout<<a<<"  "<<b<<"  "<<c<<"  "<<i<<endl;
    //cout<<xm_f[i]<<"  "<<ym_f[i]<<"  "<<xp_f[i]<<"  "<<yp_f[i]<<"  "<<i<<endl;
    //cout<<ym_f[i]<<"  "<<sm[i]<<"  "<<sm[i+1]<<"  "<<sp[i]<<"  "<<sp[i+1]<<endl;

    if(((yp_f[i]<sm[i]) && (yp_f[i]>sm[i+1])) || ((yp_f[i]>sp[i]) && (yp_f[i]<sp[i+1])) &&  (xp_f[i]<=d)  && (count == 0)){
      //if((TMath::Abs(ym_f[i])<sp[i+1]) &&  (xm_f[i]>=d)  && (count == 0)){ 
      yp_d = yp_f[i];  
      xp_d = xp_f[i];
      center = xci[i];
      //cout<<"Punti presi:  "<<ym_f[i]<<"  "<<sp[i+1]<<"  "<<xp_d<<"  "<<yp_d<<"  "<<count<<endl;
      count++;
    }
    else {
      //cout<<"Lens "<<i+1<<" not intercepted: "<<count<<endl;
    }
    
  }

  if(count > 0){
  if((m>=0) && ((yp_d -yc)>=0)){
  omega = TMath::ATan2((yp_d-yc),(xp_d-center));
  theta_i = omega-gamma;
  theta_r1 = TMath::ASin(n1*TMath::Sin(theta_i)/n);
  gamma1 = omega - theta_r1;
  m2_d1 = TMath::ATan(gamma1);
  q2_d1 = yp_d-m2_d1*xp_d;
 
  }
  if(m<0 && ((yp_d -yc)>=0)){
  omega = TMath::ATan2((yp_d-yc),(xp_d-center));
  theta_i = -gamma+omega;
  theta_r1 = TMath::ASin(n1*TMath::Sin(theta_i)/n);
  gamma1 = omega - theta_r1;
  m2_d1 = TMath::ATan(gamma1);
  q2_d1 = yp_d-m2_d1*xp_d;

  }

  if((m<=0) && ((yp_d -yc)<0)){
  gamma = -gamma;  
  omega = TMath::ATan2((yc-yp_d),(xp_d-center));
  theta_i = -gamma+omega;
  theta_r1 = TMath::ASin(n1*TMath::Sin(theta_i)/n);
  gamma1 = omega - theta_r1;
  m2_d1 = TMath::ATan(-gamma1);
  q2_d1 = yp_d-m2_d1*xp_d;
  
  }
  if(m>0 && ((yp_d -yc)<=0)){
  omega = TMath::ATan2((yc-yp_d),(xp_d-center));
  theta_i = gamma+omega;
  theta_r1 = TMath::ASin(n1*TMath::Sin(theta_i)/n);
  gamma1 = omega - theta_r1;
  m2_d1 = TMath::ATan(-gamma1);
  q2_d1 = yp_d-m2_d1*xp_d;
  //cout<<"gamma: "<<gamma*180./3.14<<" omega: "<<omega*180./3.14<<" theta_i "<<theta_r1*180./3.14<<endl;
  }
  

  //cout<<xm_d<<"  "<<ym_d<<"  "<<m2_d1<<"  "<<q2_d1<<endl;

  //cout<<"Vari angoli: "<<gamma*180./3.14<<"  "<<omega*180./3.14<<"  "<<theta_i*180./3.14<<"  "<<theta_r1*180./3.14<<"  "<<gamma1*180./3.14<<endl;
  //cout<<xm_d<<"  "<<ym_d<<"  "<<xp_d<<"  "<<yp_d<<endl;

  /*TH2F *plane = new TH2F("plane","",500.,0.,1000.,500,-1000.,1000.);
  plane->GetXaxis()->SetTitle("cm");
  plane->GetYaxis()->SetTitle("cm");

  TCanvas *cc1 = new TCanvas("cc1","cc1",700,700);
  cc1->Divide(1,1);
  cc1->cd(1);
  plane->Draw();
  
  TLine *diot1 = new TLine(xc+r-d1,sm[4],xc+r-d1,sp[4]);
  TLine *diotr = new TLine(xm_d,ym_d,center,yc);
  TArc *diot2 = new TArc(xc,yc,r,ang[1],ang[0]);
  TArc *diot3 = new TArc(xc+d1,yc,r,ang[2],ang[1]);
  TArc *diot4 = new TArc(xc+d1*2,yc,r,ang[3],ang[2]);
  TArc *diot22 = new TArc(xc,yc,r,-ang[0],-ang[1]);
  TArc *diot33 = new TArc(xc+d1,yc,r,-ang[1],-ang[2]);
  TArc *diot44 = new TArc(xc+d1*2,yc,r,-ang[2],-ang[3]);
  //TArc *diot5 = new TArc(xc-d1*3,yc,r,ang[4],ang[3]);
  TArc *diot55 = new TArc(xc+d1*3,yc,r,-ang[4],-ang[3]);
  TArc *diot5 = new TArc(xc+d1*3,yc,r,0,360);
  TLine *diot6 = new TLine(xc-r,-800,xc-r,800);
  TLine *diot7 = new TLine(xc-r+20,-800,xc-r+20,800);
  TLine *diot8 = new TLine(xc-r-20,-800,xc-r-20,800);
  diot1->Draw();
  //diotr->Draw();
  diot2->SetFillStyle(0);
  diot3->SetFillStyle(0);
  diot4->SetFillStyle(0);
  diot22->SetFillStyle(0);
  diot33->SetFillStyle(0);
  diot44->SetFillStyle(0);
  diot5->SetFillStyle(0);
  diot55->SetFillStyle(0);
  diot2->Draw("only");
  diot3->Draw("only");
  diot4->Draw("only");
  diot22->Draw("only");
  diot33->Draw("only");
  diot44->Draw("only");
  diot5->Draw("only");
  diot55->Draw("only");
  //diot6->Draw();
  //diot7->Draw();
  //diot8->Draw();

  TLine *ray1 = new TLine(-1000,m*(-1000)+q,xr,yr);
  TLine *ray2 = new TLine(xr,yr,xp_d,yp_d);
  TLine *ray3 = new TLine(xp_d,yp_d,900,m2_d1*900+q2_d1);
  ray1->SetLineColor(kRed);
  ray2->SetLineColor(kBlue);
  ray3->SetLineColor(26);
  TLine *radius = new TLine(xm_d,ym_d,center,yc);

  ray1->Draw();
  ray2->Draw();
  ray3->Draw();
  //radius->Draw();

  //cout<<xm_d<<"  "<<ym_d<<"  "<<xp_d<<"  "<<yp_d<<endl;
  //cout<<gamma<<"  "<<omega<<"  "<<theta_r<<"  "<<gamma1<<endl;*/
  }
  else{cout<<"no lens intercepted"<<endl;}

  count = 0;

}



void ray_tracing::set_detector(){

  x1_det = 250.;
  x2_det = 250.1;
  y1_det = 130.;
  y2_det = 165.;

}

void ray_tracing::trac(Double_t xc, Double_t yc, Double_t r, Double_t xc1, Double_t yc1, Double_t r1, Double_t y_sep, Double_t ang){

  //Double_t xc = 100.;
  //Double_t yc = 100.;
  //Double_t r = 300.;
  //Double_t m = 1.;
  //Double_t q = 0.;
  ray_tracing *det_pos = new ray_tracing();
  det_pos->set_detector();

  Double_t theta = ang*TMath::Pi()/180.;
  Double_t y_track = TMath::Tan(theta)*1000;

  Double_t rf[100],rf1[100],t_ang[100];
  Int_t ang_count =0;

  Double_t p_start_x[8];
  Double_t p_start_y[8];
  Double_t q[8];
  Double_t q1[8],q1_1[8],m1[8],m1_1[8];

  Double_t y_up[8], y_down[8];

  Double_t momentum =5.;
  Double_t momentumg =35.;

  Double_t theta_c =acos(sqrt(momentum*momentum+0.139570*0.139570)/momentum/1.02);
  Double_t m = TMath::Tan(theta_c+theta);
  Double_t mm = TMath::Tan(-theta_c+theta);

  Double_t theta_cg = acos(sqrt(momentumg*momentumg+0.139570*0.139570)/momentumg/1.000482);
  Double_t mg = TMath::Tan(theta_c+theta);
  Double_t mmg = TMath::Tan(-theta_c+theta);

  TH2F *plane = new TH2F("plane","",500.,0.,500.,500,0.,500.);
  plane->GetXaxis()->SetTitle("cm");
  plane->GetYaxis()->SetTitle("cm");
  TH2F *focal_plane = new TH2F("plane","",500.,0.,500.,500,0.,500.);
  focal_plane->GetXaxis()->SetTitle("cm");
  focal_plane->GetYaxis()->SetTitle("cm");
  TH2F *focal_plane_g = new TH2F("plane","",500.,0.,500.,500,0.,500.);
  focal_plane_g->GetXaxis()->SetTitle("cm");
  focal_plane_g->GetYaxis()->SetTitle("cm");
  TH2F *plane_g = new TH2F("plane","",500.,0.,500.,500,0.,500.);
  plane_g->GetXaxis()->SetTitle("cm");
  plane_g->GetYaxis()->SetTitle("cm");

  //gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  //gStyle->SetTitleW(0.4);
  //gStyle->SetTitleH(0.095);
  gStyle->SetTitleSize(0.03,"xyzt");
  //gStyle->SetTitleXOffset(.15);
  //gStyle->SetTitleYOffset(.12);
  gStyle->SetLabelSize(0.03,"XYZ");
  gStyle->SetPadTopMargin(.10);
  gStyle->SetPadLeftMargin(.14);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadBottomMargin(.15);
  gStyle->SetPalette(1);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetFrameBorderMode(0);

  TCanvas *c1 = new TCanvas("c1","c1",700,700);
  c1->Divide(1,1);
  c1->cd(1);
  plane->Draw();

  TLine *d_oriz1 = new TLine(250.,10.,415.,10.);
  d_oriz1->SetLineColor(kGreen);
  TLine *d_oriz2 = new TLine(250.,100.,415.,150.);
  d_oriz2->SetLineColor(kGreen);
  TLine *d_vert1 = new TLine(250.,10.,250.,100.);
  d_vert1->SetLineColor(kGreen);
  TLine *d_vert2 = new TLine(415.,10.,415.,150.);
  d_vert2->SetLineColor(kGreen);
  
  d_oriz1->Draw();
  d_oriz2->Draw();
  d_vert1->Draw();
  d_vert2->Draw();

  for(Double_t p=5;p<=20;p=p+5){
    
    ang = p;
    theta = ang*TMath::Pi()/180.;
    y_track = TMath::Tan(theta)*1000;
    m = TMath::Tan(theta_c+theta);
    mm = TMath::Tan(-theta_c+theta);

  for(Int_t i=0;i<8;i++){
    if(i<4){
      p_start_x[i] = 250.+ 1.*i;
      p_start_y[i] = p_start_x[i]*TMath::Tan(theta);
      q[i] = p_start_y[i]-m*p_start_x[i];
    }
    if(i>=4){
      p_start_x[i] = 250.+ 1.*(i-4);
      p_start_y[i] = p_start_x[i]*TMath::Tan(theta);
      q[i] = p_start_y[i]-mm*p_start_x[i];
    }
  }

  ray_tracing *a = new ray_tracing();
  ray_tracing *b = new ray_tracing();

  for(Int_t i=0;i<8;i++){

    if(i<4) a->spherical_mirror(xc,yc,r,m,q[i]);
    if(i>=4) a->spherical_mirror(xc,yc,r,mm,q[i]);

    if(i<4) b->spherical_mirror1(xc1,yc1,r1,m,q[i]);
    if(i>=4) b->spherical_mirror1(xc1,yc1,r1,mm,q[i]);

    //cout<<a->xp<<"  "<<a->yp<<"  "<<b->xp_1<<"  "<<b->yp_1<<endl;
    //cout<<a->xm<<"  "<<a->ym<<"  "<<a->xp<<"  "<<a->yp<<endl;
    //cout<<a->phi*180./3.14<<"  "<<a->alpha*180/3.14<<"  "<<a->alpha1*180./3.14<<endl;

    q1[i] = a->yp-TMath::Tan(a->alpha1)*a->xp;
    q1_1[i] = b->yp_1-TMath::Tan(b->alpha1_1)*b->xp_1;
    m1[i] = TMath::Tan(a->alpha1);
    m1_1[i] = TMath::Tan(b->alpha1_1);
    y_up[i] = a->yp;
    y_down[i] = b->yp_1;
    
    Double_t m_rr = (b->yp_1-yc1)/(b->xp_1-xc1); 

    Double_t ang_r = TMath::ATan(TMath::Abs((m-m_rr)/(1.+m*m_rr)))*180./3.14;
    Double_t ang_r1 = TMath::ATan(TMath::Abs((m1_1[i]-m_rr)/(1.+m1_1[i]*m_rr)))*180./3.14;

    //if(i<4) cout<<ang_r<<"  "<<ang_r1<<endl;
    //cout<<m1[i]<<"  "<<q1[i]<<endl;
    //TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
    //c1->Divide(1,1);
    //c1->cd(1);
    
    //Double_t y_sep = 100.;
    Double_t m_ang = TMath::ASin((y_sep-yc)/r)*180./TMath::Pi();
    Double_t m_ang_min = TMath::ASin((0.-yc)/r)*180./TMath::Pi();
    Double_t m_ang_max = TMath::ASin((160.-yc1)/r1)*180./TMath::Pi();
    
    TLine *x_axis = new TLine(0,0,500,0);
    TLine *y_axis = new TLine(0,0,0,500);
    TLine *detector = new TLine(det_pos->x1_det,det_pos->y1_det,det_pos->x2_det,det_pos->y2_det);
    TLine *track = new TLine(0,0,1000,y_track);
    TArc *mirror = new TArc(xc,yc,r,m_ang_min,m_ang);
    TArc *mirror1 = new TArc(xc1,yc1,r1,m_ang,m_ang_max);
    TArc *focus = new TArc(xc,yc,r/2.,-90,90);
    TArc *focus1 = new TArc(xc1,yc1,r1/2.,-90,90);
    mirror->SetFillStyle(0);
    mirror1->SetFillStyle(0);
    focus->SetFillStyle(0);
    focus->SetLineColor(kBlue);
    focus1->SetFillStyle(0);
    focus1->SetLineColor(28);
    TLine *ray1 = new TLine(p_start_x[i],p_start_y[i],a->xp,a->yp);
    TLine *ray2 = new TLine(0,q1[i],a->xp,a->yp);
    TLine *ray2u = new TLine(a->xp,a->yp,500.,m1[i]*500.+q1[i]);
    TLine *radius = new TLine(xc,yc,a->xp,a->yp);
    TLine *ray1_1 = new TLine(p_start_x[i],p_start_y[i],b->xp_1,b->yp_1);
    TLine *ray2_1 = new TLine(0,q1_1[i],b->xp_1,b->yp_1);
    TLine *ray2_1u = new TLine(b->xp_1,b->yp_1,500.,m1_1[i]*500.+q1_1[i]);
    TLine *radius_1 = new TLine(xc1,yc1,b->xp_1,b->yp_1);
    ray1->SetLineColor(kRed);
    ray2->SetLineColor(kBlue);
    ray2u->SetLineColor(kBlue);
    ray1_1->SetLineColor(kRed);
    ray2_1->SetLineColor(28);
    detector->SetLineColor(9);
    detector->SetLineWidth(3);
    radius->SetLineStyle(3);
    radius_1->SetLineStyle(3);
    x_axis->Draw();
    y_axis->Draw();
    mirror->Draw("only");
    mirror1->Draw("only");
    if(a->yp<y_sep){
      ray1->Draw();
      //ray2->Draw();
      if(q1[i]>=0)ray2->Draw();
      if(q1[i]<0)ray2u->Draw();
      radius->Draw();
      focus->Draw("only");
    }
    if(b->yp_1>=y_sep){
      ray1_1->Draw();
      //ray2_1->Draw();
      if(q1[i]>=0)ray2_1->Draw();
      if(q1[i]<0)ray2_1u->Draw();
      radius_1->Draw();
      //focus1->Draw("only");
    }
    track->Draw();
    detector->Draw();
    
  }

  Int_t c = 0;
  Double_t x_f, y_f;
  Double_t rf_sum=0.;
  Double_t rf1_sum=0.;

  while(c<4){
    for(Int_t i=0+c;i<3;i++){
      if(y_up[c]<y_sep){
	x_f = (q1[i+1]-q1[c])/(m1[c]-m1[i+1]);
        y_f = m1[c]*x_f+q1[c];
	rf_sum = rf_sum + sqrt((x_f-xc)*(x_f-xc)+(y_f-yc)*(y_f-yc));
	focal_plane->Fill(x_f,y_f);
      }
      if(y_up[c]>=y_sep){
	x_f = (q1_1[i+1]-q1_1[c])/(m1_1[c]-m1_1[i+1]);
	y_f = m1_1[c]*x_f+q1_1[c];
	focal_plane->Fill(x_f,y_f);
      }
      //cout<<rf_sum<<"  "<<x_f<<"  "<<y_f<<endl;
    }
    for(Int_t i=4+c;i<7;i++){
      if(y_up[c+4]<y_sep){
	x_f = (q1[i+1]-q1[c+4])/(m1[c+4]-m1[i+1]);
	y_f = m1[c+4]*x_f+q1[c+4];
	rf1_sum = rf1_sum + sqrt((x_f-xc)*(x_f-xc)+(y_f-yc)*(y_f-yc));
	focal_plane->Fill(x_f,y_f);
      }
      if(y_up[c+4]>=y_sep){
	x_f = (q1_1[i+1]-q1_1[c+4])/(m1_1[c+4]-m1_1[i+1]);
	y_f = m1_1[c+4]*x_f+q1_1[c+4];
	focal_plane->Fill(x_f,y_f);
      }
      //cout<<rf1_sum<<"  "<<x_f<<"  "<<y_f<<endl;
    }
    c++;
  }
  
  t_ang[ang_count] = p;
  rf[ang_count] = (rf_sum/6.)/(r/2.);
  rf1[ang_count] = (rf1_sum/6.)/(r/2.);
  ang_count++;

  rf_sum = 0.;
  rf1_sum = 0.;

  focal_plane->SetMarkerStyle(7);
  focal_plane->Draw("same");
  c1->Update(); 

  }

  TGraph *gr = new TGraph(ang_count,t_ang,rf);
  TGraph *gr1 = new TGraph(ang_count,t_ang,rf1);

  TCanvas *cf1 = new TCanvas("cf1","cf1",500,500);
  cf1->Divide(1,2);
  cf1->cd(1);
  gr->Draw("AC*");
  cf1->cd(2);
  gr1->Draw("AC*");

  //Gas generation.....................................................

  TCanvas *c2 = new TCanvas("c2","c2",700,700);
  c2->Divide(1,1);
  c2->cd(1);
  plane_g->Draw();
  
  d_oriz1->Draw();
  d_oriz2->Draw();
  d_vert1->Draw();
  d_vert2->Draw();


  for(Double_t p=5;p<=20;p=p+5){
    
    ang = p;
    theta = ang*TMath::Pi()/180.;
    y_track = TMath::Tan(theta)*1000;
    mg = TMath::Tan(theta_cg+theta);
    mmg = TMath::Tan(-theta_cg+theta);

  for(Int_t i=0;i<8;i++){
    if(i<4){
      p_start_x[i] = 250.+ 35.*i;
      p_start_y[i] = p_start_x[i]*TMath::Tan(theta);
      q[i] = p_start_y[i]-mg*p_start_x[i];
    }
    if(i>=4){
      p_start_x[i] = 250.+ 35.*(i-4);
      p_start_y[i] = p_start_x[i]*TMath::Tan(theta);
      q[i] = p_start_y[i]-mmg*p_start_x[i];
    }
  }

  ray_tracing *a = new ray_tracing();
  ray_tracing *b = new ray_tracing();

  for(Int_t i=0;i<8;i++){

    if(i<4) a->spherical_mirror(xc,yc,r,mg,q[i]);
    if(i>=4) a->spherical_mirror(xc,yc,r,mmg,q[i]);

    if(i<4) b->spherical_mirror1(xc1,yc1,r1,mg,q[i]);
    if(i>=4) b->spherical_mirror1(xc1,yc1,r1,mmg,q[i]);

    //cout<<a->xp<<"  "<<a->yp<<"  "<<b->xp_1<<"  "<<b->yp_1<<endl;
    //cout<<a->xm<<"  "<<a->ym<<"  "<<a->xp<<"  "<<a->yp<<endl;
    //cout<<a->phi*180./3.14<<"  "<<a->alpha*180/3.14<<"  "<<a->alpha1*180./3.14<<endl;

    q1[i] = a->yp-TMath::Tan(a->alpha1)*a->xp;
    q1_1[i] = b->yp_1-TMath::Tan(b->alpha1_1)*b->xp_1;
    m1[i] = TMath::Tan(a->alpha1);
    m1_1[i] = TMath::Tan(b->alpha1_1);
    y_up[i] = a->yp;
    y_down[i] = b->yp_1;
    
    Double_t m_rr = (b->yp_1-yc1)/(b->xp_1-xc1); 

    Double_t ang_r = TMath::ATan(TMath::Abs((m-m_rr)/(1.+m*m_rr)))*180./3.14;
    Double_t ang_r1 = TMath::ATan(TMath::Abs((m1_1[i]-m_rr)/(1.+m1_1[i]*m_rr)))*180./3.14;

    //if(i<4) cout<<ang_r<<"  "<<ang_r1<<endl;
    //cout<<m1[i]<<"  "<<q1[i]<<endl;
    //TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
    //c1->Divide(1,1);
    //c1->cd(1);
    
    //Double_t y_sep = 100.;
    Double_t m_ang = TMath::ASin((y_sep-yc)/r)*180./TMath::Pi();
    Double_t m_ang_min = TMath::ASin((0.-yc)/r)*180./TMath::Pi();
    Double_t m_ang_max = TMath::ASin((160.-yc1)/r1)*180./TMath::Pi();
    
    TLine *x_axis = new TLine(0,0,500,0);
    TLine *y_axis = new TLine(0,0,0,500);
    TLine *detector = new TLine(det_pos->x1_det,det_pos->y1_det,det_pos->x2_det,det_pos->y2_det);
    TLine *track = new TLine(0,0,1000,y_track);
    TArc *mirror = new TArc(xc,yc,r,m_ang_min,m_ang);
    TArc *mirror1 = new TArc(xc1,yc1,r1,m_ang,m_ang_max);
    TArc *focus = new TArc(xc,yc,r/2.,-90,90);
    TArc *focus1 = new TArc(xc1,yc1,r1/2.,-90,90);
    mirror->SetFillStyle(0);
    mirror1->SetFillStyle(0);
    focus->SetFillStyle(0);
    focus->SetLineColor(kBlue);
    focus1->SetFillStyle(0);
    focus1->SetLineColor(28);
    TLine *ray1 = new TLine(p_start_x[i],p_start_y[i],a->xp,a->yp);
    TLine *ray2 = new TLine(0,q1[i],a->xp,a->yp);
    TLine *ray2u = new TLine(a->xp,a->yp,500.,m1[i]*500.+q1[i]);
    TLine *radius = new TLine(xc,yc,a->xp,a->yp);
    TLine *ray1_1 = new TLine(p_start_x[i],p_start_y[i],b->xp_1,b->yp_1);
    TLine *ray2_1 = new TLine(0,q1_1[i],b->xp_1,b->yp_1);
    TLine *ray2_1u = new TLine(b->xp_1,b->yp_1,500.,m1_1[i]*500.+q1_1[i]);
    TLine *radius_1 = new TLine(xc1,yc1,b->xp_1,b->yp_1);
    ray1->SetLineColor(kRed);
    ray2->SetLineColor(kBlue);
    ray2u->SetLineColor(kBlue);
    ray1_1->SetLineColor(kRed);
    ray2_1->SetLineColor(28);
    detector->SetLineColor(9);
    detector->SetLineWidth(3);
    radius->SetLineStyle(3);
    radius_1->SetLineStyle(3);
    x_axis->Draw();
    y_axis->Draw();
    mirror->Draw("only");
    mirror1->Draw("only");
    if(a->yp<y_sep){
      ray1->Draw();
      if(q1[i]>=0)ray2->Draw();
      if(q1[i]<0)ray2u->Draw();
      radius->Draw();
      focus->Draw("only");
    }
    if(b->yp_1>=y_sep){
      ray1_1->Draw();
      if(q1[i]>=0)ray2_1->Draw();
      if(q1[i]<0)ray2_1u->Draw();
      radius_1->Draw();
      //focus1->Draw("only");
    }
    track->Draw();
    detector->Draw();
    
  }

  Int_t c = 0;
  Double_t x_f, y_f;

  while(c<4){
    for(Int_t i=0+c;i<3;i++){
      if(y_up[c]<y_sep){
	x_f = (q1[i+1]-q1[c])/(m1[c]-m1[i+1]);
        y_f = m1[c]*x_f+q1[c];
	focal_plane_g->Fill(x_f,y_f);
      }
      if(y_up[c]>=y_sep){
	x_f = (q1_1[i+1]-q1_1[c])/(m1_1[c]-m1_1[i+1]);
	y_f = m1_1[c]*x_f+q1_1[c];
	focal_plane_g->Fill(x_f,y_f);
      }
      //cout<<i+1<<"  "<<x_f<<"  "<<y_f<<endl;
    }
    for(Int_t i=4+c;i<7;i++){
      if(y_up[c+4]<y_sep){
	x_f = (q1[i+1]-q1[c+4])/(m1[c+4]-m1[i+1]);
	y_f = m1[c+4]*x_f+q1[c+4];
	focal_plane_g->Fill(x_f,y_f);
      }
      if(y_up[c+4]>=y_sep){
	x_f = (q1_1[i+1]-q1_1[c+4])/(m1_1[c+4]-m1_1[i+1]);
	y_f = m1_1[c+4]*x_f+q1_1[c+4];
	focal_plane_g->Fill(x_f,y_f);
      }
      //cout<<i+1<<"  "<<x_f<<"  "<<y_f<<endl;
    }
    c++;
  }

  focal_plane_g->SetMarkerStyle(7);
  focal_plane_g->Draw("same");
  c2->Update(); 

  }



}

void ray_tracing::trac_sep(Double_t xc, Double_t yc, Double_t r, Double_t xc1, Double_t yc1, Double_t r1, Double_t y_sep, Double_t ang, Int_t radiator){

  //Double_t xc = 100.;
  //Double_t yc = 100.;
  //Double_t r = 300.;
  //Double_t m = 1.;
  //Double_t q = 0.;

  ray_tracing *det_pos = new ray_tracing();
  det_pos->set_detector();

  Double_t momentum;

  Double_t theta = ang*TMath::Pi()/180.;
  Double_t y_track = TMath::Tan(theta)*1000;

  Double_t m_det = (det_pos->y2_det-det_pos->y1_det)/(det_pos->x2_det-det_pos->x1_det);
  Double_t q_det = det_pos->y1_det - m_det*det_pos->x1_det;

  Double_t p_start_x[8];
  Double_t p_start_y[8];
  Double_t q[8], mv[8];
  Double_t q1[8],q1_1[8],m1[8],m1_1[8];

  Double_t y_up[8], y_down[8];
  Double_t ayp[8],byp[8];

  Double_t theta_c;
  Double_t theta_ck;
  theta_c=acos(sqrt(momentum*momentum+0.139570*0.139570)/momentum/1.02);
  theta_ck=acos(sqrt(momentum*momentum+0.493677*0.493677)/momentum/1.02);

  Double_t spatial_sep[1000], spatial_sep1[1000], spatial_var[100], spatial_var1[100];
  Double_t angle_c[10], spatial_sep_mean[10], spatial_sep_mean1[10],momentum_points[10];

  TGraph *grup[4];
  TGraph *grdown[4];
  TGraph *grvup[4];
  TGraph *grvdown[4];  

  Double_t m = TMath::Tan(theta_c+theta);
  Double_t mm = TMath::Tan(-theta_c+theta);
  Double_t mk = TMath::Tan(theta_ck+theta);
  Double_t mmk = TMath::Tan(-theta_ck+theta);

  TH2F *plane = new TH2F("plane","",500.,0.,500.,500,0.,500.);
  plane->GetXaxis()->SetTitle("cm");
  plane->GetYaxis()->SetTitle("cm");
  TH2F *focal_plane = new TH2F("plane","",500.,0.,500.,500,0.,500.);
  focal_plane->GetXaxis()->SetTitle("cm");
  focal_plane->GetYaxis()->SetTitle("cm");

  //gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  //gStyle->SetTitleW(0.4);
  //gStyle->SetTitleH(0.095);
  gStyle->SetTitleSize(0.03,"xyzt");
  gStyle->SetTitleXOffset(1.15);
  gStyle->SetTitleYOffset(2.12);
  gStyle->SetLabelSize(0.03,"XYZ");
  gStyle->SetPadTopMargin(.10);
  gStyle->SetPadLeftMargin(.14);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadBottomMargin(.15);
  gStyle->SetPalette(1);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetFrameBorderMode(0);

  TCanvas *c1 = new TCanvas("c1","c1",700,700);
  c1->Divide(1,1);
  c1->cd(1);
  plane->Draw();

  TLine *d_oriz1 = new TLine(250.,10.,415.,10.);
  d_oriz1->SetLineColor(kGreen);
  TLine *d_oriz2 = new TLine(250.,100.,415.,150.);
  d_oriz2->SetLineColor(kGreen);
  TLine *d_vert1 = new TLine(250.,10.,250.,100.);
  d_vert1->SetLineColor(kGreen);
  TLine *d_vert2 = new TLine(415.,10.,415.,150.);
  d_vert2->SetLineColor(kGreen);
  
  d_oriz1->Draw();
  d_oriz2->Draw();
  d_vert1->Draw();
  d_vert2->Draw();

  TRandom3 rran;
  rran.SetSeed(0);

  Int_t ang_count = 0;

  for(Double_t p=5;p<=20;p=p+5){
    
    if(radiator == 0) momentum = 5;
    if(radiator == 1) momentum = 35.;

    ang = p;
    theta = ang*TMath::Pi()/180.;
    y_track = TMath::Tan(theta)*1000;
    m = TMath::Tan(theta_c+theta);
    mm = TMath::Tan(-theta_c+theta);
    mk = TMath::Tan(theta_ck+theta);
    mmk = TMath::Tan(-theta_ck+theta);

    for(Int_t mom=0;mom<10;mom++){

      TH1F *xpi_dis = new TH1F("xpi_dis","",500000.,0.,500.);
      TH1F *xpi_dis1 = new TH1F("xpi_dis1","",500000.,0.,500.);
      TH1F *xk_dis = new TH1F("xk_dis","",500000.,0.,500.);
      TH1F *xk_dis1 = new TH1F("xk_dis1","",500000.,0.,500.);

      if(radiator == 0){
      theta_c=acos(sqrt(momentum*momentum+0.139570*0.139570)/momentum/1.02);
      theta_ck=acos(sqrt(momentum*momentum+0.493677*0.493677)/momentum/1.02);
      }
      if(radiator == 1){
      theta_c=acos(sqrt(momentum*momentum+0.139570*0.139570)/momentum/1.000482);
      theta_ck=acos(sqrt(momentum*momentum+0.493677*0.493677)/momentum/1.000482);
      }
      m = TMath::Tan(theta_c+theta);
      mm = TMath::Tan(-theta_c+theta);
      mk = TMath::Tan(theta_ck+theta);
      mmk = TMath::Tan(-theta_ck+theta);

    for(Int_t ev=0;ev<20;ev++){
    
      if(radiator == 0) p_start_x[0] = 250.+rran.Uniform(0.,4.);
      if(radiator == 1) p_start_x[0] = 250.+rran.Uniform(4.,140.);
      p_start_y[0] = p_start_x[0]*TMath::Tan(theta);
      q[0] = p_start_y[0]-m*p_start_x[0];
      mv[0]=m;

      if(radiator == 0) p_start_x[1] = 250.+rran.Uniform(0.,4.);
      if(radiator == 1) p_start_x[1] = 250.+rran.Uniform(4.,140.);
      p_start_y[1] = p_start_x[1]*TMath::Tan(theta);
      q[1] = p_start_y[1]-mk*p_start_x[1];
      mv[1]=m;

      p_start_x[2] = p_start_x[0];
      p_start_y[2] = p_start_x[2]*TMath::Tan(theta);
      q[2] = p_start_y[2]-mm*p_start_x[2];
      mv[2]=mm;

      p_start_x[3] = p_start_x[1];
      p_start_y[3] = p_start_x[3]*TMath::Tan(theta);
      q[3] = p_start_y[3]-mmk*p_start_x[3];
      mv[3]=mm;


  ray_tracing *a = new ray_tracing();
  ray_tracing *b = new ray_tracing();
  ray_tracing *d = new ray_tracing();
  ray_tracing *f = new ray_tracing();


  for(Int_t i=0;i<4;i++){

    if(i==0) a->spherical_mirror(xc,yc,r,m,q[i]);
    if(i==1) a->spherical_mirror(xc,yc,r,mk,q[i]);
    if(i==2) a->spherical_mirror(xc,yc,r,mm,q[i]);
    if(i==3) a->spherical_mirror(xc,yc,r,mmk,q[i]);


    if(i==0) b->spherical_mirror1(xc1,yc1,r1,m,q[i]);
    if(i==1) b->spherical_mirror1(xc1,yc1,r1,mk,q[i]);
    if(i==2) b->spherical_mirror1(xc1,yc1,r1,mm,q[i]);
    if(i==3) b->spherical_mirror1(xc1,yc1,r1,mmk,q[i]);

    //cout<<a->xp<<"  "<<a->yp<<"  "<<b->xp_1<<"  "<<b->yp_1<<endl;
    //cout<<a->xm<<"  "<<a->ym<<"  "<<a->xp<<"  "<<a->yp<<endl;
    //cout<<a->phi*180./3.14<<"  "<<a->alpha*180/3.14<<"  "<<a->alpha1*180./3.14<<endl;

    q1[i] = a->yp-TMath::Tan(a->alpha1)*a->xp;
    q1_1[i] = b->yp_1-TMath::Tan(b->alpha1_1)*b->xp_1;
    m1[i] = TMath::Tan(a->alpha1);
    m1_1[i] = TMath::Tan(b->alpha1_1);
    y_up[i] = a->yp;
    y_down[i] = b->yp_1;
    
    Double_t m_rr = (b->yp_1-yc1)/(b->xp_1-xc1); 

    Double_t ang_r = TMath::ATan(TMath::Abs((m-m_rr)/(1.+m*m_rr)))*180./3.14;
    Double_t ang_r1 = TMath::ATan(TMath::Abs((m1_1[i]-m_rr)/(1.+m1_1[i]*m_rr)))*180./3.14;

    //if(i<4) cout<<ang_r<<"  "<<ang_r1<<endl;
    //cout<<m1[i]<<"  "<<q1[i]<<endl;
    //TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
    //c1->Divide(1,1);
    //c1->cd(1);
    
    //Double_t y_sep = 100.;
    Double_t m_ang = TMath::ASin((y_sep-yc)/r)*180./TMath::Pi();
    Double_t m_ang_min = TMath::ASin((0.-yc)/r)*180./TMath::Pi();
    Double_t m_ang_max = TMath::ASin((160.-yc1)/r1)*180./TMath::Pi();
    
    TLine *x_axis = new TLine(0,0,500,0);
    TLine *y_axis = new TLine(0,0,0,500);
    TLine *detector = new TLine(det_pos->x1_det,det_pos->y1_det,det_pos->x2_det,det_pos->y2_det);
    TLine *track = new TLine(0,0,1000,y_track);
    TArc *mirror = new TArc(xc,yc,r,m_ang_min,m_ang);
    TArc *mirror1 = new TArc(xc1,yc1,r1,m_ang,m_ang_max);
    TArc *focus = new TArc(xc,yc,r/2.,-90,90);
    TArc *focus1 = new TArc(xc1,yc1,r1/2.,-90,90);
    mirror->SetFillStyle(0);
    mirror1->SetFillStyle(0);
    focus->SetFillStyle(0);
    focus->SetLineColor(kBlue);
    focus1->SetFillStyle(0);
    focus1->SetLineColor(28);
    TLine *ray1 = new TLine(p_start_x[i],p_start_y[i],a->xp,a->yp);
    TLine *ray2 = new TLine(0,q1[i],a->xp,a->yp);
    TLine *ray2u = new TLine(a->xp,a->yp,500.,m1[i]*500.+q1[i]);
    TLine *radius = new TLine(xc,yc,a->xp,a->yp);
    TLine *ray1_1 = new TLine(p_start_x[i],p_start_y[i],b->xp_1,b->yp_1);
    TLine *ray2_1 = new TLine(0,q1_1[i],b->xp_1,b->yp_1);
    TLine *ray2_1u = new TLine(b->xp_1,b->yp_1,500.,m1_1[i]*500.+q1_1[i]);
    TLine *radius_1 = new TLine(xc1,yc1,b->xp_1,b->yp_1);
    ray1->SetLineColor(kRed);
    ray2->SetLineColor(kBlue);
    if((i==1) || (i==3)) ray2->SetLineColor(kRed);
    ray2u->SetLineColor(kBlue);
    ray1_1->SetLineColor(kRed);
    ray2_1->SetLineColor(28);
    if((i==1) || (i==3)) ray2_1->SetLineColor(kRed);
    ray2_1u->SetLineColor(28);
    detector->SetLineColor(9);
    detector->SetLineWidth(3);
    radius->SetLineStyle(3);
    radius_1->SetLineStyle(3);
    x_axis->Draw();
    y_axis->Draw();
    mirror->Draw("only");
    mirror1->Draw("only");
    if(a->yp<y_sep){
      ray1->Draw();
      if(q1[i]>=0)ray2->Draw();
      if(q1[i]<0)ray2u->Draw();
      radius->Draw();
      focus->Draw("only");
    }
    if(b->yp_1>=y_sep){
      ray1_1->Draw();
      if(q1_1[i]>=0)ray2_1->Draw();
      if(q1_1[i]<0)ray2_1u->Draw();
      radius_1->Draw();
      //focus1->Draw("only");
    }
    track->Draw();
    detector->Draw();

    ayp[i]=a->yp;
    byp[i]=b->yp_1;

  }

  Double_t xpi,xk,x1pi,x1k;
  Double_t ypi,yk,y1pi,y1k;
  
  if(ayp[0]<y_sep) xpi = (q1[0]-q_det)/(m_det-m1[0]);
  if(ayp[1]<y_sep) xk = (q1[1]-q_det)/(m_det-m1[1]);
  if(ayp[2]<y_sep) x1pi = (q1[2]-q_det)/(m_det-m1[2]);
  if(ayp[3]<y_sep) x1k = (q1[3]-q_det)/(m_det-m1[3]);

  if(ayp[0]<y_sep) ypi = m1[0]*((q1[0]-q_det)/(m_det-m1[0]))+q1[0];
  if(ayp[1]<y_sep) yk = m1[1]*((q1[1]-q_det)/(m_det-m1[1]))+q1[1];
  if(ayp[2]<y_sep) y1pi = m1[2]*((q1[2]-q_det)/(m_det-m1[2]))+q1[2];
  if(ayp[3]<y_sep) y1k = m1[3]*((q1[3]-q_det)/(m_det-m1[3]))+q1[3];
  
  
  if(byp[0]>=y_sep) xpi = (q1_1[0]-q_det)/(m_det-m1_1[0]);
  if(byp[1]>=y_sep) xk = (q1_1[1]-q_det)/(m_det-m1_1[1]);
  if(byp[2]>=y_sep) x1pi = (q1_1[2]-q_det)/(m_det-m1_1[2]);
  if(byp[3]>=y_sep) x1k = (q1_1[3]-q_det)/(m_det-m1_1[3]);

  if(byp[0]>=y_sep) ypi = m1_1[0]*((q1_1[0]-q_det)/(m_det-m1_1[0]))+q1_1[0];
  if(byp[1]>=y_sep) yk = m1_1[1]*((q1_1[1]-q_det)/(m_det-m1_1[1]))+q1_1[1];
  if(byp[2]>=y_sep) y1pi = m1_1[2]*((q1_1[2]-q_det)/(m_det-m1_1[2]))+q1_1[2];
  if(byp[3]>=y_sep) y1k = m1_1[3]*((q1_1[3]-q_det)/(m_det-m1_1[3]))+q1_1[3];
  
  spatial_sep[ev] = sqrt((xpi-xk)*(xpi-xk)+(ypi-yk)*(ypi-yk));
  spatial_sep1[ev] = sqrt((x1pi-x1k)*(x1pi-x1k)+(y1pi-y1k)*(y1pi-y1k));
  
  //if(mom==0)cout<<spatial_sep[ev]<<"  "<<spatial_sep1[ev]<<endl;
  //if(mom==0)cout<<xpi<<"  "<<xk<<endl;
  //if(mom==0)cout<<sqrt(x1pi*x1pi+y1pi*y1pi)<<"  "<<sqrt(x1k*x1k+y1k*y1k)<<endl;
  
  xpi_dis->Fill(sqrt(xpi*xpi+ypi*ypi));
  xk_dis->Fill(sqrt(xk*xk+yk*yk));
  xpi_dis1->Fill(sqrt(x1pi*x1pi+y1pi*y1pi));
  xk_dis1->Fill(sqrt(x1k*x1k+y1k*y1k));

  //momentum_points[mom] = momentum;
  //angle_c[mom] = TMath::Abs(theta_c-theta_ck)*1000.;
  //spatial_sep_mean[mom] = spatial_sep_mean[mom]+spatial_sep[ev];
  //spatial_sep_mean1[mom] = spatial_sep_mean1[mom]+spatial_sep1[ev];
    }

    /*TCanvas *p1 = new TCanvas("p1","p1",500,500);
    p1->Divide(2,1);
    p1->cd(1);
    xpi_dis->Draw();
    //xk_dis1->Draw("same");*/

    //spatial_sep_mean[mom]=spatial_sep_mean[mom]/20.;
    //spatial_sep_mean1[mom]=spatial_sep_mean1[mom]/20.;
    momentum_points[mom] = momentum;
    angle_c[mom] = TMath::Abs(theta_c-theta_ck)*1000.;
    spatial_sep_mean[mom]=TMath::Abs(xpi_dis->GetMean()-xk_dis->GetMean());
    spatial_sep_mean1[mom]=TMath::Abs(xpi_dis1->GetMean()-xk_dis1->GetMean());
    spatial_var[mom] = xk_dis->GetRMS();
    spatial_var1[mom] = xk_dis1->GetRMS();

    delete xpi_dis;
    delete xk_dis;
    delete xpi_dis1;
    delete xk_dis1;
  //cout<<spatial_sep_mean[mom]<<"  "<<spatial_sep_mean1[mom]<<endl;

    momentum++;

  }

  Int_t c = 0;
  Double_t x_f, y_f;

  while(c<4){
    for(Int_t i=0+c;i<3;i++){
      if(y_up[c]<y_sep){
	x_f = (q1[i+1]-q1[c])/(m1[c]-m1[i+1]);
        y_f = m1[c]*x_f+q1[c];
	focal_plane->Fill(x_f,y_f);
      }
      if(y_up[c]>=y_sep){
	x_f = (q1_1[i+1]-q1_1[c])/(m1_1[c]-m1_1[i+1]);
	y_f = m1_1[c]*x_f+q1_1[c];
	focal_plane->Fill(x_f,y_f);
      }
      //cout<<i+1<<"  "<<x_f<<"  "<<y_f<<endl;
    }
    for(Int_t i=4+c;i<7;i++){
      if(y_up[c+4]<y_sep){
	x_f = (q1[i+1]-q1[c+4])/(m1[c+4]-m1[i+1]);
	y_f = m1[c+4]*x_f+q1[c+4];
	focal_plane->Fill(x_f,y_f);
      }
      if(y_up[c+4]>=y_sep){
	x_f = (q1_1[i+1]-q1_1[c+4])/(m1_1[c+4]-m1_1[i+1]);
	y_f = m1_1[c+4]*x_f+q1_1[c+4];
	focal_plane->Fill(x_f,y_f);
      }
      //cout<<i+1<<"  "<<x_f<<"  "<<y_f<<endl;
    }
    c++;
  }

  focal_plane->Draw("same" "box");
  c1->Update(); 


  grup[ang_count] = new TGraph(10,momentum_points,spatial_sep_mean);
  grdown[ang_count] = new TGraph(10,momentum_points,spatial_sep_mean1);
  grvup[ang_count] = new TGraph(10,momentum_points,spatial_var);
  grvdown[ang_count] = new TGraph(10,momentum_points,spatial_var1);

  ang_count++;

  }

  TCanvas *c3 = new TCanvas("c3","c3",500,500);
  c3->Divide(2,2);
  c3->cd(1);
  grup[0]->GetXaxis()->SetTitle("momentum [GeV/c]");
  grup[0]->GetYaxis()->SetTitle("#Delta s [cm]");
  grup[0]->Draw("AC*");
  grup[1]->SetLineColor(kRed);
  grup[1]->SetMarkerColor(kRed);
  grup[1]->Draw("C*");
  grup[2]->SetLineColor(kBlue);
  grup[2]->SetMarkerColor(kBlue);
  grup[2]->Draw("C*");
  grup[3]->SetLineColor(kGreen);
  grup[3]->SetMarkerColor(kGreen);
  grup[3]->Draw("C*");
  c3->cd(2);
  grdown[0]->GetXaxis()->SetTitle("momentum [GeV/c]");
  grdown[0]->GetYaxis()->SetTitle("#Delta s [cm]");
  grdown[0]->Draw("AC*");
  grdown[1]->SetLineColor(kRed);
  grdown[1]->SetMarkerColor(kRed);
  grdown[1]->Draw("C*");
  grdown[2]->SetLineColor(kBlue);
  grdown[2]->SetMarkerColor(kBlue);
  grdown[2]->Draw("C*");
  grdown[3]->SetLineColor(kGreen);
  grdown[3]->SetMarkerColor(kGreen);
  grdown[3]->Draw("C*");
  c3->cd(3);
  grvup[0]->GetXaxis()->SetTitle("momentum [GeV/c]");
  grvup[0]->GetYaxis()->SetTitle("#sigma s [cm]");
  grvup[0]->Draw("AC*");
  grvup[1]->SetLineColor(kRed);
  grvup[1]->SetMarkerColor(kRed);
  grvup[1]->Draw("C*");
  grvup[2]->SetLineColor(kBlue);
  grvup[2]->SetMarkerColor(kBlue);
  grvup[2]->Draw("C*");
  grvup[3]->SetLineColor(kGreen);
  grvup[3]->SetMarkerColor(kGreen);
  grvup[3]->Draw("C*");
  c3->cd(4);
  grvdown[0]->GetXaxis()->SetTitle("momentum [GeV/c]");
  grvdown[0]->GetYaxis()->SetTitle("#sigma s [cm]");
  grvdown[0]->Draw("AC*");
  grvdown[1]->SetLineColor(kRed);
  grvdown[1]->SetMarkerColor(kRed);
  grvdown[1]->Draw("C*");
  grvdown[2]->SetLineColor(kBlue);
  grvdown[2]->SetMarkerColor(kBlue);
  grvdown[2]->Draw("C*");
  grvdown[3]->SetLineColor(kGreen);
  grvdown[3]->SetMarkerColor(kGreen);
  grvdown[3]->Draw("C*");


  }


void ray_tracing::trac_fres(Double_t xc, Double_t yc, Double_t r, Double_t xc1, Double_t yc1, Double_t r1, Double_t y_sep, Double_t ang){

  //Double_t xc = 100.;
  //Double_t yc = 100.;
  //Double_t r = 300.;
  //Double_t m = 1.;
  //Double_t q = 0.;

  ray_tracing *det_pos = new ray_tracing();
  det_pos->set_detector();

  Double_t theta = ang*TMath::Pi()/180.;
  Double_t y_track = TMath::Tan(theta)*1000;

  Double_t p_start_x[8];
  Double_t p_start_y[8];
  Double_t q[8], mv[8];
  Double_t q1[8],q1_1[8],m1[8],m1_1[8];

  Double_t y_up[8], y_down[8];

  Double_t momentum =5.;

  Double_t theta_c  =acos(sqrt(momentum*momentum+0.139570*0.139570)/momentum/1.02);
  Double_t m = TMath::Tan(theta_c+theta);
  Double_t mm = TMath::Tan(-theta_c+theta);

  TH2F *plane = new TH2F("plane","",500.,0.,500.,500,0.,500.);
  plane->GetXaxis()->SetTitle("cm");
  plane->GetYaxis()->SetTitle("cm");
  TH2F *focal_plane = new TH2F("plane","",500.,0.,500.,500,0.,500.);
  focal_plane->GetXaxis()->SetTitle("cm");
  focal_plane->GetYaxis()->SetTitle("cm");

  //gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  //gStyle->SetTitleW(0.4);
  //gStyle->SetTitleH(0.095);
  gStyle->SetTitleSize(0.03,"xyzt");
  //gStyle->SetTitleXOffset(.15);
  //gStyle->SetTitleYOffset(.12);
  gStyle->SetLabelSize(0.03,"XYZ");
  gStyle->SetPadTopMargin(.10);
  gStyle->SetPadLeftMargin(.14);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadBottomMargin(.15);
  gStyle->SetPalette(1);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetFrameBorderMode(0);

  TCanvas *c1 = new TCanvas("c1","c1",700,700);
  c1->Divide(1,1);
  c1->cd(1);
  plane->Draw();

  TLine *d_oriz1 = new TLine(250.,10.,415.,10.);
  d_oriz1->SetLineColor(kGreen);
  TLine *d_oriz2 = new TLine(250.,100.,415.,150.);
  d_oriz2->SetLineColor(kGreen);
  TLine *d_vert1 = new TLine(250.,10.,250.,100.);
  d_vert1->SetLineColor(kGreen);
  TLine *d_vert2 = new TLine(415.,10.,415.,150.);
  d_vert2->SetLineColor(kGreen);
  
  d_oriz1->Draw();
  d_oriz2->Draw();
  d_vert1->Draw();
  d_vert2->Draw();

  /*Double_t xxc = 300.;
  Double_t yyc = 50.;
  Double_t rr = 50.;
  Double_t n = 1.000482;
  Double_t n1 = 1.5;
  Double_t d1 = 330.;*/

  /*Double_t xxc = -25.;
  Double_t yyc = -150.;
  Double_t rr = 400.;
  Double_t n = 1.000482;
  Double_t n1 = 1.5;
  Double_t d1 = 280.;*/

  /*Double_t xxc = -45.;
  Double_t yyc = 0.;
  Double_t rr = 350.;
  Double_t n = 1.000482;
  Double_t n1 = 1.5;
  Double_t d1 = 260.;*/

  Double_t xxc = 555.;  //fresnel
  Double_t yyc = 50.;
  Double_t rr = 250.;
  Double_t n = 1.000482;
  Double_t n1 = 1.5;
  Double_t d1 = 5.;

  xxc = 256.+rr;

  /*Double_t xxc = -40.;  //fresnel1
  Double_t yyc = 55.;
  Double_t rr = 300.;
  Double_t n = 1.000482;
  Double_t n1 = 1.5;
  Double_t d1 = 5.;

  xxc = 256.-rr;*/

  /*TLine *diot1 = new TLine(d1,0,d1,100.);
  TArc *diot2 = new TArc(xxc,yyc,rr,90,270);
  diot1->Draw();
  diot2->SetFillStyle(0);
  diot2->Draw("only");*/

  for(Double_t p=5;p<=20;p=p+5){
    
    ang = p;
    theta = ang*TMath::Pi()/180.;
    y_track = TMath::Tan(theta)*1000;
    m = TMath::Tan(theta_c+theta);
    mm = TMath::Tan(-theta_c+theta);

  for(Int_t i=0;i<8;i++){
    if(i<4){
      p_start_x[i] = 250.+ 1.*i;
      p_start_y[i] = p_start_x[i]*TMath::Tan(theta);
      q[i] = p_start_y[i]-m*p_start_x[i];
      mv[i]=m;
    }
    if(i>=4){
      p_start_x[i] = 250.+ 1.*(i-4);
      p_start_y[i] = p_start_x[i]*TMath::Tan(theta);
      q[i] = p_start_y[i]-mm*p_start_x[i];
      mv[i]=mm;
    }
  }

  ray_tracing *a = new ray_tracing();
  ray_tracing *b = new ray_tracing();
  ray_tracing *d = new ray_tracing();
  ray_tracing *f = new ray_tracing();


  for(Int_t i=0;i<8;i++){

      if(i<4) f->fresnel(xxc,yyc,rr,m,q[i],n,n1,d1);                                //fresnel, reversed
      if(i>=4) f->fresnel(xxc,yyc,rr,mm,q[i],n,n1,d1);

      TLine *diot1 = new TLine(xxc-rr+d1,f->sm[4],xxc-rr+d1,f->sp[4]);             
      TArc *diot2 = new TArc(xxc,yyc,rr,f->ang[1],f->ang[0]);
      TArc *diot3 = new TArc(xxc-d1,yyc,rr,f->ang[2],f->ang[1]);
      TArc *diot4 = new TArc(xxc-d1*2,yyc,rr,f->ang[3],f->ang[2]);
      TArc *diot22 = new TArc(xxc,yyc,rr,-(f->ang[0]),-(f->ang[1]));
      TArc *diot33 = new TArc(xxc-d1,yyc,rr,-(f->ang[1]),-(f->ang[2]));
      TArc *diot44 = new TArc(xxc-d1*2,yyc,rr,-(f->ang[2]),-(f->ang[3]));
      TArc *diot5 = new TArc(xxc-d1*3,yyc,rr,f->ang[4],f->ang[3]);
      TArc *diot55 = new TArc(xxc-d1*3,yyc,rr,-(f->ang[4]),-(f->ang[3]));

    /*if(i<4) f->fresnel1(xxc,yyc,rr,m,q[i],n,n1,d1);                                //fresnel1, standard
      if(i>=4) f->fresnel1(xxc,yyc,rr,mm,q[i],n,n1,d1);
      
      TLine *diot1 = new TLine(xxc+rr-d1,f->sm[4],xxc+rr-d1,f->sp[4]);             
      TArc *diot2 = new TArc(xxc,yyc,rr,f->ang[1],f->ang[0]);
      TArc *diot3 = new TArc(xxc+d1,yyc,rr,f->ang[2],f->ang[1]);
      TArc *diot4 = new TArc(xxc+d1*2,yyc,rr,f->ang[3],f->ang[2]);
      TArc *diot22 = new TArc(xxc,yyc,rr,-(f->ang[0]),-(f->ang[1]));
      TArc *diot33 = new TArc(xxc+d1,yyc,rr,-(f->ang[1]),-(f->ang[2]));
      TArc *diot44 = new TArc(xxc+d1*2,yyc,rr,-(f->ang[2]),-(f->ang[3]));
      TArc *diot5 = new TArc(xxc+d1*3,yyc,rr,f->ang[4],f->ang[3]);
      TArc *diot55 = new TArc(xxc+d1*3,yyc,rr,-(f->ang[4]),-(f->ang[3]));*/
 

    diot1->Draw();
    diot2->SetFillStyle(0);
    diot3->SetFillStyle(0);
    diot4->SetFillStyle(0);
    diot22->SetFillStyle(0);
    diot33->SetFillStyle(0);
    diot44->SetFillStyle(0);
    diot5->SetFillStyle(0);
    diot55->SetFillStyle(0);
    diot2->Draw("only");
    diot3->Draw("only");
    diot4->Draw("only");
    diot22->Draw("only");
    diot33->Draw("only");
    diot44->Draw("only");
    diot5->Draw("only");
    diot55->Draw("only");

    //if(i<4) f->diopter1(xxc,yyc,rr,m,q[i],n,n1,d1);
    //if(i>=4) f->diopter1(xxc,yyc,rr,mm,q[i],n,n1,d1);

    //cout<<f->xm_d<<"  "<<f->ym_d<<"  "<<f->m2_d1<<"  "<<f->q2_d1<<endl;

    if(i<4) a->spherical_mirror(xc,yc,r,f->m2_d1,f->q2_d1);
    if(i>=4) a->spherical_mirror(xc,yc,r,f->m2_d1,f->q2_d1);

    if(i<4) b->spherical_mirror1(xc1,yc1,r1,f->m2_d1,f->q2_d1);
    if(i>=4) b->spherical_mirror1(xc1,yc1,r1,f->m2_d1,f->q2_d1);
    
    //if(i<4) cout<<TMath::ATan(m)*180./3.14<<"  "<<q[i]<<"  "<<TMath::ATan(d->m2_d1)*180./3.14<<"  "<<d->q2_d1<<endl;
    //cout<<a->xp<<"  "<<a->yp<<"  "<<b->xp_1<<"  "<<b->yp_1<<endl;
    //cout<<a->xm<<"  "<<a->ym<<"  "<<a->xp<<"  "<<a->yp<<endl;
    //cout<<a->phi*180./3.14<<"  "<<a->alpha*180/3.14<<"  "<<a->alpha1*180./3.14<<endl;

    q1[i] = a->yp-TMath::Tan(a->alpha1)*a->xp;
    q1_1[i] = b->yp_1-TMath::Tan(b->alpha1_1)*b->xp_1;
    m1[i] = TMath::Tan(a->alpha1);
    m1_1[i] = TMath::Tan(b->alpha1_1);
    y_up[i] = a->yp;
    y_down[i] = b->yp_1;
    
    Double_t m_rr = (a->yp-yc)/(a->xp-xc); 

    Double_t ang_r = TMath::ATan(TMath::Abs((f->m2_d1-m_rr)/(1.+f->m2_d1*m_rr)))*180./3.14;
    Double_t ang_r1 = TMath::ATan(TMath::Abs((m1[i]-m_rr)/(1.+m1[i]*m_rr)))*180./3.14;

    if(i<4) cout<<ang_r<<"  "<<ang_r1<<endl;

    //cout<<m1[i]<<"  "<<q1[i]<<endl;
    //TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
    //c1->Divide(1,1);
    //c1->cd(1);
    
    //Double_t y_sep = 100.;
    Double_t m_ang = TMath::ASin((y_sep-yc)/r)*180./TMath::Pi();
    Double_t m_ang_min = TMath::ASin((0.-yc)/r)*180./TMath::Pi();
    Double_t m_ang_max = TMath::ASin((160.-yc1)/r1)*180./TMath::Pi();
    
    TLine *x_axis = new TLine(0,0,500,0);
    TLine *y_axis = new TLine(0,0,0,500);
    TLine *detector = new TLine(det_pos->x1_det,det_pos->y1_det,det_pos->x2_det,det_pos->y2_det);
    TLine *track = new TLine(0,0,1000,y_track);
    TArc *mirror = new TArc(xc,yc,r,m_ang_min,m_ang);
    TArc *mirror1 = new TArc(xc1,yc1,r1,m_ang,m_ang_max);
    TArc *focus = new TArc(xc,yc,r/2.,-90,90);
    TArc *focus1 = new TArc(xc1,yc1,r1/2.,-90,90);
    mirror->SetFillStyle(0);
    mirror1->SetFillStyle(0);
    focus->SetFillStyle(0);
    focus->SetLineColor(kBlue);
    focus1->SetFillStyle(0);
    focus1->SetLineColor(28);
    
    TLine *ray1d = new TLine(p_start_x[i],p_start_y[i],f->xm_d,f->ym_d);                     //fresnel
    TLine *ray2d = new TLine(f->xm_d,f->ym_d,(xxc-rr+d1),f->m1_d1*(xxc-rr+d1)+f->q1_d1);
    TLine *ray1 = new TLine((xxc-rr+d1),f->m2_d1*(xxc-rr+d1)+f->q2_d1,a->xp,a->yp);
    TLine *ray1_1 = new TLine((xxc-rr+d1),f->m2_d1*(xxc-rr+d1)+f->q2_d1,b->xp_1,b->yp_1);
    TLine *fres_center = new TLine(f->xm_d,f->ym_d,f->center,yyc);

    /*TLine *ray1d = new TLine(p_start_x[i],p_start_y[i],f->xr,f->yr);                           //fresnel1
    TLine *ray2d = new TLine(f->xr,f->yr,f->xp_d,f->yp_d);
    TLine *ray1 = new TLine(f->xp_d,f->yp_d,a->xp,a->yp);
    TLine *ray1_1 = new TLine(f->xp_d,f->yp_d,b->xp_1,b->yp_1);
    TLine *fres_center = new TLine(f->center,yyc,f->xp_d,f->yp_d);*/

    /*TLine *ray1d = new TLine(p_start_x[i],p_start_y[i],f->xm_d,f->ym_d);
      TLine *ray2d = new TLine(f->xm_d,f->ym_d,d1,f->m1_d1*d1+f->q1_d1);*/
    //TLine *ray1d_1 = new TLine(p_start_x[i],p_start_y[i],d1,mv[i]*d1+q[i]);
    //TLine *ray2d_1 = new TLine(d1,mv[i]*d1+q[i],d->xp_d,d->yp_d);
    //TLine *ray1 = new TLine(d1,f->m1_d1*d1+f->q1_d1,a->xp,a->yp);
    TLine *ray2 = new TLine(0,q1[i],a->xp,a->yp);
    TLine *radius = new TLine(xc,yc,a->xp,a->yp);
    TLine *ray2_1 = new TLine(0,q1_1[i],b->xp_1,b->yp_1);
    TLine *radius_1 = new TLine(xc1,yc1,b->xp_1,b->yp_1);
    ray1d->SetLineColor(38);
    ray2d->SetLineColor(39);    
    ray1->SetLineColor(kRed);
    ray2->SetLineColor(kBlue);
    ray1_1->SetLineColor(kRed);
    ray2_1->SetLineColor(28);
    detector->SetLineColor(9);
    detector->SetLineWidth(3);
    radius->SetLineStyle(3);
    radius_1->SetLineStyle(3);
    fres_center->SetLineStyle(3);
    x_axis->Draw();
    y_axis->Draw();
    mirror->Draw("only");
    mirror1->Draw("only");
    if(a->yp<y_sep){
      ray1d->Draw();
      ray2d->Draw();
      ray1->Draw();
      ray2->Draw();
      radius->Draw();
      fres_center->Draw();
      //focus->Draw("only");
    }
    if(b->yp_1>=y_sep){
      ray1d->Draw();
      ray2d->Draw();
      ray1_1->Draw();
      ray2_1->Draw();
      radius_1->Draw();
      fres_center->Draw();
      //focus1->Draw("only");
    }
    track->Draw();
    detector->Draw();

  }

  Int_t c = 0;
  Double_t x_f, y_f;

  while(c<4){
    for(Int_t i=0+c;i<3;i++){
      if(y_up[c]<y_sep){
	x_f = (q1[i+1]-q1[c])/(m1[c]-m1[i+1]);
        y_f = m1[c]*x_f+q1[c];
	focal_plane->Fill(x_f,y_f);
      }
      if(y_up[c]>=y_sep){
	x_f = (q1_1[i+1]-q1_1[c])/(m1_1[c]-m1_1[i+1]);
	y_f = m1_1[c]*x_f+q1_1[c];
	focal_plane->Fill(x_f,y_f);
      }
      //cout<<i+1<<"  "<<x_f<<"  "<<y_f<<endl;
    }
    for(Int_t i=4+c;i<7;i++){
      if(y_up[c+4]<y_sep){
	x_f = (q1[i+1]-q1[c+4])/(m1[c+4]-m1[i+1]);
	y_f = m1[c+4]*x_f+q1[c+4];
	focal_plane->Fill(x_f,y_f);
      }
      if(y_up[c+4]>=y_sep){
	x_f = (q1_1[i+1]-q1_1[c+4])/(m1_1[c+4]-m1_1[i+1]);
	y_f = m1_1[c+4]*x_f+q1_1[c+4];
	focal_plane->Fill(x_f,y_f);
      }
      //cout<<i+1<<"  "<<x_f<<"  "<<y_f<<endl;
    }
    c++;
  }

  focal_plane->Draw("same" "box");
  c1->Update(); 

  }

  TArc *focus = new TArc(xc,yc,r/2.,-90,90);
  TArc *focus1 = new TArc(xc1,yc1,r1/2.,-90,90);

  focus->SetFillStyle(0);
  focus->SetLineColor(kBlue);
  focus1->SetFillStyle(0);
  focus1->SetLineColor(28);

   TCanvas *c2 = new TCanvas("c2","c2",700,700);
   c2->Divide(1,1);
   c2->cd(1);
   focal_plane->Draw("BOX");

   d_oriz1->Draw();
   d_oriz2->Draw();
   d_vert1->Draw();
   d_vert2->Draw();

   //focus->Draw("only");
   focus1->Draw("only");

}


void ray_tracing::trac_fres_sep(Double_t xc, Double_t yc, Double_t r, Double_t xc1, Double_t yc1, Double_t r1, Double_t y_sep, Double_t ang, Int_t radiator){

  //Double_t xc = 100.;
  //Double_t yc = 100.;
  //Double_t r = 300.;
  //Double_t m = 1.;
  //Double_t q = 0.;

  ray_tracing *det_pos = new ray_tracing();
  det_pos->set_detector();

  Double_t momentum;

  Double_t theta = ang*TMath::Pi()/180.;
  Double_t y_track = TMath::Tan(theta)*1000;

  Double_t m_det = (det_pos->y2_det-det_pos->y1_det)/(det_pos->x2_det-det_pos->x1_det);
  Double_t q_det = det_pos->y1_det - m_det*det_pos->x1_det;

  Double_t p_start_x[8];
  Double_t p_start_y[8];
  Double_t q[8], mv[8];
  Double_t q1[8],q1_1[8],m1[8],m1_1[8];

  Double_t y_up[8], y_down[8];
  Double_t ayp[8],byp[8];

  Double_t theta_c;
  Double_t theta_ck;
  theta_c=acos(sqrt(momentum*momentum+0.139570*0.139570)/momentum/1.02);
  theta_ck=acos(sqrt(momentum*momentum+0.493677*0.493677)/momentum/1.02);

  Double_t spatial_sep[1000], spatial_sep1[1000], spatial_var[1000], spatial_var1[1000];
  Double_t angle_c[10], spatial_sep_mean[10], spatial_sep_mean1[10],momentum_points[10];

  TGraph *grup[4];
  TGraph *grdown[4];
  TGraph *grvup[4];
  TGraph *grvdown[4];  

  Double_t m = TMath::Tan(theta_c+theta);
  Double_t mm = TMath::Tan(-theta_c+theta);
  Double_t mk = TMath::Tan(theta_ck+theta);
  Double_t mmk = TMath::Tan(-theta_ck+theta);

  TH2F *plane = new TH2F("plane","",500.,0.,500.,500,0.,500.);
  plane->GetXaxis()->SetTitle("cm");
  plane->GetYaxis()->SetTitle("cm");
  TH2F *focal_plane = new TH2F("plane","",500.,0.,500.,500,0.,500.);
  focal_plane->GetXaxis()->SetTitle("cm");
  focal_plane->GetYaxis()->SetTitle("cm");

  //gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  //gStyle->SetTitleW(0.4);
  //gStyle->SetTitleH(0.095);
  gStyle->SetTitleSize(0.03,"xyzt");
  gStyle->SetTitleXOffset(1.15);
  gStyle->SetTitleYOffset(2.12);
  gStyle->SetLabelSize(0.03,"XYZ");
  gStyle->SetPadTopMargin(.10);
  gStyle->SetPadLeftMargin(.14);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadBottomMargin(.15);
  gStyle->SetPalette(1);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetFrameBorderMode(0);

  TCanvas *c1 = new TCanvas("c1","c1",700,700);
  c1->Divide(1,1);
  c1->cd(1);
  plane->Draw();

  TLine *d_oriz1 = new TLine(250.,10.,415.,10.);
  d_oriz1->SetLineColor(kGreen);
  TLine *d_oriz2 = new TLine(250.,100.,415.,150.);
  d_oriz2->SetLineColor(kGreen);
  TLine *d_vert1 = new TLine(250.,10.,250.,100.);
  d_vert1->SetLineColor(kGreen);
  TLine *d_vert2 = new TLine(415.,10.,415.,150.);
  d_vert2->SetLineColor(kGreen);
  
  d_oriz1->Draw();
  d_oriz2->Draw();
  d_vert1->Draw();
  d_vert2->Draw();

  /*Double_t xxc = 300.;
  Double_t yyc = 50.;
  Double_t rr = 50.;
  Double_t n = 1.000482;
  Double_t n1 = 1.5;
  Double_t d1 = 330.;*/

  /*Double_t xxc = -25.;
  Double_t yyc = -150.;
  Double_t rr = 400.;
  Double_t n = 1.000482;
  Double_t n1 = 1.5;
  Double_t d1 = 280.;*/

  /*Double_t xxc = -45.;
  Double_t yyc = 0.;
  Double_t rr = 350.;
  Double_t n = 1.000482;
  Double_t n1 = 1.5;
  Double_t d1 = 260.;*/

  Double_t xxc = 555.;  //fresnel
  Double_t yyc = 55.;
  Double_t rr = 110.;
  Double_t n = 1.000482;
  Double_t n1 = 1.5;
  Double_t d1 = 5.;

  xxc = 256.+rr;

  /*Double_t xxc = -40.;  //fresnel1
  Double_t yyc = 55.;
  Double_t rr = 300.;
  Double_t n = 1.000482;
  Double_t n1 = 1.5;
  Double_t d1 = 5.;

  xxc = 256.-rr;*/

  /*TLine *diot1 = new TLine(d1,0,d1,100.);
  TArc *diot2 = new TArc(xxc,yyc,rr,90,270);
  diot1->Draw();
  diot2->SetFillStyle(0);
  diot2->Draw("only");*/

  TRandom3 rran;
  rran.SetSeed(0);

  Int_t ang_count =0;

  for(Double_t p=5;p<=20;p=p+5){
    
    if(radiator == 0) momentum = 5.;
    if(radiator == 1) momentum = 35.;

    ang = p;
    theta = ang*TMath::Pi()/180.;
    y_track = TMath::Tan(theta)*1000;
    m = TMath::Tan(theta_c+theta);
    mm = TMath::Tan(-theta_c+theta);
    mk = TMath::Tan(theta_ck+theta);
    mmk = TMath::Tan(-theta_ck+theta);

    for(Int_t mom=0;mom<10;mom++){

      TH1F *xpi_dis = new TH1F("xpi_dis","",500000.,0.,500.);
      TH1F *xpi_dis1 = new TH1F("xpi_dis1","",500000.,0.,500.);
      TH1F *xk_dis = new TH1F("xk_dis","",500000.,0.,500.);
      TH1F *xk_dis1 = new TH1F("xk_dis1","",500000.,0.,500.);
      if(radiator == 0){
      theta_c=acos(sqrt(momentum*momentum+0.139570*0.139570)/momentum/1.02);
      theta_ck=acos(sqrt(momentum*momentum+0.493677*0.493677)/momentum/1.02);
      }
      if(radiator == 1){
      theta_c=acos(sqrt(momentum*momentum+0.139570*0.139570)/momentum/1.000482);
      theta_ck=acos(sqrt(momentum*momentum+0.493677*0.493677)/momentum/1.000482);
      }
      m = TMath::Tan(theta_c+theta);
      mm = TMath::Tan(-theta_c+theta);
      mk = TMath::Tan(theta_ck+theta);
      mmk = TMath::Tan(-theta_ck+theta);

      for(Int_t ev=0;ev<20;ev++){
    
      if(radiator == 0) p_start_x[0] = 250.+rran.Uniform(0.,4.);
      if(radiator == 1) p_start_x[0] = 250.+rran.Uniform(4.,140.);
      p_start_y[0] = p_start_x[0]*TMath::Tan(theta);
      q[0] = p_start_y[0]-m*p_start_x[0];
      mv[0]=m;

      if(radiator == 0) p_start_x[1] = 250.+rran.Uniform(0.,4.);
      if(radiator == 1) p_start_x[1] = 250.+rran.Uniform(4.,140.);
      p_start_y[1] = p_start_x[1]*TMath::Tan(theta);
      q[1] = p_start_y[1]-mk*p_start_x[1];
      mv[1]=m;

      p_start_x[2] = p_start_x[0];
      p_start_y[2] = p_start_x[2]*TMath::Tan(theta);
      q[2] = p_start_y[2]-mm*p_start_x[2];
      mv[2]=mm;

      p_start_x[3] = p_start_x[1];
      p_start_y[3] = p_start_x[3]*TMath::Tan(theta);
      q[3] = p_start_y[3]-mmk*p_start_x[3];
      mv[3]=mm;


  ray_tracing *a = new ray_tracing();
  ray_tracing *b = new ray_tracing();
  ray_tracing *f = new ray_tracing();


  for(Int_t i=0;i<4;i++){

      if(i==0) f->fresnel(xxc,yyc,rr,m,q[i],n,n1,d1);                                //fresnel, reversed
      if(i==1) f->fresnel(xxc,yyc,rr,mk,q[i],n,n1,d1);
      if(i==2) f->fresnel(xxc,yyc,rr,mm,q[i],n,n1,d1);
      if(i==3) f->fresnel(xxc,yyc,rr,mmk,q[i],n,n1,d1);

      TLine *diot1 = new TLine(xxc-rr+d1,f->sm[4],xxc-rr+d1,f->sp[4]);             
      TArc *diot2 = new TArc(xxc,yyc,rr,f->ang[1],f->ang[0]);
      TArc *diot3 = new TArc(xxc-d1,yyc,rr,f->ang[2],f->ang[1]);
      TArc *diot4 = new TArc(xxc-d1*2,yyc,rr,f->ang[3],f->ang[2]);
      TArc *diot22 = new TArc(xxc,yyc,rr,-(f->ang[0]),-(f->ang[1]));
      TArc *diot33 = new TArc(xxc-d1,yyc,rr,-(f->ang[1]),-(f->ang[2]));
      TArc *diot44 = new TArc(xxc-d1*2,yyc,rr,-(f->ang[2]),-(f->ang[3]));
      TArc *diot5 = new TArc(xxc-d1*3,yyc,rr,f->ang[4],f->ang[3]);
      TArc *diot55 = new TArc(xxc-d1*3,yyc,rr,-(f->ang[4]),-(f->ang[3]));

    /*if(i<4) f->fresnel1(xxc,yyc,rr,m,q[i],n,n1,d1);                                //fresnel1, standard
      if(i>=4) f->fresnel1(xxc,yyc,rr,mm,q[i],n,n1,d1);
      
      TLine *diot1 = new TLine(xxc+rr-d1,f->sm[4],xxc+rr-d1,f->sp[4]);             
      TArc *diot2 = new TArc(xxc,yyc,rr,f->ang[1],f->ang[0]);
      TArc *diot3 = new TArc(xxc+d1,yyc,rr,f->ang[2],f->ang[1]);
      TArc *diot4 = new TArc(xxc+d1*2,yyc,rr,f->ang[3],f->ang[2]);
      TArc *diot22 = new TArc(xxc,yyc,rr,-(f->ang[0]),-(f->ang[1]));
      TArc *diot33 = new TArc(xxc+d1,yyc,rr,-(f->ang[1]),-(f->ang[2]));
      TArc *diot44 = new TArc(xxc+d1*2,yyc,rr,-(f->ang[2]),-(f->ang[3]));
      TArc *diot5 = new TArc(xxc+d1*3,yyc,rr,f->ang[4],f->ang[3]);
      TArc *diot55 = new TArc(xxc+d1*3,yyc,rr,-(f->ang[4]),-(f->ang[3]));*/
 

    diot1->Draw();
    diot2->SetFillStyle(0);
    diot3->SetFillStyle(0);
    diot4->SetFillStyle(0);
    diot22->SetFillStyle(0);
    diot33->SetFillStyle(0);
    diot44->SetFillStyle(0);
    diot5->SetFillStyle(0);
    diot55->SetFillStyle(0);
    diot2->Draw("only");
    diot3->Draw("only");
    diot4->Draw("only");
    diot22->Draw("only");
    diot33->Draw("only");
    diot44->Draw("only");
    diot5->Draw("only");
    diot55->Draw("only");

    //if(i<4) f->diopter1(xxc,yyc,rr,m,q[i],n,n1,d1);
    //if(i>=4) f->diopter1(xxc,yyc,rr,mm,q[i],n,n1,d1);

    //cout<<f->xm_d<<"  "<<f->ym_d<<"  "<<f->m2_d1<<"  "<<f->q2_d1<<endl;

    a->spherical_mirror(xc,yc,r,f->m2_d1,f->q2_d1);
    b->spherical_mirror1(xc1,yc1,r1,f->m2_d1,f->q2_d1);
    
    
    //if(i<4) cout<<TMath::ATan(m)*180./3.14<<"  "<<q[i]<<"  "<<TMath::ATan(d->m2_d1)*180./3.14<<"  "<<d->q2_d1<<endl;
    //cout<<a->xp<<"  "<<a->yp<<"  "<<b->xp_1<<"  "<<b->yp_1<<endl;
    //cout<<a->xm<<"  "<<a->ym<<"  "<<a->xp<<"  "<<a->yp<<endl;
    //cout<<a->phi*180./3.14<<"  "<<a->alpha*180/3.14<<"  "<<a->alpha1*180./3.14<<endl;

    q1[i] = a->yp-TMath::Tan(a->alpha1)*a->xp;
    q1_1[i] = b->yp_1-TMath::Tan(b->alpha1_1)*b->xp_1;
    m1[i] = TMath::Tan(a->alpha1);
    m1_1[i] = TMath::Tan(b->alpha1_1);
    y_up[i] = a->yp;
    y_down[i] = b->yp_1;
    
    Double_t m_rr = (a->yp-yc)/(a->xp-xc); 

    Double_t ang_r = TMath::ATan(TMath::Abs((f->m2_d1-m_rr)/(1.+f->m2_d1*m_rr)))*180./3.14;
    Double_t ang_r1 = TMath::ATan(TMath::Abs((m1[i]-m_rr)/(1.+m1[i]*m_rr)))*180./3.14;

    //if(i<4) cout<<ang_r<<"  "<<ang_r1<<endl;

    //cout<<m1[i]<<"  "<<q1[i]<<endl;
    //TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
    //c1->Divide(1,1);
    //c1->cd(1);
    
    //Double_t y_sep = 100.;
    Double_t m_ang = TMath::ASin((y_sep-yc)/r)*180./TMath::Pi();
    Double_t m_ang_min = TMath::ASin((0.-yc)/r)*180./TMath::Pi();
    Double_t m_ang_max = TMath::ASin((160.-yc1)/r1)*180./TMath::Pi();
    
    TLine *x_axis = new TLine(0,0,500,0);
    TLine *y_axis = new TLine(0,0,0,500);
    TLine *detector = new TLine(det_pos->x1_det,det_pos->y1_det,det_pos->x2_det,det_pos->y2_det);
    TLine *track = new TLine(0,0,1000,y_track);
    TArc *mirror = new TArc(xc,yc,r,m_ang_min,m_ang);
    TArc *mirror1 = new TArc(xc1,yc1,r1,m_ang,m_ang_max);
    TArc *focus = new TArc(xc,yc,r/2.,-90,90);
    TArc *focus1 = new TArc(xc1,yc1,r1/2.,-90,90);
    mirror->SetFillStyle(0);
    mirror1->SetFillStyle(0);
    focus->SetFillStyle(0);
    focus->SetLineColor(kBlue);
    focus1->SetFillStyle(0);
    focus1->SetLineColor(28);
    
    TLine *ray1d = new TLine(p_start_x[i],p_start_y[i],f->xm_d,f->ym_d);                     //fresnel
    TLine *ray2d = new TLine(f->xm_d,f->ym_d,(xxc-rr+d1),f->m1_d1*(xxc-rr+d1)+f->q1_d1);
    TLine *ray1 = new TLine((xxc-rr+d1),f->m2_d1*(xxc-rr+d1)+f->q2_d1,a->xp,a->yp);
    TLine *ray1_1 = new TLine((xxc-rr+d1),f->m2_d1*(xxc-rr+d1)+f->q2_d1,b->xp_1,b->yp_1);
    TLine *fres_center = new TLine(f->xm_d,f->ym_d,f->center,yyc);

    /*TLine *ray1d = new TLine(p_start_x[i],p_start_y[i],f->xr,f->yr);                           //fresnel1
    TLine *ray2d = new TLine(f->xr,f->yr,f->xp_d,f->yp_d);
    TLine *ray1 = new TLine(f->xp_d,f->yp_d,a->xp,a->yp);
    TLine *ray1_1 = new TLine(f->xp_d,f->yp_d,b->xp_1,b->yp_1);
    TLine *fres_center = new TLine(f->center,yyc,f->xp_d,f->yp_d);*/

    /*TLine *ray1d = new TLine(p_start_x[i],p_start_y[i],f->xm_d,f->ym_d);
      TLine *ray2d = new TLine(f->xm_d,f->ym_d,d1,f->m1_d1*d1+f->q1_d1);*/
    //TLine *ray1d_1 = new TLine(p_start_x[i],p_start_y[i],d1,mv[i]*d1+q[i]);
    //TLine *ray2d_1 = new TLine(d1,mv[i]*d1+q[i],d->xp_d,d->yp_d);
    //TLine *ray1 = new TLine(d1,f->m1_d1*d1+f->q1_d1,a->xp,a->yp);
    TLine *ray2 = new TLine(0,q1[i],a->xp,a->yp);
    TLine *radius = new TLine(xc,yc,a->xp,a->yp);
    TLine *ray2_1 = new TLine(0,q1_1[i],b->xp_1,b->yp_1);
    TLine *radius_1 = new TLine(xc1,yc1,b->xp_1,b->yp_1);
    ray1d->SetLineColor(38);
    ray2d->SetLineColor(39);    
    ray1->SetLineColor(kRed);
    ray2->SetLineColor(kBlue);
    if((i==1) || (i==3)) ray2->SetLineColor(28);
    ray1_1->SetLineColor(kRed);
    ray2_1->SetLineColor(28);
    detector->SetLineColor(9);
    detector->SetLineWidth(3);
    radius->SetLineStyle(3);
    radius_1->SetLineStyle(3);
    fres_center->SetLineStyle(3);
    x_axis->Draw();
    y_axis->Draw();
    mirror->Draw("only");
    mirror1->Draw("only");
    if(a->yp<y_sep){
      ray1d->Draw();
      ray2d->Draw();
      ray1->Draw();
      ray2->Draw();
      radius->Draw();
      fres_center->Draw();
      //focus->Draw("only");
    }
    if(b->yp_1>=y_sep){
      ray1d->Draw();
      ray2d->Draw();
      ray1_1->Draw();
      ray2_1->Draw();
      radius_1->Draw();
      fres_center->Draw();
      //focus1->Draw("only");
    }
    track->Draw();
    detector->Draw();

    ayp[i]=a->yp;
    byp[i]=b->yp_1;

  }

  Double_t xpi,xk,x1pi,x1k;
  Double_t ypi,yk,y1pi,y1k;
  
  if(ayp[0]<y_sep) xpi = (q1[0]-q_det)/(m_det-m1[0]);
  if(ayp[1]<y_sep) xk = (q1[1]-q_det)/(m_det-m1[1]);
  if(ayp[2]<y_sep) x1pi = (q1[2]-q_det)/(m_det-m1[2]);
  if(ayp[3]<y_sep) x1k = (q1[3]-q_det)/(m_det-m1[3]);

  if(ayp[0]<y_sep) ypi = m1[0]*((q1[0]-q_det)/(m_det-m1[0]))+q1[0];
  if(ayp[1]<y_sep) yk = m1[1]*((q1[1]-q_det)/(m_det-m1[1]))+q1[1];
  if(ayp[2]<y_sep) y1pi = m1[2]*((q1[2]-q_det)/(m_det-m1[2]))+q1[2];
  if(ayp[3]<y_sep) y1k = m1[3]*((q1[3]-q_det)/(m_det-m1[3]))+q1[3];
  
  
  if(byp[0]>=y_sep) xpi = (q1_1[0]-q_det)/(m_det-m1_1[0]);
  if(byp[1]>=y_sep) xk = (q1_1[1]-q_det)/(m_det-m1_1[1]);
  if(byp[2]>=y_sep) x1pi = (q1_1[2]-q_det)/(m_det-m1_1[2]);
  if(byp[3]>=y_sep) x1k = (q1_1[3]-q_det)/(m_det-m1_1[3]);

  if(byp[0]>=y_sep) ypi = m1_1[0]*((q1_1[0]-q_det)/(m_det-m1_1[0]))+q1_1[0];
  if(byp[1]>=y_sep) yk = m1_1[1]*((q1_1[1]-q_det)/(m_det-m1_1[1]))+q1_1[1];
  if(byp[2]>=y_sep) y1pi = m1_1[2]*((q1_1[2]-q_det)/(m_det-m1_1[2]))+q1_1[2];
  if(byp[3]>=y_sep) y1k = m1_1[3]*((q1_1[3]-q_det)/(m_det-m1_1[3]))+q1_1[3];
  

  spatial_sep[ev] = TMath::Abs(m1[0]*250.+q1[0]-m1[1]*250.-q1[1]);
  spatial_sep1[ev] = TMath::Abs(m1[2]*250.+q1[2]-m1[3]*250.-q1[3]);
  
  xpi_dis->Fill(sqrt(xpi*xpi+ypi*ypi));
  xk_dis->Fill(sqrt(xk*xk+yk*yk));
  xpi_dis1->Fill(sqrt(x1pi*x1pi+y1pi*y1pi));
  xk_dis1->Fill(sqrt(x1k*x1k+y1k*y1k));

  //if(mom==0)cout<<spatial_sep[ev]<<"  "<<spatial_sep1[ev]<<endl;

  }


    momentum_points[mom] = momentum;
    angle_c[mom] = TMath::Abs(theta_c-theta_ck)*1000.;
    spatial_sep_mean[mom]=TMath::Abs(xpi_dis->GetMean()-xk_dis->GetMean());
    spatial_sep_mean1[mom]=TMath::Abs(xpi_dis1->GetMean()-xk_dis1->GetMean());
    spatial_var[mom] = xk_dis->GetRMS();
    spatial_var1[mom] = xk_dis1->GetRMS();

    delete xpi_dis;
    delete xk_dis;
    delete xpi_dis1;
    delete xk_dis1;

    //cout<<spatial_sep_mean[mom]<<"  "<<spatial_sep_mean1[mom]<<endl;

    momentum++;

    }

  Int_t c = 0;
  Double_t x_f, y_f;

  while(c<4){
    for(Int_t i=0+c;i<3;i++){
      if(y_up[c]<y_sep){
	x_f = (q1[i+1]-q1[c])/(m1[c]-m1[i+1]);
        y_f = m1[c]*x_f+q1[c];
	focal_plane->Fill(x_f,y_f);
      }
      if(y_up[c]>=y_sep){
	x_f = (q1_1[i+1]-q1_1[c])/(m1_1[c]-m1_1[i+1]);
	y_f = m1_1[c]*x_f+q1_1[c];
	focal_plane->Fill(x_f,y_f);
      }
      //cout<<i+1<<"  "<<x_f<<"  "<<y_f<<endl;
    }
    for(Int_t i=4+c;i<7;i++){
      if(y_up[c+4]<y_sep){
	x_f = (q1[i+1]-q1[c+4])/(m1[c+4]-m1[i+1]);
	y_f = m1[c+4]*x_f+q1[c+4];
	focal_plane->Fill(x_f,y_f);
      }
      if(y_up[c+4]>=y_sep){
	x_f = (q1_1[i+1]-q1_1[c+4])/(m1_1[c+4]-m1_1[i+1]);
	y_f = m1_1[c+4]*x_f+q1_1[c+4];
	focal_plane->Fill(x_f,y_f);
      }
      //cout<<i+1<<"  "<<x_f<<"  "<<y_f<<endl;
    }
    c++;
  }

  focal_plane->Draw("same" "box");
  c1->Update(); 

  grup[ang_count] = new TGraph(10,momentum_points,spatial_sep_mean);
  grdown[ang_count] = new TGraph(10,momentum_points,spatial_sep_mean1);
  grvup[ang_count] = new TGraph(10,momentum_points,spatial_var);
  grvdown[ang_count] = new TGraph(10,momentum_points,spatial_var1);

  ang_count++;

  }

  TCanvas *c3 = new TCanvas("c3","c3",500,500);
  c3->Divide(2,2);
  c3->cd(1);
  grup[0]->GetXaxis()->SetTitle("momentum [GeV/c]");
  grup[0]->GetYaxis()->SetTitle("#Delta s [cm]");
  grup[0]->Draw("AC*");
  grup[1]->SetLineColor(kRed);
  grup[1]->SetMarkerColor(kRed);
  grup[1]->Draw("C*");
  grup[2]->SetLineColor(kBlue);
  grup[2]->SetMarkerColor(kBlue);
  grup[2]->Draw("C*");
  grup[3]->SetLineColor(kGreen);
  grup[3]->SetMarkerColor(kGreen);
  grup[3]->Draw("C*");
  c3->cd(2);
  grdown[0]->GetXaxis()->SetTitle("momentum [GeV/c]");
  grdown[0]->GetYaxis()->SetTitle("#Delta s [cm]");
  grdown[0]->Draw("AC*");
  grdown[1]->SetLineColor(kRed);
  grdown[1]->SetMarkerColor(kRed);
  grdown[1]->Draw("C*");
  grdown[2]->SetLineColor(kBlue);
  grdown[2]->SetMarkerColor(kBlue);
  grdown[2]->Draw("C*");
  grdown[3]->SetLineColor(kGreen);
  grdown[3]->SetMarkerColor(kGreen);
  grdown[3]->Draw("C*");
  c3->cd(3);
  grvup[0]->GetXaxis()->SetTitle("momentum [GeV/c]");
  grvup[0]->GetYaxis()->SetTitle("#sigma s [cm]");
  grvup[0]->Draw("AC*");
  grvup[1]->SetLineColor(kRed);
  grvup[1]->SetMarkerColor(kRed);
  grvup[1]->Draw("C*");
  grvup[2]->SetLineColor(kBlue);
  grvup[2]->SetMarkerColor(kBlue);
  grvup[2]->Draw("C*");
  grvup[3]->SetLineColor(kGreen);
  grvup[3]->SetMarkerColor(kGreen);
  grvup[3]->Draw("C*");
  c3->cd(4);
  grvdown[0]->GetXaxis()->SetTitle("momentum [GeV/c]");
  grvdown[0]->GetYaxis()->SetTitle("#sigma s [cm]");
  grvdown[0]->Draw("AC*");
  grvdown[1]->SetLineColor(kRed);
  grvdown[1]->SetMarkerColor(kRed);
  grvdown[1]->Draw("C*");
  grvdown[2]->SetLineColor(kBlue);
  grvdown[2]->SetMarkerColor(kBlue);
  grvdown[2]->Draw("C*");
  grvdown[3]->SetLineColor(kGreen);
  grvdown[3]->SetMarkerColor(kGreen);
  grvdown[3]->Draw("C*");
}
