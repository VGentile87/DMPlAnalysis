#include "DMPlsFuncImages.h"
#include "TMatrixDEigen.h"
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

TVectorD DMPlsFuncImages::fit_ellipse(TGraph *g)
{
  TVectorD ellipse;
  
  if (!g) return ellipse; // just a precaution
  if (g->GetN() < 6) return ellipse; // just a precaution
  
  Int_t i;
  Double_t tmp;
  
  Int_t N = g->GetN();
  Double_t xmin, xmax, ymin, ymax, X0, Y0;
  g->ComputeRange(xmin, ymin, xmax, ymax);
#if 1 /* 0 or 1 */
  X0 = (xmax + xmin) / 2.0;
  Y0 = (ymax + ymin) / 2.0;
#else /* 0 or 1 */
  X0 = Y0 = 0.0;
#endif /* 0 or 1 */

  TMatrixD D1(N, 3); // quadratic part of the design matrix
  TMatrixD D2(N, 3); // linear part of the design matrix

  for (i = 0; i < N; i++) {
    Double_t x = (g->GetX())[i] - X0;
    Double_t y = (g->GetY())[i] - Y0;
    D1[i][0] = x * x;
    D1[i][1] = x * y;
    D1[i][2] = y * y;
    D2[i][0] = x;
    D2[i][1] = y;
    D2[i][2] = 1.0;
  }

  // quadratic part of the scatter matrix
  TMatrixD S1(TMatrixD::kAtA, D1);
  // combined part of the scatter matrix
  TMatrixD S2(D1, TMatrixD::kTransposeMult, D2);
  // linear part of the scatter matrix
  TMatrixD S3(TMatrixD::kAtA, D2);
  S3.Invert(&tmp); S3 *= -1.0;
  if (tmp == 0.0) {
    std::cout << "fit_ellipse : linear part of the scatter matrix is singular!" << std::endl;
    return ellipse;
  }
  // for getting a2 from a1
  TMatrixD T(S3, TMatrixD::kMultTranspose, S2);
  // reduced scatter matrix
  TMatrixD M(S2, TMatrixD::kMult, T); M += S1;
  // premultiply by inv(C1)
  for (i = 0; i < 3; i++) {
    tmp = M[0][i] / 2.0;
    M[0][i] = M[2][i] / 2.0;
    M[2][i] = tmp;
    M[1][i] *= -1.0;
  }
  // solve eigensystem
  TMatrixDEigen eig(M); // note: eigenvectors are not normalized
  const TMatrixD &evec = eig.GetEigenVectors();
  // const TVectorD &eval = eig.GetEigenValuesRe();
  if ((eig.GetEigenValuesIm()).Norm2Sqr() != 0.0) {
    std::cout << "fit_ellipse : eigenvalues have nonzero imaginary parts!" << std::endl;
    return ellipse;
  }
  // evaluate a’Ca (in order to find the eigenvector for min. pos. eigenvalue)
  for (i = 0; i < 3; i++) {
    tmp = 4.0 * evec[0][i] * evec[2][i] - evec[1][i] * evec[1][i];
    if (tmp > 0.0) break;
  }
  if (i > 2) {
    std::cout << "fit_ellipse : no min. pos. eigenvalue found!" << std::endl;
    // i = 2;
    return ellipse;
  }
  // eigenvector for min. pos. eigenvalue
  TVectorD a1(TMatrixDColumn_const(evec, i));
  tmp = a1.Norm2Sqr();
  if (tmp > 0.0) {
    a1 *= 1.0 / std::sqrt(tmp); // normalize this eigenvector
  } else {
    std::cout << "fit_ellipse : eigenvector for min. pos. eigenvalue is NULL!" << std::endl;
    return ellipse;
  }
  TVectorD a2(T * a1);

  // ellipse coefficients
  ellipse.ResizeTo(8);
  ellipse[0] = X0; // "X0"
  ellipse[1] = Y0; // "Y0"
  ellipse[2] = a1[0]; // "A"
  ellipse[3] = a1[1]; // "B"
  ellipse[4] = a1[2]; // "C"
  ellipse[5] = a2[0]; // "D"
  ellipse[6] = a2[1]; // "E"
  ellipse[7] = a2[2]; // "F"

  return ellipse;
}


TVectorD DMPlsFuncImages::ConicToParametric(const TVectorD &conic)
{
  TVectorD ellipse;

  if (conic.GetNrows() != 8) {
    std::cout << "ConicToParametric : improper input vector length!" << std::endl;
    return ellipse;
  }

  Double_t a, b, theta;
  Double_t x0 = conic[0]; // = X0
  Double_t y0 = conic[1]; // = Y0

  // http://mathworld.wolfram.com/Ellipse.html
  Double_t A = conic[2];
  Double_t B = conic[3] / 2.0;
  Double_t C = conic[4];
  Double_t D = conic[5] / 2.0;
  Double_t F = conic[6] / 2.0;
  Double_t G = conic[7];

  Double_t J = B * B - A * C;
  Double_t Delta = A * F * F + C * D * D + J * G - 2.0 * B * D * F;
  Double_t I = - (A + C);

  // http://mathworld.wolfram.com/QuadraticCurve.html
  if (!( (Delta != 0.0) && (J < 0.0) && (I != 0.0) && (Delta / I < 0.0) )) {
    std::cout << "ConicToParametric : ellipse (real) specific constraints not met!" << std::endl;
    return ellipse;
  }

  x0 += (C * D - B * F) / J;
  y0 += (A * F - B * D) / J;

  Double_t tmp = std::sqrt((A - C) * (A - C) + 4.0 * B * B);
  a = std::sqrt(2.0 * Delta / J / (I + tmp));
  b = std::sqrt(2.0 * Delta / J / (I - tmp));

  theta = 0.0;
  if (B != 0.0) {
    tmp = (A - C) / 2.0 / B;
    theta = -45.0 * (std::atan(tmp) / TMath::PiOver2());
    if (tmp < 0.0) { theta -= 45.0; } else { theta += 45.0; }
    if (A > C) theta += 90.0;
  } else if (A > C) theta = 90.0;

  // try to keep "a" > "b"
  if (a < b) { tmp = a; a = b; b = tmp; theta -= 90.0; }
  // try to keep "theta" = -45 ... 135 degrees
  if (theta < -45.0) theta += 180.0;
  if (theta > 135.0) theta -= 180.0;

  // ellipse coefficients
  ellipse.ResizeTo(5);
  ellipse[0] = x0; // ellipse's "x" center
  ellipse[1] = y0; // ellipse's "y" center
  ellipse[2] = a; // ellipse's "semimajor" axis along "x"
  ellipse[3] = b; // ellipse's "semiminor" axis along "y"
  ellipse[4] = theta; // ellipse's axes rotation angle (in degrees)

  return ellipse;
}

TGraph * DMPlsFuncImages::TestGraphDLSF(Bool_t randomize = kFALSE) {
  Int_t i;

  // define the test ellipse
  Double_t x0 = 4; // ellipse's "x" center
  Double_t y0 = 3; // ellipse's "y" center
  Double_t a = 2; // ellipse's "semimajor" axis along "x" (> 0)
  Double_t b = 1; // ellipse's "semiminor" axis along "y" (> 0)
  Double_t theta = 100; // ellipse's axes rotation angle (-45 ... 135 degrees)

  // gRandom->SetSeed(0);
  if (randomize) {
    x0 = 10.0 - 20.0 * gRandom->Rndm();
    y0 = 10.0 - 20.0 * gRandom->Rndm();
    a = 0.5 + 4.5 * gRandom->Rndm();
    b = 0.5 + 4.5 * gRandom->Rndm();
    theta = 180.0 - 360.0 * gRandom->Rndm();
  }

  const Int_t n = 100; // number of points
  Double_t x[n], y[n];
  Double_t dt = TMath::TwoPi() / Double_t(n);
  Double_t tmp;
  theta *= TMath::PiOver2() / 90.0; // degrees -> radians
  for (i = 0; i < n; i++) {
    x[i] = a * (std::cos(dt * Double_t(i)) + 0.1 * gRandom->Rndm() - 0.05);
    y[i] = b * (std::sin(dt * Double_t(i)) + 0.1 * gRandom->Rndm() - 0.05);
    // rotate the axes
    tmp = x[i];
    x[i] = x[i] * std::cos(theta) - y[i] * std::sin(theta);
    y[i] = y[i] * std::cos(theta) + tmp * std::sin(theta);
    // shift the center
    x[i] += x0;
    y[i] += y0;
  }

  // create the test TGraph
  TGraph *g = ((TGraph *)(gROOT->FindObject("g")));
  if (g) delete g;
  g = new TGraph(n, x, y);
  g->SetNameTitle("g", "test ellipse");

  return g;
}

//
// "ROOT Script" entry point (the same name as the "filename's base")
//
void DMPlsFuncImages::fitEllipseTGraphDLSF(TGraph *g = ((TGraph *)0), double *fit_info=0)
{
  if (!g) g = TestGraphDLSF(kTRUE); // create a "random" ellipse

  // fit the TGraph
  TVectorD conic = DMPlsFuncImages::fit_ellipse(g);
  TVectorD ellipse = ConicToParametric(conic);

#if 1 // 0 or 1 
  if ( ellipse.GetNrows() == 5 ) {
    std::cout << std::endl;
    std::cout << "x0 = " << ellipse[0] << std::endl;
    std::cout << "y0 = " << ellipse[1] << std::endl;
    std::cout << "a = " << ellipse[2] << std::endl;
    std::cout << "b = " << ellipse[3] << std::endl;
    std::cout << "theta = " << ellipse[4]*TMath::Pi()/180. << std::endl;
    std::cout << std::endl;

    fit_info[0]=ellipse[3]*2*0.02752;
    fit_info[1]=ellipse[2]*2*0.02752;
    fit_info[2]=ellipse[4]*TMath::Pi()/180.;
    if(ellipse[4]>90 && ellipse[4]<270)fit_info[2]=fit_info[2]-TMath::Pi();
    if(ellipse[4]<-90 && ellipse[4]>-270)fit_info[2]=fit_info[2]+TMath::Pi();
    if(ellipse[4]>=270)fit_info[2]=fit_info[2]-2*TMath::Pi();
    if(ellipse[4]<=-270)fit_info[2]=fit_info[2]+2*TMath::Pi();
  }
#endif // 0 or 1 

#if 1 // 0 or 1 
  // draw everything
  TCanvas * c = ((TCanvas *)(gROOT->GetListOfCanvases()->FindObject("c")));
  if (c) { c->Clear(); } else { c = new TCanvas("c", "c"); }
  c->SetGrid(1, 1);
  g->Draw("A*");
  g->GetXaxis()->SetLimits(0,41);
  g->GetYaxis()->SetRangeUser(0,41);
  if ( ellipse.GetNrows() == 5 ) {
    TEllipse *e = new TEllipse(ellipse[0], ellipse[1], // "x0", "y0"
                               ellipse[2], ellipse[3], // "a", "b"
                               0, 360,
                               ellipse[4]); // "theta" (in degrees)
    e->SetFillStyle(0); // hollow
    e->SetLineColor(2);
    e->SetLineWidth(2);
    e->Draw();
  }
  c->Modified(); c->Update(); // make sure it's really drawn
#endif // 0 or 1 
  c->SaveAs("fitEll.root");
  return;
}

// end of file fitEllipseTGraphDLSF.cxx by Silesius Anonymus

TH2F * DMPlsFuncImages::merged_histo(DMRRun *aRun, DMRView *view, int ihd, int igr, int *clset)
{

  const int dr=50;
  float pixX = aRun->GetHeader()->pixX;
  float pixY = aRun->GetHeader()->pixY;
  int nx     = aRun->GetHeader()->npixX;
  int ny     = aRun->GetHeader()->npixY;
  int x0,y0;

  int tmp_max=0;
  int tmp_min=255;

  int tmp_max2=-255;
  int tmp_min2=255;

  
  const int ncp=8;


  TH2F *hsum = new TH2F();
  hsum=0;
  TH2F *h[ncp];
  DMRImage imin; imin.SetImage( 2*dr, 2*dr );
  DMRImage imax; imax.SetImage( 2*dr, 2*dr );

  for(int i=0;i<ncp;i++){
    h[i]= new TH2F("","",2*dr,0,dr,2*dr,0,dr);
    if(clset[i]!=-1){
      DMRImageCl *im    = aRun->GetCLIM(ihd,clset[i],i,dr);
      if(im){
	h[i] = im->GetHist2();
	if(h[i]->GetMaximum()>tmp_max)tmp_max = h[i]->GetMaximum();
	if(h[i]->GetMinimum()<tmp_min)tmp_min = h[i]->GetMinimum();
      }
      delete im;
    }
  }
  
  for(int i=0;i<ncp;i++){

    if(clset[i]!=-1){
      
      DMRViewHeader   *hd = view->GetHD();
      DMRCluster      *cl = view->GetCL(clset[i]);
      DMRFrame        *frcl = view->GetFR(cl->ifr);       
      DMRFrame        *fr  = view->GetFR(frcl->iz,frcl->ipol);
      
      if(i==0){
	x0 = ((cl->x+hd->x)-fr->x)/pixX + nx/2;
	y0 = ((cl->y+hd->y)-fr->y)/pixY + ny/2;
      }
    
      DMRImageCl *im    = aRun->GetCLIMBFC(ihd,clset[i],i,dr,x0,y0);
      if(im){
	h[i] = im->GetHist2();
	//cout << i << " " << h[i]->GetMaximum() << endl;
	if(!hsum){
	  hsum = (TH2F*)(h[i]->Clone("hsum"));
	  hsum->SetTitle(Form("Sum over all polarizations"));
	  hsum->SetName(Form("hsum"));
	  imin.Max(*im);
	  imax.Max(*im);
	}
	else {
	  hsum->Add(h[i]);
	  imin.Min(*im);
	  imax.Max(*im);
	  //cout << "sum "<< hsum->GetMaximum() << endl;
	}
      }
      delete im;
      //delete hd;
      //delete cl;
      //delete frcl;
      //delete fr;
    }
    delete h[i];
  }
  hsum->Scale(1./8.);
  //delete view;
  return hsum;
  delete hsum;
}

void DMPlsFuncImages::gr_spectrum(TH2F *h, double *gr_scale_result)
{
  int nbinsX;
  int nbinsY; 
  double cont=0;
  Double_t bkg_mean=0;
  Double_t sig_peak=0;

  TH1F *hspect = new TH1F("hspect","hs",255,0,255);

  nbinsX = h->GetNbinsX();
  nbinsY = h->GetNbinsY();
    
  for(int i=0;i<nbinsX;i++){
    for(int j=0;j<nbinsY;j++){
      cont = h->GetBinContent(i+1,j+1);
      if(cont>0)hspect->Fill(cont);
    }
  }
  
  if(hspect->GetEntries()!=0){
    hspect->Fit("gaus","Q");
    TF1 *gaus = hspect->GetFunction("gaus");
    bkg_mean = gaus->GetParameter(1);
    sig_peak = h->GetMaximum();
  }
  else {
    bkg_mean=-1;
    sig_peak=-1;
  }
  
  gr_scale_result[0]=bkg_mean;
  gr_scale_result[1]=sig_peak;
  
  //delete h;
  delete hspect;
}


void DMPlsFuncImages::scaling(DMRRun *aRun, int ihd, int icl, int ipol, double *scale_result)
{

  double br[6] = {20,35,50,65,80,95};
  double r0[6] = {1.076,1.157,1.185,1.184,1.215,1.229};
  double r1[6] = {1.043,1.103,1.084,1.086,1.106,1.083};
  double r2[6] = {0.996,1.013,0.962,0.947,0.973,0.973};
  double r3[6] = {0.953,0.884,0.862,0.878,0.878,0.878};
  double r4[6] = {0.947,0.886,0.869,0.872,0.889,0.889};
  double r5[6] = {0.977,0.980,0.947,0.963,0.961,0.961};
  double r6[6] = {1.024,1.077,1.075,1.082,1.109,1.082};
  double r7[6] = {1.064,1.140,1.150,1.174,1.206,1.237};

  double brM[8][6];
 
  for(int j=0;j<6;j++){
    brM[0][j]=r0[j];
    brM[1][j]=r1[j];
    brM[2][j]=r2[j];
    brM[3][j]=r3[j];
    brM[4][j]=r4[j];
    brM[5][j]=r5[j];
    brM[6][j]=r6[j];
    brM[7][j]=r7[j];
  }

  
  int nbinsX;
  int nbinsY; 
  double cont=0;
  Double_t bkg_mean=0;
  Double_t sig_peak=0;
  double cont_S=0;
  Double_t bkg_mean_S=0;
  Double_t sig_peak_S=0;
  //float w[8] = {1,0.957,0.909,0.881,0.887,0.918,0.962,0.996};
  int dr=50;
  TH2F *h=0;
  TH2F *hs=0;
  TH1F *hspect = new TH1F("hspect","hs",255,0,255);
  TH1F *hspect_S = new TH1F("hspect_S","hs_S",255,0,255);
  hspect->Reset(0);
  hspect_S->Reset(0);
  DMRImageCl *im    = aRun->GetCLIM(ihd,icl,ipol,dr);
  
  if(im){
    h = im->GetHist2();
    hs = im->GetHist2();    

    for(int k=0;k<6;k++){
      if(k==0){
	if(hs->GetMaximum()<br[1])hs->Scale(1./brM[ipol][0]);
      }
      if(k>0 && k<5){
	if(hs->GetMaximum()>=br[k] && hs->GetMaximum()<br[k+1])hs->Scale(1./brM[ipol][k]);
      }
      if(k==5){
	if(hs->GetMaximum()>=br[5])hs->Scale(1./brM[ipol][5]);
      }
    }

    nbinsX = h->GetNbinsX();
    nbinsY = h->GetNbinsY();
    
    //if(icl==168) cout <<"hello1 "<< scale_result[0] << " " << scale_result[1] << " " << scale_result[2] << " " << scale_result[3] << endl;
    for(int i=0;i<nbinsX;i++){
      for(int j=0;j<nbinsY;j++){
	cont = h->GetBinContent(i+1,j+1);
	if(cont>0)hspect->Fill(cont);
	cont_S = hs->GetBinContent(i+1,j+1);
	if(cont_S>0)hspect_S->Fill(cont_S);
      }
    }
    
    if(hspect->GetEntries()!=0){
      hspect->Fit("gaus","Q");
      TF1 *gaus = hspect->GetFunction("gaus");
      bkg_mean = gaus->GetParameter(1);
      sig_peak = h->GetMaximum();
    }
    else {
      bkg_mean=-1;
      sig_peak=-1;
    }
    if(hspect_S->GetEntries()!=0){
      hspect_S->Fit("gaus","Q");
      TF1 *gaus_S = hspect_S->GetFunction("gaus");
      bkg_mean_S = gaus_S->GetParameter(1);
      sig_peak_S = hs->GetMaximum();
    }
    else {
      bkg_mean_S=-1;
      sig_peak_S=-1;
    }
  }
  else{
    bkg_mean=-1;
    sig_peak=-1;
    bkg_mean_S=-1;
    sig_peak_S=-1;
  }

  scale_result[0]=bkg_mean;
  scale_result[1]=sig_peak;
  scale_result[2]=bkg_mean_S;
  scale_result[3]=sig_peak_S;

    
  
  delete im;  
  delete h;
  delete hspect;
  delete hs;
  delete hspect_S;
  
}

float DMPlsFuncImages::maxpeak(DMRRun *aRun, int ihd, int icl, int ipol)
{
  int dr=15;
  TH2F *h=0;
  //cout << iv << " " << icl << " " << ipol << endl;
  DMRImageCl *im    = aRun->GetCLIM(ihd,icl,ipol,dr);
  if(im){
    h = im->GetHist2();
    float maxpeak = h->GetMaximum();
    return maxpeak;
  }
  else return -1;
}



void DMPlsFuncImages::images_info(DMRRun *aRun,int ihd, int icl, int ipol, double *im_info)
{
  //cout << ihd << " " << icl <<" " << ipol << endl;
  int nbinsX=0;
  int nbinsY=0;
  double cont=0;
  Double_t bkg_mean=0;
  Double_t sig_peak=0;
  
  //TCanvas *c1 = new TCanvas();
  TH1F *hspect = new TH1F("hspect","hs",255,0,255);
  TF1 *gaus = new TF1("fit_gaus","gaus");
  gaus->SetParameter(1,23);
  gaus->SetParameter(2,1);

  int dr=15;
  TH2F *h = new TH2F();
  //cout << iv << " " << icl << " " << ipol << endl;
  DMRImageCl *im    = new DMRImageCl();
  im = aRun->GetCLIM(ihd,icl,ipol,dr);
  if(im){
    h = im->GetHist2();
    nbinsX = h->GetNbinsX();
    nbinsY = h->GetNbinsY();
    
    for(int i=0;i<nbinsX;i++){
      for(int j=0;j<nbinsY;j++){
	cont = h->GetBinContent(i+1,j+1);
	if(cont>0)hspect->Fill(cont);
      }
    }

    
  if(hspect->GetEntries()!=0){
    double max_pos = hspect->GetBinCenter(hspect->GetMaximumBin());
    hspect->Fit(gaus,"Q","",max_pos-5,max_pos+5);
    //TF1 *gaus = hspect->GetFunction("gaus");
    bkg_mean = gaus->GetParameter(1);
    sig_peak = h->GetMaximum();
  }
  else {
    bkg_mean=-1;
    sig_peak=-1;
  }

  }
  else{
    bkg_mean=-1;
    sig_peak=-1;
  }

  im_info[0]=bkg_mean;
  im_info[1]=sig_peak;
  delete hspect;
  delete im;
  delete gaus;
  delete h;
  //if (c1) { c1->Close(); gSystem->ProcessEvents(); }
  
}




int DMPlsFuncImages::otsu_method(float *histogram, long int total_pixels) {
    double probability[256], mean[256];
    double max_between, between[256];
    int threshold;

    
    //probability = class probability
    //mean = class mean
    //between = between class variance
    

    for(int i = 0; i < 256; i++) {
        probability[i] = 0.0;
        mean[i] = 0.0;
        between[i] = 0.0;
    }

    probability[0] = histogram[0];

    for(int i = 1; i < 256; i++) {
        probability[i] = probability[i - 1] + histogram[i];
        mean[i] = mean[i - 1] + i * histogram[i];
    }

    threshold = 0;
    max_between = 0.0;

    for(int i = 0; i < 255; i++) {
        if(probability[i] != 0.0 && probability[i] != 1.0)
            between[i] = pow(mean[255] * probability[i] - mean[i], 2) / (probability[i] * (1.0 - probability[i]));
    else
        between[i] = 0.0;
        if(between[i] > max_between) {
            max_between = between[i];
            threshold = i;
        }
    }

    return threshold;
}

void  DMPlsFuncImages::images_ellipse(DMRRun *aRun,DMRView* view, int ihd, int icl, int ipol, double *fit_info, TCanvas *c1)
{
  gStyle->SetOptFit(1111);
  //cout << ihd << " " << icl <<" " << ipol << endl;
  int nbinsX=0;
  int nbinsY=0;
  int cont=0;
  Double_t bkg_mean=0;
  Double_t bkg_sigma=0;
  Double_t sig_peak=0;
  int locmax,locmay,locmaz;
  int total_pixels=0;
  float histogram[256]={};
  float occurrence[256]={};
  float threshold_value=0;
  float pixX=aRun->GetHeader()->pixX;
  float pixY=aRun->GetHeader()->pixY;
  int   nx   = aRun->GetHeader()->npixX;
  int   ny   = aRun->GetHeader()->npixY;
  double fit_result[3]={};
  vector<double> dist_bar;
  TRandom3 *r = new TRandom3();

  view = aRun->GetEntry( ihd,1,1,1,1,1,1 );
  DMRViewHeader   *hd = view->GetHD(); 
  DMRCluster      *cl = view->GetCL(icl);             
  DMRFrame        *frcl = view->GetFR(cl->ifr);      
  DMRFrameRaw     *frr = aRun->GetFrameRaw( Form("ventr==%d&&iz==%d&&ipol==%d&&flags==0",ihd,frcl->iz,ipol));
  DMRFrame        *fr  = view->GetFR(frr->iz,frr->ipol);

  //TF2 *bigaus = new TF2("fit_bigaus","[0]*x*x+2*[1]*x*y+2*[2]*x+2*[3]*y+[4]");
  //TF2 *bigaus = new TF2("fit_bigaus","bigaus+[6]");
  TF2 *bigaus = new TF2("bigaus","bigaus",10,30,10,30);
  TF2 *biconst = new TF2("fit_biconst","[6]");
  TH1F *hspect = new TH1F("hspect","hs",255,0,255);
  TF1 *gaus = new TF1("fit_gaus","gaus");
  gaus->SetParameter(1,23);
  gaus->SetParameter(2,1);
  int dr=20;
  TH2F *h = new TH2F();
  TH2F *hh = new TH2F();
  TH2F *hhh = new TH2F();
  TH1F *hbar = new TH1F("hbar","hbar",60,0,30);
  TGraph *gr = new TGraph();
  //cout << iv << " " << icl << " " << ipol << endl;
  DMRImageCl *im    = new DMRImageCl();

  int x0,y0;
 
  x0 = ((cl->x+hd->x)-fr->x)/pixX + nx/2;
  y0 = ((cl->y+hd->y)-fr->y)/pixY + ny/2;
    

  //cout << "xo " << x0 << " " << y0 << endl;

  //im = aRun->GetGRIM(ihd,igr,ipol,dr);
  im = aRun->GetCLIM(ihd,icl,ipol,dr);
  if(im){
    h = im->GetHist2();
    hhh = im->GetHist2();
    nbinsX = h->GetNbinsX();
    nbinsY = h->GetNbinsY();
    //total_pixels = (nbinsX-1)*(nbinsY-1); //-1 per la riga/colonna vuota
    
    for(int i=0;i<nbinsX;i++){
      for(int j=0;j<nbinsY;j++){
	cont = h->GetBinContent(i+1,j+1);
	if(cont>0){
	  occurrence[cont] = occurrence[cont] + 1;
	  total_pixels++;
	  if(i>10 && i<30 && j>10 && j<30)hspect->Fill(cont);
	}
      }
    }
    
    
  if(hspect->GetEntries()!=0){
    double max_pos = hspect->GetBinCenter(hspect->GetMaximumBin());
    hspect->Fit(gaus,"Q","",max_pos-5,max_pos+5);
    //TF1 *gaus = hspect->GetFunction("gaus");
    bkg_mean = gaus->GetParameter(1);
    bkg_sigma =  gaus->GetParameter(2);
    sig_peak = h->GetMaximum();
  }
  else {
    bkg_mean=-1;
    sig_peak=-1;
  }
  
  }
  else{
    bkg_mean=-1;
    sig_peak=-1;
  }

  for(int i = 0; i <=255; i++) {
    occurrence[i]=0;
  }

  if(sig_peak!=-1){
    h->GetMaximumBin(locmax,locmay,locmaz);
    cout << locmax << " " << locmay << endl; 
  for(int i=0;i<nbinsX;i++){
      for(int j=0;j<nbinsY;j++){
	cont = h->GetBinContent(i+1,j+1);
	if(cont>0){
	  //hspect->Fill(cont);
	  occurrence[cont] = occurrence[cont] + 1;
	  //cout<< i << " " << j << " " << cont << " " << occurrence[cont] << endl;
	}
      }
    }

  float sum_histogram=0;

  for(int i = 0; i <=255; i++) {
    //  TAKES NUMBER OF OCCURRENCES OF A PARTICULAR PIXEL 
    //  AND DIVIDES BY THE TOTAL NUMBER OF PIXELS YIELDING 
    //  A RATIO 
    histogram[i] = (float) occurrence[i] / (float) total_pixels;
    sum_histogram += histogram[i];
    //cout << i << " " << histogram[i] << " " << sum_histogram << endl;
  }

  threshold_value = otsu_method(histogram, total_pixels);
  cout << "Threshold " << threshold_value << endl;
  //threshold_value=25;
  hh = (TH2F*)h->Clone();
  int thr_count=0;
  bool thr_is_ok=false;
  while(!thr_is_ok){
    thr_count=0;
    for(int i=dr/2;i<3*dr/2;i++){
      for(int j=dr/2;j<3*dr/2;j++){
	if(h->GetBinContent(i+1,j+1)>threshold_value)thr_count++;
      }
    }
    if(thr_count>10)threshold_value++;
    else thr_is_ok=true;
    cout << threshold_value << " " << thr_count << " " << thr_is_ok << endl;
  }
   for(int i=0;i<nbinsX;i++){
      for(int j=0;j<nbinsY;j++){
	//cout << i << " " << j << " " << threshold_value << " " << h->GetBinContent(i+1,j+1) << endl;
        if(h->GetBinContent(i+1,j+1)<threshold_value)hh->SetBinContent(i+1,j+1,0);
	if(h->GetBinContent(i+1,j+1)>=threshold_value)hh->SetBinContent(i+1,j+1,h->GetBinContent(i+1,j+1));
	//if(h->GetBinContent(i+1,j+1)<(sig_peak-2))h->SetBinContent(i+1,j+1,0);
	//else h->SetBinContent(i+1,j+1,255);
      }
    }

  
   //float average = accumulate( dist_bar.begin(), dist_bar.end(), 0.0/ dist_bar.size());
   //cout << "The average is" << average << endl;

   
   //hh->Fit(bigaus,"BQ");
}


  //hbar->Draw("");
  
  //hbar->SaveAs("hbar.root");
  //fitEllipseTGraphDLSF(gr,fit_info);
  //TFile *fit_file = TFile::Open("fitEll.root");
  //TCanvas *c = (TCanvas *)fit_file->Get("c");
  //delete fit_file;
  

  bigaus->SetParLimits(1,18,22);
  bigaus->SetParLimits(3,18,22);
  hh->Fit("bigaus","R");
  for(int i=0;i<6;i++){
    fit_info[i]=bigaus->GetParameter(i);
    bigaus->ReleaseParameter(i);
    cout << i << " " << fit_info[i] << endl;
  }

  fit_info[6] = (x0 - dr + fit_info[1] - nx/2)*pixX + fr->x -hd->x;
  fit_info[7] = (y0 - dr + fit_info[3] - ny/2)*pixY + fr->y -hd->y;
  fit_info[2] *=pixX*2;
  fit_info[4] *=-pixY*2;

  float b = (fit_info[4]-fit_info[2] + TMath::Sqrt(TMath::Power(fit_info[4]+fit_info[2],2)-4*fit_info[4]*fit_info[2]*(1-fit_info[5]*fit_info[5])))/(2*fit_info[5]*TMath::Sqrt(fit_info[4])*TMath::Sqrt(fit_info[2]));

  fit_info[8]=TMath::ACos(b);
 
  cout << "info " << nx << " " << ny << " " << pixX << " " << pixY << " " << fr->x << " " << fr->y << " " << hd->x << " " << hd->y << " " << fit_info[1] << " " << fit_info[3] << endl;

  cout << "cluster x " << cl->ID() << " " << cl->x << " " << fit_info[6] << endl;
  cout << "cluster y " << cl->ID() << " " << cl->y << " " << fit_info[7] << endl;

  cout << "cluster phi " << cl->ID() << " " << cl->phi << " " << fit_info[8] << endl;

  gROOT->SetBatch(kTRUE);
  
  //cout << fit_info[0] << " " << fit_info[1] << " " << fit_info[2] << " " << threshold_value << endl;
  c1 = new TCanvas(Form("vid_%d_gr_%d_cl_%d_pol%d",ihd,cl->igr,icl,ipol));
  c1->Divide(2,1);
  c1->cd(1);
  h->Draw("colz");
  c1->cd(2);
  hh->Draw("colz");
  //c1->cd(3);
  //hh->Draw("colz");
  //c1->cd(4);
  //c->DrawClonePad();
  c1->Modified();
  c1->Update();
  c1->SaveAs(Form("ellipse_vid_%d_gr_%d_cl_%d.root",ihd,cl->igr,icl,ipol));
  if(ipol==7){
    gSystem->Exec(Form("hadd vid_%d_gr_%d.root ellipse_vid_%d_gr_%d*",ihd,cl->igr,ihd,cl->igr));
    gSystem->Exec(Form("rm ellipse_vid_%d_gr_%d*",ihd,cl->igr));
  }
  delete bigaus;
  delete biconst;
  delete hspect;
  delete im;
  delete gaus;
  delete h;
  delete hh;
  delete hbar;
  delete c1;
}
