#include <float.h>
#include <string>
#include <sstream>
#include <iostream>
#include "TMarker.h"
#include "TText.h"
#include "TLatex.h"
#include <vector>
#include <cmath>
#include "TLegend.h"
#include <cstdio>
#include <stdlib.h>
#include <string.h>
/*
#include "/home/vale/DMDS_root6/include/DMRViewHeader.h"
#include "/home/vale/DMDS_root6/include/DMRAffine2D.h"
#include "/home/vale/DMDS_root6/include/DMRCluster.h"
#include "/home/vale/DMDS_root6/include/DMRGrain.h"
#include "/home/vale/DMDS_root6/include/DMRMicrotrack.h"
#include "/home/vale/DMDS_root6/include/DMRImage.h"
#include "/home/vale/DMDS_root6/include/DMRImage.h"
#include "/home/vale/DMDS_root6/include/DMRImage.h"
#include "/home/vale/DMDS_root6/include/DMRRun.h"
#include "/home/vale/DMDS_root6/include/DMRRunHeader.h"
*/
#include "/home/vale/Dropbox/DMDS/include/DMRViewHeader.h"
#include "/home/vale/Dropbox/DMDS/include/DMRAffine2D.h"
#include "/home/vale/Dropbox/DMDS/include/DMRCluster.h"
#include "/home/vale/Dropbox/DMDS/include/DMRGrain.h"
#include "/home/vale/Dropbox/DMDS/include/DMRMicrotrack.h"
#include "/home/vale/Dropbox/DMDS/include/DMRImage.h"
#include "/home/vale/Dropbox/DMDS/include/DMRImage.h"
#include "/home/vale/Dropbox/DMDS/include/DMRImage.h"
#include "/home/vale/Dropbox/DMDS/include/DMRRun.h"
#include "/home/vale/Dropbox/DMDS/include/DMRRunHeader.h"

using namespace std;

DMRRun  *run=0;
DMRView *v=0;

//void dgr(const char *file="/mnt/data/NEWS/C100keV/rep_nop/dm_grains.dm.root")
//void dgr(const char *file="/mnt/data/NEWS/C100keV/rep_pol/dm_grains.dm.root")
void dgr(const char *file="dm_tracks_cl.dm.root")
{
  run = new DMRRun(file);
  run->GetHeader()->Print();
  run->SetFixEncoderFaults(1);
  v = run->GetView();
  gStyle->SetOptStat("n"); 
  //grpolN2(16,296,10);
}


void drawCL( DMRCluster *cl, int igr, int icl)
{

  DMRGrain *gr = v->GetGR(cl->igr);

  Int_t cl_id = cl->ID();
  float pixX = run->GetHeader()->pixX;
  float pixY = run->GetHeader()->pixY;
  TMarker *m = new TMarker(cl->x,cl->y,8);
  if(cl->igr!=igr) m->SetMarkerColor(kGreen);
  if(cl->igr==igr) m->SetMarkerColor(kBlack);
  if(cl->ID()==icl) m->SetMarkerColor(kBlue);
  
  m->Draw();
  cout << "HELLO "<<cl->x << " " << cl->y  << " " << cl->ID() << " " << cl->igr << " " << igr << " " << icl << endl;
  printf("id:ifr:img:igr:imt %d:%d:%d:%d:%d  major/minor = %.3f/%.3f = %.3f\n",
         cl->ID(), cl->ifr, cl->img, cl->igr, cl->imt, cl->lx,cl->ly,cl->lx/cl->ly );
  float dxma = cl->lx*Cos(cl->phi)/2/pixX;
  float dyma = cl->lx*Sin(cl->phi)/2/pixY;
  TLine *linema = new TLine( cl->x-dxma,cl->y-dyma,cl->x+dxma,cl->y+dyma);
  linema->Draw();
  float dxmi = cl->ly*Cos(cl->phi+PiOver2())/2/pixX;
  float dymi = cl->ly*Sin(cl->phi+PiOver2())/2/pixY;
  TLine *linemi = new TLine( cl->x-dxmi,cl->y-dymi,cl->x+dxmi,cl->y+dymi);
  linemi->SetLineColor(kWhite);
  linemi->Draw();
  if(cl->color>0) 
  {
    TText *t = new TText();
    t->SetTextSize(0.025);
    t->DrawText(cl->x+2,cl->y+1,Form("%d %.2f",cl->igr, cl->lx/cl->ly));
  }
}


void drawMT( DMRCluster *cl, int imt, int igr)
{

  DMRGrain *gr = v->GetGR(cl->igr);
  
  float pixX = run->GetHeader()->pixX;
  float pixY = run->GetHeader()->pixY;
  TMarker *m = new TMarker(cl->x,cl->y,8);
  m->SetMarkerColor(0);
  if(cl->imt!=imt) m->SetMarkerColor(4);
  if(cl->imt==imt) m->SetMarkerColor(2);
  if(cl->igr==igr) m->SetMarkerColor(1);
  
  m->Draw();
  printf("id:ifr:img:igr:imt %d:%d:%d:%d:%d  major/minor = %.3f/%.3f = %.3f\n",
         cl->ID(), cl->ifr, cl->img, cl->igr, cl->imt, cl->lx,cl->ly,cl->lx/cl->ly );
  float dxma = cl->lx*Cos(cl->phi)/2/pixX;
  float dyma = cl->lx*Sin(cl->phi)/2/pixY;
  TLine *linema = new TLine( cl->x-dxma,cl->y-dyma,cl->x+dxma,cl->y+dyma);
  linema->Draw();
  float dxmi = cl->ly*Cos(cl->phi+PiOver2())/2/pixX;
  float dymi = cl->ly*Sin(cl->phi+PiOver2())/2/pixY;
  TLine *linemi = new TLine( cl->x-dxmi,cl->y-dymi,cl->x+dxmi,cl->y+dymi);
  linemi->SetLineColor(kWhite);
  linemi->Draw();
  if(cl->color>0) 
  {
    TText *t = new TText();
    t->SetTextSize(0.025);
    t->DrawText(cl->x+2,cl->y+1,Form("%d %.2f",cl->igr, cl->lx/cl->ly));
  }
}



TH2F* drawimcls(DMRViewHeader   *hd, DMRImageCl *im, TObjArray  *acl,int igr, int icl, bool scl)
{
  TH2F *h = 0;
  if(im)
  {
    h = im->GetHist2();
    h->Draw("colz");
    TText *t = new TText();
    t->SetTextSize(0.1);
    t->SetTextColor(1);
    t->DrawText(1,1,Form("%.2f",im->pol));	
  
    if(acl)
    {
      int ncl= acl->GetEntries();
      printf("ncl=%d\n",ncl);

	for(int ic=0; ic<ncl; ic++)
	  {
	    DMRCluster *c = (DMRCluster*)(acl->At(ic));
	    if(c->flags==0 && hd->flag==0)
	      {
		if(!scl){
		  c->x = c->x - im->x;
		  c->y = c->y - im->y;
		}
		drawCL(c,igr,icl);
	      }
	  }
    }
  }
  return h;
}



TH2F* drawimMT(DMRImageCl *im, TObjArray  *acl, int imt, int igr)
{
  TH2F *h =0;
  if(im)
  {
    h = im->GetHist2();
    h->Draw("colz");
    TText *t = new TText();
    t->SetTextSize(0.1);
    t->SetTextColor(1);
    t->DrawText(1,1,Form("%.2f",im->pol));	
  
    if(acl)
    {
      int ncl= acl->GetEntries();
      printf("ncl=%d\n",ncl);
      
      for(int ic=0; ic<ncl; ic++)
      {
        DMRCluster *c = (DMRCluster*)(acl->At(ic));
        if(c->flags==1)
        {
          c->x = c->x - im->x;
          c->y = c->y - im->y;
	  drawMT(c,imt,igr);
        }
      }      
    }
  }
  return h;
}





void bfcl(int ihd, int iv, int igr, int dr, int cl0,int cl1,int cl2,int cl3,int cl4,int cl5,int cl6,int cl7)
{

  v = run->GetEntry(ihd,1,1,1,1,1,1);
  
  TString dir = gSystem->UnixPathName(__FILE__);
  dir.ReplaceAll("dgr.C","");
  dir.ReplaceAll("/./","/");
  //cout << dir << endl;
  ifstream in;
  in.open(Form("%sweigths.txt",dir.Data()));

  
  double br[6] = {20,35,50,65,80,95};
  /*double r0[6] = {1.076,1.157,1.185,1.184,1.215,1.229};
  double r1[6] = {1.043,1.103,1.084,1.086,1.106,1.083};
  double r2[6] = {0.996,1.013,0.962,0.947,0.973,0.973};
  double r3[6] = {0.953,0.889,0.848,0.866,0.866,0.866};
  double r4[6] = {0.947,0.896,0.854,0.862,0.829,0.829};
  double r5[6] = {0.977,0.980,0.947,0.963,0.961,0.961};
  double r6[6] = {1.024,1.077,1.075,1.082,1.109,1.082};
  double r7[6] = {1.064,1.140,1.150,1.174,1.206,1.237};
*/
  //// NO BRIGHTNESS CORRECTIONS
  double r0[6] = {1,1,1,1,1,1};
  double r1[6] = {1,1,1,1,1,1};
  double r2[6] = {1,1,1,1,1,1};
  double r3[6] = {1,1,1,1,1,1};
  double r4[6] = {1,1,1,1,1,1};
  double r5[6] = {1,1,1,1,1,1};
  double r6[6] = {1,1,1,1,1,1};
  double r7[6] = {1,1,1,1,1,1};

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
  

  int filt_kern[81] = {0,-5,-4,-4,-4,-4,-4,-5, 0,-5,-4,-2, 1, 1, 1,-2,-4,-5,-4,-2, 1, 3, 5, 3, 1,-2,-4,-4, 1, 3, 8,11, 8, 3, 1,-4,-4, 1, 5,11,16,11, 5, 1,-4,-4, 1, 3, 8,11, 8, 3, 1,-4,-4,-2, 1, 3, 5, 3, 1,-2,-4,-5,-4,-2, 1, 1, 1,-2,-4,-5, 0,-5,-4,-4,-4,-4,-4,-5, 0};
  int filt_norm=128;
  int filt_thr=100;
  int filt_min=1000000;
  int filt_max=-1000000;
  //TCanvas *chist = new TCanvas("chist","chist",1200,900);
  //chist->Divide(4,3);

  float pixX = run->GetHeader()->pixX;
  float pixY = run->GetHeader()->pixY;
  int nx     = run->GetHeader()->npixX;
  int ny     = run->GetHeader()->npixY;
  int x0,y0;

  int tmp_max=0;
  int tmp_min=255;
  int tmp_maxS=0;
  int tmp_minS=1000;

  int tmp_max2=-255;
  int tmp_min2=255;

  Bool_t rmv_pnt[8]={};
  int nrmv=0;
  
  gStyle->SetOptStat(0);
  
  float xcp[8], ycp[8], pcp[8], Xcp[8], Ycp[8];
  float xcp2[8]={};
  float ycp2[8]={};
  float volcp[8], npxcp[8];
  float vol2cp[8]={};
  float npx2cp[8]={};
  float brcp[8]={};
  float br2cp[8]={};
  const int ncp=8;
  Int_t cl_id[8]={};

  float npxref=0; 
  float volref=0; 
  float npx2ref=0;
  float vol2ref=0;

  cl_id[0]=cl0;
  cl_id[1]=cl1;
  cl_id[2]=cl2;
  cl_id[3]=cl3;
  cl_id[4]=cl4;
  cl_id[5]=cl5;
  cl_id[6]=cl6;
  cl_id[7]=cl7;

  TH2F *hsum =0;
  TH2F *hsumS =0;
  TH2F *h[8]={};
  TH2F *hs[8]={};
  TH2F *h2[8]={};
  //TH2F *hfilt[8];
  DMRImage imin; imin.SetImage( 2*dr, 2*dr );
  DMRImage imax; imax.SetImage( 2*dr, 2*dr );
  Int_t im_is_ok[8]={};
  
  for(int i=0;i<8;i++){
    if(cl_id[i]!=-1){
      
      DMRImageCl *im    = run->GetCLIM(ihd,cl_id[i],i,dr);
       if(im){
	 im_is_ok[i]=1;
	 h[i] = im->GetHist2();
	 hs[i] = im->GetHist2();
	 for(int k=0;k<6;k++){
	  if(k==0){
	    if(hs[i]->GetMaximum()<br[1])hs[i]->Scale(1./brM[i][0]);
	  }
	  if(k>0 && k<5){
	    if(hs[i]->GetMaximum()>=br[k] && hs[i]->GetMaximum()<br[k+1])hs[i]->Scale(1./brM[i][k]);
	  }
	  if(k==5){
	    if(hs[i]->GetMaximum()>=br[5])hs[i]->Scale(1./brM[i][5]);
	  }
	}
	if(h[i]->GetMaximum()>tmp_max)tmp_max = h[i]->GetMaximum();
	if(h[i]->GetMinimum()<tmp_min)tmp_min = h[i]->GetMinimum();
	if(hs[i]->GetMaximum()>tmp_maxS)tmp_maxS = hs[i]->GetMaximum();
	if(hs[i]->GetMinimum()<tmp_minS)tmp_minS = hs[i]->GetMinimum();
	//cout << i << " aaaaaaaaaaaaaaaaaaaaaaaaaa " << hs[i]->GetMaximum() << " " << tmp_maxS << endl;
      }
      else {
      rmv_pnt[i]=true;
      nrmv++;
      }
      //cout << "hello "<< i << " " << rmv_pnt[i] << endl;
    }
    else {
      rmv_pnt[i]=true;
      nrmv++;
    }
    //if(im_is_ok[i]==1) cout << i << " aaaaaaaaaaaaaaaaaaaaaaaaaa " << hs[i]->GetMaximum() << " " << tmp_maxS << endl;
  }
  
  int mincont[8]={};
  int mincontS[8]={};
  
  for(int i=0;i<8;i++){
    mincont[i]=255;
    mincontS[i]=1000;
    if(!rmv_pnt[i] && im_is_ok[i]==1){
      for(int jn=0;jn<2*dr;jn++){
	for(int kn=0;kn<2*dr;kn++){
	  if(h[i]->GetBinContent(jn+1,kn+1)<mincont[i] && h[i]->GetBinContent(jn+1,kn+1)!=0 )mincont[i]=h[i]->GetBinContent(jn+1,kn+1);
	  if(hs[i]->GetBinContent(jn+1,kn+1)<mincontS[i] && hs[i]->GetBinContent(jn+1,kn+1)!=0 )mincontS[i]=hs[i]->GetBinContent(jn+1,kn+1);
	  // cout << i << " aaaaaaaaaaaaaaaaaaaaaaaaaa " << hs[i]->GetMaximum() << " " << tmp_maxS << endl;
	}
      }
    }
    //cout << mincont[i] << endl;
    //cout << i << " bbbbbbbbbbbbbbbbbbbbbbbbbb " << h[i]->GetMaximum() << " " << tmp_maxS << endl;
    //cout << i << " aaaaaaaaaaaaaaaaaaaaaaaaaa " << hs[i]->GetMaximum() << " " << tmp_maxS << endl;
  }
  
  float min_el = TMath::MinElement(8,mincont);
  float min_elS = TMath::MinElement(8,mincontS);
  bool scl = false;

  TCanvas *cnv = new TCanvas("cnv","cnv",1200,900);
  cnv->Divide(4,3);

  TCanvas *cnvS = new TCanvas("cnvS","cnvS",900,900);
  cnvS->Divide(3,3);

  TCanvas *c1S = new TCanvas("c1S","Dynamic Scaled Images",200,10,600,600);
  char gif[30] = "NULL";

  
  for(int i=0;i<8;i++){

    if(im_is_ok[i]==1)hs[i]->Reset(0);

    if(cl_id[i]!=-1){
      
    DMRViewHeader   *hd = v->GetHD();

    if(hd->flag==0){
    
    DMRCluster      *cl = v->GetCL(cl_id[i]);
    DMRFrame        *frcl = v->GetFR(cl->ifr);       
    DMRFrame        *fr  = v->GetFR(frcl->iz,frcl->ipol);

    cout << "cluster " << hd->aid << " " << hd->flag << " " << cl->igr << " " << cl->ID() << " " << cl->ipol << endl;
    
    if(i==0){
      x0 = ((cl->x+hd->x)-fr->x)/pixX + nx/2;
      y0 = ((cl->y+hd->y)-fr->y)/pixY + ny/2;
    }
    
    DMRImageCl *im    = run->GetCLIMBFC(ihd,cl_id[i],i,dr,x0,y0);
    if(im){
      TObjArray  *acl = run->GetCLCLs( ihd,cl_id[i],i,dr );  //clusters
      //cout << "ipoint " << iv << " " << cl_id[i] << " " << cl->ID() << " " << cl->ipol << " " << fr->iz << " " << cl->x << " " << cl->y << " " << cl->z << endl;
      cnv->cd(i+1)->SetGrid();
      scl=false;
      h[i] = drawimcls(hd,im,acl,igr,cl_id[i],scl);     
      h[i]->SetBinContent(2,2*dr-1,tmp_max);
      h[i]->SetBinContent(2,2,min_el);
      h[i]->SetTitle(Form("ipol %d cl %d in frame %d at xyz: %.2f %.2f %.2f",cl->ipol,cl->ID(),fr->iz,cl->x,cl->y, cl->z));
      h[i]->SetName(Form("cl%df%d",cl->ID(),cl->ifr));
      scl=true;
      for(int j=0;j<2;j++){
	if(j==0 && im_is_ok[i]==1){
	  cnvS->cd(i+1)->SetGrid();
	  hs[i] = drawimcls(hd,im,acl,igr,cl_id[i],scl);
	  for(int k=0;k<6;k++){
	    if(k==0){
	      if(hs[i]->GetMaximum()<br[1])hs[i]->Scale(1./brM[i][0]);
	    }
	    if(k>0 && k<5){
	      if(hs[i]->GetMaximum()>=br[k] && hs[i]->GetMaximum()<br[k+1])hs[i]->Scale(1./brM[i][k]);
	    }
	    if(k==5){
	      if(hs[i]->GetMaximum()>=br[5])hs[i]->Scale(1./brM[i][5]);
	    }
	  }
	  hs[i]->SetBinContent(2,2*dr-1,tmp_maxS);
	  hs[i]->SetBinContent(2,2,min_elS);
	  hs[i]->SetTitle(Form("ipol %d cl %d in frame %d at xyz: %.2f %.2f %.2f",cl->ipol,cl->ID(),fr->iz,cl->x,-cl->y, cl->z));
	  hs[i]->SetName(Form("Scl cl%df%d",cl->ID(),cl->ifr));
	    }
	
	if(j==1 && im_is_ok[i]==1){
	  c1S->cd();
	  hs[i] = drawimcls(hd,im,acl,igr,cl_id[i],scl);
	  for(int k=0;k<6;k++){
	    if(k==0){
	      if(hs[i]->GetMaximum()<br[1])hs[i]->Scale(1./brM[i][0]);
	    }
	    if(k>0 && k<5){
	      if(hs[i]->GetMaximum()>=br[k] && hs[i]->GetMaximum()<br[k+1])hs[i]->Scale(1./brM[i][k]);
	    }
	    if(k==5){
	      if(hs[i]->GetMaximum()>=br[5])hs[i]->Scale(1./brM[i][5]);
	    }
	  }
	  hs[i]->SetBinContent(2,2*dr-1,tmp_maxS);
	  hs[i]->SetBinContent(2,2,min_elS);
	  hs[i]->SetTitle(Form("ipol %d cl %d in frame %d at xyz: %.2f %.2f %.2f",cl->ipol,cl->ID(),fr->iz,cl->x,-cl->y, cl->z));
	  hs[i]->SetName(Form("Scl cl%df%d",cl->ID(),cl->ifr));
	  c1S->Modified();
	  c1S->Update();
	  sprintf(gif,"scl_pict_%.3d_%.3d_%d.gif",iv,igr,i); // add
	  c1S->SaveAs(gif);                         // add
	  //c1S->Print("file.gif+20");
	  if (gSystem->ProcessEvents()) break;
	}
	       
      }
      TText *t = new TText();
      t->SetTextSize(0.1);
      t->SetTextColor(1);
      t->DrawText(1,1,Form("%.2f",im->pol));    

      if(!hsum){
	hsum = (TH2F*)(h[i]->Clone("hsum"));
	hsum->SetTitle(Form("Sum over all polarizations"));
	hsum->SetName(Form("hsum"));
	hsumS = (TH2F*)(hs[i]->Clone("hsum"));
	hsumS->SetTitle(Form("Sum over all polarizations"));
	hsumS->SetName(Form("Scl hsum"));
	imin.Max(*im);
	imax.Max(*im);
      }
      else {
	hsum->Add(h[i]);
	hsumS->Add(hs[i]);
	imin.Min(*im);
	imax.Max(*im);
      }
      
      int ncl= acl->GetEntries();
      for(int ic=0; ic<ncl; ic++){
	DMRCluster *c = (DMRCluster*)(acl->At(ic));
	if( c && ncp<100 && cl->ifr==c->ifr && c->igr==igr){
	  //cout << "clusters " << c->x << " " << c->y << " " << c->x << " " << c->npx << " " << c->vol << " " << c->ifr << " " << cl->ifr << " " << c->pol << " " << c->igr <<  endl;	  
	  xcp2[i]+=((c->x-nx/2)*pixX + fr->x - hd->x)*c->vol;
	  ycp2[i]+=((-c->y-ny/2)*pixY + fr->y - hd->y)*c->vol;
	  npx2cp[i]+=c->npx;
	  vol2cp[i]+=c->vol;
	  //cout << "ttttttt "<<i << " " <<  xcp2[i] << " " << ycp2[i] << " " << npx2cp[i] << " " << vol2cp[i] << " " << c->x << " " << c->y << " " << c->npx << " " << c->vol << endl;
	}
      }
      
      if(vol2cp[i]!=0){
	xcp2[i]=xcp2[i]/vol2cp[i];
	ycp2[i]=ycp2[i]/vol2cp[i];
      }
      xcp[i] = cl->x;//*pixX*1000;
      ycp[i] = -cl->y;//*pixY*1000;
      pcp[i] = im->pol;
      volcp[i]=cl->vol;
      npxcp[i]=cl->npx;
 

      
      //hfilt[i] = im->GetHist2();
      //hfilt[i]->Reset();
      int nbinx = h[i]->GetNbinsX()-4;
      int nbiny = h[i]->GetNbinsY()-4;
      int min_h = 255;
      for(int ibin=0;ibin<nbinx;ibin++){
	for(int jbin=0;jbin<nbiny;jbin++){
	  if(ibin>4 && jbin>4){
	    int ii=ibin-4;
	    int ij=jbin+4;
	    int filt_value=0;
	    for(int ifilt=0;ifilt<81;ifilt++){
	      if(ifilt%9==0 && ifilt!=0){
		ij--;
		ii=ibin-4;
	      }
	      filt_value += h[i]->GetBinContent(ii+1,ij+1)*filt_kern[ifilt];
	      ii++;
	    }
	    //filt_value = (filt_value - filt_thr) / filt_norm;
	    //hfilt[i]->SetBinContent(ibin+1,jbin+1,filt_value);
	    //cout << "hello "<< " " << ibin << " " << jbin <<" " << h->GetBinContent(ibin+1,jbin+1) << " " << filt_value << endl;
	    //cout <<"hist "<< jentry << " " << kn << " " <<index_pol << endl;
	  }
	  //if(h->GetBinContent(ibin+1,jbin+1)<min_h &&  h->GetBinContent(ibin+1,jbin+1)) min_h = h->GetBinContent(ibin+1,jbin+1); 
	}
      }
      //if(hfilt[i]->GetMinimum()<filt_min)filt_min=hfilt[i]->GetMinimum();
      //if(hfilt[i]->GetMaximum()>filt_max)filt_max=hfilt[i]->GetMaximum();
      
      }
    else{
      if(i>0){
	xcp[i] = xcp[i-1];
	ycp[i] = ycp[i-1];
	pcp[i] = pcp[i-1];
	volcp[i]= volcp[i-1];
	npxcp[i]= npxcp[i-1];
	//cout << "brrr "<< xcp[i] << " " << ycp[i]  << endl;
	xcp2[i] = xcp2[i-1];
	ycp2[i] = ycp2[i-1];
      }
    }
    }
     }

   }
  
  float xlim2[2]={FLT_MAX,-FLT_MAX};
  float ylim2[2]={FLT_MAX,-FLT_MAX};
  float xlim[2]={FLT_MAX,-FLT_MAX};
  float ylim[2]={FLT_MAX,-FLT_MAX};
  for(int ic=0; ic<ncp; ic++)
    {
      if(!rmv_pnt[ic]){
	xlim[1] = Max(xcp[ic],xlim[1]);
	ylim[1] = Max(ycp[ic],ylim[1]);
	xlim[0] = Min(xcp[ic],xlim[0]);
	ylim[0] = Min(ycp[ic],ylim[0]);

	xlim2[1] = Max(xcp2[ic],xlim2[1]);
	ylim2[1] = Max(ycp2[ic],ylim2[1]);
	xlim2[0] = Min(xcp2[ic],xlim2[0]);
	ylim2[0] = Min(ycp2[ic],ylim2[0]);
      }
      
    }
  
  for(int ic=0; ic<ncp; ic++)
    {
      //cout << "2 "<<ic << " " <<  xcp2[ic] << " " << ycp2[ic] << " " << npx2cp[ic] << " " << vol2cp[ic] << endl;
      if(!rmv_pnt[ic]){
      xcp[ic] -= xlim[0];
      ycp[ic] -= ylim[0];
      xcp[ic] *= 1000;
      ycp[ic] *= 1000;

      xcp2[ic] -= xlim2[0];
      ycp2[ic] -= ylim2[0];
      xcp2[ic] *= 1000;
      ycp2[ic] *= 1000;

      if(ic==0){
	npxref 	= npxcp[0]; 
	volref 	= volcp[0]; 
	npx2ref	 = npx2cp[0];
	vol2ref	 = vol2cp[0];
      }

      npxcp[ic] = npxcp[ic]/npxref;
      volcp[ic] = volcp[ic]/volref;
      npx2cp[ic] = npx2cp[ic]/npx2ref;
      vol2cp[ic] = vol2cp[ic]/vol2ref;

      brcp[ic] = volcp[ic]/npxcp[ic];
      br2cp[ic] = vol2cp[ic]/npx2cp[ic];
      

      //cout << "3a "<<ic << " " <<  xcp[ic] << " " << ycp[ic] << " " << npxcp[ic] << " " << volcp[ic] << endl;
      //cout << "3 "<<ic << " " <<  xcp2[ic] << " " << ycp2[ic] << " " << npx2cp[ic] << " " << vol2cp[ic] << endl;
      //xscout << "3 "<<ic << " " <<  xcp[ic] << " " << ycp[ic] << " " << npxcp[ic] << " " << volcp[ic] << endl;
      }
    }

  
  Double_t x_max = TMath::MaxElement(ncp,xcp);
  Double_t y_max = TMath::MaxElement(ncp,ycp);
  Double_t xymax=0;
  if(x_max>=y_max)xymax=x_max;
  else xymax=y_max;
  
  
  cnv->cd(9)->SetGrid();

  TMultiGraph *mg = new TMultiGraph();
  TGraph *gxp2 = new TGraph(ncp,pcp,xcp2);
  for(int k=0;k<8;k++){
    //cout << k << " " << pcp[k] << " " << xcp2[k] << endl;
    if(rmv_pnt[k])gxp2->RemovePoint(k);
  }
  gxp2->SetMarkerStyle(21);
  gxp2->SetLineStyle(1);
  gxp2->SetMarkerColor(kBlack);
  //gxp2->SetTitle("x and y displacement");
  gxp2->GetXaxis()->SetRangeUser(0,180);
  gxp2->GetYaxis()->SetRangeUser(0,1.2*xymax);
  gxp2->Draw("ACP");
  mg->Add(gxp2);

  TGraph *gxp = new TGraph(ncp,pcp,xcp);
  for(int k=0;k<8;k++){
    if(rmv_pnt[k])gxp->RemovePoint(k);
  }
  gxp->SetMarkerStyle(22);
  gxp->SetLineStyle(2);
  gxp->SetMarkerColor(kGreen);
  gxp->Draw("CP");
  mg->Add(gxp);

  TGraph *gyp2 = new TGraph(ncp,pcp,ycp2);
  for(int k=0;k<8;k++){
    if(rmv_pnt[k])gyp2->RemovePoint(k);
  }
  gyp2->SetMarkerStyle(21);
  gyp2->SetLineColor(kRed);
  gyp2->SetLineStyle(1);
  gyp2->SetMarkerColor(kRed);
  gyp2->Draw("CP");
  mg->Add(gyp2);
  
  TGraph *gyp = new TGraph(ncp,pcp,ycp);
  for(int k=0;k<8;k++){
    if(rmv_pnt[k])gyp->RemovePoint(k);
  }
  gyp->SetMarkerStyle(22);
  gyp->SetLineStyle(2);
  gyp->SetLineColor(kRed);
  gyp->SetMarkerColor(kBlue);
  gyp->Draw("CP");
  mg->Add(gyp);
  mg->Draw("ACP");
  mg->SetTitle("x and y displacement");
  mg->GetYaxis()->SetTitleOffset(1.4);
  mg->GetYaxis()->SetTitle("[nm]");
  mg->GetXaxis()->SetTitle("pol angle [deg]");

  //hframe->GetXaxis()->SetRangeUser(pad->GetUxmin(),pad->GetUxmax());
  TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
  //leg3->SetHeader("");
  leg->AddEntry(gxp2,"x - bar","lp");
  leg->AddEntry(gxp,"x - bfc","lp");
  leg->AddEntry(gyp2,"y - bar","lp");
  leg->AddEntry(gyp,"y - bfc","lp");
  leg->Draw();

  cnv->cd(11)->SetGrid();

  
  TPad *pad2 = new TPad("pad2","",0,0,1,1);
  pad2->SetGrid();
  pad2->Draw();
  pad2->cd();

  // draw a frame to define the range
  //TH1F *hr2 = pad2->DrawFrame(-0.4,0,1.2,12);



  Double_t npx_max = TMath::MaxElement(ncp,npxcp);
  Double_t npx2_max = TMath::MaxElement(ncp,npx2cp);
  Double_t npxmax=0;
  if(npx_max>=npx2_max)npxmax=npx_max;
  else npxmax=npx2_max;

  Double_t vol_max = TMath::MaxElement(ncp,volcp);
  Double_t vol2_max = TMath::MaxElement(ncp,vol2cp);
  Double_t volmax=0;
  if(vol_max>=vol2_max)volmax=vol_max;
  else volmax=vol2_max;

  TMultiGraph *mg2 = new TMultiGraph(); 
  TGraph *gnpx2 = new TGraph(ncp,pcp,npx2cp);
   for(int k=0;k<8;k++){
     //cout << k << " " << rmv_pnt[k] << endl;
    if(rmv_pnt[k])gnpx2->RemovePoint(k);
  }
  gnpx2->SetMarkerStyle(21);
  gnpx2->SetTitle("npx and vol displacement");
  gnpx2->Draw("ACP");
  gnpx2->GetXaxis()->SetRangeUser(0,180);
  gnpx2->GetYaxis()->SetRangeUser(0,1.2*volmax);
  mg2->Add(gnpx2);

  
  TGraph *gnpx = new TGraph(ncp,pcp,npxcp);
   for(int k=0;k<8;k++){
    if(rmv_pnt[k])gnpx->RemovePoint(k);
  }
  gnpx->SetMarkerStyle(22);
  //gnpx->SetLineColor(kRed);
  gnpx->SetLineStyle(2);
  gnpx->SetMarkerColor(kGreen);
  gnpx->Draw("CP");
  mg2->Add(gnpx);
  
 
  TGraph *gvol2 = new TGraph(ncp,pcp,vol2cp);
  for(int k=0;k<8;k++){
    //cout << k << " " << rmv_pnt[k] << endl;
    if(rmv_pnt[k])gvol2->RemovePoint(k);
  }
  gvol2->SetMarkerStyle(21);
  gvol2->SetLineColor(kRed);
  gvol2->SetMarkerColor(kRed);
  gvol2->Draw("CP");
  mg2->Add(gvol2);
  
  TGraph *gvol = new TGraph(ncp,pcp,volcp);
   for(int k=0;k<8;k++){
    if(rmv_pnt[k])gvol->RemovePoint(k);
  }
  gvol->SetMarkerStyle(22);
  gvol->SetLineColor(kRed);
  gvol->SetLineStyle(2);
  gvol->SetMarkerColor(kBlue);
  gvol->Draw("CP");
  mg2->Add(gvol);
  //mg2->GetYaxis()->SetTitleOffset(1.4);
  
  TGraph *gbr2 = new TGraph(ncp,pcp,br2cp);
  for(int k=0;k<8;k++){
    if(rmv_pnt[k])gbr2->RemovePoint(k);
  }
  gbr2->SetMarkerStyle(21);
  gbr2->SetLineColor(kBlue);
  gbr2->SetMarkerColor(kBlue);
  gbr2->Draw("CP");
  mg2->Add(gbr2);

  TGraph *gbr = new TGraph(ncp,pcp,brcp);
  for(int k=0;k<8;k++){
    if(rmv_pnt[k])gbr->RemovePoint(k);
  }
  gbr->SetMarkerStyle(22);
  gbr->SetLineColor(kBlue);
  gbr->SetMarkerColor(kRed);
  gbr->Draw("CP");
  mg2->Add(gbr);

  TLegend *leg2 = new TLegend(0.1,0.7,0.48,0.9);
  //leg3->SetHeader("");
  leg2->AddEntry(gnpx2,"npx - bar","lp");
  leg2->AddEntry(gnpx,"npx - bfc","lp");
  leg2->AddEntry(gvol2,"vol - bar","lp");
  leg2->AddEntry(gvol,"vol - bfc","lp");
  leg2->AddEntry(gbr2,"mean vol - bar","lp");
  leg2->AddEntry(gbr,"mean vol - bfc","lp");
  leg2->Draw();

  
  cnv->cd(10)->SetGrid();

  for(int i=1;i>=0; i--)
  {
    xlim[i]-=xlim[0];
    ylim[i]-=ylim[0];
  }
  xlim[1] = ylim[1] = Max( xlim[1],ylim[1] );
  for(int i=0;i<2; i++)
  {
    xlim[i] += 0.1*(2*i-1)*xlim[1];
    ylim[i] += 0.1*(2*i-1)*ylim[1];
    xlim[i] *= 1000;
    ylim[i] *= 1000;
  }

  float maxygr=0,maxygr2, max_xygr;
  if(TMath::MaxElement(ncp,xcp2)>=TMath::MaxElement(ncp,ycp2))maxygr2=TMath::MaxElement(ncp,xcp2);
  else maxygr2 = TMath::MaxElement(ncp,ycp2);
  if(TMath::MaxElement(ncp,xcp)>=TMath::MaxElement(ncp,ycp))maxygr=TMath::MaxElement(ncp,xcp);
  else maxygr = TMath::MaxElement(ncp,ycp);
  if(maxygr2>=maxygr && maxygr2)max_xygr=maxygr2;
  else max_xygr=maxygr;

  //cout << max_xygr << endl;
  
  
  TMultiGraph *mg3 = new TMultiGraph();
  //TGraph *gxy2 = new TGraph(ncp,xcp2,ycp2);
  TGraph *gxy2 = new TGraph();
  gxy2->SetName("xy_bar");
   for(int k=0;k<8;k++){
     //cout << "points " << k << " " << xcp2[k] << " " << ycp2[k]  << endl;
     //cout << "points " << k << " " << xcp[k] << " " << ycp[k]  << endl;
     if(!rmv_pnt[k]){
       gxy2->SetPoint(k,xcp2[k],ycp2[k]);
     }
     //if(rmv_pnt[k])gxy2->RemovePoint(k);
  }
   //if(rmv_pnt[0])gxy2->SetPoint(8,xcp2[0],ycp2[0]);   // aggiunta
   
  //TGraph *gxy = new TGraph(ncp,xcp,ycp);
  TGraph *gxy = new TGraph();
  gxy->SetName("xy_bfc");
   for(int k=0;k<8;k++){
     if(!rmv_pnt[k]){
       gxy->SetPoint(k,xcp[k],ycp[k]);
     }
     //if(rmv_pnt[k])gxy->RemovePoint(k);
  }
   if(!rmv_pnt[0]){
     gxy2->SetPoint(8,xcp2[0],ycp2[0]);   // aggiunta
     gxy->SetPoint(8,xcp[0],ycp[0]);   // aggiunta
   }
   gxy2->SetMarkerColor(kBlue+1);
   gxy->SetMarkerColor(kGreen+1);
   gxy2->SetMarkerStyle(21);
   gxy2->Draw("LP");
   gxy2->GetXaxis()->SetLimits(-10,1.2*max_xygr);
   gxy2->GetYaxis()->SetLimits(-10,1.2*max_xygr);
   mg3->Add(gxy2);
   gxy->SetMarkerStyle(22);
   gxy->Draw("LP");
   mg3->Add(gxy);
   mg3->Draw("ALP");
   mg3->SetTitle("xy path");
   mg3->GetXaxis()->SetTitle("x [nm]");
   mg3->GetYaxis()->SetTitle("y [nm]");
   mg3->GetYaxis()->SetTitleOffset(1.4);
   mg3->GetXaxis()->SetLimits(-10,1.2*max_xygr);
   mg3->GetYaxis()->SetLimits(-10,1.2*max_xygr);
   TLegend *leg3 = new TLegend(0.1,0.7,0.48,0.9);
   leg3->SetName(Form("id_path_%.3d_%.3d",iv,igr));
   //leg3->SetHeader("");
   leg3->AddEntry(gxy2,"xy - bar","p");
   leg3->AddEntry(gxy,"xy - bfc","p");
   leg3->Draw();
  
  
  // Draw labels on the y axis
  //TLatex *l;// = new TLatex();
  //char *iPol[8] = {"0","22.5","45","67.5","90","112.5","135","157.5"};
  char *iPol[8] = {"0","22.5","45","67.5","90","112.5","135","157.5"};
  TText *t = new TText();
  t->SetTextAlign(32);
  t->SetTextSize(0.035);
  t->SetTextFont(72);
  t->SetTextColor(kRed);
  
  for (Int_t i=0;i<8;i++) {
    //gxy2->GetPoint(i,xcp2[i],ycp2[i]);
    //l = new TLatex(xcp2[i],ycp2[i]+1,iPol[i]);
    t->DrawText(xcp2[i],ycp2[i]+1,iPol[i]);
  }
  
  //cout << ncp << " " << nrmv << endl;
  cnv->cd(12)->SetGrid();
  hsum->Scale(1./(ncp-nrmv));
  hsum->SetContour(64);
  hsum->Draw("colz");

  cnvS->cd(9)->SetGrid();
  hsumS->Scale(1./8.);
  hsumS->Draw("colz");

  cnv->SaveAs(Form("info_%.3d_%.3d.png",iv,igr));
  cnvS->SaveAs(Form("cnv_%.3d_%.3d.png",iv,igr));

  float x_bar=0;
  float y_bar=0;
  int index_path=0;
  for(int k=0;k<8;k++){
     if(!rmv_pnt[k]){
       x_bar += xcp2[k];
       y_bar += ycp2[k];
       cout << k << " coord " << xcp2[k] << " " << ycp2[k] << endl;
       index_path++;
     }
  }
  x_bar /=index_path;
  y_bar /=index_path;

  cout << " bar " << x_bar << " " << y_bar << endl;

  float shift_x_centred = x_bar - maxygr2/2.;
  float shift_y_centred = y_bar - maxygr2/2.;
  
  //TFile *fout = new TFile(Form("xy_path_%.3d_%.3d.root",iv,igr),"RECREATE");
  //fout->cd();
  TGraph *gxy3 = new TGraph();
  gxy3->SetName(Form("xy_path_%.3d_%.3d",iv,igr));
  TCanvas *cpad = new TCanvas(Form("cnv_xy_path_%.3d_%.3d",iv,igr));   //creates new canvas
  for(int k=0;k<8;k++){
    gxy3->SetPoint(k,xcp2[k]-x_bar,ycp2[k]-y_bar);
  }
  gxy3->SetPoint(8,xcp2[0]-x_bar,ycp2[0]-y_bar);
  //TAxis *axis = gxy3->GetXaxis();
  //axis->SetLimits(-10,1.5*maxygr2);
  gxy3->GetXaxis()->SetLimits(-maxygr2,maxygr2);
  //gxy3->GetXaxis()->SetLimits(-40,40);  // per confronto
  gxy3->Draw("ALP");
  gxy3->SetMarkerColor(kBlue+1);
  gxy3->SetMarkerStyle(21);
  gxy3->SetLineWidth(2);
  gxy3->GetYaxis()->SetRangeUser(-maxygr2,maxygr2);
  //gxy3->GetYaxis()->SetRangeUser(-40,40); // per confronto
  gxy3->GetXaxis()->SetTitle("x [nm]");
  gxy3->GetYaxis()->SetTitle("y [nm]");
  for(int k=0;k<8;k++){
    t->DrawText(xcp2[k]-x_bar,ycp2[k]+1-y_bar,iPol[k]);
  }
  cpad->Draw("");
  cpad->SaveAs(Form("xy_path_%.3d_%.3d.root",iv,igr));
  delete cpad;
  

  gSystem->Exec("mv xy_path_* gif/");
  gSystem->Exec("mv cnv_* gif/");
  gSystem->Exec("mv info_* gif/");
  gSystem->Exec(Form("gifsicle --delay=30 --loop *.gif > gif/scl_anim%.3d_%.3d.gif",iv,igr));
  gSystem->Exec("rm scl_*");
  gSystem->Exec("rm file.gif");
  
  cout << "GIF saved in gif/ directory" << endl;
  
}


void mtrk(int ihd, int iv, int imt, int igr1, int igr2, int dr)
{
  float pixX = run->GetHeader()->pixX;
  float pixY = run->GetHeader()->pixY;
  int nx     = run->GetHeader()->npixX;
  int ny     = run->GetHeader()->npixY;
  int x0,y0;

  int tmp_max=0;
  int tmp_min=255;


  Bool_t rmv_pnt[8]={};
  int nrmv=0;
  
  gStyle->SetOptStat(0);
  
  float xcp[8], ycp[8], pcp[8], Xcp[8], Ycp[8];
  float xcp2[8]={};
  float ycp2[8]={};
  float volcp[8], npxcp[8];
  float vol2cp[8]={};
  float npx2cp[8]={};
  float brcp[8]={};
  float br2cp[8]={};
  const int ncp=8;
  Int_t gr_id[2]={igr1,igr2};

  float npxref=0; 
  float volref=0; 
  float npx2ref=0;
  float vol2ref=0;


  TH2F *hsum =0;
  TH2F *h[2];
  TH2F *h2[8];
  //TH2F *hfilt[8];
  DMRImage imin; imin.SetImage( 2*dr, 2*dr );
  DMRImage imax; imax.SetImage( 2*dr, 2*dr );

  TCanvas *cnv = new TCanvas("cnv","cnv",1200,600);
  cnv->Divide(2,1);

  for(int i=0;i<2;i++){
    //if(cl_id[i]!=-1){
      DMRImageCl *im    = run->GetGRIM(ihd,gr_id[i],0,dr);
      //DMRCluster *cl = v->GetCL(gr_ibfcgr_id[i]);
      if(im){
	TObjArray  *acl = run->GetGRCLs( ihd,gr_id[i],0,dr );  //clusters
	//h[i] = im->GetHist2();
	//if(h[i]->GetMaximum()>tmp_max)tmp_max = h[i]->GetMaximum();
	//if(h[i]->GetMinimum()<tmp_min)tmp_min = h[i]->GetMinimum();

	cnv->cd(i+1)->SetGrid();
	h[i] = drawimMT( im,acl,imt,gr_id[i]);
	if(h[i]->GetMaximum()>tmp_max)tmp_max = h[i]->GetMaximum();
	if(h[i]->GetMinimum()<tmp_min)tmp_min = h[i]->GetMinimum();
	//h[i]->SetTitle(Form("cl %d in frame %d at xy: %.2f %.2f",cl->ID(),cl->ifr,cl->x,-cl->y));
	//h[i]->SetName(Form("cl%df%d",cl->ID(),cl->ifr));
      }
      else {
      rmv_pnt[i]=true;
      nrmv++;
      }
      cout << "hello "<< i << endl;
      // }
  }

  char gif[30] = "NULL";
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  for(int ic=0;ic<2;ic++){
      h[ic]->Draw("colz");
      h[ic]->SetBinContent(2,2*dr-1,tmp_max);
      h[ic]->SetBinContent(2,2,tmp_min);
      c1->Modified();
      c1->Update();
      sprintf(gif,"pict_%.3d_%.3d.gif",iv,ic); // add
      c1->SaveAs(gif);                         // add
      c1->Print("file.gif+200");
      if (gSystem->ProcessEvents()) break;
  }

  gSystem->Exec(Form("gifsicle --delay=30 --loop *.gif > gif/mtrk/mtrk%d_%d.gif",iv,imt));
  gSystem->Exec("rm pict_*");
}




void fr(int ifr, int pol)
{
  float pixX = run->GetHeader()->pixX;
  float pixY = run->GetHeader()->pixY;
  int nx     = run->GetHeader()->npixX;
  int ny     = run->GetHeader()->npixY;
  DMRViewHeader   *hd = v->GetHD();
  
  DMRFrame        *fr = v->GetFR( ifr, pol );
  
  //fr->RemovePixelNoise(3);
  TH2F *hfr = fr->GetHist2();
  //hfr->Smooth();
  hfr->Draw("colz");
  /*
  int ncl = v->CL->GetEntries();
  for(int i=0; i<ncl; i++)
  {
  DMRCluster *cl = v->GetCL(i);
  if(v->GetGR(cl->igr)->ibfc==cl->ID()) cl->color=kBlack;
  if( cl && cl->ifr==ifr && v->SamePol(cl->pol,fr->pol) )
  {
  float x0 = ((cl->x+hd->x)-fr->x)/pixX + nx/2;
  float y0 = ((cl->y+hd->y)-fr->y)/pixY + ny/2;
  drawCL(x0, y0, cl);
}
}
  */
}


double scaling(int ihd, int iv, int icl, int ipol)
{

  int dr=50;
  TH2F *h;
  TH1F *hspect = new TH1F("hspect","hs",255,0,255); 
  DMRImageCl *im    = run->GetCLIM(ihd,icl,ipol,dr);
  if(im){
    h = im->GetHist2();
  }

  int nbinsX = h->GetNbinsX();
  int nbinsY = h->GetNbinsY();
  double cont=0;

  for(int i=0;i<nbinsX;i++){
    for(int j=0;j<nbinsY;j++){
      cont = h->GetBinContent(i+1,j+1);
      if(cont>0)hspect->Fill(cont);
      //cout << i+1 " " << j+1 << " " << cont << endl;
      //if()
    }
  }
  //DMRCluster *cl = v->GetCL(icl);
  //cout << cl->vol << " " << cl->npx << " " << cl->vol/cl->npx << endl;
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  h->Draw("colz");
  c1->cd(2);
  hspect->Draw("");
  hspect->Fit("gaus","Q");
  TF1 *gaus = hspect->GetFunction("gaus");
  Double_t peak_pos = gaus->GetParameter(1);
  
  return peak_pos;
}


void Nbfcl(int ifrom, int ntimes)
{

  int s;
  int ind=0;
  ifstream myfile;
  do{
    cout << "Type: " <<endl;
    cout << "1) if there are no cuts" << endl;
    cout << "2) if there is a cutlist" << endl;
    cin >> s;
    switch(s){   
    case(1):
      myfile.open("bfcl8copy.txt");
      break;
    case(2):
      char fname[100];
      cout << "Enter the filename " << endl;
      cin >> fname;
      myfile.open(fname);
      break;
    default:
      cout << "Type non correct" << endl;
    }
  }while((s!=1 && s!=2));

  
  vector<vector<string> > runDescription;
  if (myfile.is_open()){
    string line;
    // cout << line << endl;
    while (getline (myfile,line)){
      istringstream iss(line);
      string a;
      vector<string> tmpS;
      if(ind>=ifrom && ind<(ntimes+ifrom)){
	while(!iss.eof()){
	  cout << ifrom << " " << ntimes << " " << ind << endl;
	  iss>>a;
	  tmpS.push_back(a);
	  cout << a << " ";
	}
	cout << endl;
	runDescription.push_back(tmpS);
	tmpS.clear();
      }
      ind++;
      
    }
    myfile.close();
    
  }
  
  for(int i =0; i< runDescription.size();i++){
    int nstring  = runDescription[i].size();
    char row[100];
    // cout << endl;
    for(int j = 0; j<nstring; j++){
       const char *text = runDescription[i][j].c_str();
       gROOT->ProcessLine(TString::Format("%s",text));
    }
  }
  gSystem->Exec("hadd all_xy_path.root gif/xy_path*");
  gSystem->Exec("mv all_xy_path.root gif/");
  gSystem->Exec("rm AutoDict_vector_vector_string_allocator_string_____*");
  
}

void Nmtrk(int ifrom, int ntimes)
{

  int ind=0;
  ifstream myfile ("mtrk.txt");
  vector<vector<string> > runDescription;
  if (myfile.is_open()){
    string line;
    // cout << line << endl;
    while (getline (myfile,line) && ind<ntimes && ind>=ifrom){
      istringstream iss(line);
      string a;
      vector<string> tmpS;
      while(!iss.eof()){     
	iss>>a;
	tmpS.push_back(a);
	//cout << a << " ";
     }
      //cout << endl;
      runDescription.push_back(tmpS);
      tmpS.clear();
      ind++;
    }
    myfile.close();
    
  }
  
  for(int i = 0; i< runDescription.size();i++){
    int nstring  = runDescription[i].size();
    char row[100];
    // cout << endl;
    for(int j = 0; j<nstring; j++){
       const char *text = runDescription[i][j].c_str();
       gROOT->ProcessLine(TString::Format("%s",text));
    }
  }
  
}

