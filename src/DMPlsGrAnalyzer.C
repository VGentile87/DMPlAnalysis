#include "../include/DMPlsGrAnalyzer.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "TMatrixD.h"

using namespace std;

std::tuple<float,float,float,float,float,float> DMPlsGrAnalyzer::borders(int hid, double ld_area[][4], float fid_cut_par)
{

  float grOx, grOy, grRx, grRy, grNpx, fid_cut;
  
  grOx=(ld_area[hid][2]+ld_area[hid][0])/2.;
  grOy=(ld_area[hid][3]+ld_area[hid][1])/2.;
  grRx=(ld_area[hid][2]-ld_area[hid][0])/2.;
  grRy=(ld_area[hid][3]-ld_area[hid][1])/2.;
  grNpx=(ld_area[hid][2]-ld_area[hid][0])*(ld_area[hid][3]-ld_area[hid][1]);
  fid_cut=fid_cut_par*TMath::Sqrt(grNpx/TMath::Pi());//*((x_pix_size-y_pix_size)/2.);

  return std::make_tuple(grOx,grOy,grRx,grRy,grNpx,fid_cut);
}


std::tuple<bool,std::array<float,4>,float,float,float,float> DMPlsGrAnalyzer::nearby_ld_area_cut(int hid, float vx, float vy, float grOx, float grOy, float fiducial_cut, float len_view_x, float len_view_y)
{

  bool ldust;
  int  min_area_x, min_area_y, max_area_x, max_area_y;
  std::array <float,4> corners={};
  
  ldust=true;
  corners[0]=grOx+vx-fiducial_cut;
  corners[1]=grOy+vy-fiducial_cut;
  corners[2]=grOx+vx+fiducial_cut;
  corners[3]=grOy+vy+fiducial_cut;
  min_area_x = TMath::Ceil((len_view_x/2. + corners[0]-vx)*100) + 1;
  min_area_y = TMath::Ceil((len_view_y/2. + corners[1]-vy)*100) + 1;
  max_area_x = TMath::Ceil((len_view_x/2. + corners[2]-vx)*100) + 1;
  max_area_y = TMath::Ceil((len_view_y/2. + corners[3]-vy)*100) + 1;

  return std::make_tuple(ldust,corners,min_area_x,min_area_y,max_area_x,max_area_y);
    
}


std::tuple<float,Int_t,Int_t> DMPlsGrAnalyzer::bfc_fr_multiplicity(int hid, int npol, int nCopy, int ipol_gr[][8], float *fr_z, int *bfc_zfr, int *cl_ifr)
{

  ///// GAP Z (VIBRATION OR ENCODER FAULTS AND BFCFR GAP) //////////////////////
  //Double_t zbfcl[nCopy];//={};	
  Double_t *zbfcl = new Double_t[nCopy];
  Int_t iCopy=0;
  float gr_gap_z=-100;
  Int_t tmp_nfr[npol];
  Int_t tmp_clfr[npol];
  Int_t fr_sort[npol]; // fr_sort ordina in ordine decrescente per indice i molteplicità dei frame
  Int_t ifr=0;

  for(int jn=0;jn<npol;jn++){
    tmp_nfr[jn]=0;
    tmp_clfr[jn]=-1;
    fr_sort[jn]=-1;
    int polgr=ipol_gr[hid][jn];
    int clifr=cl_ifr[polgr];
    
    if(ipol_gr[hid][jn]!=-1){
      zbfcl[iCopy]=fr_z[clifr];
      iCopy++;
    }
  }
  if(nCopy!=0)gr_gap_z = TMath::MaxElement(nCopy,zbfcl) - TMath::MinElement(nCopy,zbfcl);

  ///// CORRECTION OF BFC FRAME //////////////////////////////////  (correggo i frame con diverso zeta rispetto alla popolazione dominante)
 
  
  for(int jn=0;jn<npol;jn++){
    bool match_fr=false;               // booleano
    if(ipol_gr[hid][jn]!=-1 && jn==0){  // se la polarizzazione è 0
      tmp_clfr[ifr]=bfc_zfr[jn];       // il tmp frz è il bfczfr di 0
      tmp_nfr[ifr]++;                  // è di partenza sempre 1 (numero di frame uguali)
      ifr++;                           // numero di frame tra i bfcfr (iframe diventa 1 (numero minimo di frame))
    }

    for(int kn=0;kn<ifr;kn++){         
      if(ipol_gr[hid][jn]!=-1 && jn>0 && bfc_zfr[jn]==tmp_clfr[kn]){  // solo pol diverse da 0 e con frame uguale al frame della pol jn-esima
	tmp_nfr[kn]++;    // conta numero di frame uguali al frame della polarizzazione iesima, quando è uguale match_fr diventa true
	match_fr=true;
      }
    }
    
    if(ipol_gr[hid][jn]!=-1 && jn>0 && !match_fr){   // solo se match_fr è falso, cioè se diverso dal frame 0, guardo i frame successivi
      ifr++;                                        // incrementa
      tmp_clfr[ifr]=bfc_zfr[jn];                    // nuovo frame nel tmp
      tmp_nfr[ifr]++;                               // tmp frame è 1 (molteplicità minima)
    }
  }

  TMath::Sort(npol,tmp_nfr,fr_sort);
  Int_t nbfc_same_fr = TMath::MaxElement(npol,tmp_nfr);   // Indice di massima molteplicità tra i frame
  Int_t max_fr = tmp_clfr[fr_sort[0]];                 // cerca il frame associato all'indice di massima molteplicità
  delete [] zbfcl;

  return std::make_tuple(gr_gap_z,nbfc_same_fr,max_fr);
}


float DMPlsGrAnalyzer::corr_bfc_fr_check(int hid, int npol, int nCopy, int ipol_gr[][8], float *fr_z, int *bfc_zfr, int *cl_ifr)
{
  Int_t iCopy=0;
  float gr_gap_z=-100;
  Double_t *zbfcl = new Double_t[nCopy];
	
  for(int jn=0;jn<npol;jn++){
    
    int polgr=ipol_gr[hid][jn];
    int clifr=cl_ifr[polgr];

    if(ipol_gr[hid][jn]!=-1){
      zbfcl[iCopy]=fr_z[clifr];
      iCopy++;
    }
  }

  if(nCopy!=0)gr_gap_z = TMath::MaxElement(nCopy,zbfcl) - TMath::MinElement(nCopy,zbfcl);
  else gr_gap_z=-100;
  
  delete [] zbfcl;

  return gr_gap_z;

}




std::tuple<double,double> DMPlsGrAnalyzer::cls_bigaus_corr(double mu_x, double mu_y, double sigma_x, double sigma_y, double rho, double cl_x, double cl_y){
  
  double phi = TMath::ATan(rho);
  double sigma_maj = 0;
  double sigma_min = 0;
  if(sigma_x>sigma_y){
    sigma_maj = sigma_x;///TMath::Cos(phi);
    sigma_min = sigma_y;///TMath::Cos(phi);
  }
  else{
    sigma_maj = sigma_y;///TMath::Cos(phi);
    sigma_min = sigma_x;///TMath::Cos(phi);
  }
  
  
  TMatrixD m1(2,2);
  TArrayD data1(4);
  data1[0] = TMath::Cos(phi);
  data1[1] = -TMath::Sin(phi);
  data1[2] = TMath::Sin(phi);
  data1[3] = TMath::Cos(phi);
  m1.SetMatrixArray(data1.GetArray());
  
  TMatrixD m2(2,2);
  TArrayD data2(4);
  data2[0] = 1./sigma_maj;
  data2[1] = 0;
  data2[2] = 0;
  data2[3] = 1./sigma_min;
  m2.SetMatrixArray(data2.GetArray());
  
  TMatrixD m3(2,2);
  TArrayD data3(4);
  data3[0] = TMath::Cos(phi);
  data3[1] = TMath::Sin(phi);
  data3[2] = -TMath::Sin(phi);
  data3[3] = TMath::Cos(phi);
  m3.SetMatrixArray(data3.GetArray());
  
  TMatrixD m4(2,1);
  TArrayD data4(2);
  data4[0] = cl_x - mu_x;
  data4[1] = cl_y - mu_y; 
  m4.SetMatrixArray(data4.GetArray());
  
  TMatrixD A(2,1),B(2,1),C(2,1);
  A.Mult(m3,m4);
  B.Mult(m2,A);
  C.Mult(m1,B);

  
  double x = TMath::Sqrt(sigma_maj*sigma_min)*C(0,0); 
  double y = TMath::Sqrt(sigma_maj*sigma_min)*C(1,0);  
    //double x = TMath::Sqrt(sigma_x*sigma_y)*m1*m2*m3*m4;

  //cout << x << " " << y << endl;
  return std::make_tuple(x,y);
}

