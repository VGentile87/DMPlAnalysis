// GENTILE VALERIO (2016-2019) LAST UPDATE 18/11/2019
//
#define myData_cxx
#include <TH2.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <vector>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <fstream>
#include <TGraph.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TLine.h>
#include <TGraph.h>
#include <TF1.h>
#include <numeric>
#include <tuple>
/*static constexpr Int_t kMaxcl = 8434; 
static constexpr Int_t kMaxgr = 207; 
static constexpr Int_t kMaxmt = 51; 
static constexpr Int_t kMaxim = 7759; 
static constexpr Int_t kMaxfr = 1000; 
*/
#include "TROOT.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TProfile2D.h"
#include <ROOT/RVec.hxx>

#include <cmath>
#include <iostream>

#include "/home/vale/DMPlAnalysis/include/DMPlsFuncImages.h"
#include "/home/vale/DMPlAnalysis/include/DMPlsOutput.h"
#include "/home/vale/DMPlAnalysis/include/DMPlsAnalyzer.h"
#include "/home/vale/DMPlAnalysis/include/DMPlsGrAnalyzer.h"
#include "/home/vale/DMPlAnalysis/include/DMPlsMtAnalyzer.h"
#include "DMPlsDefinitions.h"

using namespace std;
using namespace ROOT::VecOps;


myData vv;

void myData::Loop() {

  gSystem->Load("DMPlsAnalyzer_C.so");
  gSystem->Load("DMPlsGrAnalyzer_C.so");
  gSystem->Load("DMPlsMtAnalyzer_C.so");
  gSystem->Load("DMPlsFuncImages_C.so");
  gSystem->Load("DMPlsOutput_C.so");


DMPlsAnalyzer dmplAn;
DMPlsGrAnalyzer dmplGrAn;
DMPlsMtAnalyzer dmplMtAn;
DMPlsFuncImages dmplFuncIm;
DMPlsOutput dmplOut;


//---------------  SETTINGS  ARGUMENT ----------------------------------------------------------------//       
char datacard[100]="settings.mac";
cout << "You are using: " << datacard << "\n\n";
/// parameters from datacard
std::tie(x_pix_size, y_pix_size, len_view_x, len_view_y, thr_ldust_br,  thr_ldust_area,  fid_cut_par,  maxcl, cut_nofit, cut_goodzone, cut_ncl,  cut_minor,  cut_isolated, cut_npol,  cut_mtrk,  cut_bar_l,  cut_reg_u, cut_bar_u, cut_reg_l, cut_view, channel) = dmplAn.myDatacard(datacard, x_pix_size, y_pix_size,len_view_x,len_view_y,thr_ldust_br,thr_ldust_area,fid_cut_par, maxcl,cut_nofit,cut_goodzone,cut_ncl, cut_minor, cut_isolated, cut_npol, cut_mtrk, cut_bar_l, cut_reg_u, cut_bar_u, cut_reg_l, cut_view, channel);

cout << "PAY ATTENTION 1" << endl;
cout << "max number of polarization angles (excluding 180°) is "<< npol << endl;
cout << "If the number of polarizations is changed please fix in 'myData_vx.C'" << endl;
cout << "PAY ATTENTION 2" << endl;
cout << "The fiducial cut parameter has to be checked in the histo 'r_dist_ldust'" << endl;
cout << "If it is not correct please change 'fid_cut_par' in 'settings.mac'" << endl; 
cout << "End of settings (enjoy the results)" << "\n\n";
//-----------------------------------------------------------------------------------------------//

//---------------  OUTPUTS  ----------------------------------------------------------------//
 ofstream log_psf("list_psf.txt");
 //ofstream log_carbon("list_carbon.txt");
 std::tie(mybfcl, bfcl8, yandex, cut8, log_col, sig, bkg) = dmplOut.createLogs();
 dmplOut.createHistos(len_view_y,len_view_y,totXbin,totYbin);
 //dmplOut.createCanvas();  //for image study
 dmplOut.createGraphs();
 dmplOut.createFile();
 dmplOut.createTree();
 //-----------------------------------------------------------------------------------------------//
 
 //---- READ INPUT FILE --------------------------------------------------------------------------//
 Long64_t nentries = fChain->GetEntries();
 cout << "nentries "<< nentries << endl;
 bool scan_type = dmplAn.scanning_type(fChain);
 dimset = dmplAn.cl_flags(fChain,dimset);
 //----------------------------------------------------------------------------------------------//


 //---- TREE READING ----------------------------------------------------------------------------//
 for (Long64_t jentry=0; jentry<nentries;jentry++) { // loop on Views
   
   GetEntry(jentry);
   Long64_t ientry = LoadTree(jentry);
   
   if(scan_type)viewID=vid;   //  fragment scanning
   else viewID=aid;           //  area scanning  
   hID=id;
   
   grNcl->SetPoint(viewID,viewID,ncl);   // graph number of clusters vs viewId
   hcutarea->Reset();
   
   
   //---- ONLY CLEAN VIEWS ---------------------------------------------------------------------//
   if((ncl<maxcl || jentry==cut_view)){   // cut on dirty views

     
     // ----- CLEAN REGIONS -----//
     if(flag==0){
       TProfile2D *hctd = new TProfile2D("ctd","",8,-40,40,6,-30,30);
       for(int in=0; in<gr_;in++){
	 ctd[in]=0;
	 cut_ctd[in]=0;
	 for(int kn=0;kn<gr_;kn++){
	   if(kn!=in){
	     rdist=TMath::Sqrt(TMath::Power(gr_x[in]-gr_x[kn],2)+TMath::Power(gr_y[in]-gr_y[kn],2)+TMath::Power(gr_z[in]-gr_z[kn],2));
	     xydist=TMath::Sqrt(TMath::Power(gr_x[in]-gr_x[kn],2)+TMath::Power(gr_y[in]-gr_y[kn],2));
	     hdist=TMath::Abs(gr_z[in]-gr_z[kn]);
	     //hrdist->Fill(rdist);
	     // cout << xydist << " " << gr_x[in] << " " << gr_y[in] << endl;
	     if(xydist<1){                // user setting
	       ctd[in]++;
	       //ctd[kn]=true;
	     }
	   }
	 }
	 //cout << in << " " << ctd[in] << endl;
       }
       for(int in=0; in<gr_;in++){
	 hctd->Fill(gr_x[in],gr_y[in],ctd[in]);
       }
       for(int in=0; in<gr_;in++){
	 cut_ctd[in]= hctd->GetBinContent(hctd->FindBin(gr_x[in],gr_y[in]));
	 //cout << in << " " << ctd[in] << " " << hctd->FindBin(gr_x[in],gr_y[in]) << " " << cut_ctd[in] << endl;
       }
       if(viewID==10){
	 hctd->Draw("colz prof");
	 hctd->SaveAs("hctd.root");
       }
       delete hctd;
     }
     

     //------------ Views Initializations -----------------//  
     first_initializer(gr_,npol,gr_imt);  // function defined in DMPlsDefinitions.h

     cout << "view: "<<viewID << " " << ncl << " " << flag << endl;
     area_scan += len_view_x*len_view_y;    // area of clean views
     
     //----- ENCODER FAULTS CROSS-CHECK -------------------------------------------------------------//
     Double_t *cl_x2 = new Double_t[ncl];
     Double_t *cl_y2 = new Double_t[ncl];  
     for(int jn=0; jn<cl_;jn++){
       std::vector<double> cl_pos =
	 dmplAn.encoder_check(jn,fr_x[cl_ifr[jn]],fr_y[cl_ifr[jn]],fr_x[cl_ifr[0]],fr_y[cl_ifr[0]],cl_x[jn],cl_y[jn]);
       cl_x2[jn]=cl_pos[0];
       cl_y2[jn]=cl_pos[1];
       
       //if(cl_x2[jn]!=cl_pos[0])cout << cl_x2[jn] << " " << cl_pos[0] << " " << cl_x[jn] << endl;	
     }
     //--------------------------------------------------------------------------------------------//

     
     
     //------------------- PARTE 1: CARATTERIZZAZIONE DEI SINGOLI GRANI RICOSTRUITI --------------------//
     
     for(int in=0; in<gr_;in++){
       
       //----------- Grain Initializations ----------------------------------------------//
       second_initializer(in);         // function defined in DMPlsDefinitions.h
       gr_proc_definitions(in,gr_ncl);

       int ibfgr = gr_ibfc[in]; // index of bfc grain
       int ibffr = cl_ifr[ibfgr]; // frame of the bfc
       int set_start = ibffr - ibffr%dimset; // start ibffr of the set      
       //gr_isolated[in]=-1; // da inizializzare nella funzione
       gr_set_start[in]=set_start;
       //gr_npeaks_npol[in]=0;
       //gr_npeaks[in]=0;

	/// start cluster
	for(int jn=0; jn<cl_;jn++){ //loop sui clusters
	  
	  //------- selection of not merged grain clusters (with polarization) -------------------------//
	  if(cl_igr[jn]==gr_id[in]  && cl_flags[jn]==0){ 

	    // Bf polarized clusters belonging to the bfc frame
	    if(cl_ifr[jn]>=set_start && cl_ifr[jn]<(set_start+dimset)){
	      if(ipol_gr[in][cl_ipol[jn]]==-1) {
		ipol_gr[in][cl_ipol[jn]]=cl_id[jn];  // index of BFCs
		xb_frbf[in][cl_ipol[jn]]+=cl_x2[jn];//*cl_vol[jn]/cl_npx[jn];
		yb_frbf[in][cl_ipol[jn]]+=cl_y2[jn];//*cl_vol[jn]/cl_npx[jn];
		vol_frbf[in][cl_ipol[jn]]+=cl_vol[jn];
		npx_frbf[in][cl_ipol[jn]]+=cl_npx[jn];
		npk_x[in][cl_ipol[jn]].push_back(cl_x2[jn]);
		npk_y[in][cl_ipol[jn]].push_back(cl_y2[jn]);
		npk_npx[in][cl_ipol[jn]].push_back(cl_npx[jn]);
		npk_vol[in][cl_ipol[jn]].push_back(cl_vol[jn]);
		frbf_ent[in][cl_ipol[jn]]++;      // number of cl in the frame
		//if(jentry==3 && in<2)cout << gr_ibfc[in]<<" "<< ipol_gr[in][cl_ipol[jn]] << " " << cl_ipol[jn] << " " << set_start << " " << cl_ifr[in] << endl;
	      }
	      else {		
		if((double)cl_vol[jn]/cl_npx[jn]>((double)cl_vol[ipol_gr[in][cl_ipol[jn]]]/cl_npx[ipol_gr[in][cl_ipol[jn]]])){
		  ipol_gr[in][cl_ipol[jn]]=cl_id[jn];
		  xb_frbf[in][cl_ipol[jn]]+=cl_x2[jn];//*cl_vol[jn]/cl_npx[jn];
		  yb_frbf[in][cl_ipol[jn]]+=cl_y2[jn];//*cl_vol[jn]/cl_npx[jn];
		  vol_frbf[in][cl_ipol[jn]]+=cl_vol[jn];
		  npx_frbf[in][cl_ipol[jn]]+=cl_npx[jn];
		  npk_x[in][cl_ipol[jn]].push_back(cl_x2[jn]);
		  npk_y[in][cl_ipol[jn]].push_back(cl_y2[jn]);
		  npk_npx[in][cl_ipol[jn]].push_back(cl_npx[jn]);
		  npk_vol[in][cl_ipol[jn]].push_back(cl_vol[jn]);
		  frbf_ent[in][cl_ipol[jn]]++;      // number of cl in the frame		  
		}
	      }
	      //if(jentry==2 && in==8) cout <<"1 "<< jn << " " << ipol_gr[in][cl_ipol[jn]] << " " << cl_ipol[jn] << " " << cl_ifr[jn] << " " << frbf_ent[in][cl_ipol[jn]] << " " << ibfgr << " " << ibffr << " " << set_start <<  endl;
	    }	    
	    	    
	    //------------- Max brightness for each polarization --> BFC definition
	    if((double)cl_vol[jn]/cl_npx[jn]>tmp_vol[cl_ipol[jn]]){
	      tmp_vol[cl_ipol[jn]]=(double)cl_vol[jn]/cl_npx[jn];
	      ipol_maxbr[in][cl_ipol[jn]]=cl_id[jn];  // index of BFCs
	      
	      //---------- BFC frame for each polarization ------------------------------// 
		bfc_ifr[cl_ipol[jn]]=cl_ifr[ipol_maxbr[in][cl_ipol[jn]]];              // frame id
		bfc_zfr[cl_ipol[jn]]=fr_iz[cl_ifr[ipol_maxbr[in][cl_ipol[jn]]]];       // z coord

		//if(jentry==2 && in==8) cout <<"2 " << jn << " " << ipol_gr[in][cl_ipol[jn]] << " " << ipol_maxbr[in][cl_ipol[jn]] << cl_ipol[jn] << " " << cl_ifr[jn] << " "  <<  endl;
	      }
	    	   	  	  
	    //---------- Grain borders  (Loop on all linked clusters) 
	    if(index_clust==0){
	      ld_area[in][0]=cl_x2[jn];
	      ld_area[in][1]=cl_y2[jn];
	      ld_area[in][2]=cl_x2[jn];
	      ld_area[in][3]=cl_y2[jn];
	    }
	    else{
	      if(cl_x2[jn]<ld_area[in][0])ld_area[in][0]=cl_x2[jn];
	      if(cl_y2[jn]<ld_area[in][1])ld_area[in][1]=cl_y2[jn];
	      if(cl_x2[jn]>ld_area[in][2])ld_area[in][2]=cl_x2[jn];
	      if(cl_y2[jn]>ld_area[in][3])ld_area[in][3]=cl_y2[jn];
	    }
	    index_clust++;	  
	  }  /// end if-cluster grain
	
	  //----------- BFC of pol 9 (9 means 0 again) ---------------// (obsoleto, non c'è più cl.flags==2)
	  if(cl_igr[jn]==gr_id[in]  && cl_flags[jn]==2 && (double)cl_vol[jn]/cl_npx[jn]>tmp_gap ){
	    tmp_gap=(double)cl_vol[jn]/cl_npx[jn];
	    gr_180[in]=cl_id[jn];  	    
	  }
	} //--------------------- end cluster loop ----------------------------------------------------
	
	//-----------  Number of polarization of a grain ------------// 
	Int_t nCopy=npol;  //(nCopy andrebbe sostituito con gr_copy)
	gr_copy[in]=npol;
	//Bool_t npk_ipol[npol]={};
	for(int jn=0;jn<npol;jn++){
	  if(ipol_gr[in][jn]==-1){
	    nCopy--;
	    gr_copy[in]--; // number of copy in a collection
	  }
	  //if(jentry==3 && in<2) cout << in << " " << gr_copy[in] << endl;

	  //---------- NPEAKS IDENTIFICATION--------------------//
	  if(frbf_ent[in][jn]>1){  // more than 1 cluster in the bfc frame with the same polarization
	    //same_frame[in][cl_ipol[jn]]=2; // npeaks (there is another cluster)
	    gr_isolated[in]=2;
	    //cout <<"grpol "<< jentry << " " << in << " " << jn << " " << frbf_ent[in][cl_ipol[jn]] << endl;
	    //if(frbf_ent[in][jn]>1 && !npk_ipol[jn]){
	    gr_npeaks_npol[in]++;  //number of npeaks category over all the polarizations 
	    //npk_ipol[jn]=true;
	    //}
	  }
	  // ---------------------------------------------------//
	  xb_frbf[in][jn]=xb_frbf[in][jn]/frbf_ent[in][jn];//]vol_frbf[in][jn]/npx_frbf[in][jn];
	  yb_frbf[in][jn]=yb_frbf[in][jn]/frbf_ent[in][jn];//]vol_frbf[in][jn]/npx_frbf[in][jn];
	  num_peaks[in][jn]=frbf_ent[in][jn];  // tot number of peaks per polarization  (ridondante)
	  //num_tot_peaks[in] = num_tot_peaks[in] + num_peaks[in][jn];  // tot number of peaks over all the polarizations
	  gr_npeaks[in] = gr_npeaks[in] + num_peaks[in][jn]; // tot number of peaks over all the polarizations
	  //if(ipol_gr[in][jn]==ipol_maxbr[in][jn]) da inserire nel tree finale
	  
	}
	
	//------------------------- Grain extension ---------------------------------------------------------------------------------------------------//
	std::tie(grain_Ox[in], grain_Oy[in], grain_rx[in], grain_ry[in], grain_area[in],fiducial_cut) = dmplGrAn.borders(in, ld_area, fid_cut_par);
	//if(fiducial_cut!=0)cout <<"fid "<< in << " " << fiducial_cut << endl;

      
	//------------- Large dust and saturated grains, corners of dirty regions  (from merged BFC of a grain)
	if(TMath::Log10((double)cl_vol[gr_ibfc[in]]/cl_npx[gr_ibfc[in]])>thr_ldust_br || TMath::Log10(cl_npx[gr_ibfc[in]])>thr_ldust_area){  // cuts on settings.mac

	  std::array <float,4> corners={};
	  std::tie(ldust[in], corners, min_area_x, min_area_y, max_area_x, max_area_y) =
	    dmplGrAn.nearby_ld_area_cut(in, x, y, grain_Ox[in], grain_Oy[in], fiducial_cut, len_view_x,len_view_y);

	  //--------- Histogram of cutted area for each view (1um^2) -----------------------------------//
	  for(int i=min_area_x;i<(max_area_x+1);i++){
	    for(int j=min_area_y;j<(max_area_y+1);j++){
	      hcutarea->SetBinContent(i,j,1);
	    }
	  }
	  iLarge_dust++;
	}

	//---- delete objects ---//
	gr_delete_obj();
        
      
	//----------------- eGoodZone (quality cut with only grains far from large dust and saturated grains) -----------------------------------------//
	if(!ldust[in]){
	  for(int jn=0;jn<iLarge_dust;jn++){
	    if(gr_x[in]+x>=cut_area[jn][0] && gr_x[in]+x<=cut_area[jn][2] && gr_y[in]+y>=cut_area[jn][1] && gr_y[in]+y<=cut_area[jn][3])
	      grain[in]=false;   // near large dust
	  }
	} else grain[in]=false;   // is a large dust
	//--------------------------------- end eGoodZone --------------------------------------------------------------------------------------------//
	
	/*
	//--------- Frame multiplicity for BFC ----------------------------//   (NON FONDAMENTALE)
	Int_t nbfc_same_fr;   // Multiplicity of bfc with the same frame (max is 8)
	Int_t max_fr;	      // Frame index corresponding to the max multiplicity
	//cout << "gap_z "<< nbfc_same_fr << " " << max_fr << " " << gr_gap_z[in] << " " << in  << " " << nCopy << endl;
	std::tie(gr_gap_z[in],nbfc_same_fr,max_fr) =
	  dmplGrAn.bfc_fr_multiplicity(in,npol,nCopy,ipol_maxbr,fr_z,bfc_zfr,cl_ifr);


	//----------- Start correction of bfc fr for ncopy>=5 -------------------------------//
	
	if((nCopy-nbfc_same_fr)/(nCopy/2.)<1 && nCopy>=5){  
	  for(int jn=0; jn<cl_;jn++){                       // loop on all clusters
	    // look for linked clusters, where the  bfc is not in the fr of max multiplicity, the one in that frame
	    if(cl_igr[jn]==gr_id[in]  && cl_flags[jn]==0 && fr_iz[cl_ifr[jn]]==max_fr && bfc_zfr[cl_ipol[jn]]!=max_fr ){
	      ipol_maxbr[in][cl_ipol[jn]]=cl_id[jn];    // replace the previous bfc index with the new one
	      bfc_ifr[cl_ipol[jn]]=cl_ifr[ipol_maxbr[in][cl_ipol[jn]]];            // frame
	      bfc_zfr[cl_ipol[jn]]=fr_iz[cl_ifr[ipol_maxbr[in][cl_ipol[jn]]]];
	      bfc_gap[in]=1;                                          // there is a gap in the frames of bfc
	    }	  
	  } 
	}
	
	if((nCopy-nbfc_same_fr)==0 && nCopy!=0)bfc_gap[in]=0;    // no fr gap in the collection 
	gr_gap_z2[in] = dmplGrAn.corr_bfc_fr_check(in,npol,nCopy,ipol_maxbr,fr_z,bfc_zfr,cl_ifr);  // cross-check after replacement
	//--------- end bfc fr correction -----------------------------------------------------------------//
	*/
     }
     // -------------------------------------------- END OF GRAIN CHARACTERIZATION --------------------------------------//



     //OOOOOooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooo//
     //ooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooOOOOOOOOOOO//
     //OOOOOooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooo//

     
      
      //-------------- PARTE 2: MICROTRACKS CHARACTERIZATION -------------------------------------------------------------//
     
      for(int kn=0;kn<mt_;kn++){

	mt_initialization(mt_ngr[kn]);

	for(int in=0;in<gr_;in++){  // microtrack grain study
	  if(gr_imt[in]==kn){
	    br_mt[imt]=(double)cl_vol[gr_ibfc[in]]/cl_npx[gr_ibfc[in]];
	    npx_mt[imt]=cl_npx[gr_ibfc[in]];
	    x_mt[imt]=gr_x[in];
	    y_mt[imt]=gr_y[in];
	    z_mt[imt]=gr_z[in];
	    if(br_mt[imt]>tmp_br_mt){
	      if(tmp_mt_id!=-1)gr_mtrk[tmp_mt_id] = 11; // set to eleven a grain darker 
	      gr_mtrk[in] = 1; // univoce flag for microtracks counting (the brither grain in set to 1)
	      tmp_mt_id=in;    // 
	      tmp_br_mt=br_mt[imt];
	    }
	    else gr_mtrk[in] = 11;
	    imt++;
	  }
	}

	phi_mt[tmp_mt_id] = dmplMtAn.mtrk_phi(mt_ngr[kn],x_mt,y_mt);	                     // phi angle of the mtrk	
	mt_dif_npx[tmp_mt_id]=TMath::MaxElement(imt,npx_mt) - TMath::MinElement(imt,npx_mt); // gap between max and min npx of mtrk grains
	mt_dif_br[tmp_mt_id]=TMath::MaxElement(imt,br_mt) - TMath::MinElement(imt,br_mt);    // gap between max and min br of mtrk grains
	mt_dif_z[tmp_mt_id]=TMath::MaxElement(imt,z_mt) - TMath::MinElement(imt,z_mt);       // gap between max and min zcoord of mtrk grains

	mt_delete_obj();
      }
      //-------------------- END MTRK CHARACTERIZATION -------------------------------------------------------------------------------//




      //OOOOOooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooo//
      //ooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooOOOOOOOOOOO//
      //OOOOOooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooo//

    
      
      /////////////////// PARTE 3: COLLEZIONI DI GRANI  ////////////////////////////////////////////////////////////////////////////////////////////////
      
      for(int in = 0;in<gr_; in++){

	// phi angle (raw data) is set in the range from -pi/2 to pi/2
	if(gr_phi[in]<=(TMath::Pi()/2.))eGrainPhi=gr_phi[in];
	if(gr_phi[in]>(TMath::Pi()/2.) && gr_phi[in]<=(3*TMath::Pi()/2.))eGrainPhi=gr_phi[in] - TMath::Pi();
	if(gr_phi[in]>(3*TMath::Pi()/2.) && gr_phi[in]<=(2*TMath::Pi()))eGrainPhi=gr_phi[in] - 2*TMath::Pi();
	/////
	
	mybfcl  << "Id candidate: "<< gr_id[in] << endl;
	if(gr_imt[in]!=-1)mybfcl << "It's a microtrack " << endl;

	///////////////// ELEMENTI IN UNA COLLEZIONE
	//gr_copy[in]=npol;     
	//for(int jpol=0;jpol<npol;jpol++){
	// if(ipol_gr[in][jpol]==-1)gr_copy[in]--;
	//}
	///////////////////////////////////////////

	//////////////// INIZIALIZATIONS
	index_pol=0;  // index of polarization
	index_2pol=0; // index for mutual comparison
	match=0;
	Int_t max_br_pol=0;
	Int_t min_br_pol=255;
	Int_t ndelta= (gr_copy[in]*(gr_copy[in]-1))/2; // number of mutual comparison
	Double_t *gr_phi_pol = new Double_t[gr_copy[in]];
        Double_t *gr_min_pol = new Double_t[gr_copy[in]];
	Double_t *gr_maj_pol = new Double_t[gr_copy[in]];
	Double_t *gr_ell_pol = new Double_t[gr_copy[in]];
	//if(ndelta<2)ndelta=2;
	Int_t idelta=0;
	//Double_t *delta_phi_pol = new Double_t[gr_copy[in]];
        Double_t *delta_phi_pol = new Double_t[ndelta];
	Double_t *delta_phi_fit = new Double_t[ndelta];
	Double_t *gr_x_pol = new Double_t[gr_copy[in]];
	Double_t *gr_y_pol = new Double_t[gr_copy[in]];
	Double_t *gr_x_pol_bar = new Double_t[gr_copy[in]];
	Double_t *gr_y_pol_bar = new Double_t[gr_copy[in]];
	Double_t *gr_mx_pol_bar = new Double_t[gr_copy[in]];
	Double_t *gr_my_pol_bar = new Double_t[gr_copy[in]];
	Double_t *gr_npx_pol_bar = new Double_t[gr_copy[in]];
	Double_t *gr_vol_pol_bar = new Double_t[gr_copy[in]];
	Double_t *gr_bright_pol_bar = new Double_t[gr_copy[in]];
	Double_t *gr_z_pol = new Double_t[gr_copy[in]];
	Double_t *gr_angpol = new Double_t[gr_copy[in]];
	Double_t *gr_npx_pol = new Double_t[gr_copy[in]];
	Double_t *gr_pol_id = new Double_t[gr_copy[in]];
	Double_t *gr_amp = new Double_t[gr_copy[in]];
	Double_t *gr_step = new Double_t[npol];

	Double_t *bkg_mean = new Double_t[gr_copy[in]];
	Double_t *sig_peak = new Double_t[gr_copy[in]];
	Double_t *br_amp = new Double_t[gr_copy[in]];
	Double_t *gr_xy_dist_bar = new Double_t[gr_copy[in]];

        Double_t clx_rel_dist=0;
	Double_t cly_rel_dist=0;

	q_line[in]=0;
	m_line[in]=0;
	gr_max_amp[in]=0;
	gr_max_dist_bar[in]=0;
	gr_max_dist_fit[in]=0;
	gr_max_dist[in]=0;
	gr_path[in]=0;
	gr_fr[in]=0;
	//gr_isolated[in]=-1;
	gr_same_pos[in]=0;
	//gr_npeaks_npol[in]=0;
	//gr_npeaks[in]=0;
      
	Double_t sort_gr[npol]={-1,-1,-1,-1,-1,-1,-1,-1};

	////////////////////////////////////////////////////

	
	/////////// CARATTERIZZAZIONE ELEMENTI DI UNA COLLEZIONE
	for(int jpol=0;jpol<npol;jpol++){
	  if(ipol_gr[in][jpol]!=-1){
	    gr_phi_pol[index_pol]=cl_phi[ipol_gr[in][jpol]];	    
	    gr_x_pol_bar[index_pol]=xb_frbf[in][jpol];
	    gr_y_pol_bar[index_pol]=yb_frbf[in][jpol];
	    gr_npx_pol_bar[index_pol]=npx_frbf[in][jpol];
	    gr_vol_pol_bar[index_pol]=vol_frbf[in][jpol];
	    gr_bright_pol_bar[index_pol]=vol_frbf[in][jpol]/npx_frbf[in][jpol];
	    gr_x_pol[index_pol]=cl_x2[ipol_gr[in][jpol]];
	    gr_y_pol[index_pol]=cl_y2[ipol_gr[in][jpol]];
	    gr_z_pol[index_pol]=cl_z[ipol_gr[in][jpol]];

	    gr_maj_pol[index_pol]=cl_lx[ipol_gr[in][jpol]];
	    gr_min_pol[index_pol]=cl_ly[ipol_gr[in][jpol]];
	    gr_ell_pol[index_pol]=cl_lx[ipol_gr[in][jpol]]/cl_ly[ipol_gr[in][jpol]];
	   
	    
	    gr_npx_pol[index_pol]=cl_npx[ipol_gr[in][jpol]];
	    gr_pol_id[index_pol]=cl_ipol[ipol_gr[in][jpol]];
	    gr_angpol[index_pol]=cl_pol[ipol_gr[in][jpol]];
	    if(cl_ipol[ipol_gr[in][jpol]]==0)gr_vol_ref[in]=vol_frbf[in][jpol];	  
	    //if(same_frame[in][jpol]==2){
	      //gr_isolated[in]=2;
	      //gr_npeaks_npol[in]++;  //number of npeaks over all the polarizations 
	      //gr_npeaks[in] = gr_npeaks[in] + num_peaks[in][jpol]-1; // -1 per non avere molteplicità 2
	    //}

	    mybfcl  << "Bfc candidates: "<< ipol_gr[in][jpol] << " " << gr_pol_id[index_pol] << " " << gr_x_pol[index_pol] << " " << gr_y_pol[index_pol] << " " <<  gr_z_pol[index_pol] << endl;
	    sort_gr[cl_ipol[ipol_gr[in][jpol]]]=ipol_gr[in][jpol];	    
	    index_pol++;
	  }
	  
	}

	
	///// SET STATISTICS INFO /////////////////////
	if(index_pol>1){ // at least two clusters in a collection
	  gr_max_br_amp[in] = TMath::Mean(index_pol,br_amp);
	  gr_x_rms[in] = TMath::RMS(index_pol,gr_x_pol);
	  gr_y_rms[in] = TMath::RMS(index_pol,gr_y_pol);
	  gr_z_rms[in] = TMath::RMS(index_pol,gr_z_pol);
	  gr_npx_rms_bar[in] = TMath::RMS(index_pol,gr_npx_pol_bar);
	  gr_vol_rms_bar[in] = TMath::RMS(index_pol,gr_vol_pol_bar);
	  gr_vol_mean_bar[in] = TMath::Mean(index_pol,gr_bright_pol_bar);   // luminosità (vol_mean nome errato)   
	  gr_x_mean[in] = TMath::Mean(index_pol,gr_x_pol);
	  gr_y_mean[in] = TMath::Mean(index_pol,gr_y_pol);
	  //cout << gr_phi_rms[in] << endl;
	}


	//// CORREZIONE FIT BIGAUS PER CLUSTERS
	for(int iset = 0;iset<index_pol;iset++){
	  std::tie(gr_x_pol_bar[iset], gr_y_pol_bar[iset]) = dmplGrAn.cls_bigaus_corr(mu_x[ePolID],mu_y[ePolID],sigma_x[ePolID],sigma_y[ePolID],rho[ePolID],(gr_x_pol_bar[iset]-gr_x_mean[in])*1000,(gr_y_pol_bar[iset]-gr_y_mean[in])*1000);
	  gr_x_pol_bar[iset] /= 1000; //um
	  gr_y_pol_bar[iset] /= 1000; //um
	  gr_x_pol_bar[iset] += gr_x_mean[in];
	  gr_y_pol_bar[iset] += gr_y_mean[in];
	  xb_frbf_corr[in][iset] = gr_x_pol_bar[iset];
	  yb_frbf_corr[in][iset] = gr_y_pol_bar[iset];
	}
      
	///////// MAX DISTANCE E PHI PER BFCL E BARYCENTER IN BFCLFR
		
	for(int iset = 0;iset<index_pol;iset++){
	  /////
	  for(int jset = 0;jset<index_pol;jset++){
	    
	    //// BFCLFR
	    if(TMath::Sqrt(TMath::Power(gr_x_pol_bar[iset]-gr_x_pol_bar[jset],2)+TMath::Power(gr_y_pol_bar[iset]-gr_y_pol_bar[jset],2))>gr_max_dist_bar[in]){
	      gr_max_dist_bar[in]=TMath::Sqrt(TMath::Power(gr_x_pol_bar[iset]-gr_x_pol_bar[jset],2)+TMath::Power(gr_y_pol_bar[iset]-gr_y_pol_bar[jset],2));
	      phi_set_bar[in] = TMath::ATan((gr_y_pol_bar[iset]-gr_y_pol_bar[jset])/(gr_x_pol_bar[iset]-gr_x_pol_bar[jset]));
	      the_set_bar[in] = TMath::ATan((TMath::Power((gr_y_pol_bar[iset]-gr_y_pol_bar[jset]),2)+TMath::Power((gr_x_pol_bar[iset]-gr_x_pol_bar[jset]),2))/(gr_z_pol[iset]-gr_z_pol[jset]));
	      m_line[in] = TMath::Tan(phi_set_bar[in]);                             // servono per l'ampiezza  (coefficiente angolare della retta)
	      q_line[in] = gr_y_pol_bar[jset] - m_line[in]*gr_x_pol_bar[jset];      // servono per l'ampiezza  (intercetta)
	      //cout << in << " " << iset << " " << cl_x_pos[in][iset] << " " << cl_y_pos[in][iset] << " " << cl_z_pos[in][iset] <<  endl;
	      //cout << in << " " << jset << " " << cl_x_pos[in][jset] << " " << cl_y_pos[in][jset] << " " << cl_z_pos[in][jset] <<  endl;
	      //cout << the_set_bar[in] << " " << gr_z_pol[iset]-gr_z_pol[jset] << endl;

	      // index of polarizations with the farther distance
	      if(iset>jset){
		gr_maxpol1[in]=gr_angpol[iset];
		gr_maxpol2[in]=gr_angpol[jset];
	      }
	      else{
		gr_maxpol1[in]=gr_angpol[jset];
		gr_maxpol2[in]=gr_angpol[iset];
	      }
	    }

	    /// RELATIVE DISTANCE BETWEEN CLUSTER OF SAME SET
	    if( jset>iset){
	      clx_rel_dist = (gr_x_pol_bar[iset]-gr_x_pol_bar[jset])*1000.; //nm
	      cly_rel_dist = (gr_y_pol_bar[iset]-gr_y_pol_bar[jset])*1000.; //nm
	      if(abs(clx_rel_dist)==0 && abs(cly_rel_dist)==0) gr_same_pos[in]++; 
		//cout << iset << " " << jset << " " << gr_x_pol_bar[iset] << " " << gr_x_pol_bar[jset] << endl;
	      if(gr_copy[in]==8 && (clx_rel_dist!=0 && cly_rel_dist!=0)){
		hclx->Fill(clx_rel_dist);
		hcly->Fill(cly_rel_dist);
	      }

	      //// VARIAZIONE DELL'ANGOLO PHI DAL FIT GAUSSIANO
	      delta_phi_pol[index_2pol] = TMath::Abs(gr_phi_pol[iset]-gr_phi_pol[jset]);
	      if(delta_phi_pol[index_2pol]>(TMath::Pi()/2.)) delta_phi_pol[index_2pol] = TMath::Pi()- TMath::Abs(delta_phi_pol[index_2pol]);
	      index_2pol++;
	    }
	    
	  }
	  
	}

	
	///////////////////////////////////////////////////////
	
	////////////////////// MAX AMPLITUDE    (distanza massima di un bfcl dalla retta che collega i bfc che danno origine alla massima distanza)
	for(int iset = 0;iset<index_pol;iset++){
	  gr_amp[iset]=TMath::Abs((gr_y_pol_bar[iset]) - (m_line[in]*gr_x_pol_bar[iset] + q_line[in]))/(TMath::Sqrt(1+m_line[in]*m_line[in]));
	  if(TMath::Abs((gr_y_pol_bar[iset]) - (m_line[in]*gr_x_pol_bar[iset] + q_line[in]))/(TMath::Sqrt(1+m_line[in]*m_line[in]))>gr_max_amp[in]) gr_max_amp[in]=TMath::Abs((gr_y_pol_bar[iset]) - (m_line[in]*gr_x_pol_bar[iset] + q_line[in]))/(TMath::Sqrt(1+m_line[in]*m_line[in]));
	  //cout << in << " " << gr_maxpol1[in] << " " << gr_maxpol2[in] << " " << gr_pol_id[iset] << " " << q_line[in] << " " << m_line[in] << " " << gr_x_pol_bar[iset] << " " << gr_y_pol_bar[iset] << " " << gr_max_amp[in] << endl; 
	}
	///////////////////////////////
	
	///////// LUNGHEZZA DEL PERCORSO (LINEA CHIUSA) CALCOLATO A PARTIRE DA POL 0 (ha senso solo per numero di copie uguale a 8)
	for(int pid = 0;pid<index_pol; pid++){
	  gr_step[pid]=-1;
	  for(int iset = 0;iset<index_pol;iset++){
	    if(gr_pol_id[iset]==pid && match==0){   // solo per match=0 (ovvero il punto di origine che può essere anche diverso da pol=0°)
	      pre_step_x=gr_x_pol_bar[iset];
	      pre_step_y=gr_y_pol_bar[iset];
	      //cout << pid << " " << gr_pol_id[iset] << " " << gr_path[in] << endl;
	      first_step_x=gr_x_pol_bar[iset];
	      first_step_y=gr_y_pol_bar[iset];
	      match++;
	    }
	    if(gr_pol_id[iset]==pid && match>0 && pid>0){ /// punti successivi
	      gr_path[in] += TMath::Sqrt(TMath::Power(gr_x_pol_bar[iset]-pre_step_x,2)+TMath::Power(gr_y_pol_bar[iset]-pre_step_y,2));
	      if(gr_copy[in]==npol)gr_step[pid] = TMath::Sqrt(TMath::Power(gr_x_pol_bar[iset]-pre_step_x,2)+TMath::Power(gr_y_pol_bar[iset]-pre_step_y,2));
	      pre_step_x=gr_x_pol_bar[iset];
	      pre_step_y=gr_y_pol_bar[iset];
	      //cout << pid << " " << gr_pol_id[iset] << " " << gr_path[in] << " " << match <<  endl;
	      match++;
	    }
	  }
	  //cout << pid <<  " " << gr_path[in] << " " << match <<  endl;
	  if(match==index_pol){ // IL PATH HA SENSO SOLO PER COLLEZIONI COMPLETE
	    gr_path[in]+=TMath::Sqrt(TMath::Power(first_step_x-pre_step_x,2)+TMath::Power(first_step_y-pre_step_y,2));  // chiudo la linea
	    if(gr_copy[in]==npol)gr_step[0] = TMath::Sqrt(TMath::Power(first_step_x-pre_step_x,2)+TMath::Power(first_step_y-pre_step_y,2));
	    //cout << pid <<  " " << gr_path[in] << " " << match <<  " " << gr_step[0] << endl;
	  }
	  //cout << pid <<  " " << gr_path[in] << " " << match <<  " " << gr_step[pid] << endl;
	}
	//////////////////////////////////////////// end 

	
	gr_phi_mean[in] = TMath::Mean(ndelta,delta_phi_pol);
	gr_phi_rms[in] = TMath::RMS(ndelta,delta_phi_pol);
	gr_phi_max[in] = TMath::MaxElement(ndelta,delta_phi_pol);
	gr_mean_path[in] = TMath::Mean(npol,gr_step);
	gr_rms_path[in] = TMath::RMS(npol,gr_step);
	gr_max_path[in] = TMath::MaxElement(npol,gr_step);
	gr_rms_amp[in] = TMath::RMS(index_pol,gr_amp);
	gr_mean_amp[in] = TMath::Mean(index_pol,gr_amp);


	delete [] gr_amp;
	delete [] gr_step;
	delete [] delta_phi_pol;
	delete [] delta_phi_fit;
	delete [] gr_phi_pol;
	delete [] gr_x_pol;
	delete [] gr_y_pol;
	delete [] gr_z_pol;
	delete [] gr_min_pol;
	delete [] gr_maj_pol;
	delete [] gr_ell_pol;
	delete [] gr_npx_pol;
	delete [] gr_x_pol_bar;
	delete [] gr_y_pol_bar;
	delete [] gr_npx_pol_bar;
	delete [] gr_vol_pol_bar;
	delete [] gr_bright_pol_bar;
	delete [] gr_angpol;
	delete [] gr_pol_id;
	delete [] bkg_mean;
	delete [] sig_peak;
	delete [] br_amp;
	delete [] gr_xy_dist_bar;
	
	/////// OUTPUT LOGS ////////////////
	//mybfcl << "View: "<< viewID << " - Grain ID: " << in << " - # Elements: " << index_pol << endl;
	//mybfcl << "bfcl("<< hID << "," << viewID << "," << in << ",40,";
	if(gr_copy[in]==npol)bfcl8 << "bfcl("<< hID << "," << viewID << "," << in << ",40,";
	
	//if(gr_copy[in]==cut_npol && gr_lx[in]/gr_ly[in]!=cut_nofit && grain[in]==cut_goodzone && gr_ly[in]>=cut_minor && gr_ncl[in]<cut_ncl && gr_imt[in]==cut_mtrk && gr_isolated[in]==cut_isolated)log_col << hID << " " << viewID << " " << in << " " << gr_max_dist_bar[in] << " " << phi_set_bar[in] << endl; 
	
	//if(gr_copy[in]==cut_npol && gr_lx[in]/gr_ly[in]!=cut_nofit && grain[in]==cut_goodzone && gr_ly[in]>=cut_minor && gr_ncl[in]<cut_ncl && gr_imt[in]==cut_mtrk && gr_isolated[in]==cut_isolated && gr_max_dist_bar[in]>cut_bar_l && gr_max_dist_bar[in]<cut_bar_u &&  gr_path[in]/gr_max_dist_bar[in]>cut_reg_l && gr_path[in]/gr_max_dist_bar[in]<cut_reg_u)cut8 << "bfcl("<< hID << "," << viewID << "," << in << ",40,";
	
	
	for(int k=0;k<npol;k++){
	  if(k<7){
	    //mybfcl <<sort_gr[k] << ",";
	    if(gr_copy[in]==npol)bfcl8 <<sort_gr[k] << ",";
	    //if(gr_copy[in]==cut_npol && gr_lx[in]/gr_ly[in]!=cut_nofit && grain[in]==cut_goodzone && gr_ly[in]>=cut_minor && gr_ncl[in]<cut_ncl && gr_imt[in]==cut_mtrk && gr_isolated[in]==cut_isolated && gr_max_dist_bar[in]>cut_bar_l && gr_max_dist_bar[in]<cut_bar_u &&  gr_path[in]/gr_max_dist_bar[in]>cut_reg_l && gr_path[in]/gr_max_dist_bar[in]<cut_reg_u)cut8 <<sort_gr[k] << ",";
	  }
	  if(k==7){
	    //mybfcl <<sort_gr[k] << ")";
	    if(gr_copy[in]==npol)bfcl8 <<sort_gr[k] << ")";
	    //if(gr_copy[in]==cut_npol && gr_lx[in]/gr_ly[in]!=cut_nofit && grain[in]==cut_goodzone && gr_ly[in]>=cut_minor && gr_ncl[in]<cut_ncl && gr_imt[in]==cut_mtrk && gr_isolated[in]==cut_isolated && gr_max_dist_bar[in]>cut_bar_l && gr_max_dist_bar[in]<cut_bar_u &&  gr_path[in]/gr_max_dist_bar[in]>cut_reg_l && gr_path[in]/gr_max_dist_bar[in]<cut_reg_u)cut8 <<sort_gr[k] << ")";
	  }
	}
	//mybfcl << endl;
	if(gr_copy[in]==npol)bfcl8 << endl;
        //if(gr_copy[in]==cut_npol && gr_lx[in]/gr_ly[in]!=cut_nofit && grain[in]==cut_goodzone && gr_ly[in]>=cut_minor && gr_ncl[in]<cut_ncl && gr_imt[in]==cut_mtrk && gr_isolated[in]==cut_isolated && gr_max_dist_bar[in]>cut_bar_l && gr_max_dist_bar[in]<cut_bar_u &&  gr_path[in]/gr_max_dist_bar[in]>cut_reg_l && gr_path[in]/gr_max_dist_bar[in]<cut_reg_u)cut8 << endl;
	///////////////////////////////////
	
      }
      
      ////////////////// END 3
     
      /////////////////// PARTE 4: CARATTERIZZAZIONE DEI SINGOLI GRANI RICOSTRUITI (CONTROL HISTROGRAMS) //////
      
      if(flag==0){
      for(int in=0;in<gr_;in++){
		
	//// INIZIALIZATIONS
	tmp_gr_dist=0;
	tmp_rdist=100;
	tmp_rdist2=100;
	tmp_rdist_ldust=100;
	gr_npeaks_dist_rms[in]=-1;
	gr_npeaks_max_phi[in]=-5;
	gr_npeaks_mean_phi[in]=-5;
      
	if(grain[in] && gr_lx[in]!=gr_ly[in] && gr_ncl[in]<=30 && gr_ly[in]>=0.100  && gr_isolated[in]==-1){ // SOLO GRANI ISOLATI CHE PASSANO I CUT
	  /// MINIMA DISTANZA 
	  for(int kn=0;kn<gr_;kn++){
	    if(grain[kn] && gr_lx[kn]!=gr_ly[kn] && gr_ncl[kn]<=30 && gr_ly[kn]>=0.100 &&  gr_isolated[kn]==-1 && kn!=in){
	      // SOLO GRANI ISOLATI CHE PASSANO I CUT
	      bfcl_rdist=TMath::Sqrt(TMath::Power(gr_x[kn]-gr_x[in],2)+TMath::Power(gr_y[kn]-gr_y[in],2)+TMath::Power(gr_z[kn]-gr_z[in],2));
	      if(bfcl_rdist<tmp_rdist){
		tmp_rdist=bfcl_rdist;
	      }	    
	    }	    	    
	  }	  
	  hrdist->Fill(tmp_rdist);	  
	}
	// DEBUGGING ZONE
	
	if( gr_isolated[in]==2){ // NPEAKS
	 

	  int npk_index=0;
	  for(int jn=0;jn<npol;jn++){
	    npk_index = npk_index + (npk_x[in][jn].size()*(npk_x[in][jn].size()-1))/2.;
	  }

	  Int_t npeaks_index=0;
	  Int_t npeaks_index2=0;
	  RVec<float> gr_npeaks_dist;
	  RVec<float> gr_npeaks_phi;
	  Int_t inpeaks=0; // index npeaks

	  /// DISTANZA TRA I DUE PICCHI

	  for(int jn=0;jn<npol;jn++){
	    tmp_rdist2=100;
	    //cout << npk_x[in][jn].size() << endl;
	    for(int kn=0;kn<npk_x[in][jn].size();kn++){
	      //cout << npk_x[in][jn].size() << endl;
	      if(npk_x[in][jn].size()>1){
		for(int ln=(kn+1);ln<npk_x[in][jn].size();ln++){
		  bfcl_rdist2=TMath::Sqrt(TMath::Power(npk_x[in][jn][ln]-npk_x[in][jn][kn],2)+TMath::Power(npk_y[in][jn][ln]-npk_y[in][jn][kn],2));
		  // cout << jentry << " " << in << " " <<jn << " " <<  kn << " " << ln << " " << npk_x[in][jn].size() << " " << bfcl_rdist2 << " " << inpeaks << " " << npk_index << endl;
		  gr_npeaks_dist.push_back(bfcl_rdist2);
		  //gr_npeaks_dist[inpeaks]=bfcl_rdist2;
		  gr_npeaks_phi.push_back(TMath::ATan((npk_y[in][jn][ln]-npk_y[in][jn][kn])/(npk_x[in][jn][ln]-npk_x[in][jn][kn])));
		  if(bfcl_rdist2<tmp_rdist2){
		    tmp_rdist2=bfcl_rdist2;
		    //phiang2 = TMath::ATan((npk_y[in][jn][ln]-npk_y[in][jn][kn])/(npk_x[in][jn][ln]-npk_x[in][jn][kn]));
		    
		    twopeak_dist[in]=bfcl_rdist2;
		    twopeak_phi[in]=gr_npeaks_phi[inpeaks];		  		  
		    //theang2 = TMath::ATan((cl_y2[ipol_gr[in][cl_ipol[jn]]]-cl_y2[jn])/(cl_x2[ipol_gr[in][cl_ipol[jn]]]-cl_x2[jn]));
		    twopeak_dvol[in] = (npk_vol[in][jn][ln]-npk_vol[in][jn][kn]);   // differenza volume tra i due picchi
		    twopeak_dnpx[in] = (npk_npx[in][jn][ln]-npk_npx[in][jn][kn]);   // differenza area tra i due picchi
		    twopeak_dbri[in] = ((double)npk_vol[in][jn][ln]/npk_npx[in][jn][ln] - (double)npk_vol[in][jn][kn]/npk_npx[in][jn][kn]);  // differenza luminosità tra i due picchi
		  }
		  //cout << "debug in progress "<< npeaks_index << " " << inpeaks <<  endl;
		  
		  inpeaks++;
		  
		}
	      }
	    }	     
	  }
	  
	  npeaks_index = (inpeaks*(inpeaks-1))/2.;
	  RVec<float> delta_npeaks_phi;
	  /////// Distanza angolare su npeaks elements
	  for(int kn=0;kn<inpeaks;kn++){
	    for(int ln=0;ln<inpeaks;ln++){
	      if(ln > kn){
		delta_npeaks_phi.push_back(TMath::Abs(gr_npeaks_phi[kn] - gr_npeaks_phi[ln]));
		if(delta_npeaks_phi[npeaks_index2]>(TMath::Pi()/2.)) delta_npeaks_phi[npeaks_index2] = TMath::Pi()- TMath::Abs(delta_npeaks_phi[npeaks_index2]);
		npeaks_index2++;
	      }
	    }
	  }
	  gr_npeaks_dist_rms[in] = StdDev(gr_npeaks_dist);
	  if(delta_npeaks_phi.size()==0)delta_npeaks_phi.push_back(-5);
	  if(delta_npeaks_phi.size()==1)delta_npeaks_phi.push_back(0);
	  gr_npeaks_max_phi[in] = Max(delta_npeaks_phi);
	  gr_npeaks_mean_phi[in] = Mean(delta_npeaks_phi);
	  
	  hrdist2->Fill(bfcl_rdist2);
	  
	}
	//END OF DEBUGGING ZONE
	////////////////// end
	
	//////  MIN DIST BETWEEN LARGE DUST AND NEARBY LARGE DUST
	if(ldust[in] && gr_lx[in]!=gr_ly[in]){  /// prendo solo i large dust
	  ///  DISTANZA DAL LARGE DUST  (espressa in unità di raggio medio del grano)
	  for(int jn=0;jn<gr_;jn++){
	    if(!ldust[jn] &&  gr_lx[jn]!=gr_ly[jn]) { ///// tutti i grani tranne i large dust
	      bfcl_rdist_ldust=(TMath::Sqrt(TMath::Power(gr_x[jn]-grain_Ox[in],2)+TMath::Power(gr_y[jn]-grain_Oy[in],2)))/(TMath::Sqrt(grain_area[in]/TMath::Pi()));
	      hrdist_ldust->Fill(bfcl_rdist_ldust);
	      // cout << jentry << " " << ent << endl;
	    }
	  }
	}
	//////////////////// end
	
      } /// relative distance
      
      // } 
	
      ///////////////////////////// END 4

      /////////////////// PARTE 5: DATASET (costruzione del Tree)  ////////////////////////////////////////////////////////////////////////////////////////////////
      
      for(int in = 0;in<gr_; in++){   /// indice sui grani
	//cout << "grains " << in << endl;
	  int clset[npol]={};
	  for(int jn=0; jn<npol;jn++){     /// indice sulle polarizzazioni
	    
	    eChannel=channel;             // for smart comparison
	    eHeaderID=id;
	    eEvent=event&0x3F;
	    eViewID=viewID;
	    eFlag=flag;
	    eGrainID=in;
	    eNgr=ngr;
	    ePolID=jn;
	    ePuls=gr_ncl[in];
	    eBfcPolID=ipol_gr[in][jn];
	    eBfcID=gr_ibfc[in];
	    eBfcGap=bfc_gap[in];
	    eNcl=ncl;
	    eNclFr=fr_ncl[cl_ifr[gr_ibfc[in]]];
	    eGrainx=gr_x[in];
	    eGrainy=gr_y[in];
	    eGrainz=gr_z[in];
	    eClustx=xb_frbf[in][jn];
	    eClusty=yb_frbf[in][jn];
	    eMClustx = xb_frbf_corr[in][jn];
	    eMClusty = yb_frbf_corr[in][jn];
	    eClustMin=cl_ly[ipol_gr[in][jn]];
	    eClustMaj=cl_lx[ipol_gr[in][jn]];
	    eClustEll=cl_lx[ipol_gr[in][jn]]/cl_ly[ipol_gr[in][jn]];
	    eClustPhi=cl_phi[ipol_gr[in][jn]];
	    //eClustMeanBkg=cl_bkg_mean[in][jn]; // da verificare
	    //eClustMaxPeak=cl_sig_peak[in][jn]; // da verificare
	    //eSetBrAmp=gr_max_br_amp[in];
	    //eSetBkgAmp=gr_mean_bkg_amp[in];
	    // eSetPeakAmp=gr_max_peak_amp[in];
	    eSetBrMaxPol=gr_max_br_pol[in];
	    eSetBrMinPol=gr_min_br_pol[in];
	    eGrainMin=gr_ly[in];
	    eGrainMaj=gr_lx[in];
	    eGrainEll=gr_lx[in]/gr_ly[in];
	    if(gr_phi[in]<=(TMath::Pi()/2.))eGrainPhi=gr_phi[in];
	    if(gr_phi[in]>(TMath::Pi()/2.) && gr_phi[in]<=(3*TMath::Pi()/2.))eGrainPhi=gr_phi[in] - TMath::Pi();
	    if(gr_phi[in]>(3*TMath::Pi()/2.) && gr_phi[in]<=(2*TMath::Pi()))eGrainPhi=gr_phi[in] - 2*TMath::Pi();
	    eGrainTheta=gr_theta[in];
	    eVolume=gr_vol[in];
	    eBfcVolume=cl_vol[gr_ibfc[in]];
	    eFrBfVolume=vol_frbf[in][jn];
	    eArea=gr_npx[in];
	    eBfcArea=cl_npx[gr_ibfc[in]];
	    eFrBfArea=npx_frbf[in][jn];
	    eGrainArea=grain_area[in];
	    eXView=x;
	    eYView=y;
	    eZ=(zmin+zmax)*0.5;
	    eZlen=gr_lz[in];
	    eIsolated=gr_isolated[in];
	    eNFrame=gr_frLast[in]-gr_frFirst[in]+1;
	    eSetFrame=cl_ifr[gr_ibfc[in]];
	    eSetNCopy=gr_copy[in];
	    eSetNStatic=gr_same_pos[in];
	    eSetGapZ=gr_gap_z[in];
	    eSetGapZ2=gr_gap_z2[in];
	    
	    if(ipol_gr[in][jn]!=-1)eSetGap=TMath::Sqrt(TMath::Power(cl_y2[gr_180[in]] - cl_y2[ipol_gr[in][jn]],2) + TMath::Power(cl_x2[gr_180[in]] - cl_x2[ipol_gr[in][jn]],2));
	    
	    eSetPath=gr_path[in];
	    eSetMeanPath=gr_mean_path[in];
	    eSetRmsPath=gr_rms_path[in];
	    eSetMaxPath=gr_max_path[in];
	    eSetMaxDist=gr_max_dist[in];
	    
	    eSetNpeaksDist=twopeak_dist[in];
	    eSetNpeaksDistRms=gr_npeaks_dist_rms[in];
	    eSetNpeaksPhi=twopeak_phi[in];
	    eSetNpeaksDVol=twopeak_dvol[in];
	    eSetNpeaksDNpx=twopeak_dnpx[in];
	    eSetNpeaksDBri=twopeak_dbri[in];
	    //cout << "hey " << gr_npeaks[in] << " " << gr_copy[in] << endl;
	    if(gr_copy[in]!=0)eSetNpeaks=gr_npeaks[in]/gr_copy[in];
	    else eSetNpeaks=-1;
	    eSetNpeaksNcopy=gr_npeaks_npol[in];
	    eSetNpeaksMaxPhiAmp = gr_npeaks_max_phi[in];
	    eSetNpeaksMeanPhi = gr_npeaks_mean_phi[in];
	    eSetNpeaksMax = TMath::MaxElement(npol,num_peaks[in]);
	    eSetMaxPol1=gr_maxpol1[in];
	    eSetMaxPol2=gr_maxpol2[in];
	    eSetMaxBar=gr_max_dist_bar[in];
	    eSetMaxAmp=gr_max_amp[in];
	    eSetRmsAmp=gr_rms_amp[in];
	    eSetMeanAmp=gr_mean_amp[in];
	    eSetXRms=gr_x_rms[in];
	    eSetYRms=gr_y_rms[in];
	    eSetXBar=gr_x_mean[in];
	    eSetYBar=gr_y_mean[in];
	    eSetXMaxBar=gr_x_maxbar[in];
	    eSetYMaxBar=gr_y_maxbar[in];
	    eSetXMinBar=gr_x_minbar[in];
	    eSetYMinBar=gr_y_minbar[in];
	    eSetPhiRms=gr_phi_rms[in];
	    eSetPhiMean=gr_phi_mean[in];
	    eSetPhiMaxAmp=gr_phi_max[in];
	    eSetPhi=phi_set[in];
	    eSetPhiBar=phi_set_bar[in];
	    eSetTheBar=the_set_bar[in];
	    eSetMeanBright=gr_vol_mean_bar[in];
	    eSetNpxRms=gr_npx_rms_bar[in];
	    eSetVolRms=gr_vol_rms_bar[in];

	    
	    if(gr_vol_ref[in]!=0)eSetVolRatio=vol_frbf[in][jn]/gr_vol_ref[in];
	    
	    
	    if((eSetPhiBar>=0 && eGrainPhi>=0) || (eSetPhiBar<0 && eGrainPhi<0))eDeltaPhi = TMath::Abs(eSetPhiBar-eGrainPhi);
	    if(eSetPhiBar>=0 && eGrainPhi<0 && (eSetPhiBar - eGrainPhi)>(TMath::Pi()/2.)) eDeltaPhi = TMath::Pi() - TMath::Abs(eSetPhiBar-eGrainPhi);
	    if(eSetPhiBar<0 && eGrainPhi>=0 && (eGrainPhi - eSetPhiBar)>(TMath::Pi()/2.)) eDeltaPhi = TMath::Pi() - TMath::Abs(eSetPhiBar-eGrainPhi);
	    
	    if(TMath::Abs(eSetMaxPol1 - eSetMaxPol2)>90) eDeltaPol = 180 - TMath::Abs(eSetMaxPol1-eSetMaxPol2);
	    else eDeltaPol = TMath::Abs(eSetMaxPol1-eSetMaxPol2);
	    
	    
	    eLargeDust=ldust[in];
	    eGoodZone=grain[in];
	    eShadow=shadow[in];
	    eClustTypeDust=ctd[in];
	    eCleanPar=cut_ctd[in];

	    if(TMath::Abs(eGrainMaj*TMath::Cos(eGrainPhi))>TMath::Abs(eGrainMin*TMath::Sin(eGrainPhi)))eEllPrjX = TMath::Abs(eGrainMaj*TMath::Cos(eGrainPhi));
	    else eEllPrjX = TMath::Abs(eGrainMin*TMath::Sin(eGrainPhi));
	    if(TMath::Abs(eGrainMin*TMath::Cos(eGrainPhi))>TMath::Abs(eGrainMaj*TMath::Sin(eGrainPhi)))eEllPrjY = TMath::Abs(eGrainMin*TMath::Cos(eGrainPhi));
	    else eEllPrjY = TMath::Abs(eGrainMaj*TMath::Sin(eGrainPhi));


	    //std::tie(eMClustx, eMClusty) = dmplGrAn.cls_bigaus_corr(mu_x[ePolID],mu_y[ePolID],sigma_x[ePolID],sigma_y[ePolID],rho[ePolID],eClustx-eSetXBar,eClusty-eSetYBar);
	    //eMClustx += eSetXBar;
	    //eMClusty += eSetYBar;
	    
	    
	    //////// MICROTRACK //////////////////////
	    eMTrk=-1;
	    eMTGr=-1;
	    eMTNfr=-1;
	    eMTPhi=phi_mt[in];
	    eMTThe=0;
	    eMTLen=-1;
	    eMTChi2=-1;
	    eMTBrDif=mt_dif_br[in];
	    eMTNpxDif=mt_dif_npx[in];
	    eMTZDif=mt_dif_z[in];
	    eMTID=gr_imt[in];
	    
	    if(gr_imt[in]!=-1){
	      eMTrk=gr_mtrk[in];
	      eMTGr=mt_ngr[gr_imt[in]];
	      eMTNfr=mt_frLast[gr_imt[in]];	  	  
	      eMTThe=mt_theta[gr_imt[in]];
	      eMTLen=mt_len[gr_imt[in]];
	      eMTChi2=mt_chi2[gr_imt[in]];

	      if((eMTGr-1)/eMTLen>2.45 && eMTThe>0.4)eIsolated=3;
	      else {
		if(gr_isolated[in]==-1)eIsolated=-11;
		if(gr_isolated[in]==2)eIsolated=22;
	      }
	      
	    }
	    
	  // YANDEX LOG FILE
	  /*if(eBfcArea>10 && eGoodZone==true && eFlag==0){
	    if(jn==0) yandex << eHeaderID << ","  << eViewID << "," << eGrainID << "," << eBfcPolID << ",";
	    if(jn>0) yandex << eBfcPolID << ",";
	    if(jn==7) yandex << eIsolated << "," << eSetNCopy << endl;

	    if(jn==0 && eSetMaxBar>0.04 && eSetMaxBar<0.100 && eIsolated<0) sig << eHeaderID << "," << eViewID << "," << eGrainID << "," << eSetMaxBar << "," << eSetPhiBar << endl;
	    if(jn==0 && eSetMaxBar<0.02 && eSetMaxBar>0 && eIsolated<0) bkg << eHeaderID << "," << eViewID << "," << eGrainID << "," << eSetMaxBar << "," << eSetPhiBar << endl;
	  }*/

	  // PSF LOG FILE
	  if(eBfcArea>10 && eGoodZone==true && eFlag==0 && ePuls>1 && (eIsolated==-1 || (eIsolated==-2 && eSetNpeaksPhi==10)) && eSetNCopy==8 && eSetMaxBar<0.02 && eCleanPar==0 && eGrainMin>0.13 && eGrainMin<0.18 && eGrainEll<1.25 ){

	    if(jn==0) log_psf << eHeaderID << " "  << eViewID << " " << eGrainID << " " << eBfcPolID << " ";
	    if(jn>0) log_psf << eBfcPolID << " ";
	    if(jn==7) log_psf << eIsolated << endl;
	    
	  }
	  
	  Tree_out->Fill();
	  if(jentry%100==0)Tree_out->AutoSave("SaveSelf");
	  }
	  }
      }//end flag
      ////////////////// END 5
      
      delete [] cl_x2;
      delete [] cl_y2;
      
    } //end crash

    //// CUT AREA /////
    for(int i=0;i<hcutarea->GetNbinsX();i++){
      for(int j=0;j<hcutarea->GetNbinsY();j++){
	if(hcutarea->GetBinContent(i+1,j+1)!=0)area_cut++;
      }
    }
    //cout << "area_tagliata "<<eViewID << " " << area_cut/10000. << endl;
  } // end su tutte le View


 // comment from here for only gelatin sample
  
  ofstream log_run("log_run.txt");
  log_run << "area_scannata [um]^2: " << area_scan << endl;
  log_run << "area_tagliata [um]^2: " << area_cut/10000. << endl;
  log_run << "area_effettiva [um^2]: "<< area_scan - area_cut/10000. << endl;
  //log_run.close();

  f_out->cd();
  
  grNcl->Write();
  hcutarea->Write();
  hrdist->Write();
  hrdist2->Write();
  hrdist_ldust->Write();
  //hmtrk->Write();
  hclx->Write();
  hcly->Write();

  Tree_out->Write();
  mybfcl.close();
  bfcl8.close();
  yandex.close();
  sig.close();
  bkg.close();
  cut8.close();
  log_col.close();
  log_psf.close();
  f_out->Close();
  
} // end Loop



int myrun(){
  vv.Loop();  
}




