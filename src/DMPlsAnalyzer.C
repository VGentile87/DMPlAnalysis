#include "../include/DMPlsAnalyzer.h"
#include "TMatrixDEigen.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

std::vector<double> DMPlsAnalyzer::encoder_check(int ipol, double fr_cl_x, double fr_cl_y, double fr_cl_x_pol0, double fr_cl_y_pol0, double cl_x, double cl_y)
{
  vector <double> cl_pos(2);  
  if(ipol==0){
    cl_pos[0]=cl_x;
    cl_pos[1]=cl_y;
  }
  if(ipol!=0){
    if(fr_cl_x!=fr_cl_x_pol0) cl_pos[0]=cl_x - (fr_cl_x - fr_cl_x_pol0);	  	  
    else cl_pos[0]=cl_x;
    if(fr_cl_y!=fr_cl_y_pol0) cl_pos[1]=cl_y - (fr_cl_y - fr_cl_y_pol0);	  	  
    else cl_pos[1]=cl_y;
    }
  
  return cl_pos;
}


std::tuple <double, double, double, double, double, double, double, int, int, int, int, double, int, int, int, double, double, double, double, int, int> DMPlsAnalyzer::myDatacard (char*datacard, double x_pix_size, double y_pix_size, double len_view_x, double len_view_y, double thr_ldust_br, double thr_ldust_area, double fid_cut_par, int maxcl, int cut_nofit, int cut_goodzone, int cut_ncl, double cut_minor, int cut_isolated, int cut_npol, int cut_mtrk, double cut_bar_l, double cut_reg_u, double cut_bar_u, double cut_reg_l, int cut_view, int channel)
/*
void DMPlsAnalyzer::myDatacard(char* datacard, double x_pix_size, double y_pix_size, double len_view_x, double len_view_y, double thr_ldust_br, double thr_ldust_area, double fid_cut_par, int maxcl, int cut_nofit, int cut_goodzone, int cut_ncl, double cut_minor, int cut_isolated, int cut_npol, int cut_mtrk, double cut_bar_l, double cut_reg_u, double cut_bar_u, double cut_reg_l, int cut_view, int channel)
*/
{
  ifstream myfile(datacard);
  vector<vector<string> > runDescription;
  if (myfile.is_open()){
    string line;
    // cout << line << endl;
    while (getline (myfile,line)){
      istringstream iss(line);
      string a;
      vector<string> tmpS;
      while(!iss.eof()){     
	iss>>a;
	tmpS.push_back(a);
	//cout << a << " ";
      }
      runDescription.push_back(tmpS);
      tmpS.clear();
    }
    myfile.close();
    
  }
  else cout << "Unable to open file: "<< datacard << endl; 

  cout << "SETTINGS from " << datacard << endl;
  for(int i = 0; i< runDescription.size();i++){
    int nstring  = runDescription[i].size();
    // cout << endl;
    for(int j = 0; j<nstring; j++){
      //cout << runDescription[i][j]<<endl;
      if(runDescription[i][j] == "pix_size_x") {
	x_pix_size = atof( runDescription[i][2].c_str());
	cout << "pixel size x [um] is " << x_pix_size << endl;
      }
      if(runDescription[i][j] == "pix_size_y") {
	y_pix_size = atof( runDescription[i][2].c_str());
	cout << "pixel size y [um] is "<< y_pix_size << endl;
      }
      if(runDescription[i][j] == "len_view_x") {
	len_view_x = atof( runDescription[i][2].c_str());
	cout << "view size x [um] is " << len_view_x << endl;
      }
      if(runDescription[i][j] == "len_view_y") {
	len_view_y = atof( runDescription[i][2].c_str());
	cout << "view size y [um] is " << len_view_y << endl;
      }
      if(runDescription[i][j] == "thr_ldust_br") {
	thr_ldust_br = atof( runDescription[i][2].c_str());
	cout << "threshold brigthness for large dust  (log10 scale) is " << thr_ldust_br << endl;
      }
      if(runDescription[i][j] == "thr_ldust_area") {
	thr_ldust_area = atof( runDescription[i][2].c_str());
	cout << "threshold number of pixel for large dust (log 10 scale) is " << thr_ldust_area << endl;
      }
      if(runDescription[i][j] == "fid_cut_par") {
	fid_cut_par = atof( runDescription[i][2].c_str());
	cout << "fiducial cut parameter for nearby large dust (unit of large dust area) is "<< fid_cut_par << endl;
      }
      if(runDescription[i][j] == "maxcl") {
	maxcl = atof( runDescription[i][2].c_str());
	cout << "max cluster number per view to start the  analysis (no dirty view control) is "<< maxcl << endl;
      }
      if(runDescription[i][j] == "cut_nofit") {
	cut_nofit = atof( runDescription[i][2].c_str());
	cout << "No fit cut for ellipticity equal to "<< cut_nofit << endl;
      }
      if(runDescription[i][j] == "cut_goodzone") {
	cut_goodzone = atof( runDescription[i][2].c_str());
	cout << "Boolean Good zone cut is "<< cut_goodzone << endl;
      }
      if(runDescription[i][j] == "cut_minor") {
	cut_minor = atof( runDescription[i][2].c_str());
	cout << "Lower cut on minor length is "<< cut_minor << " [um]" <<  endl;
      }
      if(runDescription[i][j] == "cut_ncl") {
	cut_ncl = atof( runDescription[i][2].c_str());
	cout << "Cut on number of merged cluster liked to a grain is "<< cut_ncl << endl;
      }
      if(runDescription[i][j] == "cut_mtrk") {
	cut_mtrk = atof( runDescription[i][2].c_str());
	cout << "Mtrk index is set to "<< cut_mtrk << endl;
      }
      if(runDescription[i][j] == "cut_npol") {
	cut_npol = atof( runDescription[i][2].c_str());
	cout << "Lower cut of number of polarization is "<< cut_npol << endl;
      } 
      if(runDescription[i][j] == "cut_isolated") {
	cut_isolated = atof( runDescription[i][2].c_str());
	cout << "Isolated index is set to "<< cut_isolated << endl;
      }
      if(runDescription[i][j] == "cut_bar_l") {
	cut_bar_l = atof( runDescription[i][2].c_str());
	cout << "Lower cut on barshift is "<< cut_bar_l << " [um]" << endl;
      }
      if(runDescription[i][j] == "cut_reg_l") {
	cut_reg_l = atof( runDescription[i][2].c_str());
	cout << "Lower cut on regularity index is "<< cut_reg_l << endl;
      }
      if(runDescription[i][j] == "cut_bar_u") {
	cut_bar_u = atof( runDescription[i][2].c_str());
	cout << "Upper cut on barshift is "<< cut_bar_u << " [um]" << endl;
      }
      if(runDescription[i][j] == "cut_reg_u") {
	cut_reg_u = atof( runDescription[i][2].c_str());
	cout << "Upper cut on regularity index is "<< cut_reg_u << endl;
      }
      if(runDescription[i][j] == "cut_view") {
	cut_view = atof( runDescription[i][2].c_str());
	cout << "Views removed "<< cut_view << endl;
      }

      if(runDescription[i][j] == "channel") {
	channel = atof( runDescription[i][2].c_str());
	cout << "channel "<< channel << endl;
      } 
    }
  }
  return std::make_tuple(x_pix_size, y_pix_size, len_view_x, len_view_y, thr_ldust_br, thr_ldust_area,  fid_cut_par, maxcl, cut_nofit, cut_goodzone, cut_ncl, cut_minor, cut_isolated, cut_npol, cut_mtrk, cut_bar_l, cut_reg_u, cut_bar_u, cut_reg_l, cut_view, channel);
}


bool DMPlsAnalyzer::scanning_type(TTree *fChain)
{
 
 fChain->Draw("vid","","goff");
 Float_t vid_mean = TMath::Mean(fChain->GetSelectedRows(),fChain->GetV1());
 fChain->Draw("aid","","goff");
 Float_t aid_mean = TMath::Mean(fChain->GetSelectedRows(),fChain->GetV1());
 //cout << aid_mean << " " << vid_mean << " " << viewID << endl;

 if(vid_mean>aid_mean) return true;
 else return false;
}
