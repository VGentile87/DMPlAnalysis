// point2d.h
#ifndef DMPlsGrAnalyzer_H
#define DMPlsGrAnalyzer_H

#include <TVectorD.h>
#include <TGraph.h>
#include <TTree.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <array>

using namespace std;

/*namespace DMPlsAnalyzer {
	class Plasmon;
}*/
class DMPlsGrAnalyzer
{
 public:
  
  std::tuple<float,float,float,float,float,float> borders(int hid, double ld_area[][4], float fid_cut_par);
  std::tuple<bool,std::array<float,4>,float,float,float,float> nearby_ld_area_cut(int hid, float vx, float vy, float grOx, float grOy, float fiducial_cut, float len_view_x, float len_view_y);

  std::tuple<float,Int_t,Int_t> bfc_fr_multiplicity(int hid, int npol, int nCopy, int ipol_gr[][8], float *fr_z, int *bfc_zfr, int *cl_ifr);

  float corr_bfc_fr_check(int hid, int npol, int nCopy, int ipol_gr[][8], float *fr_z, int *bfc_zfr, int *cl_ifr);
  std::tuple<double,double> cls_bigaus_corr(double mu_x, double mu_y, double sigma_x, double sigma_y, double rho, double cl_x, double cl_y);

private:

};
#endif // DMPlsGrAnalyzer_H


