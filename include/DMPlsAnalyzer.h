// point2d.h
#ifndef DMPlsAnalyzer_H
#define DMPlsAnalyzer_H

#include <TVectorD.h>
#include <TGraph.h>
#include <TTree.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>

using namespace std;

/*namespace DMPlsAnalyzer {
	class Plasmon;
}*/
class DMPlsAnalyzer
{
 public:
  std::vector <double> encoder_check(int ipol, double fr_x, double fr_y, double fr_x0, double fr_y0, double cl_x, double cl_y);

   std::tuple <double, double, double, double, double, double, double, int, int, int, int, double, int, int, int, double, double, double, double, int, int> myDatacard (char*datacard, double x_pix_size, double y_pix_size, double len_view_x, double len_view_y, double thr_ldust_br, double thr_ldust_area, double fid_cut_par, int maxcl, int cut_nofit, int cut_goodzone, int cut_ncl, double cut_minor, int cut_isolated, int cut_npol, int cut_mtrk, double cut_bar_l, double cut_reg_u, double cut_bar_u, double cut_reg_l, int cut_view, int channel);

   bool scanning_type(TTree *fChain);

     private:
   double x;
   double y;
};
#endif // DMPlsAnalyzer_H
