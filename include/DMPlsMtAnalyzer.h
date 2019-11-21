// point2d.h
#ifndef DMPlsMtAnalyzer_H
#define DMPlsMtAnalyzer_H

#include <TVectorD.h>
#include <TGraph.h>
#include <TTree.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <array>
#include <TF1.h>

using namespace std;

/*namespace DMPlsAnalyzer {
	class Plasmon;
}*/
class DMPlsMtAnalyzer
{
 public:

  double mtrk_phi(int mt_ngr, double *x_mt, double *y_mt);  
  
private:

};
#endif // DMPlsMtAnalyzer_H


