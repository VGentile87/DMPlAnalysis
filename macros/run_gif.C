#include "TSystem.h"
#include "TROOT.h"
#include "Riostream.h"
#include <utility>
#include <vector>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <fstream>
#include <stdio.h>
#include "TString.h"
#include "TUnixSystem.h"
#include "TTree.h"
#include <map>


using namespace std;

void run_gif(){

  gSystem->Exec("mkdir gif");
  gROOT->ProcessLine("gSystem->Load(\"libDMRoot\")");
  gROOT->ProcessLine(".x dgr.C");
  cout << "Copy a line from bfcl8.txt" << endl;
  gSystem->Exec("emacs bfcl8copy.txt &");
}
