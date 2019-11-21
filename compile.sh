#!bin/bash
#
#
cd src
root -l  <<EOC
gROOT->ProcessLine(".L DMPlsAnalyzer.C+"); 
gROOT->ProcessLine(".L DMPlsFuncImages.C+"); 
gROOT->ProcessLine(".L DMPlsGrAnalyzer.C+"); 
gROOT->ProcessLine(".L DMPlsMtAnalyzer.C+"); 
gROOT->ProcessLine(".L DMPlsOutput.C+"); 
.q
EOC
mv *.so *.pcm ../lib
cd ..
