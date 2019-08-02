{
   // load library
   gROOT->ProcessLine(".L AnalysisFunctions.C");

   // pull histograms from file
   TH1D* data = (TH1D*) gFile->Get("ac_2013_586");
   TH1D* Z0hst = (TH1D*) gFile->Get("Z0hst_noW14");
   TH1D* Z2hst = (TH1D*) gFile->Get("Z2hst_noW14");
   TH1D* Z4hst = (TH1D*) gFile->Get("Z4hst_noW14");

   // check to make sure all histograms were found
   if (!data) { cout <<"Data histogram not found. Exiting..." <<endl; break; }
   if (!Z0hst) { cout <<"Z0 histogram not found. Exiting..." <<endl; break; }
   if (!Z2hst) { cout <<"Z2 histogram not found. Exiting..." <<endl; break; }
   if (!Z4hst) { cout <<"Z4 histogram not found. Exiting..." <<endl; break; }
   // if all found, continue...

   Method2a2a4(data,Z0hst,Z2hst,Z4hst,false,false,true); // boolean order is: folded, fixa4, visualization
//   Method2mixing(data,Z0hst,Z2hst,Z4hst,"results2.txt",1,3,7);
}
