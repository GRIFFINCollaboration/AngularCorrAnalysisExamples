{
   // load library
   gROOT->ProcessLine(".L AnalysisFunctions.C");

   // pull data histogram from file
   TH1D* data = (TH1D*) gFile->Get("ac_2013_586");

   // convert data histogram to graph
   SetupAC();
   TGraphAsymmErrors* graph = ac->CreateGraphFromHst(ac_2013_586,false,false);

   // 0.9560 and 0.849 are the 16-clover parametrized beta and gamma for 586-2013
   // 0.0008 and 0.002 are the 16-clover parametrized beta and gamma uncertainties for 586-2013
   Method4a2a4(graph,0.9560,0.0008,0.849,0.002,false,true);
   //Method4mixing(graph,0.9560,0.0008,0.849,0.002,"results4.txt",1,3,7);
}
