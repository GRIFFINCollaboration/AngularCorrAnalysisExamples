{
   gROOT->ProcessLine(".L Method4_a2a4.C");
   //gROOT->ProcessLine(".L Method4_mixing.C");
   TH1D* ac_2013_586 = (TH1D*) gFile->Get("ac_2013_586");
   TH1D* ac_2013_565 = (TH1D*) gFile->Get("ac_2013_565");
   SetupAC();
   TGraphAsymmErrors* graph = ac->CreateGraphFromHst(ac_2013_586,false,false);
   Method4a2a4(graph,0.9560,0.0008,0.849,0.002,true); // these values are the 16-clover parametrized beta and gamma for 586-2013
   //Method4mixing(graph,0.9560,0.0008,0.849,0.002,"results4.txt"); // these values are the 16-clover parametrized beta and gamma for 586-2013
}
