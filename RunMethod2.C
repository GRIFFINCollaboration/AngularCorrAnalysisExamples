{
//   gROOT->ProcessLine(".L Method2_a2a4.C");
   gROOT->ProcessLine(".L Method2_mixing.C");
   TH1D* ac_2013_586 = (TH1D*) gFile->Get("ac_2013_586");
   TH1D* ac_2013_565 = (TH1D*) gFile->Get("ac_2013_565");
   TH1D* Z0hst_noW14 = (TH1D*) gFile->Get("Z0hst_noW14");
   TH1D* Z2hst_noW14 = (TH1D*) gFile->Get("Z2hst_noW14");
   TH1D* Z4hst_noW14 = (TH1D*) gFile->Get("Z4hst_noW14");
//   Method2a2a4(ac_2013_586,Z0hst_noW14,Z2hst_noW14,Z4hst_noW14,false);
   Method2mixing(ac_2013_586,Z0hst_noW14,Z2hst_noW14,Z4hst_noW14,"results2.txt");
}
