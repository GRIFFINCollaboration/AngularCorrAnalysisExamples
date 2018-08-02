// ------------- Instructions -----------------------//
// In SetupAC, users should fill the vector arraynums with the crystal numbers that were present in the experiment.
// In SetupAC, users should fill the vector distances with the distances of the crystals that were present in the experiment.
// The Method4a2a4 function needs TGraphAsymmErrors of the data AC, with x-axis of cos(theta).
// The Method4a2a4 function needs beta/gamma values and associated errors. These will be specific to the energies
// of the particular cascade.

 
// --------- TAngularCorrelation setup -----------//
TAngularCorrelation *ac;

void SetupAC(){
   // create vectors for my detector configuration
   vector<Int_t> arraynums, distances;
   for (int i=1;i<=64;i++)
   {
      if ( i==4 || i==8 || i==12 || i==16 || (i<53 && i>48) ) continue;
      arraynums.push_back(i);
      distances.push_back(110); // distance of crystal #i in mm
   }
   ac = new TAngularCorrelation();
   ac->GenerateMaps(arraynums,distances);
}

// --------- Formatting setup -----------//
int labelsize, titlesize,yNdivisions;
void StandardFormatting() {
	gStyle->SetOptStat(0);

	gStyle->SetTitleFont(132,"xyz");
	gStyle->SetLabelFont(132,"xyz");
	gStyle->SetTextFont(132);
	gStyle->SetLegendFont(132);

	labelsize = 20; // should be on the order of 12 or larger
	titlesize = 20; // should be on the order of 12 or larger
	yNdivisions = 503; // for the non-residual plots
}

// --------- Fitting functions --------- //
TH1D *Z0, *Z2, *Z4, *hD;

Double_t Zfit(Double_t *x, Double_t *p){
	Double_t i = x[0];

	Int_t bin_0 = Z0->GetXaxis()->FindBin(i);
	Double_t z_0 = (1-p[1]-p[2])*Z0->GetBinContent(bin_0);

	Int_t bin_2 = Z2->GetXaxis()->FindBin(i);
	Double_t z_2 = p[1]*Z2->GetBinContent(bin_2);

	Int_t bin_4 = Z4->GetXaxis()->FindBin(i);
	Double_t z_4 = p[2]*Z4->GetBinContent(bin_4);

	Double_t out = p[0]*(z_0 + z_2 + z_4);
	return out;
}

int npfits; // number of points that are used in fitting
void ZFcn(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
	TAxis *xaxis1  = hD->GetXaxis();

	int nbinX1 = hD->GetNbinsX();
	double chi2 = 0;
	double content,data_error,func_val,error,tmp;
	double x[1];
	double z_0error,z_2error,z_4error,Zerror;
   // scaling factors for each Z-component
   double dZ0 = p[0]*(1-p[1]-p[2]);
	double dZ2 = p[0]*p[1];
	double dZ4 = p[0]*p[2];
   npfits = 0;
	for (int bin = 1; bin <= nbinX1; ++bin) {
      // get data values
		content = hD->GetBinContent(bin);
      if (content<=0) continue;
		data_error = hD->GetBinError(bin);
		if (data_error<=0) continue;

      // get simulation values
		x[0] = xaxis1->GetBinCenter(bin);
		func_val = Zfit(x,p);
		z_0error = Z0->GetBinError(bin);
		z_2error = Z2->GetBinError(bin);
		z_4error = Z4->GetBinError(bin);

      // calculate errors
      // simulation error
		Zerror = sqrt(pow(dZ0*z_0error,2)+pow(dZ2*z_2error,2)+pow(dZ4*z_4error,2));
      // total error
		error = sqrt(pow(data_error,2) + pow(Zerror,2));

      // calculate chi^2
		tmp = (content-func_val)/error;
		chi2 += tmp*tmp;
      ++npfits;
	}
	fval = chi2; // final value
}
// ------------------------------------- //

void Method4a2a4(TGraphAsymmErrors* graph, double beta, double betaerr, double gamma, double gammaerr, bool visualization=false) {

   // debug flag - set to true for more output
   bool debug = false;

   if (debug) {cout <<"In Method4 function\n" <<endl;}

   // -------------------------------------------------------------------//
   //                       Fitting
   // -------------------------------------------------------------------//

   SetupAC();

   // the Legendre polynomial function is defined on a domain of -1 to 1,
   TF1 *Method4Fit = new TF1("Method4Fit",TGRSIFunctions::LegendrePolynomial,-1,1,3);
   Method4Fit->SetParNames("A_{0}","c_{2}","c_{4}");
   Method4Fit->SetParameters(1,0.5,0.5);
   TFitResultPtr result;
   result = graph->Fit(Method4Fit,"EMRS","",-1,1);

   // -------------------------------------------------------------------//
   //                  Result extraction and printing
   // -------------------------------------------------------------------//
	double minParams[3];
	double parErrors[3];
	for (int i = 0; i < 3; ++i) {
		minParams[i] = Method4Fit->GetParameter(i);
		parErrors[i] = Method4Fit->GetParError(i);
	}
   // most of these are dummy variables that GetStats needs
	double Zchi2 = Method4Fit->GetChisquare();
	int ndf = Method4Fit->GetNDF(); // degrees of freedom = number of points used in fitting - number of parameters used
   // the next three lines print results
   cout <<"Chi^2/NDF = " <<Zchi2 <<" / " <<ndf <<" = " <<Zchi2/ndf <<endl;
	if(debug == true){std::cout << "Zchi2 is : " << Zchi2  << std::endl;}

   // the covariance matrix tells us how related our fitted parameters are.
   // for example, perhaps a_2 can vary greatly, but the best fit for a larger
   // a_2 value also requires a larger a_4 value. this would be a positive correlation.
   // the full covariance matrix deals with A_0, a_2, and a_4,
   // but we only care about a_2 and a_4, so we only extract some elements.
	double *covarmatrixZ = result->GetCovarianceMatrix();
	TMatrixD matrixZ(2,2);
	matrixZ[0][0] = covarmatrixZ[4];
	matrixZ[0][1] = covarmatrixZ[5];
	matrixZ[1][0] = covarmatrixZ[7];
	matrixZ[1][1] = covarmatrixZ[8];
   cout <<"---------- c2/c4 covariance matrix ----------" <<endl;
   matrixZ.Print();
   cout <<"---------------------------------------" <<endl;
	if(debug){std::cout << "I MADE IT HERE"  << std::endl;}

   // -------------------------------------------------------------------//
   //                   Convert from c2/c4 --> a2/a4
   //                      and do error analysis
   // -------------------------------------------------------------------//
   // We'll do this with the user-supplied beta and gamma values.
   // Note: This is not complete!!! What more do we need?
   // -JKS, 2 August 2018
   // -------------------------------------------------------------------//

   double a2 = minParams[1]/beta;
   double a4 = minParams[2]/gamma;
   double a2err = parErrors[1]/beta;
   double a4err = parErrors[2]/gamma;
   cout <<"a2: " <<a2 <<" +/- " <<a2err <<endl;
   cout <<"a4: " <<a4 <<" +/- " <<a4err <<endl;

   // -------------------------------------------------------------------//
   //                            Visualization
   // -------------------------------------------------------------------//
   // This section is optional and will make a 2-paned figure with a
   // comparison of the data and the fit, the fit statistics, and the
   // residual. You might use this to check the fit visually - you might
   // use it as a starting point for a paper/thesis figure.
   // -------------------------------------------------------------------//
   if (visualization) {
      // set up some sizes and formatting for output
      StandardFormatting();

      // The first pane will have the graph (data) and the fit (which is already attached to the TGraphAsymmErrors*)
   	graph->SetMarkerStyle(8);
   	graph->SetTitle(
   			";;Normalized Counts;"
   			);
   
      // This text box will display the fit statistics
   	TPaveText* pt_Z = new TPaveText(0.37,0.5,0.8,0.8);
   	pt_Z->SetTextFont(133);
   	pt_Z->SetTextSize(20);
   	pt_Z->SetFillStyle(0);
   	pt_Z->SetBorderSize(0);
   	pt_Z->AddText(Form("a_{2} = %f #pm %f",minParams[1],parErrors[1]));
   	pt_Z->AddText(Form("a_{4} = %f #pm %f",minParams[2],parErrors[2]));
   	pt_Z->AddText(Form("#chi^{2}/NDF = %.2f",(TMath::Nint(Zchi2/ndf * 100) / 100.0)));
   
      // making the residual plot
      TGraphAsymmErrors* Zres = new TGraphAsymmErrors();
      double x,y,yerrhigh,yerrlow,value,diff;
   	for (int i=0;i<graph->GetN();i++)
   	{
         graph->GetPoint(i,x,y);
         value = Method4Fit->Eval(x);
         yerrhigh = graph->GetErrorYhigh(i);
         yerrlow = graph->GetErrorYlow(i);
         diff = y-value;
         int point = Zres->GetN();
         Zres->SetPoint(point,x,diff);
         Zres->SetPointError(point,0,0,yerrlow,yerrhigh);
   	}
   	Zres->SetMarkerStyle(8);
   	Zres->RemovePoint(0);
   	Zres->SetTitle(
   			""
   			"cos(#theta);"
   			"Normalized Counts;"
   		      );
   
      // this canvas is prepared specially for two differently
      // sized pads: a larger one (pads[1]) for the data and sim
      // and a smaller one (pads[0]) below that for the residual.
   	TCanvas* c1b = new TCanvas("c2","Z-distribution fit",500,500);
   	pads[0] = new TPad("pad0","pad0",0,0,1,0.3);
   	pads[1] = new TPad("pad1","pad1",0,0.3,1,1);
   	pads[1]->SetBottomMargin(0);
   	pads[1]->SetTopMargin(0.01);
   	pads[1]->SetLeftMargin(0.17);
   	pads[1]->SetRightMargin(0.02);
   	pads[0]->SetBottomMargin(0.22);
   	pads[0]->SetTopMargin(0);
   	pads[0]->SetLeftMargin(0.17);
   	pads[0]->SetRightMargin(0.02);
   	pads[0]->Draw();
   	pads[1]->Draw();
   	pads[0]->SetTickx(1);
   	pads[0]->SetTicky(1);
   	pads[1]->SetTickx(1);
   	pads[1]->SetTicky(1);
   
   	pads[1]->cd();
      graph->Draw("ap");
   	pt_Z->Paint("NDC");
      c1b->Modified();
      c1b->Update();
   	graph->GetXaxis()->SetLabelFont(133);
   	graph->GetYaxis()->SetTitleFont(133);
   	graph->GetYaxis()->SetLabelFont(133);
   	graph->GetXaxis()->SetLabelSize(0);
   	graph->GetYaxis()->SetLabelSize(labelsize);
   	graph->GetYaxis()->SetTitleSize(titlesize);
   	graph->GetYaxis()->CenterTitle(kTRUE);
   	graph->GetYaxis()->SetTitleOffset(1.7);
   	pt_Z->Draw();
   
   	pads[0]->cd();
   	Zres->Draw("ap");
   	Zres->SetMarkerSize(0.8);
   	Zres->GetYaxis()->SetNdivisions(305);
   	Zres->SetTitle(";cos(#theta);Residual");
   	Zres->Draw("pa");
   	Zres->GetXaxis()->SetTitleOffset(2.4);
   	Zres->GetYaxis()->SetTitleOffset(1.7);
   	Zres->GetXaxis()->CenterTitle(kTRUE);
   	Zres->GetYaxis()->CenterTitle(kTRUE);
   	Zres->GetXaxis()->SetLabelFont(133);
   	Zres->GetYaxis()->SetTitleFont(133);
   	Zres->GetYaxis()->SetLabelFont(133);
   	Zres->GetYaxis()->SetLabelSize(labelsize);
   	Zres->GetYaxis()->SetTitleSize(titlesize);
   	Zres->GetYaxis()->CenterTitle(kTRUE);
   	Zres->GetXaxis()->SetTitleFont(133);
   	Zres->GetXaxis()->SetLabelSize(labelsize);
   	Zres->GetXaxis()->SetTitleSize(titlesize);
   	Zres->GetXaxis()->CenterTitle(kTRUE);
   	Zres->GetXaxis()->SetTickLength(0.07);
   	pads[0]->SetTicks(1,1);
   	TLine* line = new TLine(-1.1,0,1.1,0);
   	line->Draw("same");
   }
}
