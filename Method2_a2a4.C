// ------------- Instructions -----------------------//
// In SetupAC, users should fill the vector arraynums with the crystal numbers that were present in the experiment.
// In SetupAC, users should fill the vector distances with the distances of the crystals that were present in the experiment.
// The Method2a2a4 function needs histograms of:
//    - data AC
//    - Z0 simulated AC
//    - Z2 simulated AC
//    - Z4 simulated AC

 
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

void Method2a2a4(TH1* datahst, TH1* Z0hst, TH1* Z2hst, TH1* Z4hst, bool visualization) {

   // debug flag - set to true for more output
   bool debug = false;

   if (debug) {cout <<"In Method2 function\n" <<endl;}

   // assigning the Z-distribution histograms to the global pointers
   Z0 = (TH1D*) Z0hst;
   Z2 = (TH1D*) Z2hst;
   Z4 = (TH1D*) Z4hst;
   hD = (TH1D*) datahst;

   // -------------------------------------------------------------------//
   //                       Preliminary fitting
   // -------------------------------------------------------------------//
   // this fit is quite sensitive to initial conditions,
   // so we'll do a quick fit to a bare Legendre polynomial
   // to get some initial estimates of a2 and a4.
   // -------------------------------------------------------------------//
   SetupAC();
   int nfitpars = 3; // number of parameters used in the fitting

   // the Legendre polynomial function is defined on a domain of -1 to 1,
   // so we need to create a graph real quick.
   TGraphAsymmErrors *acD = ac->CreateGraphFromHst(hD,false,false);

   TF1 *Method2Fit = new TF1("Method2Fit",TGRSIFunctions::LegendrePolynomial,-1,1,nfitpars);
   Method2Fit->SetParNames("A_{0}","a_{2}","a_{4}");
   Method2Fit->SetParameters(1,0.5,0.5);
   acD->Fit(Method2Fit,"QN0","",-1,1); // if you want to see the results printed, remove the Q option
   Method2Fit->SetParameter(0,1);

   // -------------------------------------------------------------------//
   //                                Fitting
   // -------------------------------------------------------------------//
   // We're going to include the uncertainty of the Monte Carlo simulation
   // in this chi^2 calculation. ROOT assumes all functions are perfectly
   // known with no error, so that means we need to dictate our own
   // minimization function. This takes a bit of time.
   // -------------------------------------------------------------------//
	TVirtualFitter *minuitZ = TVirtualFitter::Fitter(0,nfitpars);
   // initialize the parameters based on the preliminary fitting 
	for (int i = 0; i < nfitpars; ++i) {
		minuitZ->SetParameter(i, Method2Fit->GetParName(i), Method2Fit->GetParameter(i), 0.0001, -10, 10); // arguments are parameter number, parameter name, parameter initial value, value error, value minimum, and value maximum
	}
	minuitZ->SetFCN(ZFcn);
	double arglist[100];
	arglist[0] = 0;
	minuitZ->ExecuteCommand("SET PRINT",arglist,2);
	arglist[0] = 5000; // number of function calls
	arglist[1] = 0.0001; // tolerance
	minuitZ->ExecuteCommand("MIGRAD",arglist,2); // arguments here are minimization type (MIGRAD, SIMPLEX, MINIMIZE, SCAN, SEEK, MINOS), argument list, and number of arguments
	minuitZ->ExecuteCommand("MINOS",arglist,2); // arguments here are minimization type (MIGRAD, SIMPLEX, MINIMIZE, SCAN, SEEK, MINOS), argument list, and number of arguments
   // this additional call seems necessary to make the chi^2 value correct.
	minuitZ->ExecuteCommand("MIGRAD",arglist,2); // arguments here are minimization type (MIGRAD, SIMPLEX, MINIMIZE, SCAN, SEEK, MINOS), argument list, and number of arguments

   // -------------------------------------------------------------------//
   //                  Result extraction and printing
   // -------------------------------------------------------------------//
   // now the fitting is done and we can go about extracting the results
   // -------------------------------------------------------------------//
	double minParamsZ[3];
	double parErrorsZ[3];
	for (int i = 0; i < nfitpars; ++i) {
		minParamsZ[i] = minuitZ->GetParameter(i);
		parErrorsZ[i] = minuitZ->GetParError(i);
	}
   // most of these are dummy variables that GetStats needs
	double Zchi2, edm, errdef;
	int nvpar, nparx;
	minuitZ->GetStats(Zchi2,edm,errdef,nvpar,nparx);
	int ndf = npfits-nvpar; // degrees of freedom = number of points used in fitting - number of parameters used
   // the next three lines print results
	minuitZ->PrintResults(1,Zchi2);
   cout <<"NDF = " <<ndf <<endl;
	if(debug == true){std::cout << "Zchi2 is : " << Zchi2  << std::endl;}

   // the covariance matrix tells us how related our fitted parameters are.
   // for example, perhaps a_2 can vary greatly, but the best fit for a larger
   // a_2 value also requires a larger a_4 value. this would be a positive correlation.
   // the full covariance matrix deals with A_0, a_2, and a_4,
   // but we only care about a_2 and a_4, so we only extract some elements.
	TMatrixD matrixZ(2,2);
   for (int i=0;i<2;i++) {
      for (int j=0;j<2;j++) {
         matrixZ[i][j] = minuitZ->GetCovarianceMatrixElement(i+1,j+1);
      }
   }
   cout <<"---------- Covariance matrix ----------" <<endl;
   matrixZ.Print();
   cout <<"---------------------------------------" <<endl;
	if(debug){std::cout << "I MADE IT HERE"  << std::endl;}

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

      // clone Z distributions so we can manipulate them without changing
      // the originals
   	TH1D* Z = (TH1D*) Z0->Clone();
   	TH1D* Z_2 = (TH1D*) Z2->Clone();
   	TH1D* Z_4 = (TH1D*) Z4->Clone();
   
      // extract fitted scaling parameters and scale Z distributions
   	double A0 = minParamsZ[0];
   	double a2 = minParamsZ[1];
   	double a4 = minParamsZ[2];
      if (debug) cout <<"A0: " <<A0 <<"\ta2: " <<a2 <<"\ta4: " <<a4 <<endl;
   	Z->Scale(1.0-a2-a4);
   	Z->Add(Z_2,a2);
   	Z->Add(Z_4,a4);
   	Z->Scale(A0);// This is now the fitted simulated distribution.
   
      // In the next loop, we'll calculate the X2 again (as a check) and
      // extract bin-by-bin errors so we can use them for the residual
      // plot (where the error bars will include the data and simulated
      // uncertainties, added in quadrature).
      double X2 = 0;
   	double yerror[52];
   	double diff[52];
      for (int i=0;i<52;i++) {diff[i]=0;yerror[i]=0;} // initialization of arrays
   	if(debug){std::cout << "//////////////////////////////////Z-DISTRIBUTION CHI SQUARES/////////////////////////////"  << std::endl;}
   	for(int i=0;i<52;i++){
         int bin=i+1;
   		double data = hD->GetBinContent(bin);
         if (data<=0) continue;
   		double data_error2 = hD->GetBinError(bin);
         if (data_error2<=0) continue;
   		double z = Z->GetBinContent(bin);
   		double z_error = Z->GetBinError(bin);
   
         // this is used for the residual plot later
   		yerror[i] = TMath::Sqrt((data_error2*data_error2) + (z_error*z_error));
         if (yerror[i]==0) continue;
   		diff[i] = data - z;
   		if(debug == true){std::cout <<i <<  "\t" << diff[i] <<"\t" <<yerror[i] << std::endl;}
   		if(debug == true){std::cout <<i <<  " Chi square is = " << pow(diff[i]/yerror[i],2) << std::endl;}
   		X2 += pow(diff[i]/yerror[i],2);
   	}
   	if (debug) {std::cout << "/////////////////////////////////////////////////////////////////////////////////////////"  << std::endl;}
   	if (debug) cout <<"My calculated chi^2 = " <<X2 <<". This should be equal to FCN above, which should be " <<Zchi2 <<". Reduced is = " << X2/ndf <<endl;
  
      // multi_Zplot will have the acD (data) and Zac (simulated) TGraphs
   	TMultiGraph* multi_Zplot = new TMultiGraph();

   	TGraphAsymmErrors *Zac = ac->CreateGraphFromHst(Z,false,false);
   	Zac->SetMarkerStyle(8);
   	Zac->SetMarkerColor(kRed);
   	Zac->SetLineColor(kRed);
   	Zac->SetFillColor(kRed);
   	Zac->GetXaxis()->SetRangeUser(-1.1,1.1);
   	multi_Zplot->Add(Zac,"l3"); // the l3 option paints the histogram as a filled contour between the upper and lower error bars

   	acD->SetMarkerStyle(8);
   	multi_Zplot->Add(acD,"p");
   	multi_Zplot->SetTitle(
   			";;Normalized Counts;"
   			);
   
      // This text box will display the fit statistics
   	TPaveText* pt_Z = new TPaveText(0.37,0.5,0.8,0.8);
   	pt_Z->SetTextFont(133);
   	pt_Z->SetTextSize(20);
   	pt_Z->SetFillStyle(0);
   	pt_Z->SetBorderSize(0);
   	pt_Z->AddText(Form("a_{2} = %f #pm %f",a2,parErrorsZ[1]));
   	pt_Z->AddText(Form("a_{4} = %f #pm %f",a4,parErrorsZ[2]));
   	pt_Z->AddText(Form("#chi^{2}/NDF = %.2f",(TMath::Nint(Zchi2/ndf * 100) / 100.0)));
   
      // making the residual plot
      TGraphErrors* Zres = new TGraphErrors();
   	for (int i=0;i<52;i++)
   	{
         if (diff[i]==0) continue;
   		double angle = ac->GetAngleFromIndex(i);
   		double cosine = TMath::Cos(angle); // no deg->rad needed
         int point = Zres->GetN();
         Zres->SetPoint(point,cosine,diff[i]);
         Zres->SetPointError(point,0,yerror[i]);
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
      TPad* pads[2];
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
      multi_Zplot->Draw("ap");
   	pt_Z->Paint("NDC");
      c1b->Modified();
      c1b->Update();
   	multi_Zplot->GetXaxis()->SetLabelFont(133);
   	multi_Zplot->GetYaxis()->SetTitleFont(133);
   	multi_Zplot->GetYaxis()->SetLabelFont(133);
   	multi_Zplot->GetXaxis()->SetLabelSize(0);
   	multi_Zplot->GetYaxis()->SetLabelSize(labelsize);
   	multi_Zplot->GetYaxis()->SetTitleSize(titlesize);
   	multi_Zplot->GetYaxis()->CenterTitle(kTRUE);
   	multi_Zplot->GetYaxis()->SetTitleOffset(1.7);
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
