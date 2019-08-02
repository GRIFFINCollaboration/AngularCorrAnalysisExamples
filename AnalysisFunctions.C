// ------------- Instructions -----------------------//
// - In SetupAC, users should fill the vector arraynums with the crystal numbers that were present in the experiment.
// - In SetupAC, users should fill the vector distances with the distances of the crystals that were present in the experiment.
// - The Method2a2a4 function needs histograms of:
//    - data AC
//    - Z0 simulated AC
//    - Z2 simulated AC
//    - Z4 simulated AC
// - In Method2mixing, users should include whatever physics information they have about the cascade.
//   This could include spins, angular momenta, or mixing ratios. This function will change significantly
//   depending on the physics of a given cascade.
// - The Method2mixing function needs histograms of:
//    - data AC
//    - Z0 simulated AC
//    - Z2 simulated AC
//    - Z4 simulated AC
// - The Method4a2a4 function needs TGraphAsymmErrors of the data AC, with x-axis of cos(theta).
// - The Method4a2a4 function needs beta/gamma values and associated errors. These will be specific to the energies
//   of the particular cascade.
// - In Method4mixing, users should include whatever physics information they have about the cascade.
//   This could include spins, angular momenta, or mixing ratios. This function will change significantly
//   depending on the physics of a given cascade.
// - The Method4mixing function needs TGraphAsymmErrors of the experimental AC with an x-axis of cos(theta).
// - The Method4mixing function needs beta/gamma values and associated errors. These will be specific to the energies
//   of the particular cascade.

 
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

// --------- Physics functions --------- //
double ClebschGordan(double j1,double m1,double j2,double m2,double j,double m){

    double term=0., cg=0., term1=0., sum=0.;
    int k;
    // Conditions check
    if( 2*j1 != floor(2*j1) || 
        2*j2 != floor(2*j2) || 
        2*j != floor(2*j) || 
        2*m1 != floor(2*m1) || 
        2*m2 != floor(2*m2) || 
        2*m != floor(2*m) ){
//	printf("All arguments must be integers or half-integers.\n");
        return 0;
    }
    if((m1+m2) != m){
        //printf("m1 + m2 must equal m.\n");
        return 0;
    }
    if((j1-m1) != floor(j1-m1)){
        //printf("2*j1 and 2*m1 must have the same parity");
        return 0;
    }
    if((j2-m2) != floor(j2-m2)){
        //printf("2*j2 and 2*m2 must have the same parity");
        return 0;
    }
    if( j - m != floor(j-m) ){
        //printf("2*j and 2*m must have the same parity");
        return 0;
    }
    if(j>(j1+j2) || j < abs(j1-j2)){
        //printf("j is out of bounds.");
        return 0;
    }
    if(abs(m1) > j1){
        //printf("m1 is out of bounds.");
        return 0;
    }
    if(abs(m2) > j2){
        //printf("m2 is out of bounds.");
        return 0;
    }
    if(abs(m) > j){
        //printf("m is out of bounds.\n");
        return 0 ;
    }
    //
    term1 = pow((((2*j+1)/TMath::Factorial(j1+j2+j+1))*TMath::Factorial(j2+j-j1)*TMath::Factorial(j+j1-j2)*TMath::Factorial(j1+j2-j)*TMath::Factorial(j1+m1)*TMath::Factorial(j1-m1)*TMath::Factorial(j2+m2)*TMath::Factorial(j2-m2)*TMath::Factorial(j+m)*TMath::Factorial(j-m)),(0.5));
    sum = 0;
    
    for(k = 0 ; k <= 99 ; k++ ){
        if( (j1+j2-j-k < 0) || (j1-m1-k < 0) || (j2+m2-k < 0) ){
            //no further terms will contribute to sum, exit loop
            break;
        } else if( (j-j1-m2+k < 0) || (j-j2+m1+k < 0)  ){
            //jump ahead to next term that will contribute
	    const Int_t a1 = (j-j1-m2);
	    const Int_t a2 = (j-j2+m1);
            k = TMath::Max(-TMath::Min(a1, a2) - 1, k);
        }else{
            term = TMath::Factorial(j1+j2-j-k)*TMath::Factorial(j-j1-m2+k)*TMath::Factorial(j-j2+m1+k)*TMath::Factorial(j1-m1-k)*TMath::Factorial(j2+m2-k)*TMath::Factorial(k);
            if((k%2) == 1){
                term = -1*term;
            }
            sum = sum + 1.0/term;
        }
    }
    cg = term1*sum;
    return cg;
    // Reference: An Effective Algorithm for Calculation of the C.G.
    // Coefficients Liang Zuo, et. al.
    // J. Appl. Cryst. (1993). 26, 302-304
}

double RacahW(double a, double b,double c,double d, double e, double f){
    return TMath::Power((-1),int(a+b+d+c))*ROOT::Math::wigner_6j(int(2*a),int(2*b),int(2*e),int(2*d),int(2*c),int(2*f));
}

double F(double k,double jf,double L1,double L2,double ji){
    double W;
    double CG = ClebschGordan(L1,1,L2,-1,k,0);

    if(CG == 0){
        return 0;
    }
    W = RacahW(ji,ji,L1,L2,k,jf);
    if(W == 0){
        return 0;
    }
    return TMath::Power((-1),(jf-ji-1))*(TMath::Power((2*L1+1)*(2*L2+1)*(2*ji+1),(1.0/2.0)))*CG*W;
    // Reference: Tables of coefficients for angular distribution of gamma rays from aligned nuclei
    // T. Yamazaki. Nuclear Data A, 3(1):1?23, 1967.
}

double A(double k,double ji,double jf,double L1,double L2,double delta){
    double f1 = F(k,ji,L1,L1,jf);
    double f2 = F(k,ji,L1,L2,jf);
    double f3 = F(k,ji,L2,L2,jf);

    return (1/(1+TMath::Power(delta,2)))*(f1+2*delta*f2+delta*delta*f3);
}

double B(double k,double ji,double jf,double L1,double L2,double delta){
    double f1 = F(k,jf,L1,L1,ji);
    double f2 = F(k,jf,L1,L2,ji);
    double f3 = F(k,jf,L2,L2,ji);
    return (1/(1+TMath::Power(delta,2)))*(f1+(TMath::Power((-1),((L1+L2))))*2*delta*f2+delta*delta*f3);
}

double calculate_a2(double j1,double j2,double j3,double l1a,double l1b,double l2a,double l2b,double delta1,double delta2){
    return B(2,j2,j1,l1a,l1b,delta1)*A(2,j3,j2,l2a,l2b,delta2);
}

double calculate_a4(double j1,double j2,double j3,double l1a,double l1b,double l2a,double l2b,double delta1,double delta2){
    return B(4,j2,j1,l1a,l1b,delta1)*A(4,j3,j2,l2a,l2b,delta2);
}
// ------------------------------------- //

void Method2a2a4(TH1* datahst, TH1* Z0hst, TH1* Z2hst, TH1* Z4hst, bool folded, bool fixa4, bool visualization) {

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
   TGraphAsymmErrors *acD = ac->CreateGraphFromHst(hD,folded,false);

   TF1 *Method2Fit = new TF1("Method2Fit",TGRSIFunctions::LegendrePolynomial,-1,1,nfitpars);
   Method2Fit->SetParNames("A_{0}","a_{2}","a_{4}");
   Method2Fit->SetParameters(1,0.5,0.5);
   if (fixa4) Method2Fit->FixParameter(2,0);
   acD->Fit(Method2Fit,"QN0","",-1,1); // if you want to see the results printed, remove the Q option
   Method2Fit->SetParameter(0,hD->GetMaximum()/Z0->GetMaximum()); // this sets the scale appropriately initially

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
   if (fixa4) { minuitZ->FixParameter(2); }
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
   cout <<"Reduced chi^2 is " <<Zchi2/ndf <<endl;
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

   	TGraphAsymmErrors *Zac = ac->CreateGraphFromHst(Z,folded,false);
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
// ------------------------------------- //

void Method2mixing(TH1* datahst, TH1* Z0hst, TH1* Z2hst, TH1* Z4hst, const char* outputfile, int twoJhigh, int twoJmid, int twoJlow) {

   // debug flag - set to true for more output
   bool debug = false;

   if (debug) {cout <<"In Method2 function\n" <<endl;}

   // assigning the Z-distribution histograms to the global pointers
   Z0 = (TH1D*) Z0hst;
   Z2 = (TH1D*) Z2hst;
   Z4 = (TH1D*) Z4hst;
   hD = (TH1D*) datahst;

   SetupAC();
   int nfitpars = 3; // number of parameters used in the fitting

   // -------------------------------------------------------------------//
   //                         Fitting setup
   // -------------------------------------------------------------------//
   // We're going to include the uncertainty of the Monte Carlo simulation
   // in this chi^2 calculation. ROOT assumes all functions are perfectly
   // known with no error, so that means we need to dictate our own
   // minimization function. This takes a bit of time.
   // -------------------------------------------------------------------//
	TVirtualFitter *minuitZ = TVirtualFitter::Fitter(0,nfitpars);
	double arglist[100];
	arglist[0] = -1;
	minuitZ->ExecuteCommand("SET PRINT",arglist,2);
   // initialize the parameters with arbitrary values 
	minuitZ->SetParameter(0, "A_{0}",0.5, 0.0001, -10, 10); // arguments are parameter number, parameter name, parameter initial value, value error, value minimum, and value maximum
	minuitZ->SetParameter(1, "a_{2}", 0.5, 0.0001, -10, 10); // arguments are parameter number, parameter name, parameter initial value, value error, value minimum, and value maximum
	minuitZ->SetParameter(2, "a_{4}", 0.5, 0.0001, -10, 10); // arguments are parameter number, parameter name, parameter initial value, value error, value minimum, and value maximum
	minuitZ->SetFCN(ZFcn);
	arglist[0] = 5000; // number of function calls
	arglist[1] = 0.0001; // tolerance


   // -------------------------------------------------------------------//
   //                         Physics setup
   // -------------------------------------------------------------------//
   // Here's where you include the relevant physics for your situation.
   // In this example, I know the three spins, but am varying the two
   // mixing ratios.
   // -------------------------------------------------------------------//
   
   // j1 is the spin of the highest level
   // j2 is the spin of the middle level
   // j3 is the spin of the bottom level
   double j1 = 0.5*twoJhigh;
   double j2 = 0.5*twoJmid;
   double j3 = 0.5*twoJlow;
   // l1 is the transition between j1 and j2
   // a is the lowest allowed spin
   // b is the mixing spin
   int l1a = TMath::Abs(twoJhigh-twoJmid)/2;
   if (l1a == 0) l1a = 1;
   int l1b = l1a + 1;
   // l2 is the transition between j2 and j3
   // a is the lowest allowed spin
   // b is the mixing spin
   int l2a = TMath::Abs(twoJmid-twoJlow)/2;
   if (l2a == 0) l2a = 1;
   int l2b = l2a + 1;

   double delta1, delta2; // mixing ratios that I will vary
   double a2,a4; // a2 and a4 values

   // -------------------------------------------------------------------//
   //                       Constrained fitting
   // -------------------------------------------------------------------//
   // Here's where you'll iterate over the physical quantities you don't
   // know, whether that be mixing ratios or spins or both. Here's an
   // example of two unknown mixing ratios but known spins.
   //
   // The basic idea is to select a particular set of physical quantities,
   // calculate a2/a4, fix the minuitZ parameters, and fit the scaling
   // factor A0. Then output the specifications for that set of physical
   // quantities and the chi^2 for further analysis.
   // -------------------------------------------------------------------//

   // delta runs from -infinity to infinity (unless constrained by known physics)
   // in this case, it then makes more sense to sample evenly from tan^{-1}(delta)
   // The next few lines are where you may want to include limits to significantly speed up calculations
   // mixing for the high-middle transition
   double mixanglemin1 = -TMath::Pi()/2;
   double mixanglemax1 = TMath::Pi()/2;
   int steps1 = 100;
   double stepsize1 = (mixanglemax1-mixanglemin1)/double(steps1);
   // mixing for the middle-low transition
   // to constrain this to zero, set mixanglemin2 to 0, mixanglemax2 to 1; and steps2 to 1
   double mixanglemin2 = -TMath::Pi()/2;
   double mixanglemax2 = TMath::Pi()/2;
   int steps2 = 100;
   double stepsize2 = (mixanglemax2-mixanglemin2)/double(steps2);

	double Zchi2, edm, errdef;
	int nvpar, nparx;
   ofstream outfile;
   outfile.open(outputfile);
   for (int i=0;i<steps1;i++) {
      double mixangle1 = mixanglemin1 + i*stepsize1;
      double delta1 = TMath::Tan(mixangle1);
      if (i%10==0) cout <<i <<"/" <<steps1 <<endl;
      if (debug) cout <<mixangle1 <<"\t" <<delta1 <<endl;
      for (int j=0;j<steps2;j++) {
         double mixangle2 = mixanglemin2 + j*stepsize2;
         double delta2 = TMath::Tan(mixangle2);
         if (debug) cout <<"\t" <<mixangle2 <<"\t" <<delta2 <<endl;
         // calculate a2
         a2 = calculate_a2(j1,j2,j3,l1a,l1b,l2a,l2b,delta1,delta2);
         // fix a2 in minuitZ
         minuitZ->ReleaseParameter(1);      
	      minuitZ->SetParameter(1, "a_{2}", a2, 0.0001, -10, 10); // arguments are parameter number, parameter name, parameter initial value, value error, value minimum, and value maximum
         minuitZ->FixParameter(1);
         // calculate a4
         a4 = calculate_a4(j1,j2,j3,l1a,l1b,l2a,l2b,delta1,delta2);
         // fix a4 in minuitZ
         minuitZ->ReleaseParameter(2);      
	      minuitZ->SetParameter(2, "a_{4}", a4, 0.0001, -10, 10); // arguments are parameter number, parameter name, parameter initial value, value error, value minimum, and value maximum
         minuitZ->FixParameter(2);
         // fit scaling factor to data
	      minuitZ->ExecuteCommand("MIGRAD",arglist,2); // arguments here are minimization type (MIGRAD, SIMPLEX, MINIMIZE, SCAN, SEEK, MINOS), argument list, and number of arguments
         // extract chi^2
	      minuitZ->GetStats(Zchi2,edm,errdef,nvpar,nparx);
         // output to file
         outfile <<i <<"\t" <<mixangle1 <<"\t" <<j <<"\t" <<mixangle2 <<"\t" <<Zchi2 <<endl;
      }
   }
   outfile.close();
}
// ------------------------------------- //

void Method4a2a4(TGraphAsymmErrors* graph, double beta, double betaerr, double gamma, double gammaerr, bool fixa4, bool visualization=false) {

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
  if (fixa4) Method4Fit->FixParameter(2,0);
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
  // for example, perhaps c_2 can vary greatly, but the best fit for a larger
  // c_2 value also requires a larger c_4 value. this would be a positive correlation.
  // the full covariance matrix deals with A_0, c_2, and c_4,
  // but we only care about c_2 and c_4, so we only extract some elements.
  TMatrixD matrixZ(2,2);
  for (int i=0;i<2;i++) {
    for (int j=0;j<2;j++) {
      matrixZ[i][j] = result->GetCovarianceMatrix()[i+1][j+1];
    }
  }
  cout <<"---------- c2/c4 covariance matrix ----------" <<endl;
  matrixZ.Print();
  cout <<"---------------------------------------" <<endl;
  if(debug){std::cout << "I MADE IT HERE"  << std::endl;}

  // -------------------------------------------------------------------//
  //                   Convert from c2/c4 --> a2/a4
  //                      and do error analysis
  // -------------------------------------------------------------------//
  // Error on a2,a4 parameters are given by the sqrt of the diagonal 
  // elements of the coviarance matrix defined in Eq. 27 of the NIM paper.
  // - ASC 15 August 2018
  // -------------------------------------------------------------------//

  // par[1] is c2, par[2] is c4
  double a2 = minParams[1]/beta;
  double a4 = minParams[2]/gamma;
  double a2err = TMath::Sqrt( TMath::Power(beta,-4.)*TMath::Power(minParams[1],2.)*TMath::Power(betaerr,2.) 
			      + TMath::Power(beta,-2.)*TMath::Power(parErrors[1],2.) );
  double a4err = TMath::Sqrt( TMath::Power(gamma,-4.)*TMath::Power(minParams[2],2.)*TMath::Power(gammaerr,2.) 
			      + TMath::Power(gamma,-2.)*TMath::Power(parErrors[2],2.) );
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
    pt_Z->AddText(Form("a_{2} = %f #pm %f",a2,a2err));
    pt_Z->AddText(Form("a_{4} = %f #pm %f",a4,a4err));
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
    TPad *pads[2];
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
// ------------------------------------- //

void Method4mixing(TGraphAsymmErrors* graph, double beta, double betaerr, double gamma, double gammaerr, const char* outputfile, int twoJhigh, int twoJmid, int twoJlow) {

   // debug flag - set to true for more output
   bool debug = false;

   if (debug) {cout <<"In Method4 function\n" <<endl;}

   SetupAC();
   int nfitpars = 3; // number of parameters used in the fitting

   // -------------------------------------------------------------------//
   //                         Fitting setup
   // -------------------------------------------------------------------//
   // the Legendre polynomial function is defined on a domain of -1 to 1,
   TF1 *Method4Fit = new TF1("Method4Fit",TGRSIFunctions::LegendrePolynomial,-1,1,3);
   Method4Fit->SetParNames("A_{0}","c_{2}","c_{4}");
   Method4Fit->SetParameters(1,0.5,0.5);

   // -------------------------------------------------------------------//
   //                         Physics setup
   // -------------------------------------------------------------------//
   // Here's where you include the relevant physics for your situation.
   // In this example, I know the three spins, but am varying the two
   // mixing ratios.
   // -------------------------------------------------------------------//
   
   // j1 is the spin of the highest level
   // j2 is the spin of the middle level
   // j3 is the spin of the bottom level
   double j1 = 0.5*twoJhigh;
   double j2 = 0.5*twoJmid;
   double j3 = 0.5*twoJlow;
   // l1 is the transition between j1 and j2
   // a is the lowest allowed spin
   // b is the mixing spin
   int l1a = TMath::Abs(twoJhigh-twoJmid)/2;
   if (l1a == 0) l1a = 1;
   int l1b = l1a + 1;
   // l2 is the transition between j2 and j3
   // a is the lowest allowed spin
   // b is the mixing spin
   int l2a = TMath::Abs(twoJmid-twoJlow)/2;
   if (l2a == 0) l2a = 1;
   int l2b = l2a + 1;

   double delta1, delta2; // mixing ratios that I will vary
   double a2,a4,a2err,a4err; // a2 and a4 values
   double c2,c4,c2err,c4err; // c2 and c4 values

   // -------------------------------------------------------------------//
   //                       Constrained fitting
   // -------------------------------------------------------------------//
   // Here's where you'll iterate over the physical quantities you don't
   // know, whether that be mixing ratios or spins or both. Here's an
   // example of two unknown mixing ratios but known spins.
   //
   // The basic idea is to select a particular set of physical quantities,
   // calculate a2/a4, fix the Method4Fit parameters, and fit the scaling
   // factor A0. Then output the specifications for that set of physical
   // quantities and the chi^2 for further analysis.
   // -------------------------------------------------------------------//

   // delta runs from -infinity to infinity (unless constrained by known physics)
   // in this case, it then makes more sense to sample evenly from tan^{-1}(delta)
   // The next few lines are where you may want to include limits to significantly speed up calculations
   // mixing for the high-middle transition
   double mixanglemin1 = -TMath::Pi()/2;
   double mixanglemax1 = TMath::Pi()/2;
   int steps1 = 100;
   double stepsize1 = (mixanglemax1-mixanglemin1)/double(steps1);
   // mixing for the middle-low transition
   // to constrain this to zero, set mixanglemin2 to 0, mixanglemax2 to 1; and steps2 to 1
   double mixanglemin2 = -TMath::Pi()/2;
   double mixanglemax2 = TMath::Pi()/2;
   int steps2 = 100;
   double stepsize2 = (mixanglemax2-mixanglemin2)/double(steps2);

	double Zchi2;
   ofstream outfile;
   outfile.open(outputfile);
   for (int i=0;i<steps1;i++) {
      double mixangle1 = mixanglemin1 + i*stepsize1;
      double delta1 = TMath::Tan(mixangle1);
      if (debug) cout <<mixangle1 <<"\t" <<delta1 <<endl;
      for (int j=0;j<steps2;j++) {
         double mixangle2 = mixanglemin2 + j*stepsize2;
         double delta2 = TMath::Tan(mixangle2);
         if (debug) cout <<"\t" <<mixangle2 <<"\t" <<delta2 <<endl;
         // calculate a2
         a2 = calculate_a2(j1,j2,j3,l1a,l1b,l2a,l2b,delta1,delta2);
         c2 = a2*beta;
         // fix a2 in Method4Fit
         Method4Fit->FixParameter(1,c2);
         // calculate a4
         a4 = calculate_a4(j1,j2,j3,l1a,l1b,l2a,l2b,delta1,delta2);
         c4 = a4*gamma;
         // fix a2 in Method4Fit
         Method4Fit->FixParameter(2,c4);
         // fit scaling factor to data
         graph->Fit(Method4Fit,"QN0ERS","",-1,1);
         // extract chi^2
         Zchi2 = Method4Fit->GetChisquare();
         // output to file
         outfile <<i <<"\t" <<mixangle1 <<"\t" <<j <<"\t" <<mixangle2 <<"\t" <<Zchi2 <<endl;
      }
   }
   outfile.close();
}
