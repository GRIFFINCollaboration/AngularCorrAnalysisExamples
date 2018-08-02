// ------------- Instructions -----------------------//
// - In SetupAC, users should fill the vector arraynums with the crystal numbers that were present in the experiment.
// - In SetupAC, users should fill the vector distances with the distances of the crystals that were present in the experiment.
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

void Method4mixing(TGraphAsymmErrors* graph, double beta, double betaerr, double gamma, double gammaerr, const char* outputfile) {

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
   TFitResultPtr result;

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
   double j1 = 0.5;
   double j2 = 1.5;
   double j3 = 3.5;
   // l1 is the transition between j1 and j2
   // a is the lowest allowed spin
   // b is the mixing spin
   double l1a = 1;
   double l1b = 2;
   // l2 is the transition between j2 and j3
   // a is the lowest allowed spin
   // b is the mixing spin
   double l2a = 2;
   double l2b = 3;

   double delta1, delta2; // mixing ratios that I will vary
   double a2,a4,a2err,a4err; // a2 and a4 values
   double c2,c4,c2err,c4err; // c2 and c4 values

   // in this example, I happen to know that a4 is always zero because j2<2
   // the following lines will fix that value
   a4 = 0;
   a4err = 0;
   c4 = a4*gamma;
   c4err = a4err*gammaerr;
	Method4Fit->FixParameter(2,a4);

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
   double mixanglemin = -TMath::Pi()/2;
   double mixanglemax = TMath::Pi()/2;
   int steps = 100;
   double stepsize = (mixanglemax-mixanglemin)/double(steps);
	double Zchi2;
   ofstream outfile;
   outfile.open(outputfile);
   for (int i=0;i<steps;i++) {
      double mixangle1 = mixanglemin + i*stepsize;
      double delta1 = TMath::Tan(mixangle1);
      if (debug) cout <<mixangle1 <<"\t" <<delta1 <<endl;
      for (int j=0;j<steps;j++) {
         double mixangle2 = mixanglemin + j*stepsize;
         double delta2 = TMath::Tan(mixangle2);
         if (debug) cout <<"\t" <<mixangle2 <<"\t" <<delta2 <<endl;
         // calculate a2
         a2 = calculate_a2(j1,j2,j3,l1a,l1b,l2a,l2b,delta1,delta2);
         c2 = a2*beta;
         // fix a2 in Method4Fit
         Method4Fit->FixParameter(1,c2);
         // fit scaling factor to data
         result = graph->Fit(Method4Fit,"QN0EMRS","",-1,1);
         // extract chi^2
         Zchi2 = Method4Fit->GetChisquare();
         // output to file
         outfile <<i <<"\t" <<mixangle1 <<"\t" <<j <<"\t" <<mixangle2 <<"\t" <<Zchi2 <<endl;
      }
   }
   outfile.close();
}
