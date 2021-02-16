#define NueSterileFit_C

#include "NueAna/MultiBinAna/NueSterileFit.h"
#include <iomanip>

MultiVarGauss::MultiVarGauss(TH1 *hMean, TH2 *hCov)
    : _generator(0)
{
  _generator = new TRandom3();
  _generator->SetSeed(0);

  if ((hMean->GetXaxis()->GetNbins() != hCov->GetXaxis()->GetNbins()) || (hCov->GetXaxis()->GetNbins() != hCov->GetYaxis()->GetNbins()))
  {
    cout << "hMean and hCov don't have the same number of bins! Aborting." << endl;
    exit(0);
  }

  _npar = hMean->GetXaxis()->GetNbins();

  _bestValues.ResizeTo(_npar);
  _cov.ResizeTo(_npar, _npar);

  double max_abs(-1000);

  for (int ipar = 0; ipar < _npar; ++ipar)
  {
    _bestValues(ipar) = hMean->GetBinContent(ipar + 1);
    for (int jpar = 0; jpar < _npar; ++jpar)
    {
      _cov(ipar, jpar) = hCov->GetBinContent(ipar + 1, jpar + 1);
      if (Abs(_cov(ipar, jpar)) > max_abs)
        max_abs = Abs(_cov(ipar, jpar));
    }
  }

  std::cout << "Cov matrix" << std::endl;
  _cov.Print();

  _hVaried = (TH1 *)hMean->Clone("VariedHisto");
  _hVaried->SetDirectory(0);

  DecomposeCov();

  return;
}

MultiVarGauss::MultiVarGauss(TH1 *hMean, TMatrixD *hCov)
    : _generator(0)
{
  _generator = new TRandom3();
  _generator->SetSeed(0);

  if ((hMean->GetXaxis()->GetNbins() != hCov->GetNrows()) || (hCov->GetNrows() != hCov->GetNcols()))
  {
    cout << "hMean and hCov don't have the same number of bins! Aborting." << endl;
    exit(0);
  }

  _npar = hMean->GetXaxis()->GetNbins();

  _bestValues.ResizeTo(_npar);
  _cov.ResizeTo(_npar, _npar);

  double max_abs(-1000);

  for (int ipar = 0; ipar < _npar; ++ipar)
  {
    _bestValues(ipar) = hMean->GetBinContent(ipar + 1);
    for (int jpar = 0; jpar < _npar; ++jpar)
    {
      _cov(ipar, jpar) = (*hCov)(ipar, jpar);
      if (Abs(_cov(ipar, jpar)) > max_abs)
        max_abs = Abs(_cov(ipar, jpar));
    }
  }

  std::cout << "Cov matrix" << std::endl;
  _cov.Print();

  _hVaried = (TH1 *)hMean->Clone("VariedHisto");
  _hVaried->SetDirectory(0);

  DecomposeCov();

  return;
}

void MultiVarGauss::DecomposeCov()
{
  TDecompSVD svd(_cov);
  bool ok = svd.Decompose();
  if (!ok)
  {
    cout << "SVD Decomposition failed!" << endl;
    exit(1);
  }
  TMatrixD U = svd.GetU();
  TMatrixD TU(TMatrixD::kTransposed, U);
  TVectorD Sig = svd.GetSig();
  TMatrixD SqrtS(_npar, _npar);
  for (int ipar = 0; ipar < _npar; ++ipar)
  {
    for (int jpar = 0; jpar < _npar; ++jpar)
    {
      if (ipar != jpar)
        SqrtS(ipar, jpar) = 0;
      else
        SqrtS(ipar, jpar) = sqrt(Sig(ipar));
    }
  }
  TMatrixD QR = SqrtS * TU;
  TDecompQRH qr(QR);
  ok = qr.Decompose();
  if (!ok)
  {
    cout << "QR Decomposition failed!" << endl;
    exit(1);
  }
  TMatrixD R = qr.GetR();
  TMatrixD TR(TMatrixD::kTransposed, R);
  TR.Print();
  _UT.ResizeTo(TR);
  _UT = TR;

  /*
      TDecompChol tdc(_cov);
      double tol = tdc.GetTol();
      //cout << "Original tolerance = " << tol << endl;
      //cout << "Max_Abs = " << max_abs << endl;
      //const Double_t scale = TMath::Min(max_abs, 1.);
      //tdc.SetTol(tol*scale);
      //cout << "New tolerance = " << tdc.GetTol() << endl;

      tdc.Decompose();
      TMatrixD U = tdc.GetU();
      //TMatrixD TU(U);
      TMatrixD TU(TMatrixD::kTransposed, U);
      std::cout << "TU matrix" << std::endl;
      TU.Print();

      _UT.ResizeTo(TU);
      _UT = TU;

  */
}

MultiVarGauss::~MultiVarGauss()
{
  delete _generator;
  delete _hVaried;
  return;
}

TH1 *MultiVarGauss::Generate()
{
  TVectorD nextValues;
  nextValues.ResizeTo(_npar);

  for (int ipar = 0; ipar < _npar; ++ipar)
  {
    nextValues(ipar) = _generator->Gaus(0.0, 1.0);
  }
  nextValues *= _UT;
  //for (int ipar = 0; ipar < _npar; ++ipar)
  //{
  //  cout << "bin " << ipar << " shift = " << nextValues(ipar) << endl;
  //}
  nextValues += _bestValues;

  for (int ipar = 0; ipar < _npar; ++ipar)
  {
    _hVaried->SetBinContent(ipar + 1, nextValues(ipar));
  }

  return _hVaried;
}

NueSterileFit::NueSterileFit() {
  NObs = 0;
  ErrCalc = 0;
  ErrorMatrix = 0;
  InvErrorMatrix = 0;
  ExternalErrorMatrix = 0;

  SetOutputFile();
  SetGridNorm();
  SetNExperiments();

  nBins = 0;
}

NueSterileFit::~NueSterileFit() {
}

void NueSterileFit::AddError(ErrorCalc *Err) {
  // References the ErrorCalc object
  ErrCalc = Err;
  return;
}

void NueSterileFit::AddExtrap(Extrapolate2D *E) {
  Extrap.push_back(E);
  return;
}

void NueSterileFit::SetNObs(TH1D *n) {
  NObs = (TH1D *)n->Clone("NObs");
  nBins = NObs->GetNbinsX();
  return;
}

void NueSterileFit::SetSinSqTh14Range(double l, double h) {
  if (l < 0 || l > 1) {
    std::cout << "Unphysical value of SinSqTh14.  Setting low value to 0." << std::endl;
    l = 0;
  }
  if (h < 0 || h > 1) {
    std::cout << "Unphysical value of SinSqTh14.  Setting high value to 1." << std::endl;
    h = 1;
  }

  SinSqTh14Low = l;
  SinSqTh14High = h;

  return;
}

void NueSterileFit::SetSinSqTh24Range(double l, double h) {
  if (l < 0 || l > 1) {
    std::cout << "Unphysical value of SinSqTh24.  Setting low value to 0." << std::endl;
    l = 0;
  }
  if (h < 0 || h > 1) {
    std::cout << "Unphysical value of SinSqTh24.  Setting high value to 1." << std::endl;
    h = 1;
  }

  SinSqTh24Low = l;
  SinSqTh24High = h;

  return;
}

void NueSterileFit::SetOutputFile(string s) {
  outFileName = s;
  return;
}

void NueSterileFit::DefineStdDlnLMinuit() { // Function sets the maximum size of Minuit.
  // Clear things first:
  // minuit->mncler();

  int npar = 0;

  // Number of systematics inputted:
  if (FracErr_Bkgd_List.size() != 0) {
    npar = FracErr_Bkgd_List.size();
  }

  int nb = 0; // Number of bins (for HOOHE):
  if (nBins != 0) {
    nb = nBins;
  } else if (NObs != 0) {
    nb = NObs->GetNbinsX();
  } else if (Extrap.size() != 0) {
    nb = Extrap[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  } else {
    std::cout << "ERROR: Add extrapolation or NObs before initializing Minuit size." << std::endl;
    return;
  }

  npar += nb;

  // make new minuit
  minuit = new TMinuit(npar + nb);
  minuit->SetPrintLevel(-1);
}

void NueSterileFit::DefineBinDlnLMinuit() { // Function sets the maximum size of Minuit.
  // Clear things first:
  // minuit->mncler();

  int npar = 0;
  if (nBins != 0) {
    npar = nBins;
  } else if (NObs != 0) {
    npar = NObs->GetNbinsX();
  } else if (Extrap.size() != 0) {
    npar = Extrap[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  } else {
    std::cout << "ERROR: Add extrapolation or NObs before initializing Minuit size." << std::endl;
    return;
  }

  // Make new minuit:
  minuit = new TMinuit(npar);
  minuit->SetPrintLevel(-1);

  double arglist[1];
  int ierflg = 0;

  arglist[0] = 1.0E-3;
  minuit->mnexcm("SET EPS", arglist, 1, ierflg);
}

double NueSterileFit::PoissonChi2(TH1D *nexp) {
  double chi2 = 0.;
  double no, ne;

  for (unsigned int i = 0; i < nBins; i++) {
    ne = nexp->GetBinContent(i + 1);
    no = NObs->GetBinContent(i + 1);
    if (no > 0) {
      chi2 += (2 * (ne - no + no * TMath::Log(no / ne)));
    } else if (no == 0) {
      chi2 += (2 * ne);
    }
    // chi2 is undefined if ne is 0
  }

  return chi2;
}

double NueSterileFit::ScaledChi2(TH1D *nexp_bkgd, TH1D *nexp_signal) {
  double chi2 = 0;
  double nb, ns, nobs, eb, es, nexp, errscale;

  for (unsigned int i = 0; i < nBins; i++) {
    nb = nexp_bkgd->GetBinContent(i + 1);
    ns = nexp_signal->GetBinContent(i + 1);
    nexp = nb + ns;
    nobs = NObs->GetBinContent(i + 1);
    eb = FracErr_Bkgd->GetBinContent(i + 1);
    es = FracErr_Sig->GetBinContent(i + 1);
    errscale = nexp / (eb * eb * nb * nb + es * es * ns * ns + nexp);
    if (nobs > 0) {
      chi2 += (2 * (nexp - nobs + nobs * TMath::Log(nobs / nexp)) * errscale);
    } else if (nobs == 0) {
      chi2 += (2 * nexp * errscale);
    }
    // chi2 is undefined if nexp is 0
  }

  return chi2;
}

double NueSterileFit::StandardChi2(TH1D *nexp) {
  double chi2 = 0;

  CalculateErrorMatrixInv(nexp);
  for (unsigned int i = 0; i < nBins; i++) {
    for (unsigned int j = 0; j < nBins; j++) {
      chi2 += ((NObs->GetBinContent(i + 1) - nexp->GetBinContent(i + 1)) *
               (NObs->GetBinContent(j + 1) - nexp->GetBinContent(j + 1)) *
               InvErrorMatrix->GetBinContent(i + 1, j + 1));
    }
  }

  return chi2;
}

void NueSterileFit::CalculateErrorMatrixInv(TH1D *nexp) {
  ErrorMatrix->Reset();
  InvErrorMatrix->Reset();

  if (ErrCalc != 0) { // if ErrorCalc object has been added, use it to calculate covariance matrix
    ErrCalc->CalculateSystErrorMatrix();
    ErrorMatrix->Add(ErrCalc->CovMatrix);
    ErrCalc->CalculateHOOError();
    ErrorMatrix->Add(ErrCalc->CovMatrix_Decomp);
  } else if (ExternalErrorMatrix != 0) { // if setting a systematic matrix externally, use it
    ErrorMatrix->Add(ExternalErrorMatrix);
  }

  // otherwise, it'll be statistics only
  double ele;
  for (unsigned int j = 0; j < nBins; j++) {
    ele = ErrorMatrix->GetBinContent(j + 1, j + 1);
    ele += nexp->GetBinContent(j + 1);
    ErrorMatrix->SetBinContent(j + 1, j + 1, ele);
  }

  // take CovMatrix and make array appropiate for TMatrix
  double *Varray = new double[nBins * nBins];
  for (unsigned int i = 0; i < nBins; i++) {
    for (unsigned int j = 0; j < nBins; j++) {
      unsigned int k = i * nBins + j;
      Varray[k] = ErrorMatrix->GetBinContent(i + 1, j + 1);
    }
  }

  // make TMatrix
  TMatrixD *Vmatrix = new TMatrixD(nBins, nBins, Varray);

  // invert
  TMatrixD Vinvmatrix = Vmatrix->Invert();

  // make array out of inverse
  Double_t *Vinvarray = Vinvmatrix.GetMatrixArray();

  // make Vinv out of array
  for (unsigned int i = 0; i < nBins; i++) {
    for (unsigned int j = 0; j < nBins; j++) {
      InvErrorMatrix->SetBinContent(i + 1, j + 1, Vinvarray[i * nBins + j]);
    }
  }

  delete[] Varray;

  return;
}

static void dlnlFunction(int &npar, double * /*gin*/, double &result, double *par, int /*iflag*/) {
  std::vector<double> params;
  NueSterileFit *nuesf = (NueSterileFit*)gMinuit->GetObjectFit();

  for (int i = 0; i < npar; i++) {
    params.push_back(par[i]);
  }
  result = nuesf->StdLikeComparison(params);
}

double NueSterileFit::StdLikeComparison(std::vector<double> npar) { // Actual likelihood comparison happens here.  This gets minimized.

  double dlnl = 0.0;
  for (unsigned int i = 0; i < nBins; i++) {
    // Signal and background predictions:
    double sig = Sig->GetBinContent(i + 1);
    double bkg = Bkgd->GetBinContent(i + 1);

    // Loop through systematics and add contributions
    for (unsigned int j = 0; j < FracErr_Bkgd_List.size(); j++) {
      // Get the nuisance parameter:
      double f = npar.at(j);

      // Total shift on background:
      double sigma_bkg = FracErr_Bkgd_List[j]->GetBinContent(i + 1);
      bkg += f * sigma_bkg * bkg;

      // Total shift on signal:
      double sigma_sig = FracErr_Sig_List[j]->GetBinContent(i + 1);
      sig += f * sigma_sig * sig;
    }

    // "Data" distribution:
    double nobs = NObs->GetBinContent(i + 1);

    // Expected distribution (with appropriate nuisance shifts)
    double nexp = sig + bkg;
    double temp = 0;

    // Likelihood comparison:
    if (nobs > 0 && nexp > 0) {
      // Regular (both distributions positive):
      temp = nexp - nobs + nobs * TMath::Log(nobs) - nobs * TMath::Log(nexp);
    } else if (nobs == 0 && nexp > 0) {
      // No data was seen in this bin:
      temp = nexp;
    } else if (nexp == 0 && nobs == 0) {
      // Nothing was seen, nothing expected:
      temp = 0;
    } else {
      // Something weird happened:
      return 1.0e10;
    }

    // Turn it into a chi2:
    dlnl += 2.0 * temp;
  }

  // Add penalty term for nuisance parameters:
  for (unsigned int j = 0; j < FracErr_Bkgd_List.size(); j++) {
    dlnl += npar.at(j) * npar.at(j);
  }

  // Return the log likelihood:
  return dlnl;
}

double NueSterileFit::DoStdMinParam() {
  // Size of systematics.  The +nBins will be for HOOHE.
  int sys_size = FracErr_Bkgd_List.size() + nBins;

  if (sys_size > minuit->fMaxpar) {
    std::cout << "WARNING: WRONG MINUIT SIZE!" << std::endl;
  }

  Double_t *vstrt = new Double_t[sys_size];
  Double_t *stp   = new Double_t[sys_size];
  Double_t *bmin  = new Double_t[sys_size];
  Double_t *bmax  = new Double_t[sys_size];
  Int_t ierflg = 0;

  // Pass on the object and set the static function:
  gMinuit = minuit;
  minuit->SetObjectFit(this);
  minuit->SetFCN(dlnlFunction);

  // Set the parameters
  for (int i = 0; i < sys_size; i++) {
    vstrt[i] = 0.0;
    stp[i] = 1.0;
    bmin[i] = 0.0;
    bmax[i] = 0.0;

    if (i >= (int) FracErr_Bkgd_List.size()) {
      vstrt[i] = 1.0;
    }
    minuit->mnparm(i, Form("f%i", i), vstrt[i], stp[i], bmin[i], bmax[i], ierflg);
    if (i >= (int) FracErr_Bkgd_List.size()) {
      minuit->FixParameter(i);
    }
  }

  // Resets function value and errors to UNDEFINED:
  minuit->mnrset(1);

  // 2*dlnl style error definition
  minuit->SetErrorDef(1.0);

  // Max iterations:
  minuit->SetMaxIterations(500);

  // Go minimize!
  minuit->Migrad();

  double minpoint = minuit->fAmin;
  return minpoint;
}

double NueSterileFit::StandardLikelihood() {
  // Call the likelihood equation either with or without systematics.
  double dlnl = 1.0E10;

  // empty vector of fit parameters
  std::vector<double> npar; npar.clear();

  unsigned int sys_size = FracErr_Bkgd_List.size();
  // If we are using systematics, do a minimization
  if (sys_size > 0) {
    dlnl = DoStdMinParam();

    // Check if the fit converged.
    Double_t fmin, fedm, errdef;
    Int_t npari, nparx, istat;
    minuit->mnstat(fmin, fedm, errdef, npari, nparx, istat);
    if (istat != 3 && fedm > 0.01) {
      cout << "Fit failed with status " << istat << endl;
      cout << "---> FMIN=" << fmin << " FEDM=" << fedm << " ERRDEF=" << errdef
           << " NPARI=" << npari << " NPARX=" << nparx << endl;
    }

  } else { // This will bypass the minimization and systematics.
    dlnl = StdLikeComparison(npar);
  }

  return dlnl;
}

void NueSterileFit::CalculateDlnLMatrix(TH2D *SysMatrix = 0, TH2D *HOOHEMatrix = 0) {
  // Calculate the covariance matrix for the "bin by bin" method

  // Make new inverted matrix if not done already
  if (InvErrorMatrix == 0) {
    if (SysMatrix) {
      InvErrorMatrix = (TH2D*)SysMatrix->Clone("InvErrorMatrix");
    } else if (HOOHEMatrix) {
      InvErrorMatrix = (TH2D*)HOOHEMatrix->Clone("InvErrorMatrix");
    } else {
      std::cout << "Didn't give me any matrices to invert!" << std::endl;
      return;
    }
  }

  // Reset everything to zero:
  InvErrorMatrix->Reset();
  unsigned int totbins = nBins;
  // Add together error elements into a single array
  double *Varray = new double[totbins * totbins];
  for (unsigned int i = 0; i < totbins; i++) {
    for (unsigned int j = 0; j < totbins; j++) {
      unsigned int k = i * totbins + j;
      Varray[k] = 0;
      if (SysMatrix != 0) {
        Varray[k] += SysMatrix->GetBinContent(i + 1, j + 1);
      }
      if (HOOHEMatrix != 0) {
        Varray[k] += HOOHEMatrix->GetBinContent(i + 1, j + 1);
      }
    }
  }

  // Hand elements to a TMatrix
  TMatrixD *Vmatrix = new TMatrixD(totbins, totbins, Varray);

  // Insert it (determ found for debugging purposes)!
  double determ;
  TMatrixD Vinvmatrix = Vmatrix->Invert(&determ);

  // Extract the array of the inverted matrix:
  Double_t *Vinvarray = Vinvmatrix.GetMatrixArray();

  // Turn it into the inverted covariance matrix
  // make Vinv out of array
  for (unsigned int i = 0; i < totbins; i++) {
    for (unsigned int j = 0; j < totbins; j++) {
      InvErrorMatrix->SetBinContent(i + 1, j + 1, Vinvarray[i * totbins + j]);
    }
  }

  delete[] Varray;
  return;
}

double NueSterileFit::BinLikeComparison(std::vector<double> npar) {
  double dlnl = 0.0;

  // Loop over fit bins:
  for (unsigned int i = 0; i < nBins; i++) {
    // Get Signal and background for prediction:
    double sig = Sig->GetBinContent(i + 1);
    double bkg = Bkgd->GetBinContent(i + 1);

    // Prediction:
    double nexp = sig + bkg;

    // Add nuisance parameter shift if using matrix:
    if (npar.size() > 0 && ErrCalc != 0 && InvErrorMatrix != 0) {
      double f = npar.at(i);
      nexp += f;
    }

    // Observed distribution:
    double nobs = NObs->GetBinContent(i + 1);
    double temp = 0;
    // Likelihood comparison:
    if (nobs > 0 && nexp > 0) {
      // Regular (both distributions positive):
      temp = nexp - nobs + nobs * TMath::Log(nobs) - nobs * TMath::Log(nexp);
    } else if (nobs == 0 && nexp > 0) {
      // No data was seen in this bin:
      temp = nexp;
    } else if (nexp == 0 && nobs == 0) {
      // Nothing was seen, nothing expected:
      temp = 0;
    } else {
      // Something weird happened:
      return 1.0E10;
    }

    dlnl += 2.0 * temp;
  }

  if (npar.size() > 0 && InvErrorMatrix != 0) {
    for (unsigned int i = 0; i < npar.size(); i++) {
      for (unsigned int j = 0; j < npar.size(); j++) {
        // Covariance terms
        dlnl += npar.at(i) * InvErrorMatrix->GetBinContent(i + 1, j + 1) * npar.at(j);
      }
    }
  }

  return dlnl;
}

static void binFunction(int &npar, double * /*gin*/, double &result, double *par, int /*iflag*/) {
  std::vector<double> params;
  NueSterileFit *nuesf = (NueSterileFit*)gMinuit->GetObjectFit();

  for (int i = 0; i < npar; i++) {
    params.push_back(par[i]);
  }
  result = nuesf->BinLikeComparison(params);
}

double NueSterileFit::DoBinMinParam() {
  int sys_size = nBins;
  if (sys_size > minuit->fMaxpar) {
    std::cout << "WARNING: WRONG MINUIT SIZE!" << std::endl;
  }
  Double_t *vstrt = new Double_t[sys_size];
  Double_t *stp   = new Double_t[sys_size];
  Double_t *bmin  = new Double_t[sys_size];
  Double_t *bmax  = new Double_t[sys_size];
  Int_t ierflg = 0;
  // Pass on the object and set the static function:
  gMinuit = minuit;
  minuit->SetObjectFit(this);
  minuit->SetFCN(binFunction);

  // Set the parameters
  for (int i = 0; i < sys_size; i++) {
    vstrt[i] = 0.0;
    stp[i]   = 1.0;
    bmin[i]  = 0.0;
    bmax[i]  = 0.0;
    minuit->mnparm(i, Form("f%i", i), vstrt[i], stp[i], bmin[i], bmax[i], ierflg);
  }

  // Resets function value and errors to UNDEFINED:
  minuit->mnrset(1);

  // 1 std dev for dlnl:
  minuit->SetErrorDef(1.0);

  // Max iterations:
  minuit->SetMaxIterations(500);
  // Go minimize!
  minuit->SetPrintLevel(-1);
  // cout << "Migrading..." << endl;
  minuit->Migrad();
  std::cout << "Done Migrading..." << std::endl;

  // Get the minimum for the function
  double minpoint = minuit->fAmin;
  return minpoint;
}

double NueSterileFit::BinLikelihood() {
  // Define the covariance matrix
  if (ErrCalc != 0) {
    std::cout << "Running CalculateSystErrorMatrix() ..." << std::endl;
    // Calculate the systematic error covariance matrix:
    ErrCalc->CalculateSystErrorMatrix();

    // Combine HOOHE and Syst into matrix, and invert
    CalculateDlnLMatrix(ErrCalc->CovMatrix, ErrCalc->CovMatrix_Decomp);
  }

  double dlnl = 1.0E10;
  std::vector<double> npar;
  if (ErrCalc == 0) {
    // Just use the empty parameter array if no error matrix set:
    dlnl = BinLikeComparison(npar);
  } else {
    // Otherwise, do a Minuit minimization
    dlnl = DoBinMinParam();

    Double_t fmin, fedm, errdef;
    Int_t npari, nparx, istat;
    minuit->mnstat(fmin, fedm, errdef, npari, nparx, istat);
    if (istat != 3 && fedm > 0.01) {
      cout << "Fit failed with status " << istat << endl;
      cout << "---> FMIN=" << fmin << " FEDM=" << fedm << " ERRDEF=" << errdef << " NPARI=" << npari << " NPARX=" << nparx << endl;
    }
  }

  return dlnl;
}

void NueSterileFit::Run2DSterileSlice() {
  // Given a value of Dm41, produces a 2D chisquared graph of both sin(th14)^2
  // vs sin(th24)^2 and |Ue4|^2 vs |Umu4|^2

  if (NObs == 0) {
    std::cout << "NObs not set.  Quitting..." << std::endl;
    return;
  }

  if (FitMethod == 1 && (FracErr_Bkgd == 0 || FracErr_Sig == 0)) {
    std::cout << "FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2. Quitting..." << std::endl;
    return;
  }

  if (Extrap.size() == 0) {
    std::cout << "No Extrapolate2D input. Quitting..." << std::endl;
    return;
  }

  if (FitMethod == 3) {
    DefineStdDlnLMinuit();
  }

  if (FitMethod == 4) {
    DefineBinDlnLMinuit();
  }

  Bkgd = (TH1D *)NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D *)NObs->Clone("Sig");
  Sig->Reset();
  NExp = (TH1D *)NObs->Clone("NExp");
  NExp->Reset();

  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
  }

  if (ErrorMatrix == 0) {
    ErrorMatrix = new TH2D("ErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }
  if (InvErrorMatrix == 0) {
    InvErrorMatrix = new TH2D("InvErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }

  int arraysize = (nSinSqTh14Steps + 1) * (nSinSqTh24Steps + 1);
  std::cout << "Array size : " <<  arraysize << std::endl;

  double Ue42[10201];
  double Um42[10201];
  double Ssqth24[10201];
  double chi[10201];
  double dchi[10201];

  Double_t Th14increment = 0;
  if (nSinSqTh14Steps > 0) {
    Th14increment = (SinSqTh14High - SinSqTh14Low) / (nSinSqTh14Steps);
  }

  Double_t Th24increment = 0;
  if (nSinSqTh24Steps > 0) {
    Th24increment = (SinSqTh24High - SinSqTh24Low) / (nSinSqTh24Steps);
  }

  Int_t idx = 0; // Array indexer
  for (int i = 0; i < nSinSqTh14Steps + 1; i++) {
    Double_t ssth14 = i * Th14increment + SinSqTh14Low;
    Double_t Theta14 = TMath::ASin(TMath::Sqrt(ssth14));
    Double_t csth14 = TMath::Cos(Theta14) * TMath::Cos(Theta14);

    for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
      Extrap[ie]->SetSinSqTh14(ssth14);
    }

    for (int j = 0; j < nSinSqTh24Steps + 1; j++) {
      Double_t ssth24 = j * Th24increment + SinSqTh24Low;

      for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
        Extrap[ie]->SetSinSqTh24(ssth24);
        Extrap[ie]->OscillatePrediction();
      }

      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();

      for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
        Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
        Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
      }
      NExp->Add(Bkgd);
      NExp->Add(Sig);

      if (FitMethod == 0) {
        Ue42[idx] = ssth14;
        Um42[idx] = csth14 * ssth24;
        Ssqth24[idx] = ssth24;
        chi[idx] = PoissonChi2(NExp);
      } else if (FitMethod == 1) {
        Ue42[idx] = ssth14;
        Um42[idx] = csth14 * ssth24;
        Ssqth24[idx] = ssth24;
        chi[idx] = ScaledChi2(Bkgd, Sig);
      } else if (FitMethod == 2) {
        Ue42[idx] = ssth14;
        Um42[idx] = csth14 * ssth24;
        Ssqth24[idx] = ssth24;
        chi[idx] = StandardChi2(NExp);
      } else if (FitMethod == 3) {
        // Likelihood: "Standard" (N syst, N nuisance)
        // Calculate the likelihood (x2 for chi)
        Ue42[idx] = ssth14;
        Um42[idx] = csth14 * ssth24;
        Ssqth24[idx] = ssth24;
        chi[idx] = StandardLikelihood();
      } else if (FitMethod == 4) {
        // Likelihood: Bin by Bin Calculation of Systematics
        // Calculate the likelihood (x2 for chi)
        Ue42[idx] = ssth14;
        Um42[idx] = csth14 * ssth24;
        Ssqth24[idx] = ssth24;
        chi[idx] = BinLikelihood();
      } else {
        std::cout << "Error in Run2DSterileSlice(): Unknown 'FitMethod'." << std::endl;
      }
      idx++;
    } // TH24 loop
  } // TH14 loop

  double Ue42best[1];
  double Um42best[1];
  double Ssqth24best[1];

  double mc2 = 1E9;
  for (int k = 0; k < arraysize; k++) {
    if (chi[k] < mc2) {
      mc2 = chi[k];
      Ue42best[0] = Ue42[k];
      Um42best[0] = Um42[k];
      Ssqth24best[0] = Ssqth24[k];
    }
  }

  for (int s = 0; s < arraysize; s++) {
    dchi[s] = chi[s] - mc2;
    if (abs(dchi[s]) < 0.0001) {
      std::cout << "I have a zero!" << std::endl;
    }
  }

  // make a TGraph
  TFile *w = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "RECREATE");
  Int_t d = 1;

  TGraph *bpss = new TGraph(d, Ue42best, Ssqth24best);
  bpss->SetName("bpss");
  bpss->Write();

  TGraph *bpUU = new TGraph(d, Ue42best, Um42best);
  bpUU->SetName("bpUU");
  bpUU->Write();

  TGraph2D *ss = new TGraph2D(arraysize, Ue42, Ssqth24, chi);
  ss->SetName("ss");
  ss->GetXaxis()->SetTitle("sin^{2}#theta_{14}");
  ss->GetYaxis()->SetTitle("sin^{2}#theta_{24}");
  ss->GetZaxis()->SetTitle("-2lnL");
  ss->Write();

  TGraph2D *dss = new TGraph2D(arraysize, Ue42, Ssqth24, dchi);
  dss->SetName("dss");
  dss->GetXaxis()->SetTitle("sin^{2}#theta_{14}");
  dss->GetYaxis()->SetTitle("sin^{2}#theta_{24}");
  dss->GetZaxis()->SetTitle("-2#DeltalnL");
  dss->Write();

  TGraph2D *UU = new TGraph2D(arraysize, Ue42, Um42, chi);
  UU->SetName("UU");
  UU->GetXaxis()->SetTitle("sin^{2}#theta_{14}");
  UU->GetYaxis()->SetTitle("cos^{2}#theta_{14}sin^{2}#theta_{24}");
  UU->GetZaxis()->SetTitle("-2lnL");
  UU->Write();

  TGraph2D *dUU = new TGraph2D(arraysize, Ue42, Um42, dchi);
  dUU->SetName("dUU");
  dUU->GetXaxis()->SetTitle("sin^{2}#theta_{14}");
  dUU->GetYaxis()->SetTitle("cos^{2}#theta_{14}sin^{2}#theta_{24}");
  dUU->GetZaxis()->SetTitle("-2#DeltalnL");
  dUU->Write();

  w->Close();
  return;

} // 2D SterileSlice

bool NueSterileFit::VariateParameter(std::vector<OscPar::EOscPar> ParameterEnumToVariates, std::vector<std::pair<double, double> > ParameterValueToVariates, std::vector<unsigned int> numberOfGridPoint, int idx) {
  if (numberOfGridPoint.size() != ParameterEnumToVariates.size() || (numberOfGridPoint.size() != ParameterValueToVariates.size())) {
    std::cout << "Mismatch between number of parameters and size of grid vector!" << std::endl;
    return false;
  }

  int TotalNumberOfGridPoints = 1;
  for (unsigned int i = 0; i < numberOfGridPoint.size(); i++) {
    TotalNumberOfGridPoints *= numberOfGridPoint.at(i);
  }

  if ((idx + 1 > TotalNumberOfGridPoints) || (idx < 0))  {
    return false;
  }

  std::vector<unsigned int> idx_component;
  idx_component.reserve(numberOfGridPoint.size());
  for (unsigned int i = 0; i < numberOfGridPoint.size(); i++) {
    int comp_index = idx % numberOfGridPoint.at(i);
    idx = (idx - comp_index) / numberOfGridPoint.at(i);
    idx_component.push_back(comp_index);
  }

  for (unsigned int i = 0; i < ParameterValueToVariates.size(); i++) {
    double param_high = ParameterValueToVariates.at(i).first + ParameterValueToVariates.at(i).second;
    double param_low  = ParameterValueToVariates.at(i).first - ParameterValueToVariates.at(i).second;
    double step_param = (param_high - param_low) / (numberOfGridPoint.at(i) - 1);
    double value = param_low + idx_component.at(i) * step_param;
    for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
      Extrap[ie]->SetOscPar(ParameterEnumToVariates.at(i), value);
      Extrap[ie]->OscillatePrediction();
    }
  }

  return true;
}

void NueSterileFit::Run3FlavorDeltaCPFit() {
  if (NObs == 0) {
    std::cout << "NObs not set.  Quitting..." << std::endl;
    return;
  }

  if (FitMethod == 1 && (FracErr_Bkgd == 0 || FracErr_Sig == 0)) {
    std::cout << "FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2. Quitting..." << std::endl;
    return;
  }

  if (Extrap.size() == 0) {
    std::cout << "No Extrapolate2D input. Quitting..." << std::endl;
    return;
  }

  if (FitMethod == 3) {
    DefineStdDlnLMinuit();
  }

  if (FitMethod == 4) {
    DefineBinDlnLMinuit();
  }

  Bkgd = (TH1D *)NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D *)NObs->Clone("Sig");
  Sig->Reset();
  NExp = (TH1D *)NObs->Clone("NExp");
  NExp->Reset();

  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
  }

  if (ErrorMatrix == 0) {
    ErrorMatrix = new TH2D("ErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }
  if (InvErrorMatrix == 0) {
    InvErrorMatrix = new TH2D("InvErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }

  int arraysize = nDeltaCP13Steps + 1;
  std::cout << "Array size : " <<  arraysize << std::endl;
  double* dcp_array  = (double*) malloc(sizeof(double) * arraysize);
  double* chi        = (double*) malloc(sizeof(double) * arraysize);
  double* dchi       = (double*) malloc(sizeof(double) * arraysize);

  double deltacp13Increment = 0.;
  if (nDeltaCP13Steps > 0) {
    deltacp13Increment = (DeltaCP13High - DeltaCP13Low) / (nDeltaCP13Steps);
  }

  for (int i = 0; i < nDeltaCP13Steps + 1; i++) {
    double deltacp13 = ((double) i) * deltacp13Increment + DeltaCP13Low;
    *(dcp_array + i) = deltacp13;

    for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
      Extrap[ie]->SetOscPar(OscPar::kDelta, deltacp13);
      Extrap[ie]->OscillatePrediction();
    }

    // Need a loop to variate the dm32 and theta23 and theta13
    int idx = 0;
    std::vector<unsigned int> numberOfGridPoint; numberOfGridPoint.clear();
    numberOfGridPoint.push_back(5);
    numberOfGridPoint.push_back(5);
    numberOfGridPoint.push_back(5);
    std::vector<OscPar::EOscPar> ParameterEnumToVariates; ParameterEnumToVariates.clear();
    ParameterEnumToVariates.push_back(OscPar::kTh13);
    ParameterEnumToVariates.push_back(OscPar::kTh23);
    ParameterEnumToVariates.push_back(OscPar::kDeltaM23);
    std::vector<std::pair<double, double> > ParameterValueToVariates; ParameterValueToVariates.clear();
    ParameterValueToVariates.push_back(std::make_pair(0.1482, 0.0110));  //theta_13
    ParameterValueToVariates.push_back(std::make_pair(0.8214, 0.06));    //theta_23
    ParameterValueToVariates.push_back(std::make_pair(2.55E-3, 0.10E-3)); //dmsq_32
    double profiling_minchi2 = 9.0E9;
    while (VariateParameter(ParameterEnumToVariates, ParameterValueToVariates, numberOfGridPoint, idx)) {
      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();

      for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
        Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
        Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
      }
      NExp->Add(Bkgd);
      NExp->Add(Sig);

      double temp_chi = 0.;
      if (FitMethod == 0) {
        temp_chi = PoissonChi2(NExp);
      } else if (FitMethod == 1) {
        temp_chi = ScaledChi2(Bkgd, Sig);
      } else if (FitMethod == 2) {
        temp_chi = StandardChi2(NExp);
      } else if (FitMethod == 3) {
        // Likelihood: "Standard" (N syst, N nuisance)
        // Calculate the likelihood (x2 for chi)
        temp_chi = StandardLikelihood();
      } else if (FitMethod == 4) {
        // Likelihood: Bin by Bin Calculation of Systematics
        // Calculate the likelihood (x2 for chi)
        temp_chi = BinLikelihood();
      } else {
        std::cout << "Error in Run2DSterileSlice(): Unknown 'FitMethod'." << std::endl;
      }
      profiling_minchi2 = profiling_minchi2 > temp_chi ? temp_chi : profiling_minchi2;
    }
    *(chi + i) = profiling_minchi2;
    idx++;
  }

  double mc2 = 1.0E9;
  for (int k = 0; k < arraysize; k++) {
    if (chi[k] < mc2) {
      mc2 = chi[k];
    }
  }

  for (int s = 0; s < arraysize; s++) {
    dchi[s] = chi[s] - mc2;
    if (abs(dchi[s]) < 0.0001) {
      std::cout << "I have a zero!" << std::endl;
    }
  }

  // make a TGraph
  TFile *w = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "RECREATE");

  TGraph *bchi = new TGraph(nDeltaCP13Steps + 1, dcp_array, chi);
  bchi->SetName("bchi");
  bchi->Write();

  TGraph *bdchi = new TGraph(nDeltaCP13Steps + 1, dcp_array, dchi);
  bdchi->SetName("bdchi");
  bdchi->Write();

  w->Close();
  return;
} // Run3FlavorDeltaCPFit

void NueSterileFit::SetGridFiles(std::string snorm="Grid_Ana_Norm.root", std::string sinvt="Grid_Ana_Invt.root") {
  GridFileName_Normal = snorm;
  GridFileName_Inverted = sinvt;
  return;
}

void NueSterileFit::ReadGridFiles() {

  double farpot;
  for (unsigned int i = 0; i < nBins; i++) {
    grid_bin_oscparerr.push_back(0);
  }

  // Processing Normal Hierarchy file
  GridTree_Normal.clear();
  GridTree_2_Normal.clear();
  nPts_Normal = 0;
  if (gSystem->AccessPathName(gSystem->ExpandPathName(GridFileName_Normal.c_str()))) {
    std::cout << "Grid file (normal hierarchy) doesn't exist." << std::endl;
    return;
  } else {
    TFile* fnorm = new TFile(gSystem->ExpandPathName(GridFileName_Normal.c_str()), "READ");
    for (unsigned int i = 0; i < nBins; i++) {
      TTree* temp = (TTree*)fnorm->Get(Form("Bin_%i", i));
      if (i == 0) nPts_Normal = temp->GetEntries();
      temp->SetName(Form("Bin_%i_Normal", i));
      temp->SetBranchAddress("Background", &grid_background);
      temp->SetBranchAddress("Signal",     &grid_signal);
      temp->SetBranchAddress("SinSqTh14",  &grid_sinsqth14);
      temp->SetBranchAddress("SinSqTh24",  &grid_sinsqth24);
      temp->SetBranchAddress("Dmsq41",     &grid_dmsq41);
      temp->SetBranchAddress("Delta13",    &grid_delta13);
      temp->SetBranchAddress("Delta14",    &grid_delta14);
      temp->SetBranchAddress("Delta24",    &grid_delta24);
      temp->SetBranchAddress("Theta34",    &grid_th34);
      if (IncludeOscParErrs) {
        temp->SetBranchAddress("DNExp_DOscPars", &grid_bin_oscparerr[i]);
      }
      GridTree_Normal.push_back(temp);
      GridTree_2_Normal.push_back(vector<TTree*>());
      for (unsigned int j = 0; j < Extrap.size(); j++) {
        TTree* temp2 = (TTree*) fnorm->Get(Form("Bin_%i_Run_%i", i, j));
        temp2->SetName(Form("Bin_%i_Run_%i_Normal", i, j));
        temp2->SetBranchAddress("NC",         &grid_nc);
        temp2->SetBranchAddress("NuMuCC",     &grid_numucc);
        temp2->SetBranchAddress("BNueCC",     &grid_bnuecc);
        temp2->SetBranchAddress("NuTauCC",    &grid_nutaucc);
        temp2->SetBranchAddress("Signal",     &grid_nue);
        temp2->SetBranchAddress("SinSqTh14",  &grid_sinsqth14);
        temp2->SetBranchAddress("SinSqTh24",  &grid_sinsqth24);
        temp2->SetBranchAddress("Dmsq41",     &grid_dmsq41);
        temp2->SetBranchAddress("Delta13",    &grid_delta13);
        temp2->SetBranchAddress("Delta14",    &grid_delta14);
        temp2->SetBranchAddress("Delta24",    &grid_delta24);
        temp2->SetBranchAddress("Theta34",    &grid_th34);
        GridTree_2_Normal.at(i).push_back(temp2);
      }
    }
    std::cout << "nPts_Normal = " << nPts_Normal << std::endl;

    paramtree_Normal = (TTree *)fnorm->Get("paramtree");
    paramtree_Normal->SetName("paramtree_Normal");
    paramtree_Normal->SetBranchAddress("farPOT",     &farpot);
    paramtree_Normal->SetBranchAddress("Theta12",    &grid_normal_th12);
    paramtree_Normal->SetBranchAddress("Theta13",    &grid_normal_th13);
    paramtree_Normal->SetBranchAddress("Theta23",    &grid_normal_th23);
    paramtree_Normal->SetBranchAddress("DeltaMSq23", &grid_normal_dm2_32);
    paramtree_Normal->SetBranchAddress("DeltaMSq12", &grid_normal_dm2_21);
    paramtree_Normal->SetBranchAddress("DeltaMSq41", &grid_dmsq41);
    paramtree_Normal->SetBranchAddress("SinSqTh14",  &grid_sinsqth14);
    paramtree_Normal->SetBranchAddress("SinSqTh24",  &grid_sinsqth24);
    paramtree_Normal->GetEntry(0);

    if (GridNorm > 0) {
      GridScale_Normal = GridNorm / farpot;
    }
  }

  return;
}

void NueSterileFit::PrintMatrix(TMatrixD& matrix) {
  unsigned int nRows = matrix.GetNrows();
  unsigned int nCols = matrix.GetNcols();

  for (unsigned int i = 0; i < nRows; i++) {
    for (unsigned int j = 0; j < nCols; j++) {
      std::cout << std::setw(16) << Form("%8.5f ", matrix[i][j]);
    }
    std::cout << std::endl;
  }
}

std::vector<TMatrixD> NueSterileFit::GenerateDiagonalizedCovarianceMatrix(TH1D *nexp, TH2D *err) {
  nBins = nexp->GetNbinsX();

  // error matrix
  double *Varray = new double[nBins * nBins];
  for (unsigned int i = 0; i < nBins; i++) {
    for (unsigned int j = 0; j < nBins; j++) {
      int k = i * nBins + j;
      Varray[k] = err->GetBinContent(i + 1, j + 1);
    }
  }
  const TMatrixD NonDiagonalCovarianceMatrix(nBins, nBins, Varray);

  // get eigenvectors of M
  TMatrixDEigen CovarianceEigenMatrix(NonDiagonalCovarianceMatrix);
  // construct matrix of eigenvectors
  TMatrixD TransposeBasisTransformationMatrix = CovarianceEigenMatrix.GetEigenVectors();


  // transpose to get transformation matrix A
  TMatrixD BasisTransformationMatrix = TransposeBasisTransformationMatrix;
  BasisTransformationMatrix.Transpose(BasisTransformationMatrix);

  // diagonalize M
  TMatrixD DiagonalCovarianceMatrix(nBins, nBins);
  DiagonalCovarianceMatrix = BasisTransformationMatrix * NonDiagonalCovarianceMatrix * TransposeBasisTransformationMatrix;

  std::vector<TMatrixD> Matrices;
  Matrices.push_back(DiagonalCovarianceMatrix);
  Matrices.push_back(BasisTransformationMatrix);
  Matrices.push_back(TransposeBasisTransformationMatrix);

  return Matrices;
}

void NueSterileFit::GenerateOneCorrelatedExp(TH1D *nexp, std::vector<TMatrixD>& Matrices) {
  NObs = (TH1D*) nexp->Clone("NObs");
  NObs->Reset();
  nBins = NObs->GetNbinsX();

  // transform prediction to diagonal basis using A
  vector<double> nnew;
  for (unsigned int i = 0; i < nBins; i++) {
    double t = 0;
    for (unsigned int j = 0; j < nBins; j++) {
      t += (Matrices.at(1))[i][j] * nexp->GetBinContent(j + 1);
    }
    // std::cout << "Bin #" << i << ": " << t << std::endl;
    nnew.push_back(t);
  }

  // randomize prediction in diagonal basis
  vector<double> nrand;
  for (unsigned int i = 0; i < nBins; i++) {
    if ((Matrices.at(0))[i][i] < 0) {
      std::cout << "Warning in NueFit2D::GenerateOneCorrelatedExp(): Negative element in diagonalized systematic error  matrix ... setting to 0." << std::endl;
      (Matrices.at(0))[i][i] = 0;
    }
    nrand.push_back(gRandom->Gaus(nnew[i], sqrt((Matrices.at(0))[i][i])));
    // std::cout << "NRAND Bin #" << i << ": " << nrand[i] << std::endl;
  }

  // Transform randomized prediction back to original basis
  // Then perform digitization with Poisson
  for (unsigned int i = 0; i < nBins; i++) {
    double t = 0;
    for (unsigned int j = 0; j < nBins; j++) {
      t += (Matrices.at(2))[i][j] * nrand[j];
    }
    // std::cout << "AT.NRAND Bin #" << i << ": " << t << std::endl;
    t = gRandom->Poisson(t);
    NObs->SetBinContent(i + 1, t);
  }

  return;
}


double NueSterileFit::FindBestFitForPseudoExperiment() {
  if (ErrCalc != 0) {
    std::cout << "Running CalculateSystErrorMatrix() ..." << std::endl;
    // Calculate the systematic error covariance matrix:
    ErrCalc->CalculateSystErrorMatrix();

    // Combine HOOHE and Syst into matrix, and invert
    CalculateDlnLMatrix(ErrCalc->CovMatrix, ErrCalc->CovMatrix_Decomp);
  }

  double dlnl = 1.0E10;
  std::vector<double> npar;

  if (ErrCalc == 0) {
    std::cout << "Error matrix is not set. Quitting..." << std::endl;
    return 0;
  }

  Sig->Reset();
  Bkgd->Reset();
  unsigned int nPoints = 10;
  double max_theta14  = +TMath::Pi() / 4.;
  double min_theta14  = -TMath::Pi() / 4.;
  double step_theta14 = (max_theta14 - min_theta14) / (nPoints - 1);
  double max_theta24  = +TMath::Pi() / 4.;
  double min_theta24  = -TMath::Pi() / 4.;
  double step_theta24 = (max_theta24 - min_theta24) / (nPoints - 1);
  double chisq_val_min = 1.0E9;
  if (!FixTheta14Theta24SterileNeutrinoFitFeldmanCousins) {
    for (unsigned i14 = 0; i14 < nPoints; i14++) {
      double current_theta14 = min_theta14 + i14 * step_theta14;
      for (unsigned int i24 = 0; i24 < nPoints; i24++) {
        double current_theta24 = min_theta24 + i24 * step_theta24;
        Sig->Reset();
        Bkgd->Reset();
        for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
          Extrap[ie]->SetOscPar(OscPar::kTh14, current_theta14);
          Extrap[ie]->SetOscPar(OscPar::kTh24, current_theta24);
          Extrap[ie]->OscillatePrediction();
          Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
          Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
        }
        double chisq_val = DoBinMinParam();
        chisq_val_min = chisq_val < chisq_val_min ? chisq_val : chisq_val_min;
      }
    }
  } else {
    Sig->Reset();
    Bkgd->Reset();
    double grid_theta14 = TMath::ASin(TMath::Sqrt(grid_sinsqth14));
    double grid_theta24 = TMath::ASin(TMath::Sqrt(grid_sinsqth24));
    for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
      Extrap[ie]->SetOscPar(OscPar::kTh14, grid_theta14);
      Extrap[ie]->SetOscPar(OscPar::kTh24, grid_theta24);
      Extrap[ie]->OscillatePrediction();
      Bkgd->Add(Extrap[ie]->Pred_TotalBkgd_VsBinNumber);
      Sig->Add(Extrap[ie]->Pred_Signal_VsBinNumber);
    }
    chisq_val_min = DoBinMinParam();
  }

  dlnl = chisq_val_min;
  Double_t fmin, fedm, errdef;
  Int_t npari, nparx, istat;
  minuit->mnstat(fmin, fedm, errdef, npari, nparx, istat);
  if (istat != 3 && fedm > 0.01) {
    std::cout << "Fit failed with status " << istat << std::endl;
    std::cout << "---> FMIN=" << fmin << " FEDM=" << fedm << " ERRDEF=" << errdef << " NPARI=" << npari << " NPARX=" << nparx << std::endl;
  }

  return dlnl;
}

void NueSterileFit::RunMultiBinPseudoExpts(bool Print) {

  if (FitMethod == 1 && (FracErr_Bkgd == 0 || FracErr_Sig == 0)) {
    std::cout << "FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2. Quitting..." << std::endl;
    return;
  }

  if (Extrap.size() == 0) {
    std::cout << "No Extrapolate2D input.  Quitting..." << std::endl;
    return;
  }

  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    Extrap[ie]->GetPrediction();
  }

  if (ErrCalc == 0) {
    std::cout << "Need to set ErrorCalc object!  Quitting..." << std::endl;
    return;
  }

  if (FitMethod == 3) {
    DefineStdDlnLMinuit();
  } else if (FitMethod == 4) {
    DefineBinDlnLMinuit();
  }

  nBins = Extrap[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  if (ErrorMatrix == 0) {
    ErrorMatrix = new TH2D("ErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }
  if (InvErrorMatrix == 0) {
    InvErrorMatrix = new TH2D("InvErrorMatrix", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);
  }

  TH2D *Error4Expts = new TH2D("Error4Expts", "", nBins, -0.5, nBins - 0.5, nBins, -0.5, nBins - 0.5);

  ReadGridFiles();

  if (nPts_Normal == 0 || nPts_Inverted == 0) {
    return;
  }

  gRandom->SetSeed(time(NULL));

  TH1D *nexp_bkgd = new TH1D("nexp_bkgd", "", nBins, -0.5, nBins - 0.5);
  TH1D *nexp_signal = new TH1D("nexp_signal", "", nBins, -0.5, nBins - 0.5);
  TH1D *nexp = new TH1D("nexp", "", nBins, -0.5, nBins - 0.5);
  double delchi2, chi2min;
  TH1D *dchi2hist = new TH1D("dchi2hist", "", 110000, -10, 100);
  int noff;

  std::vector<std::vector<double> > nc, numucc, bnuecc, nutaucc, sig;
  for (unsigned int j = 0; j < nBins; j++) {
    nc.push_back(std::vector<double>());
    numucc.push_back(std::vector<double>());
    bnuecc.push_back(std::vector<double>());
    nutaucc.push_back(std::vector<double>());
    sig.push_back(std::vector<double>());
    for (unsigned int k = 0; k < Extrap.size(); k++) {
      nc.at(j).push_back(0);
      numucc.at(j).push_back(0);
      bnuecc.at(j).push_back(0);
      nutaucc.at(j).push_back(0);
      sig.at(j).push_back(0);
    }
  }

  Bkgd = (TH1D*) nexp->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D*) nexp->Clone("Sig");
  Sig->Reset();

  ofstream myfile;
  string file, ofile;

  // normal hierarchy

  TFile *rootfile = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "RECREATE");

  if (Print) {
    ofile = gSystem->ExpandPathName(outFileName.c_str());
    file = ofile.substr(0, ofile.length() - 5) + "_Normal.dat";
    myfile.open(gSystem->ExpandPathName(file.c_str()));
  }

  for (int i = 0; i < nPts_Normal; i++) {
  // for (int i = 0; i < 1; i++) {
    std::cout << "point " << (i + 1) << "/" << nPts_Normal << " (normal hierarchy)" << std::endl;

    nexp_bkgd->Reset();
    nexp_signal->Reset();
    nexp->Reset();

    for (unsigned int j = 0; j < nBins; j++) {
      GridTree_Normal.at(j)->GetEntry(i);
      nexp_bkgd->SetBinContent(j + 1, grid_background);
      nexp_signal->SetBinContent(j + 1, grid_signal);
      nexp->SetBinContent(j + 1, grid_background + grid_signal);

      for (unsigned int k = 0; k < Extrap.size(); k++) {
        GridTree_2_Normal.at(j).at(k)->GetEntry(i);
        nc.at(j).at(k) = grid_nc;
        numucc.at(j).at(k) = grid_numucc;
        bnuecc.at(j).at(k) = grid_bnuecc;
        nutaucc.at(j).at(k) = grid_nutaucc;
        sig.at(j).at(k) = grid_nue;
      }
    }


    // for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    //   Extrap.at(ie)->OscillatePrediction();
    //   for (unsigned int i = 0; i < nBins; i++) {
    //     int ir = (int) (i / Extrap[ie]->GetNPID());
    //     int ip = i % Extrap[ie]->GetNPID();

    //     double tmp_bkgd       = Extrap.at(ie)->Pred_TotalBkgd_VsBinNumber->GetBinContent(i + 1);
    //     double tmp_sig        = Extrap.at(ie)->Pred_Signal_VsBinNumber->GetBinContent(i + 1);
    //     nc.at(i).at(ie)       = Extrap.at(ie)->Pred[Background::kNC]      ->GetBinContent(ip + 1, ir + 1);
    //     numucc.at(i).at(ie)   = Extrap.at(ie)->Pred[Background::kNuMuCC]  ->GetBinContent(ip + 1, ir + 1);
    //     bnuecc.at(i).at(ie)   = Extrap.at(ie)->Pred[Background::kBNueCC]  ->GetBinContent(ip + 1, ir + 1);
    //     nutaucc.at(i).at(ie)  = Extrap.at(ie)->Pred[Background::kNuTauCC] ->GetBinContent(ip + 1, ir + 1);
    //     sig.at(i).at(ie)      = Extrap.at(ie)->Pred[Background::kNueCC]   ->GetBinContent(ip + 1, ir + 1);
    //     nexp->SetBinContent(i + 1, tmp_bkgd + tmp_sig);
    //     nexp_bkgd->SetBinContent(i + 1, tmp_bkgd);
    //     nexp_signal->SetBinContent(i + 1, tmp_sig);
    //   }
    // }

    ErrCalc->SetGridPred(nBins, nc, numucc, bnuecc, nutaucc, sig);
    Error4Expts->Reset();
    ErrCalc->SetUseGrid(true);
    ErrCalc->CalculateSystErrorMatrix();
    Error4Expts->Add(ErrCalc->CovMatrix);
    if (IncludeOscParErrs) {
      noff = 0;
      for (unsigned int j = 0; j < nBins; j++) {
        double ele = Error4Expts->GetBinContent(j + 1, j + 1);
        ele += (grid_bin_oscparerr[j] * grid_bin_oscparerr[j] * nexp->GetBinContent(j + 1) * nexp->GetBinContent(j + 1));
        Error4Expts->SetBinContent(j + 1, j + 1, ele);
        for (unsigned int k = j + 1; k < nBins; k++) {
          ele = Error4Expts->GetBinContent(j + 1, k + 1);
          ele += (grid_bin_oscparerr[j] * grid_bin_oscparerr[k] * nexp->GetBinContent(j + 1) * nexp->GetBinContent(k + 1));
          Error4Expts->SetBinContent(j + 1, k + 1, ele);
          ele = Error4Expts->GetBinContent(k + 1, j + 1);
          ele += (grid_bin_oscparerr[j] * grid_bin_oscparerr[k] * nexp->GetBinContent(j + 1) * nexp->GetBinContent(k + 1));
          Error4Expts->SetBinContent(k + 1, j + 1, ele);
          noff++;
        }
      }
    }

    // Printing the non-diagonal error covariance matrix
    std::cout << "Covariance matrix from ErrorCalc: " << std::endl;
    for (unsigned int ix = 0; ix < nBins; ix++) {
      for (unsigned int iy = 0; iy < nBins; iy++) {
        std::cout << std::setw(15) << Form("%8.5f", Error4Expts->GetBinContent(ix + 1, iy + 1));
      }
      std::cout << std::endl;
    }

    dchi2hist->Reset();
    dchi2hist->SetName(Form("DeltaChi2Hist_Normal_%i", i));
    std::vector<TMatrixD> ErrorMatrices = GenerateDiagonalizedCovarianceMatrix(nexp, Error4Expts);
    std::cout << "Error Matrix: " << std::endl;
    PrintMatrix(ErrorMatrices.at(0));

    MultiVarGauss* mvg = new MultiVarGauss(nexp, Error4Expts);

    for (unsigned int u = 0; u < NumExpts; u++) {
      std::cout << "expt " << (u + 1) << "/" << NumExpts << std::endl;

      // GenerateOneCorrelatedExp(nexp, ErrorMatrices); // Old method by Lisa Whitehead
      // MultiVarGauss method by Adam Aurisano
      NObs = (TH1D*) nexp->Clone("NObs");
      NObs->Reset();
      NObs = (TH1D*) mvg->Generate();

      FixTheta14Theta24SterileNeutrinoFitFeldmanCousins = false;
      std::cout << "Varying Theta14 and Theta24 to find the best fit point." << std::endl;
      double chisq_bf = FindBestFitForPseudoExperiment();
      FixTheta14Theta24SterileNeutrinoFitFeldmanCousins = true;
      std::cout << "Fixing Theta14 and Theta24." << std::endl;
      double chisq    = FindBestFitForPseudoExperiment();
      double dchisq   = chisq - chisq_bf;
      if (dchisq < 0) {
        u--;
      } else {
        dchi2hist->Fill(dchisq);
      }

      if (Print) {
        for (unsigned int j = 0; j < nBins; j++) {
          myfile << std::setw(8) << Form("%4i", (int) NObs->GetBinContent(j + 1)) << " ";
        }
        myfile << std::setw(8) << Form("%4.2f", chisq_bf) << " " << Form("%4.2f", chisq) << std::endl;
      }

    } // End pseudo-experiment loop
    for (unsigned int j = 0; j < nBins; j++) {
      // myfile << std::setw(8) << Form("%4i: %4.2f", (int) nexp->GetBinContent(j + 1), (ErrorMatrices.at(0))[j][j]) << " ";
      myfile << std::setw(8) << Form("%4i", (int) nexp->GetBinContent(j + 1)) << " ";
    }
    myfile << std::endl;
    rootfile->cd();
    dchi2hist->Write();
  } // End nPtsNormal loop

  if (Print) {
    myfile.close();
  }

  rootfile->Close();

  return;
}