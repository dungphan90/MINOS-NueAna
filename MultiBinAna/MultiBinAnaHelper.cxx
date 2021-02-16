#define MultiBinAnaHelper_C

#include "NueAna/MultiBinAna/MultiBinAnaHelper.h"

MultiBinAnaHelper::MultiBinAnaHelper()
{
}
MultiBinAnaHelper::~MultiBinAnaHelper()
{
}

void MultiBinAnaHelper::Rebin2DHist(TH2D *h,Int_t nx,Double_t *x,Int_t ny,Double_t *y)
{
  //Be careful using this function!! It only works correctly if you give it correct input.  It is not smart enough to make sure that your new bins can really be constructed out of the old ones.
  
  Bool_t rebinx=true,rebiny=true;
  if(!x || nx==0)
  {
    rebinx = false;
  }
  if(!y || ny==0)
  {
    rebiny = false;
  }
  
  if(!rebinx && !rebiny)
  {
    //do nothing
    return;
  }
  
  Int_t i,ix,iy,jx,jy;
  
  Int_t noldx = h->GetNbinsX();
  Int_t noldy = h->GetNbinsY();
  double *oldxb = new double[noldx+1];
  double *oldyb = new double[noldy+1];
  for(ix=0;ix<noldx+1;ix++)
  {
    oldxb[ix] = h->GetXaxis()->GetBinLowEdge(ix+1);
  }
  for(iy=0;iy<noldy+1;iy++)
  {
    oldyb[iy] = h->GetYaxis()->GetBinLowEdge(iy+1);
  }
  
  Int_t nnewx = nx;
  Int_t nnewy = ny;
  if(!rebinx)
  {
    nnewx = noldx;
  }
  if(!rebiny)
  {
    nnewy = noldy;
  }
  
  double *xb = new double[nnewx+1];
  double *yb = new double[nnewy+1];
  for(i=0;i<nnewx+1;i++)
  {
    if(rebinx)
    {
      xb[i] = x[i];
    }
    else
    {
      xb[i] = h->GetXaxis()->GetBinLowEdge(i+1);
    }
  }
  for(i=0;i<nnewy+1;i++)
  {
    if(rebiny)
    {
      yb[i] = y[i];
    }
    else
    {
      yb[i] = h->GetYaxis()->GetBinLowEdge(i+1);
    }
  }
  
  if(nnewx>noldx)
  {
    cout<<"Error in MultiBinAnaHelper::Rebin2DHist(): number of new x bins is greater than the number of original x bins!  Quitting..."<<endl;
    return;
  }
  if(nnewy>noldy)
  {
    cout<<"Error in MultiBinAnaHelper::Rebin2DHist(): number of new y bins is greater than the number of original y bins!  Quitting..."<<endl;
    return;
  }
  
  double dx = h->GetXaxis()->GetBinWidth(1)*1e-4;
  double dy = h->GetYaxis()->GetBinWidth(1)*1e-4;
  
  string name;
  
  TH2D *old = (TH2D*)h->Clone();
  name = h->GetName();
  name += "_old";
  old->SetNameTitle(name.c_str(),"");
//   cout<<old->Integral()<<endl;
  
  h->Reset();
  h->SetBins(nnewx,xb,nnewy,yb);
  
  name = h->GetName();
  name += "_mid";
  TH2D *mid = new TH2D(name.c_str(),"",nnewx,xb,noldy,oldyb);
  mid->Sumw2();
  
  TH1D *proj;
  Double_t s,ds;
  
  if(rebinx)
  {
    for(iy=0;iy<noldy;iy++)
    {
      name = old->GetName();
      name += "_px";
      proj = old->ProjectionX(name.c_str(),iy+1,iy+1);
      //should be able to do a proj->Rebin(x) here but it doesn't work
      for(ix=0;ix<nnewx;ix++)
      {
        s=0;
        ds=0;
        for(jx=0;jx<noldx;jx++)
        {
          if(TMath::Abs(xb[ix]-oldxb[jx])<dx || (oldxb[jx]>xb[ix] && oldxb[jx]<xb[ix+1]))//low edges of old and new bin matches or low edge of old bin is in between new bin edges
          {
            s+=proj->GetBinContent(jx+1);
            ds+=proj->GetBinError(jx+1)*proj->GetBinError(jx+1);
          }
        }
        ds = TMath::Sqrt(ds);
        
        mid->SetBinContent(ix+1,iy+1,s);
        mid->SetBinError(ix+1,iy+1,ds);
      }
      proj->Reset();
    }
  }
  else
  {
    mid->Add(old);
  }
  
//   cout<<mid->Integral()<<endl;
  
  double tot = 0;
  
  if(rebiny)
  {
    for(ix=0;ix<nnewx;ix++)
    {
      name = old->GetName();
      name += "_py";
      proj = mid->ProjectionY(name.c_str(),ix+1,ix+1);
      //should be able to do a proj->Rebin(y) here but it doesn't work
      for(iy=0;iy<nnewy;iy++)
      {
//       cout<<"new bin "<<iy<<endl;
        s=0;
        ds=0;
        for(jy=0;jy<noldy;jy++)
        {
          if(TMath::Abs(yb[iy]-oldyb[jy])<dy || (oldyb[jy]>yb[iy] && oldyb[jy]<yb[iy+1]))//low edges of old and new bin matches or low edge of old bin is in between new bin edges
          {
            s+=proj->GetBinContent(jy+1);
            ds+=proj->GetBinError(jy+1)*proj->GetBinError(jy+1);
//          cout<<"old bin "<<jy<<endl;
          }
        }
        ds = TMath::Sqrt(ds);
        
        h->SetBinContent(ix+1,iy+1,s);
        h->SetBinError(ix+1,iy+1,ds);
        tot+=s;
      }
      proj->Reset();
    }
  }
  else
  {
    h->Add(mid);
  }
  
  mid->Delete();
  
  delete [] oldxb;
  delete [] oldyb;
  delete [] xb;
  delete [] yb;
  
//   cout<<h->Integral()<<endl;
//   cout<<tot<<endl;
  return;
}

void MultiBinAnaHelper::Rebin3DHist(TH3D*& h,Int_t nx,Double_t *x,Int_t ny,Double_t *y,Int_t nz,Double_t *z)
{
  //Be careful using this function!! It only works correctly if you give it correct input.  It is not smart enough to make sure that your new bins can really be constructed out of the old ones.
  
  Bool_t rebinx=true,rebiny=true,rebinz=true;
  if(!x || nx==0) rebinx = false;
  if(!y || ny==0) rebiny = false;
  if(!z || nz==0) rebinz = false;

  if(!rebinx && !rebiny && !rebinz) return; //do nothing
  
  Int_t i,ix,iy,iz;
  
  Int_t noldx = h->GetNbinsX();
  Int_t noldy = h->GetNbinsY();
  Int_t noldz = h->GetNbinsZ();
  double *oldxb = new double[noldx+1];
  double *oldyb = new double[noldy+1];
  double *oldzb = new double[noldz+1];
  for(ix = 0; ix < noldx + 1; ix++) oldxb[ix] = h->GetXaxis()->GetBinLowEdge(ix+1);
  for(iy = 0; iy < noldy + 1; iy++) oldyb[iy] = h->GetYaxis()->GetBinLowEdge(iy+1);
  for(iz = 0; iz < noldz + 1; iz++) oldzb[iz] = h->GetZaxis()->GetBinLowEdge(iz+1);
  
  Int_t nnewx = nx;
  Int_t nnewy = ny;
  Int_t nnewz = nz;
  if(!rebinx) nnewx = noldx;
  if(!rebiny) nnewy = noldy;
  if(!rebinz) nnewz = noldz;
  
  double *xb = new double[nnewx+1];
  double *yb = new double[nnewy+1];
  double *zb = new double[nnewz+1];
  for(i = 0; i < nnewx + 1; i++) {
    if(rebinx) xb[i] = x[i];
    else xb[i] = h->GetXaxis()->GetBinLowEdge(i+1);
  }
  for(i = 0; i < nnewy + 1; i++) {
    if(rebiny) yb[i] = y[i];
    else yb[i] = h->GetYaxis()->GetBinLowEdge(i+1);
  }
  for(i = 0; i < nnewz + 1; i++) {
    if(rebinz) zb[i] = z[i];
    else zb[i] = h->GetZaxis()->GetBinLowEdge(i+1);
  }
  
  if(nnewx > noldx) {
    cout<<"Error in MultiBinAnaHelper::Rebin3DHist(): number of new x bins is greater than the number of original x bins!  Quitting..."<<endl;
    return;
  }
  if(nnewy > noldy) {
    cout<<"Error in MultiBinAnaHelper::Rebin3DHist(): number of new y bins is greater than the number of original y bins!  Quitting..."<<endl;
    return;
  }
  if(nnewz > noldz) {
    cout<<"Error in MultiBinAnaHelper::Rebin3DHist(): number of new z bins is greater than the number of original z bins!  Quitting..."<<endl;
    return;
  }
  
  string origname = h->GetName();
  bool issumw2 = h->GetSumw2N();

  string name;

  TH3D *old = (TH3D*)h->Clone();
  name = origname + "_old";
  old->SetNameTitle(name.c_str(),"");

  int ik;
  double jx,jy,jz;
  double ds;
  double bincont,binerr;

  delete h;
  h = new TH3D(origname.c_str(),"",nnewx,xb,nnewy,yb,nnewz,zb);
  if (issumw2) h->Sumw2();

  for (iz = 1; iz <= noldz; iz++) {
    for (iy = 1; iy <= noldy; iy++) {
      for (ix = 1; ix <= noldx; ix++) {

	jx = old->GetXaxis()->GetBinCenter(ix);
	jy = old->GetYaxis()->GetBinCenter(iy);
	jz = old->GetZaxis()->GetBinCenter(iz);
	bincont = old->GetBinContent(ix,iy,iz);
	binerr = old->GetBinError(ix,iy,iz);

	ik = h->FindBin(jx,jy,jz);
	ds = h->GetBinError(ik);
	ds *= ds;
	ds += binerr*binerr;
	ds = TMath::Sqrt(ds);
	h->AddBinContent(ik,bincont);
	h->SetBinError(ik,ds);

      }
    }
  }
 
  //  cout << "MBH says integral was " << old->Integral() << ", now is " << h->Integral() << endl;
  delete old;
  
  delete [] oldxb;
  delete [] oldyb;
  delete [] oldzb;
  delete [] xb;
  delete [] yb;
  delete [] zb;
  
  return;
}

void MultiBinAnaHelper::CopyYtoZ(TH2D* h2, TH3D*& h3) {

  //Takes in filled h2 and empty h3
  //Copies X and Y axes of h2 to h3 AND Y axis of h2 to Z axis of h3
  //Fills h3 with contents of h2 such that YZ proj of h3 is diagonal

  if (!h2) {
    cout << "The TH2D pointer is null! Quitting..." << endl;
    return;
  }

  int nx = h2->GetNbinsX();
  int ny = h2->GetNbinsY();
  int ix,iy;
  double *x = new double[nx+1];
  double *y = new double[ny+1];

  for(ix = 0; ix < nx + 1; ix++) x[ix] = h2->GetXaxis()->GetBinLowEdge(ix+1);
  for(iy = 0; iy < ny + 1; iy++) y[iy] = h2->GetYaxis()->GetBinLowEdge(iy+1);
  
  string histname = string(h2->GetName());
  histname += "_CopyYtoZ";

  h3 = new TH3D(histname.c_str(),"",nx,x,ny,y,ny,y);
  double bincont,binerr;

  for (ix = 1; ix <= nx; ix++) {
    for (iy = 1; iy <= ny; iy++) {
      bincont = h2->GetBinContent(ix,iy);
      binerr = h2->GetBinError(ix,iy);
      h3->SetBinContent(ix,iy,iy,bincont);
      h3->SetBinError(ix,iy,iy,binerr);
    }
  }

  delete [] x;
  delete [] y;

}
