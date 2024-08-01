/// \file
/// \ingroup Analysis
/// This program:
///       - reads charge measurements using TCT setup with a Red Laser (~600 nm)
///       - fits the charge profiles using s-curves
///       - use the fitted curve to measure the inter-strip distance between two DC strips
/// \author Ashish Bisht (abisht@fbk.eu)
/// ______________________________________________

#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TStyle.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>


using namespace std;

void SetStyle(Bool_t threeDimension = kFALSE);
void readTextFile(const TString& filename, Int_t& dataPoints, vector<Float_t>& column1, vector<Float_t>& column2, vector<Float_t>& column3, vector<Float_t>& column4);

int interStripDistance_Analysis()
{
  SetStyle();
  
  TCanvas *canvas = new TCanvas("canvas","",1200,1000);  
  TLegend *leg = new TLegend(0.6196674,0.7419355,0.8293565,0.9330025,NULL,"brNDC");

  // Input data files
  TString leftStripDataFile = "data/leftStripData.txt";
  TString rightStripDataFile = "data/rightStripData.txt";


  // Define a fit function for the left strip small curve
  TF1 *erf_left_small_curve = new TF1("erf_left_small_curve","(TMath::Erf(-(x-[0])/[1])+1)*[2]+[3]",155,200);
  erf_left_small_curve->SetLineColor(kRed);
  erf_left_small_curve->SetLineWidth(3);
   
  // Define a fit function for the left strip large curve 
  TF1 *erf_left_large_curve = new TF1("erf_left_large_curve","(TMath::Erf(-(x-[0])/[1])+1)*[2]+[3]",110,150);
  erf_left_large_curve->SetLineColor(kGreen);
  erf_left_large_curve->SetLineWidth(3);
  
  // Define a fit function for the right strip small curve 
  TF1 *erf_right_small_curve = new TF1("erf_right_small_curve","(TMath::Erf((x-[0])/[1])+1)*[2]+[3]",170,215);
  erf_right_small_curve->SetLineColor(kRed);
  erf_right_small_curve->SetLineWidth(3);
  
  // Define a fit function for the right strip large curve 
  TF1 *erf_right_large_curve = new TF1("erf_right_large_curve","(TMath::Erf((x-[0])/[1])+1)*[2]+[3]",220,260);
  erf_right_large_curve->SetLineColor(kGreen);
  erf_right_large_curve->SetLineWidth(3);
   
  
  Int_t nPointsLeft=0, nPointsRight=0;
  vector <Float_t> scanningDistance_left;
  vector <Float_t> errX_left; 
  vector <Float_t> normalizedCharge_left;
  vector <Float_t> errCharge_left;
  readTextFile(leftStripDataFile, nPointsLeft, scanningDistance_left, errX_left, normalizedCharge_left, errCharge_left);
  
  vector <Float_t> scanningDistance_right;
  vector <Float_t> errX_right; 
  vector <Float_t> normalizedCharge_right;
  vector <Float_t> errCharge_right;
  readTextFile(rightStripDataFile, nPointsRight, scanningDistance_right, errX_right, normalizedCharge_right, errCharge_right);
  
  // Convert vectors to arrays of doubles
  auto toDoubleArray = [](const vector<Float_t>& vec) -> double*
  {
    double* arr = new double[vec.size()];
    for (size_t i = 0; i < vec.size(); ++i)
      {
	arr[i] = static_cast<double>(vec[i]);
      }
    return arr;
  };

  double* scanningDistance_left_arr = toDoubleArray(scanningDistance_left);
  double* errX_left_arr = toDoubleArray(errX_left);
  double* normalizedCharge_left_arr = toDoubleArray(normalizedCharge_left);
  double* errCharge_left_arr = toDoubleArray(errCharge_left);

  double* scanningDistance_right_arr = toDoubleArray(scanningDistance_right);
  double* errX_right_arr = toDoubleArray(errX_right);
  double* normalizedCharge_right_arr = toDoubleArray(normalizedCharge_right);
  double* errCharge_right_arr = toDoubleArray(errCharge_right);
  
  TGraphErrors *gr_leftStrip = new TGraphErrors(nPointsLeft, scanningDistance_left_arr, normalizedCharge_left_arr, errX_left_arr, errCharge_left_arr);
  gr_leftStrip->SetMarkerStyle(29);
  gr_leftStrip->SetMarkerSize(1);
  gr_leftStrip->SetMarkerColor(kBlack);
  gr_leftStrip->SetLineColor(kBlack);
  gr_leftStrip->GetYaxis()->SetRangeUser(-20, 400);
  gr_leftStrip->GetXaxis()->SetRangeUser(90, 280);      

  erf_left_small_curve->SetParameters(130, 10, 35, 0);
  gr_leftStrip->Fit("erf_left_small_curve","RMS");
  erf_left_large_curve->SetParameters(110, 10, 125);
  erf_left_large_curve->FixParameter(3, erf_left_small_curve->GetParameter(3));
  gr_leftStrip->Fit("erf_left_large_curve","RMS+");
  gr_leftStrip->Draw("ap");
  gr_leftStrip->GetXaxis()->SetTitle("Scanning Distance (#mum)");
  gr_leftStrip->GetYaxis()->SetTitle("Norm. Charge (arb.)");
  erf_left_large_curve->SetLineStyle(3);
  erf_left_large_curve->SetRange(90, 220);
  erf_left_large_curve->Draw("same");
  Double_t chi_left = erf_left_large_curve->GetChisquare()/erf_left_large_curve->GetNDF();
  leg->AddEntry(gr_leftStrip, Form("#frac{#chi^{2}}{ndf} = %2.2lf ", chi_left), "pl");

  TGraphErrors *gr_rightStrip = new TGraphErrors(nPointsRight, scanningDistance_right_arr, normalizedCharge_right_arr, errX_right_arr, errCharge_right_arr);
  gr_rightStrip->SetMarkerStyle(29);
  gr_rightStrip->SetMarkerSize(1);
  gr_rightStrip->SetMarkerColor(kBlack);
  gr_rightStrip->SetLineColor(kBlack);
  gr_rightStrip->GetYaxis()->SetRangeUser(-20, 400);
  gr_rightStrip->GetXaxis()->SetRangeUser(90, 280);      

  erf_right_small_curve->SetParameters(190, 10, 35, 0);
  gr_rightStrip->Fit("erf_right_small_curve","RMS");
  erf_right_large_curve->SetParameters(220, 10, 125);
  erf_right_large_curve->FixParameter(3, erf_right_small_curve->GetParameter(3));
  gr_rightStrip->Fit("erf_right_large_curve","RMS+");
  gr_rightStrip->Draw("psame");
  gr_rightStrip->GetXaxis()->SetTitle("Scanning Distance (#mum)");
  gr_rightStrip->GetYaxis()->SetTitle("Norm. Charge (arb.)");
  erf_right_large_curve->SetLineStyle(3);
  erf_right_large_curve->SetRange(150, 280);
  erf_right_large_curve->Draw("same");
  Double_t chi_right = erf_right_large_curve->GetChisquare()/erf_right_large_curve->GetNDF();
  leg->AddEntry(gr_rightStrip, Form("#frac{#chi^{2}}{ndf} = %2.2lf ", chi_right), "pl");

  Double_t interStripDistance = erf_right_large_curve->GetParameter(0)-erf_left_large_curve->GetParameter(0);
  Double_t errISD = TMath::Sqrt(TMath::Power(erf_right_large_curve->GetParError(0),2)+TMath::Power(erf_left_large_curve->GetParError(0),2));

  TPaveText *pt = new TPaveText(0.163606,0.7981557,0.4941569,0.9241803,"nbNDC");
  pt->SetFillColor(kWhite);
  pt->AddText(Form("Inter-strip Distance|_{Red} = %2.2lf #pm %2.2lf #mum", interStripDistance, errISD))->SetTextColor(kBlue);
  pt->Draw();
  leg->Draw();
  canvas->SaveAs("figures/interStripDistance.png","png");
  
  return 0;
}

// Function to read text files
void readTextFile(const TString& filename, Int_t& dataPoints, vector<Float_t>& column1, vector<Float_t>& column2, vector<Float_t>& column3, vector<Float_t>& column4)
{
  ifstream file(filename);
  if(!file.is_open())
    {
      cerr<<"Error opening file\n";
      return;
    }

  string line;
  Bool_t isFirstLine = kTRUE;

  while(getline(file, line))
    {
      if(isFirstLine)
	{
	  isFirstLine = kFALSE;
	  continue; //Skip the first line
	}
      dataPoints++; //Increase for each data line
      stringstream ss(line);
      string item;
      Float_t value;
      Int_t colIndex = 0;

      while(getline(ss, item, ','))
	{
	  stringstream itemStream(item);
	  itemStream>>value;

	  switch(colIndex)
	    {
	    case 0:
	      column1.push_back(value);
	      break;
	    case 1:
	      column2.push_back(value);
	      break;
	    case 2:
	      column3.push_back(value);
	      break;
	    case 3:
	      column4.push_back(value);
	      break;
	    default:
	      cerr << "Unexpected column index" << std::endl;
	      return;
	    }
	  colIndex++;
	}
    }
  file.close();
}

// Function to increase the plot beauty
void SetStyle(Bool_t threeDimension = kFALSE)
{
  gErrorIgnoreLevel=kError; //Removes annoying Potential memory leak warnings
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTopMargin(0.04854282);
  //gStyle->SetPadRightMargin(0.03025943);
  gStyle->SetPadBottomMargin(0.1353861);//0.015
  gStyle->SetPadLeftMargin(0.1418293);//0.21
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kBlue);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.04,"xyz");
  //gStyle->SetLabelSize(0.06,"x");
  gStyle->SetLabelOffset(0.005,"y");
  gStyle->SetLabelOffset(0.006,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  //  gStyle->SetTitleOffset(1.14,"y");//1.81
  if(threeDimension)
    {
      gStyle->SetPadRightMargin(0.18);
      gStyle->SetTitleOffset(1.15,"y");//1.81
    }
  else
    gStyle->SetTitleOffset(1.10,"y");//1.81
  gStyle->SetTitleOffset(0.95,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetTickLength(0.03,"X");
  gStyle->SetTickLength(0.03,"Y"); 
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendTextSize(0.04);
  gStyle->SetNdivisions(505,"xy");
}
