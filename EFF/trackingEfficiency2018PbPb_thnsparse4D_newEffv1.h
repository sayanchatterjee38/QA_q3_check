#ifndef TRKEFF2018PBPB
#define TRKEFF2018PBPB

#include "TFile.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TMath.h"
#include <iostream>
#include <string>
#include "THnSparse.h"


class TrkEff2018PbPb{
public:

  TrkEff2018PbPb( std::string collectionName = "general", bool isQuiet_ = false ,std::string filePath_pt = "", std::string filePath_mb = "", std::string filePath_plus = "", std::string filePath_minus = "", std::string filePath_pix = "");


  ~TrkEff2018PbPb();
  /*
  float getCorrection(float pt, float eta, int hiBin);
  float getEfficiency( float pt, float eta, int hiBin, bool passesCheck = false);
  float getFake( float pt, float eta, int hiBin, bool passesCheck = false);
  float getMultiple( float pt, float eta, int hiBin, bool passesCheck = false);
  float getSecondary( float pt, float eta, int hiBin, bool passesCheck = false);
  */

  double getCorrection(double pt, double eta, double phi, double hiBin);
  double getEfficiency( double pt, double eta, double phi, double hiBin, bool passesCheck = false);
  double getFake( double pt, double eta, double phi, double hiBin, bool passesCheck = false);
  double getMultiple( double pt, double eta, double phi, double hiBin, bool passesCheck = false);
  double getSecondary( double pt, double eta, double phi, double hiBin, bool passesCheck = false);
  
private:

  inline bool checkBounds(double pt, double eta, double phi, double hiBin);

  std::string mode;
  bool isQuiet;

  TFile * trkEff;
  TFile * trkFake;
  THnSparseD * eff;
  THnSparseD * fake;

  //new
  TFile * trkMul;
  TFile * trkSec;
  THnSparseD * mul;
  THnSparseD * sec;
  
};

inline bool TrkEff2018PbPb::checkBounds(double pt, double eta, double phi, double hiBin){
  if( TMath::Abs(eta) > 2.4 ){
    if( ! isQuiet) std::cout << "TrkEff2018PbPb: track outside |eta|<2.4, please apply this cut!  I am returning a correction factor of 0 for this track for now." << std::endl;
    return false;
  }
  
  if( hiBin <0. || hiBin > 199.){
    if( ! isQuiet) std::cout << "TrkEff2018PbPb: hiBin is outside the range [0,199].  Please fix this issue!  I am returning a correction factor of 0 for this track for now." << std::endl;
    return false;
  }
  
  if( pt< 0.5 || pt > 500.0 ){
    if( ! isQuiet) std::cout << "TrkEff2018PbPb: pT is outside the range [0,500].  I am returning a correction factor of 0 for this track for now." << std::endl;
    return false;
  }

  if( phi < -3.14159265 || phi > 3.14159265 ){
    if( ! isQuiet) std::cout << "TrkEff2018PbPb: phi is outside the range [-pi,pi].  I am returning a correction factor of 0 for this track for now." << std::endl;
    return false;
  }

  
  return true;
}

double TrkEff2018PbPb::getCorrection(double pt, double eta, double phi, double hiBin){
  if( !checkBounds(pt, eta, phi, hiBin) ) return 0;

  double efficiency = getEfficiency(pt, eta, phi, hiBin, true);
  double fake = getFake(pt, eta, phi, hiBin, true);

  double multiple = getMultiple(pt, eta, phi, hiBin, true);
  double secon = getSecondary(pt, eta, phi, hiBin, true);

  //protect against dividing by 0
  if(efficiency > 0.001){
    return (1-fake)*(1-secon)/(efficiency*(1+multiple));
  } else {
    if( ! isQuiet ) std::cout << "TrkEff2018PbPb: Warning! Tracking efficiency is very low for this track (close to dividing by 0).  Returning correction factor of 0 for this track for now." << std::endl;
    return 0;
  }
}

double TrkEff2018PbPb::getEfficiency( double pt, double eta, double phi, double hiBin, bool passesCheck){
  if( !passesCheck){
    if(  !checkBounds(pt, eta, phi, hiBin) ) return 0;
  }

  double eff4d[4] = {eta, pt, phi, hiBin};
  
  if( mode.compare("general") == 0 || mode.compare("generalMB-") == 0 || mode.compare("generalMB+") == 0){
    return eff->GetBinContent( eff->GetBin(eff4d) );
  }
  
  return 0;
}

double TrkEff2018PbPb::getFake( double pt, double eta, double phi, double hiBin, bool passesCheck){
  if( !passesCheck){
    if(  !checkBounds(pt, eta, phi, hiBin) ) return 0;
  }

  double fake4d[4] = {eta, pt, phi, hiBin};
  
  if( mode.compare("general") == 0 || mode.compare("generalMB-") == 0 || mode.compare("generalMB+") == 0){
    return fake->GetBinContent( fake->GetBin(fake4d) );
  }

  return 0;
}

double TrkEff2018PbPb::getMultiple( double pt, double eta, double phi, double hiBin, bool passesCheck){
  if( !passesCheck){
    if(  !checkBounds(pt, eta, phi, hiBin) ) return 0;
  }

  double mul4d[4] = {eta, pt, phi, hiBin};
  
  if( mode.compare("general") == 0 || mode.compare("generalMB-") == 0 || mode.compare("generalMB+") == 0){
    return mul->GetBinContent( mul->GetBin(mul4d) );
  }
  
  return 0;
}


double TrkEff2018PbPb::getSecondary( double pt, double eta, double phi, double hiBin, bool passesCheck){
  if( !passesCheck){
    if(  !checkBounds(pt, eta, phi, hiBin) ) return 0;
  }

  double sec4d[4] = {eta, pt, phi, hiBin};
  
  if( mode.compare("general") == 0 || mode.compare("generalMB-") == 0 || mode.compare("generalMB+") == 0){
    return sec->GetBinContent( sec->GetBin(sec4d) );
  }
  
  return 0;
}




TrkEff2018PbPb::TrkEff2018PbPb(std::string collectionName, bool isQuiet_, std::string filePath_pt, std::string filePath_mb, std::string filePath_plus, std::string filePath_minus,  std::string filePath_pix)

{
  isQuiet = isQuiet_;
  mode = collectionName;
  
  if( collectionName.compare("general") == 0 ) 
    {
      trkEff = new TFile((filePath_mb).c_str(),"READ");
      eff = (THnSparseD*) trkEff->Get("Eff4D");

      trkFake = new TFile((filePath_mb).c_str(),"READ");
      fake = (THnSparseD*) trkFake->Get("Fak4D");

      trkMul = new TFile((filePath_mb).c_str(),"READ");
      mul = (THnSparseD*) trkMul->Get("Mul4D");

      trkSec = new TFile((filePath_mb).c_str(),"READ");
      sec = (THnSparseD*) trkSec->Get("Sec4D");
    } 
  
  else if( collectionName.compare("generalMB+") == 0 )
    {
      trkEff = new TFile((filePath_plus).c_str(),"READ");
      eff = (THnSparseD*) trkEff->Get("Eff4D");
    
      trkFake = new TFile((filePath_plus).c_str(),"READ");
      fake = (THnSparseD*) trkFake->Get("Fak4D");
    
      trkMul = new TFile((filePath_plus).c_str(),"READ");
      mul = (THnSparseD*) trkMul->Get("Mul4D");

      trkSec = new TFile((filePath_plus).c_str(),"READ");
      sec = (THnSparseD*) trkSec->Get("Sec4D");
    } 
  
  else if( collectionName.compare("generalMB-") == 0 )
    {
      trkEff = new TFile((filePath_minus).c_str(),"READ");
      eff = (THnSparseD*) trkEff->Get("Eff4D");
      
      
      trkFake = new TFile((filePath_minus).c_str(),"READ");
      fake = (THnSparseD*) trkFake->Get("Fak4D");

      trkMul = new TFile((filePath_minus).c_str(),"READ");
      mul = (THnSparseD*) trkMul->Get("Mul4D");
      
      trkSec = new TFile((filePath_minus).c_str(),"READ");
      sec = (THnSparseD*) trkSec->Get("Sec4D");
    } 

  else 
   {
      std::cout << "Error! incorrect collectionName parameter in TrkEff2018PbPb constructor.  Please try 'general', 'generalMB+', 'generalMB-', or 'pixel.'" << std::endl;
   } 
}

TrkEff2018PbPb::~TrkEff2018PbPb(){
  trkEff->Close();
  trkFake->Close();
  trkMul->Close();
  trkSec->Close();
}

#endif
