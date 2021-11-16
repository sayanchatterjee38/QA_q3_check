// -*- C++ -*-
//   
//         Author: Sayan Chatterjee
//         Class: SayanCMW 
//         Checked version 24/06/2020
//         Last Modified on 12/11/2021

// This is a fresh cleaned version of code of QA. I just added mean pt vs eta in the track loop on 16/02/2021.



// CMSSW include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// user include files

#include "Analyzers/SayanCMW/interface/SayanCMW.h"
#include <string>

//  constructors and destructor

SayanCMW::SayanCMW(const edm::ParameterSet& iConfig) :  //Parametrized Constructor
  //******TRACKED PARAMETER********
  
  //tracks & vertex 
  trackTags_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
  trackTagsgen_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("tracksgen"))),
  mvaSrc_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("mvaSrc"))),

  vtxTags_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex"))),

  //centrality bin
  cent_(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centrality"))),
  cent_bin_(consumes<int>(iConfig.getParameter<edm::InputTag>("centbin"))),

  //******UNTRACKED PARAMETER****************
  
  //Event classifier
  centmin_(iConfig.getUntrackedParameter<int>("centmin")),
  centmax_(iConfig.getUntrackedParameter<int>("centmax")),

  //vertex selection
  zminVtx_(iConfig.getUntrackedParameter<double>("zminVtx")),
  zmaxVtx_(iConfig.getUntrackedParameter<double>("zmaxVtx")),
  rhomaxVtx_(iConfig.getUntrackedParameter<double>("rhomaxVtx")),

  //EFF Correction
  fpt_(iConfig.getUntrackedParameter<edm::InputTag>("fpt")),
  fmb_(iConfig.getUntrackedParameter<edm::InputTag>("fmb")),
  fplus_(iConfig.getUntrackedParameter<edm::InputTag>("fplus")),
  fminus_(iConfig.getUntrackedParameter<edm::InputTag>("fminus")),
  fpix_(iConfig.getUntrackedParameter<edm::InputTag>("fpix")),

  //track selection
  pTmin_(iConfig.getUntrackedParameter< std::vector< double > >("pTminTrk")),
  pTmax_(iConfig.getUntrackedParameter< std::vector< double > >("pTmaxTrk")),
  pTmin_trg_(iConfig.getUntrackedParameter< std::vector< double > >("pTminTrk_trg")),
  pTmax_trg_(iConfig.getUntrackedParameter< std::vector< double > >("pTmaxTrk_trg")),
  pTmin_ass_(iConfig.getUntrackedParameter< std::vector< double > >("pTminTrk_ass")),
  pTmax_ass_(iConfig.getUntrackedParameter< std::vector< double > >("pTmaxTrk_ass")),
  ptmin_(iConfig.getUntrackedParameter<double>("ptmin")),
  ptmax_(iConfig.getUntrackedParameter<double>("ptmax")),
  etamin_(iConfig.getUntrackedParameter<double>("etamin")),
  etamax_(iConfig.getUntrackedParameter<double>("etamax")),
  //~20.10.2021
  algoParameters_(iConfig.getParameter<std::vector<int> >("algoParameters")),

  //bool expression
  ifMcgen_(iConfig.getUntrackedParameter<bool>("ifMcgen"))
{
  //---Dihedron header file----  
  evt_ = new DiHadronCorrelationEvt(pTmin_.size(), pTmin_trg_.size(), pTmin_ass_.size());
  
  //*****Defining Histograms & Profile Histograms********************
   
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  TProfile::SetDefaultSumw2();

  TString f_PT(fpt_.label().c_str());
  edm::FileInPath f1(Form("Analyzers/SayanCMW/data/EFF/effHpT/%s",f_PT.Data()));
 
  TString f_MB(fmb_.label().c_str());
  edm::FileInPath f2(Form("Analyzers/SayanCMW/data/EFF/effMB/%s",f_MB.Data()));

  TString f_Plus(fplus_.label().c_str());
  edm::FileInPath f3(Form("Analyzers/SayanCMW/data/EFF/plus/%s",f_Plus.Data()));

  TString f_Minus(fminus_.label().c_str());
  edm::FileInPath f4(Form("Analyzers/SayanCMW/data/EFF/minus/%s",f_Minus.Data()));

  TString f_Pix(fpix_.label().c_str());
  edm::FileInPath f5(Form("Analyzers/SayanCMW/data/EFF/effPix/%s",f_Pix.Data()));

   
  TrkEff = new TrkEff2018PbPb(  "general", false, f1.fullPath(), f2.fullPath(), f3.fullPath(), f4.fullPath(), f5.fullPath());
  TrkEff1 = new TrkEff2018PbPb(  "generalMB+", false, f1.fullPath(), f2.fullPath(), f3.fullPath(), f4.fullPath(), f5.fullPath());
  TrkEff2 = new TrkEff2018PbPb(  "generalMB-", false, f1.fullPath(), f2.fullPath(), f3.fullPath(), f4.fullPath(), f5.fullPath());   
  
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  
  TFileDirectory fGlobalHist = fs->mkdir("QAplots");
 
  //const Int_t  nptbin              = 26;
  //Double_t  ptbining[nptbin+1]     = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.10, 1.15, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.0, 2.5, 3.0, 4.0, 5.0};

  const Int_t netaBin = 18;
  Double_t etabining[netaBin+1] = {-2.4, -2.0, -1.6, -1.4, -1.3, -1.2, -1.0, -0.8, -0.4, 0.0, 0.4, 0.8, 1.0, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4};
   
  
  hpt          = fGlobalHist.make<TH1F>("hpt", "", 45, 0.5, 5.0 );
  heta         = fGlobalHist.make<TH1F>("heta", "", 48, -2.4, 2.4 );
  heta_nbin    = fGlobalHist.make<TH1F>("heta_nbin", "", netaBin, etabining );
  hptP         = fGlobalHist.make<TH1F>("hptP", "", 45, 0.5, 5.0 );
  hetaP        = fGlobalHist.make<TH1F>("hetaP", "", 48, -2.4, 2.4 );
  hetaP_nbin   = fGlobalHist.make<TH1F>("hetaP_nbin", "", netaBin, etabining );
  hptN         = fGlobalHist.make<TH1F>("hptN", "", 45, 0.5, 5.0 );
  hetaN        = fGlobalHist.make<TH1F>("hetaN", "", 48, -2.4, 2.4 );
  hetaN_nbin   = fGlobalHist.make<TH1F>("hetaN_nbin", "", netaBin, etabining );
  hphi         = fGlobalHist.make<TH1F>("hphi", "", 63, -3.15, 3.15 );
  hphiP        = fGlobalHist.make<TH1F>("hphiP", "", 63, -3.15, 3.15 );
  hphiN        = fGlobalHist.make<TH1F>("hphiN", "", 63, -3.15, 3.15 );


  hZBestVtx    = fGlobalHist.make<TH1F>("hZvtx", "", 600, -30.0, 30.0);
  hcent_bin    = fGlobalHist.make<TH1I>("hcent", "", 100, 0, 100);
  
  th2d_etaphi    = fGlobalHist.make<TH2D>("th2d_etaphi", "", netaBin, etabining, 63, -3.15, 3.15 );
  th2d_etaphi_pos    = fGlobalHist.make<TH2D>("th2d_etaphi_pos", "", netaBin, etabining, 63, -3.15, 3.15 );
  th2d_etaphi_neg    = fGlobalHist.make<TH2D>("th2d_etaphi_neg", "", netaBin, etabining, 63, -3.15, 3.15 );

  //th2d_algomva    = fGlobalHist.make<TH2F>("th2d_algomva", "", 40, -20, 20, 3000, -1.5, 1.5);
  //th1d_algo    = fGlobalHist.make<TH1D>("th1d_algo", "", 40, -20, 20);
  //  th1d_mva    = fGlobalHist.make<TH1D>("th1d_mva", "", 300, -1.5, 1.5);

  tp1d_mpteta_nbin   = fGlobalHist.make<TProfile>("tp1d_mpt_nbin", "", netaBin, etabining, -1e10, 1e10);
  tp1d_mptetaP_nbin   = fGlobalHist.make<TProfile>("tp1d_mptP_nbin", "", netaBin, etabining, -1e10, 1e10);
  tp1d_mptetaN_nbin   = fGlobalHist.make<TProfile>("tp1d_mptN_nbin", "", netaBin, etabining, -1e10, 1e10);

  tp1d_mptphi   = fGlobalHist.make<TProfile>("tp1d_mptphi", "", 63, -3.15, 3.15, -1e10, 1e10);
  tp1d_mptphiP   = fGlobalHist.make<TProfile>("tp1d_mptphiP", "", 63, -3.15, 3.15, -1e10, 1e10);
  tp1d_mptphiN   = fGlobalHist.make<TProfile>("tp1d_mptphiN", "", 63, -3.15, 3.15, -1e10, 1e10);

  tp1d_charge_eta   = fGlobalHist.make<TProfile>("tp1d_charge_eta", "", netaBin, etabining, -1e10, 1e10);
  tp1d_charge_phi   = fGlobalHist.make<TProfile>("tp1d_charge_phi", "", 63, -3.15, 3.15, -1e10, 1e10);
  
  //thnsparse4d histogram
  const Int_t ndims =4;    //eta, pT, phi, centBin                                                                                                    
  Int_t bins[ndims] = {18, 45, 63, 8};
  Double_t xmin[ndims]={0., 0.5, -3.15, 0.};
  Double_t xmax[ndims] = {10., 5.0, 3.15, 10.};

  const int nVarBins0 = 18;
  Double_t varBins0[nVarBins0+1] = { -2.4, -2.0, -1.6, -1.4, -1.3, -1.2, -1.0, -0.8, -0.4, 0.0, 0.4, 0.8, 1.0, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4 };

  //  const int nVarBins1 = 50;
  //Double_t varBins1[nVarBins1+1] = {  0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
  //				      1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0, 12.0, 15.0, 
  //				      20.0, 25.0, 30.0, 45.0, 60.0, 90.0, 120.0, 180.0, 300.0, 500.0 };

  const int nVarBins3 = 8;
  Double_t varBins3[nVarBins3+1] = { 0.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 140.0, 200.0 };

  hist4D = fGlobalHist.make<THnSparseD>("hist4D",";#eta;p_{T};#phi;occ",  ndims, bins, xmin, xmax);
  hist4DP = fGlobalHist.make<THnSparseD>("hist4DP",";#eta;p_{T};#phi;occ",  ndims, bins, xmin, xmax);
  hist4DN = fGlobalHist.make<THnSparseD>("hist4DN",";#eta;p_{T};#phi;occ",  ndims, bins, xmin, xmax);

  hist4D->GetAxis(0)->Set(nVarBins0, varBins0);
  //hist4D->GetAxis(1)->Set(nVarBins1, varBins1);
  hist4D->GetAxis(3)->Set(nVarBins3, varBins3);

  hist4DP->GetAxis(0)->Set(nVarBins0, varBins0);
  //hist4DP->GetAxis(1)->Set(nVarBins1, varBins1);
  hist4DP->GetAxis(3)->Set(nVarBins3, varBins3);

  hist4DN->GetAxis(0)->Set(nVarBins0, varBins0);
  //hist4DN->GetAxis(1)->Set(nVarBins1, varBins1);
  hist4DN->GetAxis(3)->Set(nVarBins3, varBins3);

  hist4D->Sumw2();
  hist4DP->Sumw2();
  hist4DN->Sumw2();

  TFileDirectory f_effw = fs->mkdir("Effeciency_weight");

  hpt_w          = f_effw.make<TH1F>("hpt_w", "", 45, 0.5, 5.0 );
  heta_w         = f_effw.make<TH1F>("heta_w", "", 48, -2.4, 2.4 );
  heta_nbin_w    = f_effw.make<TH1F>("heta_nbin_w", "", netaBin, etabining );
  hptP_w         = f_effw.make<TH1F>("hptP_w", "", 45, 0.5, 5.0 );
  hetaP_w        = f_effw.make<TH1F>("hetaP_w", "", 48, -2.4, 2.4 );
  hetaP_nbin_w   = f_effw.make<TH1F>("hetaP_nbin_w", "", netaBin, etabining );
  hptN_w         = f_effw.make<TH1F>("hptN_w", "", 45, 0.5, 5.0 );
  hetaN_w        = f_effw.make<TH1F>("hetaN_w", "", 48, -2.4, 2.4 );
  hetaN_nbin_w   = f_effw.make<TH1F>("hetaN_nbin_w", "", netaBin, etabining );
  hphi_w         = f_effw.make<TH1F>("hphi_w", "", 63, -3.15, 3.15 );
  hphiP_w        = f_effw.make<TH1F>("hphiP_w", "", 63, -3.15, 3.15 );
  hphiN_w        = f_effw.make<TH1F>("hphiN_w", "", 63, -3.15, 3.15 );

  th2d_etaphi_w    = f_effw.make<TH2D>("th2d_etaphi_w", "", netaBin, etabining, 63, -3.15, 3.15 );
  th2d_etaphi_pos_w    = f_effw.make<TH2D>("th2d_etaphi_pos_w", "", netaBin, etabining, 63, -3.15, 3.15 );
  th2d_etaphi_neg_w    = f_effw.make<TH2D>("th2d_etaphi_neg_w", "", netaBin, etabining, 63, -3.15, 3.15 );

  tp1d_mpteta_nbin_w   = f_effw.make<TProfile>("tp1d_mpt_nbin_w", "", netaBin, etabining, -1e10, 1e10);
  tp1d_mptetaP_nbin_w   = f_effw.make<TProfile>("tp1d_mptP_nbin_w", "", netaBin, etabining, -1e10, 1e10);
  tp1d_mptetaN_nbin_w   = f_effw.make<TProfile>("tp1d_mptN_nbin_w", "", netaBin, etabining, -1e10, 1e10);
 
  tp1d_mptphi_w   = f_effw.make<TProfile>("tp1d_mptphi_w", "", 63, -3.15, 3.15, -1e10, 1e10);
  tp1d_mptphiP_w   = f_effw.make<TProfile>("tp1d_mptphiP_w", "", 63, -3.15, 3.15, -1e10, 1e10);
  tp1d_mptphiN_w   = f_effw.make<TProfile>("tp1d_mptphiN_w", "", 63, -3.15, 3.15, -1e10, 1e10);

  tp1d_charge_eta_w   = f_effw.make<TProfile>("tp1d_charge_eta_w", "", netaBin, etabining, -1e10, 1e10);
  tp1d_charge_phi_w   = f_effw.make<TProfile>("tp1d_charge_phi_w", "", 63, -3.15, 3.15, -1e10, 1e10);
  
  //corrected 4D
  hist4D_w = f_effw.make<THnSparseD>("hist4D_w",";#eta;p_{T};#phi;occ",  ndims, bins, xmin, xmax);
  hist4DP_w = f_effw.make<THnSparseD>("hist4DP_w",";#eta;p_{T};#phi;occ",  ndims, bins, xmin, xmax);
  hist4DN_w = f_effw.make<THnSparseD>("hist4DN_w",";#eta;p_{T};#phi;occ",  ndims, bins, xmin, xmax);
  
  hist4D_w->GetAxis(0)->Set(nVarBins0, varBins0);
  //  hist4D_w->GetAxis(1)->Set(nVarBins1, varBins1);
  hist4D_w->GetAxis(3)->Set(nVarBins3, varBins3);

  hist4DP_w->GetAxis(0)->Set(nVarBins0, varBins0);
  //hist4DP_w->GetAxis(1)->Set(nVarBins1, varBins1);
  hist4DP_w->GetAxis(3)->Set(nVarBins3, varBins3);

  hist4DN_w->GetAxis(0)->Set(nVarBins0, varBins0);
  //hist4DN_w->GetAxis(1)->Set(nVarBins1, varBins1);
  hist4DN_w->GetAxis(3)->Set(nVarBins3, varBins3);

  hist4D_w->Sumw2();
  hist4DP_w->Sumw2();
  hist4DN_w->Sumw2();
 }

SayanCMW::~SayanCMW() // Destructor 
{
  delete evt_;
}

//----------------method called once each job just before starting event loop---------------------------

void
SayanCMW::beginJob()
{
}

//---------------method called for each event-------------------------------------------------------

void
SayanCMW::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;

  // ----------------- centrality selection -------------------------------
  edm::Handle< reco::Centrality > centrality;
  iEvent.getByToken(cent_, centrality);

  edm::Handle< int > cbin;
  iEvent.getByToken(cent_bin_, cbin);
  centBin = *cbin;
  if(centBin < 0)
    {
      edm::LogWarning ("Invalid value") <<"Invalid centrality value";
    }

  //----Vertex selection-------
  edm::Handle< reco::VertexCollection > vertices;
  iEvent.getByToken(vtxTags_ , vertices);
  if( !vertices->size() )
    {
      edm::LogWarning ("Missing Collection") <<"Invalid or empty vertex collection";
      return;
    }
  
  nvtx_ = 0;  
  xBestVtx_ = -999.;
  yBestVtx_ = -999.;
  zBestVtx_ = -999.;
  rhoBestVtx_ = -999.;
  xBestVtxError_ = -999.;
  yBestVtxError_ = -999.;
  zBestVtxError_ = -999.;

  LoopCMWVertices(iEvent, iSetup);
  
  //  if ( nvtx_ <= 0) return;
  if ( zBestVtx_ < zminVtx_ || zBestVtx_ > zmaxVtx_ ) return; 
  if ( rhoBestVtx_ > rhomaxVtx_ ) return;
  if ( centBin < centmin_*2 || centBin >= centmax_*2 ) return;

  hcent_bin -> Fill(centBin/2.);
  hZBestVtx -> Fill(zBestVtx_); 
  //------Track selection---------
  if(ifMcgen_)  LoopCMWTracksgen(iEvent, iSetup);
  else  LoopCMWTracks(iEvent, iSetup);


  //========Fill and push event===================
  evt_->run = iEvent.id().run();
  evt_->event = iEvent.id().event();
  evt_->zvtx = zBestVtx_;
  evt_->cent = (centBin/2.);

  evtVec_.push_back(*evt_);

  //reset evt container
  evt_->reset();
  

}

//-------------method called once each job just before ending the event loop--------------------------------------------

void
SayanCMW::endJob()
{

  std::cout<< "Start sorting the events!" << std::endl;
  std::sort(evtVec_.begin(),evtVec_.end());
  std::cout<< "Finish sorting the events!" << std::endl;

  std::cout<< "Total of " << evtVec_.size() << " events are selected!" << std::endl;
      
}
//==============================================================================================

void
SayanCMW::fillDescriptions(edm::ConfigurationDescriptions&  descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//=============================================================================

void
SayanCMW::LoopCMWVertices( const edm::Event& iEvent, const edm::EventSetup& iSetup)

{

  using namespace edm;
  using namespace std;

  edm::Handle< reco::VertexCollection > vertices;
  iEvent.getByToken(vtxTags_, vertices);
    

  if(!vertices->size())
    {
      std::cout<<"Invalid or empty vertex collection!"<<std::endl;
      return;
    }
  
  reco::VertexCollection recoVertices = *vertices;

  std::sort( recoVertices.begin(), recoVertices.end(),
	     [](const reco::Vertex &a, const reco::Vertex &b)
	     {
	       if ( a.tracksSize() == b.tracksSize() ) return a.chi2() < b.chi2();
	       return a.tracksSize() > b.tracksSize();
	     }
	     );

  //*****************Vertex loop started**************
  
  for ( reco::VertexCollection::const_iterator itVtx = recoVertices.begin(); itVtx != recoVertices.end(); ++itVtx ) // start event loop

    {
     
      if ( !itVtx->isFake() && itVtx->tracksSize() >= 2 ) // Drop fake vertex 
	{
	  //x,y,z vertex position
	  double xVtx = itVtx->x();
	  double yVtx = itVtx->y();
	  double zVtx = itVtx->z();

	  //x,y,z vertex position error
	  double xVtxError = itVtx->xError();
	  double yVtxError = itVtx->yError();
	  double zVtxError = itVtx->zError();

	  //radial vertex position
	  double rho = sqrt(xVtx*xVtx + yVtx*yVtx);
	  
	  ++nvtx_;

	  //Get the first vertex as the best one
	  if( itVtx == recoVertices.begin() )
	    {
	      xBestVtx_ = xVtx;
	      yBestVtx_ = yVtx;
	      zBestVtx_ = zVtx;
	      xBestVtxError_ = xVtxError;
	      yBestVtxError_ = yVtxError;
	      zBestVtxError_ = zVtxError;
	      rhoBestVtx_ = rho;

	    }
	}
    } 
}


void
SayanCMW::LoopCMWTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<std::vector<float>> mvaoutput;
  iEvent.getByToken(mvaSrc_, mvaoutput);

  //********* start track loop *********

  edm::Handle< reco::TrackCollection > tracks;
  iEvent.getByToken(trackTags_, tracks);

  if( !tracks->size() )
    {
      edm::LogWarning ("Missing Collection") <<"Invalid or empty collection!";
      return;
    }

  // Loop over tracks
  int it = 0;

  for ( reco::TrackCollection::const_iterator itTrk = tracks->begin(); itTrk != tracks->end(); ++itTrk )
    {
      
      // Select tracks based on proximity to best vertex
      math::XYZPoint bestvtx(xBestVtx_, yBestVtx_,zBestVtx_);
      double dzvtx = itTrk->dz(bestvtx);
      double dxyvtx = itTrk->dxy(bestvtx);
      double dzerror = sqrt(itTrk->dzError()*itTrk->dzError() + zBestVtxError_*zBestVtxError_);
      double dxyerror = sqrt(itTrk->d0Error()*itTrk->d0Error() + xBestVtxError_*yBestVtxError_);
      double pterror = itTrk->ptError();

      // Get eta, pt, phi and charge of the track
      double pt = itTrk->pt();
      double eta = itTrk->eta();
      double phi = itTrk->phi();
      int charge = itTrk->charge();

      //HI specific cuts
      double chi2n   = itTrk->normalizedChi2();
      double nlayers = itTrk->hitPattern().trackerLayersWithMeasurement();
      chi2n = chi2n/nlayers;
      int nHits = itTrk->numberOfValidHits();
      int algo  = itTrk->originalAlgo();
      float trkMVA = (*mvaoutput)[it];    
      it++;
      
      
      //selected tracks
      if( !itTrk->quality(reco::TrackBase::highPurity) ) continue;
      if( charge == 0 ) continue;
      if( fabs(dzvtx / dzerror) >= 3.0 ) continue;
      if( fabs(dxyvtx / dxyerror) >= 3.0  ) continue;
      if( fabs(pterror) / pt >= 0.1 ) continue;
      if( nHits < 11 ) continue;
      if( chi2n >= 0.18 ) continue;
      if( algo==6 && trkMVA < 0.98 ) continue;
      
      if( pt<0.5 || pt > 5.0 ) continue;
      if( eta < -2.4 || eta > 2.4 ) continue;


      //for photon conversion                                                                                                                        
      //const reco::HitPattern& hit_pattern = itTrk->hitPattern();
      //if(hit_pattern.pixelLayersWithMeasurement()==0)continue;

      if(algo==6 && trkMVA < 0.98) 
	{
	  std::cout<<"...................................................Particle found..............................................."<<std::endl;
	}

      //th2d_algomva -> Fill(algo, trkMVA);  //algo and MVA checks
      // th1d_algo -> Fill(algo);  //algo and MVA checks
      //th1d_mva -> Fill(trkMVA);  //algo and MVA checks

      int index = GetpTbin(pt);
      if(index == -1) continue;


      //double weight, weight_p, weight_n;
          
      double weight = TrkEff->getCorrection(pt, eta, phi, centBin);
      
      //if (charge > 0) { weight_p = TrkEff1->getCorrection(pt, eta, phi, centBin); }
      //else { weight_n = TrkEff2->getCorrection(pt, eta, phi, centBin); }
      
      //AssignpTbins(pt, eta, phi, weight, weight_p, weight_n, charge, index);
      
      //=========Filling histograms======================
      
      double x[4] = {eta, pt, phi, centBin};

      hpt->Fill(pt);
      heta->Fill(eta);
      heta_nbin->Fill(eta);
      hphi->Fill(phi);
      
      th2d_etaphi->Fill(eta, phi);
      
      tp1d_mpteta_nbin ->Fill(eta, pt);
      tp1d_mptphi ->Fill(phi, pt);
      tp1d_charge_eta ->Fill(eta, charge);
      tp1d_charge_phi ->Fill(phi, charge);
      
      hist4D->Fill(x);

      //corrected
      hpt_w ->Fill(pt, weight);
      heta_w ->Fill(eta, weight);
      heta_nbin_w ->Fill(eta, weight);
      hphi_w->Fill(phi, weight);

      th2d_etaphi_w->Fill(eta, phi, weight);

      tp1d_mpteta_nbin_w ->Fill(eta, pt, weight);
      tp1d_mptphi_w ->Fill(phi, pt, weight);
        
      //4D
      hist4D_w->Fill(x, weight);
    
      if (charge > 0)
	{
	  double xp[4] = {eta, pt, phi, centBin};
	  double weight_p = TrkEff1->getCorrection(pt, eta, phi, centBin);

	  hptP ->Fill(pt);
	  hetaP ->Fill(eta);
	  hetaP_nbin ->Fill(eta);
	  hphiP->Fill(phi);

	  th2d_etaphi_pos->Fill(eta, phi);

	  tp1d_mptetaP_nbin ->Fill(eta, pt);
	  tp1d_mptphiP ->Fill(phi, pt);
	  
	  hist4DP->Fill(xp);

	  //corrected
	  hptP_w ->Fill(pt, weight_p);
	  hetaP_w ->Fill(eta, weight_p);
	  hetaP_nbin_w ->Fill(eta, weight_p);
	  hphiP_w->Fill(phi, weight_p);
	  
	  th2d_etaphi_pos_w->Fill(eta, phi, weight_p);

	  tp1d_mptetaP_nbin_w ->Fill(eta, pt, weight_p);
	  tp1d_mptphiP_w ->Fill(phi, pt, weight_p);
	  
	  tp1d_charge_eta_w ->Fill(eta, charge, weight_p);
	  tp1d_charge_phi_w ->Fill(phi, charge, weight_p);
	  
	  hist4DP_w->Fill(xp, weight_p);
	}

      if (charge < 0)
	{
	  double xn[4] = {eta, pt, phi, centBin};
	  double weight_n = TrkEff2->getCorrection(pt, eta, phi, centBin);
	  
	  hptN ->Fill(pt);
	  hetaN ->Fill(eta);
	  hetaN_nbin ->Fill(eta);
	  hphiN->Fill(phi);

	  th2d_etaphi_neg->Fill(eta, phi);

	  tp1d_mptetaN_nbin ->Fill(eta, pt);
	  tp1d_mptphiN ->Fill(phi, pt);

	  hist4DN->Fill(xn);

	  //corrected
	  hptN_w ->Fill(pt, weight_n);
	  hetaN_w ->Fill(eta, weight_n);
	  hetaN_nbin_w ->Fill(eta, weight_n);
	  hphiN_w->Fill(phi, weight_n);

	  th2d_etaphi_neg_w->Fill(eta, phi, weight_n);

	  tp1d_mptetaN_nbin_w ->Fill(eta, pt, weight_n);
	  tp1d_mptphiN_w ->Fill(phi, pt, weight_n);
	  
	  tp1d_charge_eta_w ->Fill(eta, charge, weight_n);
	  tp1d_charge_phi_w ->Fill(phi, charge, weight_n);
	
	  hist4DN_w->Fill(xn, weight_n);
	}
    } 

}//end of LoopCMWTracks

// loop for Mcgen
void SayanCMW::LoopCMWTracksgen(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  edm::Handle< reco::GenParticleCollection > tracksgen;
  iEvent.getByToken(trackTagsgen_, tracksgen);
  if( !tracksgen->size() )
    {
      edm::LogWarning ("Missing MC Gen Collection") <<"Invalid or empty MC-Gen track collection!";
      return;
    }
  
  // Loop over tracks
  for( reco::GenParticleCollection::const_iterator itTrk = tracksgen->begin();itTrk != tracksgen->end(); ++itTrk )
    {
      
      //status and charge of the track
      if( itTrk->status() != 1 || itTrk->charge() == 0 ) continue;
      
      // Get eta, pt, phi of the track
      double eta      = itTrk->eta();
      double pt       = itTrk->pt();
      double phi      = itTrk->phi();
      int charge      = itTrk->charge();

      float weight = 1.0, weight_p = 1.0, weight_n = 1.0;
      if( pt < ptmin_ || pt > ptmax_ ) continue;
      if( eta < etamin_ || eta > etamax_ ) continue;

      int index = GetpTbin(pt);
      if(index == -1) continue;

      hpt->Fill(pt);
      heta->Fill(eta);
      heta_nbin->Fill(eta);
      hphi->Fill(phi);
      
      tp1d_mpteta_nbin ->Fill(eta, pt);
      tp1d_mptphi ->Fill(phi, pt);
      tp1d_charge_eta ->Fill(eta, charge);
      tp1d_charge_phi ->Fill(phi, charge);

      double x[4] = {eta, pt, phi, centBin};
      hist4D->Fill(x);

      if (charge > 0)
	{
	  double xp[4] = {eta, pt, phi, centBin};

	  hptP ->Fill(pt);
	  hetaP ->Fill(eta);
	  hetaP_nbin ->Fill(eta);
	  hphiP->Fill(phi);

	  tp1d_mptetaP_nbin ->Fill(eta, pt);
	  tp1d_mptphiP ->Fill(phi, pt);

	  hist4DP->Fill(xp);
	}
	    
      if (charge < 0)
	{
	  double xn[4] = {eta, pt, phi, centBin};

	  hptN ->Fill(pt);
	  hetaN ->Fill(eta);
	  hetaN_nbin ->Fill(eta);
	  hphiN->Fill(phi);

	  tp1d_mptetaN_nbin ->Fill(eta, pt);
	  tp1d_mptphiN ->Fill(phi, pt);

	  hist4DN->Fill(xn);
	} 
	
      // AssignpTbins(pt, eta, phi, weight, weight_p, weight_n, charge, index);

    }
}



//=============================== Getptbin================================
int SayanCMW::GetpTbin(double pt)
{
  int idx = -1;
  
  for(unsigned int i = 0; i<pTmin_.size(); ++i)
    {
      if( pt >= pTmin_[i] && pt <= pTmax_[i]) idx = i;
    }
  
  return idx;
}




//=========================== AssignpTbins =================================
  
void
SayanCMW::AssignpTbins(double pt, double eta, double phi, double weight, double weight_p, double weight_n, int charge, int idx)
{

  TLorentzVector pvector;
  pvector.SetPtEtaPhiM(pt, eta, phi, 0.140);
    
  (evt_->pVect[idx]).push_back(pvector);
  (evt_->chgVect[idx]).push_back(charge);
  (evt_->weightVect[idx]).push_back(weight);
  
  if ( charge > 0)
    {
      (evt_->pVect_trg[idx]).push_back(pvector);
      (evt_->chgVect_trg[idx]).push_back(charge);
      (evt_->weightVect_trg[idx]).push_back(weight_p);
    }
  else
    {
      (evt_->pVect_ass[idx]).push_back(pvector);
      (evt_->chgVect_ass[idx]).push_back(charge);
      (evt_->weightVect_ass[idx]).push_back(weight_n);
    }
}

DEFINE_FWK_MODULE(SayanCMW);
