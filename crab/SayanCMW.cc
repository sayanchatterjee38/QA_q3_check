// -*- C++ -*-
//   
//         Author: Sayan Chatterjee
//         Class: SayanCMW 
//         Checked version 24/06/2020
//         Last Modified on 27/10/2021

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
 
  const Int_t  nptbin              = 26;
  Double_t  ptbining[nptbin+1]     = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.10, 1.15, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.0, 2.5, 3.0, 4.0, 5.0};

  const Int_t netaBin = 18;
  Double_t etabining[netaBin+1] = {-2.4, -2.0, -1.6, -1.4, -1.3, -1.2, -1.0, -0.8, -0.4, 0.0, 0.4, 0.8, 1.0, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4};
   
  //  const Int_t netaBin = 12;
  //Double_t etabining[netaBin+1] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};

  hpt          = fGlobalHist.make<TH1F>("hpt", "", 45, 0.5, 5.0 );
  hpt_nbin     = fGlobalHist.make<TH1F>("hpt_nbin", "", nptbin, ptbining );
  heta         = fGlobalHist.make<TH1F>("heta", "", 48, -2.4, 2.4 );
  heta_nbin    = fGlobalHist.make<TH1F>("heta_nbin", "", netaBin, etabining );
  hptP         = fGlobalHist.make<TH1F>("hptP", "", 45, 0.5, 5.0 );
  hptP_nbin    = fGlobalHist.make<TH1F>("hptP_nbin", "", nptbin, ptbining );
  hetaP        = fGlobalHist.make<TH1F>("hetaP", "", 48, -2.4, 2.4 );
  hetaP_nbin   = fGlobalHist.make<TH1F>("hetaP_nbin", "", netaBin, etabining );
  hptN         = fGlobalHist.make<TH1F>("hptN", "", 45, 0.5, 5.0 );
  hptN_nbin    = fGlobalHist.make<TH1F>("hptN_nbin", "", nptbin, ptbining );
  hetaN        = fGlobalHist.make<TH1F>("hetaN", "", 48, -2.4, 2.4 );
  hetaN_nbin   = fGlobalHist.make<TH1F>("hetaN_nbin", "", netaBin, etabining );
  hphi         = fGlobalHist.make<TH1F>("hphi", "", 64, -3.2, 3.2 );
  hZBestVtx    = fGlobalHist.make<TH1F>("hZvtx", "", 600, -30.0, 30.0);
  hcent_bin    = fGlobalHist.make<TH1I>("hcent", "", 100, 0, 100);
  
  th2d_etaphi    = fGlobalHist.make<TH2D>("th2d_etaphi", "", netaBin, etabining, 64, -3.2, 3.2 );
  th2d_etaphi_pos    = fGlobalHist.make<TH2D>("th2d_etaphi_pos", "", netaBin, etabining, 64, -3.2, 3.2 );
  th2d_etaphi_neg    = fGlobalHist.make<TH2D>("th2d_etaphi_neg", "", netaBin, etabining, 64, -3.2, 3.2 );

  th2d_algomva    = fGlobalHist.make<TH2F>("th2d_algomva", "", 40, -20, 20, 3000, -1.5, 1.5);
  th1d_algo    = fGlobalHist.make<TH1D>("th1d_algo", "", 40, -20, 20);
  th1d_mva    = fGlobalHist.make<TH1D>("th1d_mva", "", 300, -1.5, 1.5);

  tp1d_mpteta   = fGlobalHist.make<TProfile>("tp1d_mpt", "", 48, -2.4, 2.4, -1e10, 1e10);
  tp1d_mpteta_0p3   = fGlobalHist.make<TProfile>("tp1d_mpt_0p3", "", 16, -2.4, 2.4, -1e10, 1e10);
  tp1d_mpteta_nbin   = fGlobalHist.make<TProfile>("tp1d_mpt_nbin", "", netaBin, etabining, -1e10, 1e10);
  tp1d_mptetaP   = fGlobalHist.make<TProfile>("tp1d_mptP", "", 48, -2.4, 2.4, -1e10, 1e10);
  tp1d_mptetaP_0p3   = fGlobalHist.make<TProfile>("tp1d_mptP_0p3", "", 16, -2.4, 2.4, -1e10, 1e10);
  tp1d_mptetaP_nbin   = fGlobalHist.make<TProfile>("tp1d_mptP_nbin", "", netaBin, etabining, -1e10, 1e10);
  tp1d_mptetaN   = fGlobalHist.make<TProfile>("tp1d_mptN", "", 48, -2.4, 2.4, -1e10, 1e10);
  tp1d_mptetaN_0p3   = fGlobalHist.make<TProfile>("tp1d_mptN_0p3", "", 16, -2.4, 2.4, -1e10, 1e10);
  tp1d_mptetaN_nbin   = fGlobalHist.make<TProfile>("tp1d_mptN_nbin", "", netaBin, etabining, -1e10, 1e10);

  tp1d_charge_eta   = fGlobalHist.make<TProfile>("tp1d_charge_eta", "", netaBin, etabining, -1e10, 1e10);
  tp1d_charge_etaP   = fGlobalHist.make<TProfile>("tp1d_charge_etaP", "", netaBin, etabining, -1e10, 1e10);
  tp1d_charge_etaN   = fGlobalHist.make<TProfile>("tp1d_charge_etaN", "", netaBin, etabining, -1e10, 1e10);

  TFileDirectory f_effw = fs->mkdir("Effeciency_weight");

  hpt_w          = f_effw.make<TH1F>("hpt_w", "", 45, 0.5, 5.0 );
  hpt_nbin_w     = f_effw.make<TH1F>("hpt_nbin_w", "", nptbin, ptbining );
  heta_w         = f_effw.make<TH1F>("heta_w", "", 48, -2.4, 2.4 );
  heta_nbin_w    = f_effw.make<TH1F>("heta_nbin_w", "", netaBin, etabining );
  hptP_w         = f_effw.make<TH1F>("hptP_w", "", 45, 0.5, 5.0 );
  hptP_nbin_w    = f_effw.make<TH1F>("hptP_nbin_w", "", nptbin, ptbining );
  hetaP_w        = f_effw.make<TH1F>("hetaP_w", "", 48, -2.4, 2.4 );
  hetaP_nbin_w   = f_effw.make<TH1F>("hetaP_nbin_w", "", netaBin, etabining );
  hptN_w         = f_effw.make<TH1F>("hptN_w", "", 45, 0.5, 5.0 );
  hptN_nbin_w    = f_effw.make<TH1F>("hptN_nbin_w", "", nptbin, ptbining );
  hetaN_w        = f_effw.make<TH1F>("hetaN_w", "", 48, -2.4, 2.4 );
  hetaN_nbin_w   = f_effw.make<TH1F>("hetaN_nbin_w", "", netaBin, etabining );

  tp1d_mpteta_w   = f_effw.make<TProfile>("tp1d_mpt_w", "", 48, -2.4, 2.4, -1e10, 1e10);
  tp1d_mpteta_0p3_w   = f_effw.make<TProfile>("tp1d_mpt_0p3_w", "", 16, -2.4, 2.4, -1e10, 1e10);
  tp1d_mpteta_nbin_w   = f_effw.make<TProfile>("tp1d_mpt_nbin_w", "", netaBin, etabining, -1e10, 1e10);
  tp1d_mptetaP_w   = f_effw.make<TProfile>("tp1d_mptP_w", "", 48, -2.4, 2.4, -1e10, 1e10);
  tp1d_mptetaP_0p3_w   = f_effw.make<TProfile>("tp1d_mptP_0p3_w", "", 16, -2.4, 2.4, -1e10, 1e10);
  tp1d_mptetaP_nbin_w   = f_effw.make<TProfile>("tp1d_mptP_nbin_w", "", netaBin, etabining, -1e10, 1e10);
  tp1d_mptetaN_w   = f_effw.make<TProfile>("tp1d_mptN_w", "", 48, -2.4, 2.4, -1e10, 1e10);
  tp1d_mptetaN_0p3_w   = f_effw.make<TProfile>("tp1d_mptN_0p3_w", "", 16, -2.4, 2.4, -1e10, 1e10);
  tp1d_mptetaN_nbin_w   = f_effw.make<TProfile>("tp1d_mptN_nbin_w", "", netaBin, etabining, -1e10, 1e10);
 
  tp1d_charge_eta_w   = f_effw.make<TProfile>("tp1d_charge_eta_w", "", netaBin, etabining, -1e10, 1e10);
  tp1d_charge_etaP_w   = f_effw.make<TProfile>("tp1d_charge_etaP_w", "", netaBin, etabining, -1e10, 1e10);
  tp1d_charge_etaN_w   = f_effw.make<TProfile>("tp1d_charge_etaN_w", "", netaBin, etabining, -1e10, 1e10);

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

      th2d_algomva -> Fill(algo, trkMVA);  //algo and MVA checks
      th1d_algo -> Fill(algo);  //algo and MVA checks
      th1d_mva -> Fill(trkMVA);  //algo and MVA checks

      int index = GetpTbin(pt);
      if(index == -1) continue;


      double weight, weight_p, weight_n;
          
      weight = TrkEff->getCorrection(pt, eta, centBin);
      
      if (charge > 0) { weight_p = TrkEff1->getCorrection(pt, eta, centBin); }
      else { weight_n = TrkEff2->getCorrection(pt, eta, centBin); }
      
      AssignpTbins(pt, eta, phi, weight, weight_p, weight_n, charge, index);
      
      //=========Filling histograms======================
      
      hpt->Fill(pt);
      hpt_nbin->Fill(pt);
      heta->Fill(eta);
      heta_nbin->Fill(eta);
      hphi->Fill(phi);
      hpt_w ->Fill(pt, weight);
      hpt_nbin_w ->Fill(pt, weight);
      heta_w ->Fill(eta, weight);
      heta_nbin_w ->Fill(eta, weight);

      tp1d_mpteta ->Fill(eta, pt);
      tp1d_mpteta_0p3 ->Fill(eta, pt);
      tp1d_mpteta_nbin ->Fill(eta, pt);

      tp1d_charge_eta ->Fill(eta, charge);

      tp1d_mpteta_w ->Fill(eta, pt, weight);
      tp1d_mpteta_0p3_w ->Fill(eta, pt, weight);
      tp1d_mpteta_nbin_w ->Fill(eta, pt, weight);
      
      th2d_etaphi->Fill(eta, phi);
      
      if (charge > 0)
	{
	  hptP ->Fill(pt);
	  hptP_nbin ->Fill(pt);
	  hetaP ->Fill(eta);
	  hetaP_nbin ->Fill(eta);
	  tp1d_mptetaP ->Fill(eta, pt);
	  tp1d_mptetaP_0p3 ->Fill(eta, pt);
	  tp1d_mptetaP_nbin ->Fill(eta, pt);
	  hptP_w ->Fill(pt, weight_p);
	  hptP_nbin_w ->Fill(pt, weight_p);
	  hetaP_w ->Fill(eta, weight_p);
	  hetaP_nbin_w ->Fill(eta, weight_p);
	  tp1d_mptetaP_w ->Fill(eta, pt, weight_p);
	  tp1d_mptetaP_0p3_w ->Fill(eta, pt, weight_p);
	  tp1d_mptetaP_nbin_w ->Fill(eta, pt, weight_p);
	  
	  tp1d_charge_etaP ->Fill(eta, charge);
	  tp1d_charge_etaP_w ->Fill(eta, charge, weight_p);
	  tp1d_charge_eta_w ->Fill(eta, charge, weight_p);

	  th2d_etaphi_pos->Fill(eta, phi);
	}

      if (charge < 0)
	{
	  hptN ->Fill(pt);
	  hptN_nbin ->Fill(pt);
	  hetaN ->Fill(eta);
	  hetaN_nbin ->Fill(eta);
	  tp1d_mptetaN ->Fill(eta, pt);
	  tp1d_mptetaN_0p3 ->Fill(eta, pt);
	  tp1d_mptetaN_nbin ->Fill(eta, pt);
	  hptN_w ->Fill(pt, weight_n);
	  hptN_nbin_w ->Fill(pt, weight_n);
	  hetaN_w ->Fill(eta, weight_n);
	  hetaN_nbin_w ->Fill(eta, weight_n);
	  tp1d_mptetaN_w ->Fill(eta, pt, weight_n);
	  tp1d_mptetaN_0p3_w ->Fill(eta, pt, weight_n);
	  tp1d_mptetaN_nbin_w ->Fill(eta, pt, weight_n);

	  tp1d_charge_etaN ->Fill(eta, charge);
	  tp1d_charge_etaN_w ->Fill(eta, charge, weight_n);
	  tp1d_charge_eta_w ->Fill(eta, charge, weight_n);
	
	  th2d_etaphi_neg->Fill(eta, phi);
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
      hpt_nbin->Fill(pt);
      heta->Fill(eta);
      heta_nbin->Fill(eta);
      hphi->Fill(phi);
      
      tp1d_mpteta ->Fill(eta, pt);
      tp1d_mpteta_0p3 ->Fill(eta, pt);
      tp1d_mpteta_nbin ->Fill(eta, pt);
            
      tp1d_charge_eta ->Fill(eta, charge);

      if (charge > 0)
	{
	  hptP ->Fill(pt);
	  hptP_nbin->Fill(pt);
	  hetaP ->Fill(eta);
	  hetaP_nbin ->Fill(eta);
	  tp1d_mptetaP ->Fill(eta, pt);
	  tp1d_mptetaP_0p3 ->Fill(eta, pt);
	  tp1d_mptetaP_nbin ->Fill(eta, pt);
	}
	    
      if (charge < 0)
	{
	  hptN ->Fill(pt);
	  hptN_nbin->Fill(pt);
	  hetaN ->Fill(eta);
	  hetaN_nbin ->Fill(eta);
	  tp1d_mptetaN ->Fill(eta, pt);
	  tp1d_mptetaN_0p3 ->Fill(eta, pt);
	  tp1d_mptetaN_nbin ->Fill(eta, pt);
	} 
	
      AssignpTbins(pt, eta, phi, weight, weight_p, weight_n, charge, index);

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
