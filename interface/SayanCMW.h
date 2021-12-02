// -*- Header -*-
//
// Package:       Analyzers/SayanCMW
// Class:         SayanCMW
//
//
// Author:        Sayan Chatterjee
// Last modified: 27/10/2021




// system include files

#include <memory>

//CMSSW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "FWCore/Utilities/interface/InputTag.h"
//#include "Analyzers/SayanCMW/data/EFF/trackingEfficiency2018PbPb.h"
#include "Analyzers/SayanCMW/data/EFF/trackingEfficiency2018PbPb_newEffv1.h"
#include "Analyzers/SayanCMW/interface/DiHadronCorrelationEvt.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

//user include files

#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TTree.h"
#include "TVector3.h"
#include <string>
#include "TProfile.h"


class SayanCMW : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
 public:
     explicit SayanCMW(const edm::ParameterSet&);
     ~SayanCMW();

     static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
     
     void LoopCMWVertices(const edm::Event& iEvent, const edm::EventSetup& iSetup);
     void LoopCMWTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup);
     void LoopCMWTracksgen(const edm::Event& iEvent, const edm::EventSetup& iSetup);
     void AssignpTbins(double pt, double eta, double phi, double weight, double weight_p, double weight_n, int charge, int idx);
     int GetpTbin(double pt);
          
 private:
     virtual void beginJob() override;
     virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
     virtual void endJob() override;

     //****************member data*****************************************
     //## tracks ##
     // used to select what tracks to read from configuration file
     edm::EDGetTokenT<reco::TrackCollection> trackTags_;
     edm::EDGetTokenT<reco::GenParticleCollection> trackTagsgen_;
     edm::EDGetTokenT<std::vector<float>> mvaSrc_;

     //## vertex ##     
     // used to select what vertex to read from configuration file
     edm::EDGetTokenT<reco::VertexCollection> vtxTags_;

     //## centrality ##
     // used to select what centrality collection to read from configuration file
     edm::EDGetTokenT<reco::Centrality> cent_;

     // used to access centrality bins
     edm::EDGetTokenT<int> cent_bin_;
     // ## Event selection ##
     int centmin_;
     int centmax_;
     float centBin;

     //Event Id
     int ev_id;
     int Ntracks;
     
     // ## Vertex variables ##
     double xBestVtx_;
     double yBestVtx_;
     double zBestVtx_;
     double rhoBestVtx_;
     double xBestVtxError_;
     double yBestVtxError_;
     double zBestVtxError_;
     double zminVtx_; 
     double zmaxVtx_;
     double rhomaxVtx_;
     double nvtx_;

     //efficiency
     edm::InputTag fpt_;
     edm::InputTag fmb_;
     edm::InputTag fplus_;
     edm::InputTag fminus_;
     edm::InputTag fpix_;


     // ## track variables ##

     std::vector< double > pTmin_;
     std::vector< double > pTmax_;
     std::vector< double > pTmin_trg_;
     std::vector< double > pTmax_trg_;
     std::vector< double > pTmin_ass_;
     std::vector< double > pTmax_ass_;

     double ptmin_;
     double ptmax_;
     double etamin_;
     double etamax_;

     double dzdzerror_;
     double d0d0error_;
     double pterrorpt_;
     
     std::vector<int> algoParameters_;

     bool ifMcgen_;

     // ## Dihadron corr events ##
     DiHadronCorrelationEvt* evt_;
     std::vector< DiHadronCorrelationEvt > evtVec_;
    
     TrkEff2018PbPb* TrkEff;
     TrkEff2018PbPb* TrkEff1;
     TrkEff2018PbPb* TrkEff2;
     
     // ## histograms ##
     // QA_plots

     TH1F* hpt;
     TH1F* hpt_nbin;
     TH1F* heta;
     TH1F* heta_nbin;
     TH1F* hptP;
     TH1F* hptP_nbin;
     TH1F* hetaP;
     TH1F* hetaP_nbin;
     TH1F* hptN;
     TH1F* hptN_nbin;
     TH1F* hetaN;
     TH1F* hetaN_nbin;
     TH1F* hphi;
     TH1F* hZBestVtx;
     TH1I* hcent_bin;
     TH1I* hcentBin;

     TH2D* th2d_etaphi;
     TH2D* th2d_etaphi_pos;
     TH2D* th2d_etaphi_neg;

     TH2F* th2d_algomva;
     TH1D* th1d_algo;
     TH1D* th1d_mva;

     TProfile* tp1d_mpteta;
     TProfile* tp1d_mpteta_0p3;
     TProfile* tp1d_mpteta_nbin;
     TProfile* tp1d_mptetaP;
     TProfile* tp1d_mptetaP_0p3;
     TProfile* tp1d_mptetaP_nbin;
     TProfile* tp1d_mptetaN;
     TProfile* tp1d_mptetaN_0p3;
     TProfile* tp1d_mptetaN_nbin;

     TProfile* tp1d_charge_eta;
     TProfile* tp1d_charge_etaP;
     TProfile* tp1d_charge_etaN;

     //efficiency weight
     TH1F* hpt_w;
     TH1F* hpt_nbin_w;
     TH1F* heta_w;
     TH1F* heta_nbin_w;
     TH1F* hptP_w;
     TH1F* hptP_nbin_w;
     TH1F* hetaP_w;
     TH1F* hetaP_nbin_w;
     TH1F* hptN_w;
     TH1F* hptN_nbin_w;
     TH1F* hetaN_w;
     TH1F* hetaN_nbin_w;

     TProfile* tp1d_mpteta_w;
     TProfile* tp1d_mpteta_0p3_w;
     TProfile* tp1d_mpteta_nbin_w;
     TProfile* tp1d_mptetaP_w;
     TProfile* tp1d_mptetaP_0p3_w;
     TProfile* tp1d_mptetaP_nbin_w;
     TProfile* tp1d_mptetaN_w;
     TProfile* tp1d_mptetaN_0p3_w;
     TProfile* tp1d_mptetaN_nbin_w;

     TProfile* tp1d_charge_eta_w;
     TProfile* tp1d_charge_etaP_w;
     TProfile* tp1d_charge_etaN_w;

};
