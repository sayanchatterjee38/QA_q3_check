import FWCore.ParameterSet.Config as cms

defaultCPDC = cms.EDAnalyzer('SayanCMW', #Analyzer named: Correspond to the class name in 'plugin' folder

                             tracks = cms.InputTag("generalTracks"),
                             tracksgen = cms.InputTag("genParticles"),
                             mvaSrc = cms.InputTag('generalTracks','MVAValues'),
                             vertex = cms.InputTag('offlinePrimaryVertices'),
                             centrality = cms.InputTag("hiCentrality"),
                             centbin = cms.InputTag("centralityBin", "HFtowers"),
                             # Event selection
                             centmin = cms.untracked.int32(30),
                             centmax = cms.untracked.int32(40),

                             # Vertex selection
                             zminVtx = cms.untracked.double(-15.0),
                             zmaxVtx = cms.untracked.double(15.0),
                             rhomaxVtx = cms.untracked.double(0.2),
                             
                             # Efficiency 
                             fpt = cms.untracked.InputTag("2018PbPb_Efficiency_GeneralTracks_highPt.root"),
                             #fmb = cms.untracked.InputTag("2018PbPb_Efficiency_GeneralTracks_MB.root"),
                             #fplus = cms.untracked.InputTag("2018PbPb_Efficiency_GeneralTracks_MB_ChargePlus.root"),
                             #fminus = cms.untracked.InputTag("2018PbPb_Efficiency_GeneralTracks_MB_ChargeMinus.root"),
                             fmb = cms.untracked.InputTag("newEffv1_ptetamvacuts_noalgoparameters_alltrk_Oct29_2021.root"),
                             fplus = cms.untracked.InputTag("newEffv1_ptetamvacuts_noalgoparameters_postrk_Oct29_2021.root"),
                             fminus = cms.untracked.InputTag("newEffv1_ptetamvacuts_noalgoparameters_negtrk_Oct29_2021.root"),
                             
                             #fmb = cms.untracked.InputTag("newEffv1_ptetamvacuts_noalgoparameters_eta04_alltrk_Nov03_2021.root"),
                             #fplus = cms.untracked.InputTag("newEffv1_ptetamvacuts_noalgoparameters_eta04_postrk_Nov03_2021.root"),
                             #fminus = cms.untracked.InputTag("newEffv1_ptetamvacuts_noalgoparameters_eta04_negtrk_Nov03_2021.root"),
                             
                             fpix = cms.untracked.InputTag("2018PbPb_Efficiency_PixelTracks.root"),

                             # Track selection
                             pTminTrk = cms.untracked.vdouble(0.5),
                             pTmaxTrk = cms.untracked.vdouble(5.0),
                             pTminTrk_trg = cms.untracked.vdouble(0.5),
                             pTmaxTrk_trg = cms.untracked.vdouble(5.0),
                             pTminTrk_ass = cms.untracked.vdouble(0.5),
                             pTmaxTrk_ass = cms.untracked.vdouble(5.0),
                             ptmin = cms.untracked.double(0.5),
                             ptmax = cms.untracked.double(5.0),
                             etamin = cms.untracked.double(-2.4),
                             etamax = cms.untracked.double(2.4),
                             
                             algoParameters = cms.vint32(3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),

                             #ifMcgen = cms.untracked.bool(True), #for gen level study
                             ifMcgen = cms.untracked.bool(False), #for reco and data 
                             


)
