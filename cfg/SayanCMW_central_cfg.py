import FWCore.ParameterSet.Config as cms

process = cms.Process("SayanCMW")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("CondCore.CondDB.CondDB_cfi")
process.load('Configuration.EventContent.EventContent_cff')
# Add PbPb collision event selection
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesRecovery_cfi")
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('HeavyIonsAnalysis.Configuration.hfCoincFilter_cff')

process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


# __________________ General _________________

# Configure the logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500

# Configure the number of maximum event the analyser run on in interactive mode
# -1 == ALL
process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1) 
    #input = cms.untracked.int32(1000) 
    )


# __________________ I/O files _________________

# Define the input file to run on in interactive mode
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # '/store/user/clindsey/MinBias_Hydjet_Drum5F_2018_5p02TeV/RECODEBUG_20190625/190626_194626/0000/step2_RAW2DIGI_L1Reco_RECO_1.root' #for RecoDebug
        '/store/hidata/HIRun2018A/HIMinimumBias11/AOD/04Apr2019-v1/610000/54E9C8CE-C886-194D-8B7E-5A2017D83DF8.root'  #for data
        #'/store/hidata/HIRun2018A/HIMinimumBias11/AOD/04Apr2019-v1/00000/D9709818-F9EA-194F-A0AF-909F37F2D510.root'
 )
)

# Define output file name
import os
process.TFileService = cms.Service("TFileService",
          #fileName = cms.string('SayanCMW_PbPb2018_QA_q3_Oct29newEffv1_ptetamvacuts_noalgoparameters_efffakesecmul_pt2p0_5p0_cent3040_Nov02_2021_0.root')
                                   fileName = cms.string('SayanCMW_PbPb2018_QA_q3_etaphi2d_Oct29newEffv1_ptetamvacuts_noalgoparameters_efffakesecmul_pt0p5_5p0_cent3040_Nov08_2021_0.root')
          #fileName = cms.string('SayanCMW_MCDEBUGreco_QA_q3_pt0p5_8p0_cent0100_Sep23_2021_0.root')
)

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltSelect = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltSelect.HLTPaths = [
 "HLT_HIMinimumBias_*",
 ]
process.hltSelect.andOr = cms.bool(True)
process.hltSelect.throw = cms.bool(False)



# __________________ Detector conditions _________________
# Set the GlobalTag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '103X_dataRun2_Prompt_v2', '')   #for dat
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_dataRun2_v6', '')   #for data
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic_hi', '')  #for MC
#process.HiForest.GlobalTagLabel = process.GlobalTag.globaltag


# __________________ Event selection _________________
# Load centrality producer for centrality calculation
# Add PbPb centrality
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.hiCentrality.produceHFhits = False
process.hiCentrality.produceHFtowers = False
process.hiCentrality.produceEcalhits = False
process.hiCentrality.produceZDChits = False
process.hiCentrality.produceETmidRapidity = False
process.hiCentrality.producePixelhits = False
process.hiCentrality.produceTracks = False
process.hiCentrality.producePixelTracks = False
process.hiCentrality.reUseCentrality = True
process.hiCentrality.srcReUse = cms.InputTag("hiCentrality","","RECO")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.newCentralityBin = process.centralityBin.clone()
process.centralityBin.nonDefaultGlauberModel = cms.string("")
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
   cms.PSet(record = cms.string("HeavyIonRcd"),
            tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run2v1033p1x01_offline"),  #for data
            #tag = cms.string("CentralityTable_HFtowers200_HydjetDrum5F_v1032x02_mc"),
            connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
            label = cms.untracked.string("HFtowers")
   ),
])
process.p = cms.Path(process.hiCentrality * process.centralityBin)

process.ana_step = cms.Path(
process.offlinePrimaryVerticesRecovery )

# Load HI event selection modules

# __________________ Analyse Sequence _________________


# __________________ Analyse Sequence _________________                                                                                               
#====================                                                                                                                                 
# ZDC                                                                                                                                                 
#====================                                                                                                                                 

process.load('QWAna.QWZDC2018RecHit.QWZDC2018Producer_cfi')
process.load('QWAna.QWZDC2018RecHit.QWZDC2018RecHit_cfi')

process.newCent = cms.EDProducer("CentralityProducer",
                                 produceHFhits = cms.bool(False),
                                 produceHFtowers = cms.bool(False),
                                 produceEcalhits = cms.bool(False),
                                 produceZDChits = cms.bool(True),
                                 produceETmidRapidity = cms.bool(False),
                                 producePixelhits = cms.bool(False),
                                 produceTracks = cms.bool(False),
                                 producePixelTracks = cms.bool(False),

                                 reUseCentrality = cms.bool(True),
                                 srcZDChits = cms.InputTag("QWzdcreco"),
                                 lowGainZDC = cms.bool(True),

                                 doPixelCut = cms.bool(True),
                                 UseQuality = cms.bool(True),
                                 TrackQuality = cms.string('highPurity'),
                                 trackEtaCut = cms.double(2),
                                 trackPtCut = cms.double(1),
                                 hfEtaCut = cms.double(4), #hf above the absolute value of this cut is used                                         
  
                                 midRapidityRange = cms.double(1),
                                 srcReUse = cms.InputTag("hiCentrality")
                             )
#=====================================                                                                                                              
#   ZDC PU filter                                                                                                                                   
#   mode                                                                                                                                            
#   0: pos AND neg                                                                                                                                  
#   1: pos                                                                                                                                          
#   2: neg                                                                                                                                          
#   3: pos OR neg (default)                                                                                                                          
#======================================                                                                                                             
  
process.zdcPUfilter = cms.EDFilter('QWZDCPUFilter',
                                   src = cms.untracked.InputTag("newCent"),
                                   posPars = cms.untracked.vdouble(1098914.0, -174.),
                                   negPars = cms.untracked.vdouble(1322071.0, -193.),
                                   mode = cms.untracked.int32(3)
                               )

# Load you analyzer with initial configuration

process.load("Analyzers.SayanCMW.SayanCMW_cff")

#process.defaultAnalysis_05     = process.CPDC05.clone()
'''
process.defaultAnalysis_510    = process.CPDC510.clone()
process.defaultAnalysis_1020   = process.CPDC1020.clone()
'''
#process.defaultAnalysis_2030   = process.CPDC2030.clone()
process.defaultAnalysis_3040   = process.CPDC3040.clone()
#process.defaultAnalysis_4050   = process.CPDC4050.clone()
#process.defaultAnalysis_5070   = process.CPDC5070.clone()
#process.defaultAnalysis_70100  = process.CPDC70100.clone()

process.p = cms.Path(#process.offlinePrimaryVerticesRecovery *
                     process.hltSelect *
                     process.hfCoincFilter2Th4 *             # Requier HF coincidence with 3 GeV  
                     process.primaryVertexFilter *           # Clean up on vertices
                     process.clusterCompatibilityFilter *    # Clean up on pileup
                     process.centralityBin *                 # Compute centrality
                     process.zdcdigi *
                     process.QWzdcreco *
                     process.newCent *
                     process.zdcPUfilter *
#                     process.defaultAnalysis_05 )
                    # process.defaultAnalysis_510 *
                    # process.defaultAnalysis_1020 *
 #                    process.defaultAnalysis_2030 *
                     process.defaultAnalysis_3040 )
  #                   process.defaultAnalysis_4050)
#process.defaultAnalysis_5070 *
 #                    process.defaultAnalysis_70100 )


#process.schedule = cms.Schedule(process.p)

# peripheral pv recovery 
from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
process = MassReplaceInputTag(process,"offlinePrimaryVertices","offlinePrimaryVerticesRecovery")
process.offlinePrimaryVerticesRecovery.oldVertexLabel = "offlinePrimaryVertices"
