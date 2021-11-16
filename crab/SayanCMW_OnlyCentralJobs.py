if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    from CRABClient.UserUtilities import config # , getUsernameFromSiteDB
    config = config()
  
    config.General.workArea = 'PbPb2018_QA_q3_etaphi2d_cent3040_efffakesecmul_thnsparse4d_pt0p5_5p0_Nov10_2021_0' 
    config.General.transferOutputs = True
    config.General.transferLogs = False
    config.JobType.pluginName = 'Analysis'
    config.Data.unitsPerJob = 15          #LumiBased  #40 is good
    config.Data.totalUnits = -1
    config.Data.inputDBS = 'global'
    config.Data.splitting = 'LumiBased'
    config.Data.outLFNDirBase = '/store/user/sayan/'
    config.Data.publication = False
    config.Site.storageSite = 'T3_CH_CERNBOX'
    
    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    ################################0#############################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    ###############################
    #    Standard analysis        #
    ###############################
    
    ##### Central events 0-10% ####
    ##config.Data.unitsPerJob = 46
   
    config.General.requestName = 'PbPb2018_QA_q3_etaphi2d_cent3040_efffakesecmul_thnsparse4d_pt0p5_5p0_Nov10_2021_0' 
    config.JobType.psetName = '../cfg/SayanCMW_central_cfg.py'
    config.Data.inputDataset = '/HIMinimumBias11/HIRun2018A-04Apr2019-v1/AOD'
    config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/PromptReco/Cert_326381-327564_HI_PromptReco_Collisions18_JSON.txt'
    config.Data.outputDatasetTag = 'PbPb2018_QA_q3_etaphi2d_cent3040_efffakesecmul_thnsparse4d_pt0p5_5p0_Nov10_2021_0'
    submit(config)
            

    
    
    


