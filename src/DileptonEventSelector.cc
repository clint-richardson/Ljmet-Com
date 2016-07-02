// -*- C++ -*-
//
// FWLite PAT analyzer-selector for dilepton analyses
//
// Adapted from StopEventSelector
// Aram Avetisyan, September 2012 - TEST
//
#ifndef LJMet_Com_interface_DileptonEventSelector_h
#define LJMet_Com_interface_DileptonEventSelector_h



#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector3.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h"
//#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
//#include "FWCore/ParameterSet/interface/ProcessDesc.h"
//#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
//#include "PhysicsTools/SelectorUtils/interface/PFElectronSelector.h"
#include "LJMet/Com/interface/TopElectronSelector.h"
#include "LJMet/Com/interface/MiniIsolation.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFMuonSelector.h"
#include "PhysicsTools/SelectorUtils/interface/PVSelector.h"
#include "LJMet/Com/interface/BaseEventSelector.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"


using trigger::TriggerObject;



class DileptonEventSelector : public BaseEventSelector {
    
public:
    
    
    DileptonEventSelector();
    ~DileptonEventSelector();
    
    
    // executes before loop over events
    virtual void BeginJob(std::map<std::string, edm::ParameterSet const > par);
    
    // main method where the cuts are applied
    virtual bool operator()( edm::EventBase const & event, pat::strbitset & ret);
    
    // executes after loop over events
    virtual void EndJob(){}
    
    
    virtual void AnalyzeEvent( edm::EventBase const & event, LjmetEventContent & ec );
    
    
    boost::shared_ptr<PFJetIDSelectionFunctor> const & jetSel()        const { return jetSel_;}
    boost::shared_ptr<PVSelector>              const & pvSel()         const { return pvSel_;}
    
    
protected:
    
    std::string legend;
    bool bFirstEntry;
    
    
    boost::shared_ptr<PFJetIDSelectionFunctor> jetSel_;
    boost::shared_ptr<PVSelector>              pvSel_;
    
    edm::Handle<edm::TriggerResults >           mhEdmTriggerResults;
    edm::Handle<std::vector<pat::Jet> >         mhJets;
    edm::Handle<std::vector<pat::Muon> >        mhMuons;
    edm::Handle<std::vector<pat::Electron> >    mhElectrons;
    edm::Handle<std::vector<pat::MET> >         mhMet;
    edm::Handle<double>                         h_rho;
    edm::Handle<std::vector<reco::Vertex> >     h_primVtx;
    
    std::vector<edm::Ptr<reco::Vertex> >  good_pvs_;
    
private:
    
    void initialize(std::map<std::string, edm::ParameterSet const> par);
    
    
};


static int reg = LjmetFactory::GetInstance()->Register(new DileptonEventSelector(), "DileptonSelector");


DileptonEventSelector::DileptonEventSelector(){
}


DileptonEventSelector::~DileptonEventSelector(){
}


void DileptonEventSelector::BeginJob( std::map<std::string, edm::ParameterSet const> par){
    
    BaseEventSelector::BeginJob(par);
    
    std::string _key;
    
    _key = "pfJetIDSelector";
    if ( par.find(_key)!=par.end() ){
        jetSel_ = boost::shared_ptr<PFJetIDSelectionFunctor>( new PFJetIDSelectionFunctor(par[_key]) );
    }
    else {
        std::cout << mLegend << "jet ID selector not configured, exiting"
        << std::endl;
        std::exit(-1);
    }
    
    _key = "pvSelector";
    if ( par.find(_key)!=par.end() ){
        pvSel_ = boost::shared_ptr<PVSelector>( new PVSelector(par[_key]) );
    }
    else {
        std::cout << mLegend << "PV selector not configured, exiting"
        << std::endl;
        std::exit(-1);
    }
    
    
    _key = "event_selector";
    if ( par.find(_key)!=par.end() ){
        mbPar["trigger_cut"]              = par[_key].getParameter<bool>         ("trigger_cut");
        mbPar["dump_trigger"]             = par[_key].getParameter<bool>         ("dump_trigger");
        mvsPar["trigger_path_ee"]         = par[_key].getParameter<std::vector<std::string> >  ("trigger_path_ee");
        mvsPar["trigger_path_em"]         = par[_key].getParameter<std::vector<std::string> >  ("trigger_path_em");
        mvsPar["trigger_path_mm"]         = par[_key].getParameter<std::vector<std::string> >  ("trigger_path_mm");
        
        mbPar["pv_cut"]                   = par[_key].getParameter<bool>         ("pv_cut");
        mbPar["hbhe_cut"]                 = par[_key].getParameter<bool>         ("hbhe_cut");
        msPar["hbhe_cut_value"]           = par[_key].getParameter<std::string>  ("hbhe_cut_value");
        mbPar["hbheiso_cut"]              = par[_key].getParameter<bool>         ("hbheiso_cut");
	mbPar["cscHalo_cut"]              = par[_key].getParameter<bool>         ("cscHalo_cut");
	mbPar["eesc_cut"]                 = par[_key].getParameter<bool>         ("eesc_cut");
	mtPar["flag_tag"]                 = par[_key].getParameter<edm::InputTag>("flag_tag");
        mbPar["jet_cuts"]                 = par[_key].getParameter<bool>         ("jet_cuts");
        mdPar["jet_minpt"]                = par[_key].getParameter<double>       ("jet_minpt");
        mdPar["jet_maxeta"]               = par[_key].getParameter<double>       ("jet_maxeta");
        miPar["min_jet"]                  = par[_key].getParameter<int>          ("min_jet");
        miPar["max_jet"]                  = par[_key].getParameter<int>          ("max_jet");
        
        mbPar["muon_cuts"]                = par[_key].getParameter<bool>         ("muon_cuts");
        miPar["min_muon"]                 = par[_key].getParameter<int>          ("min_muon");
        miPar["max_muon"]                 = par[_key].getParameter<int>          ("max_muon");
        mdPar["muon_minpt"]               = par[_key].getParameter<double>       ("muon_minpt");
        mdPar["muon_maxeta"]              = par[_key].getParameter<double>       ("muon_maxeta");
        
        mbPar["electron_cuts"]            = par[_key].getParameter<bool>         ("electron_cuts");
        miPar["min_electron"]             = par[_key].getParameter<int>          ("min_electron");
        miPar["max_electron"]             = par[_key].getParameter<int>          ("max_electron");
        mdPar["electron_minpt"]           = par[_key].getParameter<double>       ("electron_minpt");
        mdPar["electron_maxeta"]          = par[_key].getParameter<double>       ("electron_maxeta");
        
        miPar["min_lepton"]               = par[_key].getParameter<int>          ("min_lepton");
        
        mbPar["met_cuts"]                 = par[_key].getParameter<bool>         ("met_cuts");
        
        mtPar["trigger_collection"]       = par[_key].getParameter<edm::InputTag>("trigger_collection");
        mtPar["pv_collection"]            = par[_key].getParameter<edm::InputTag>("pv_collection");
        mtPar["jet_collection"]           = par[_key].getParameter<edm::InputTag>("jet_collection");
        mtPar["muon_collection"]          = par[_key].getParameter<edm::InputTag>("muon_collection");
        mtPar["electron_collection"]      = par[_key].getParameter<edm::InputTag>("electron_collection");
        mtPar["met_collection"]           = par[_key].getParameter<edm::InputTag>("met_collection");
	mbPar["doLepJetCleaning"]         = par[_key].getParameter<bool>         ("doLepJetCleaning");
        mbPar["doNewJEC"]                 = par[_key].getParameter<bool>         ("doNewJEC");
	mbPar["isMc"]                     = par[_key].getParameter<bool>         ("isMc");

	//mva value
	mbPar["UseElMVA"]                 = par[_key].getParameter<bool>         ("UseElMVA");

    }
    else {
        std::cout << mLegend << "event selector not configured, exiting"
        << std::endl;
        std::exit(-1);
    }
    
    std::cout << mLegend << "initializing dilepton selection" << std::endl;
    
    
    bFirstEntry = true;
    
    push_back("No selection");
    set("No selection");
    
    push_back("Trigger");
    push_back("Primary vertex");
    push_back("HBHE noise and scraping filter");
    push_back("HBHE Iso noise filter");
    push_back("CSC Tight Halo filter");
    push_back("EE Bad SC filter");
    push_back("Min muon");
    push_back("Max muon");
    push_back("Min electron");
    push_back("Max electron");
    push_back("Min lepton");
    push_back("One jet or more");
    push_back("Two jets or more");
    push_back("Three jets or more");
    push_back("Min jet multiplicity");
    push_back("Max jet multiplicity");
    push_back("Min MET");
    push_back("All cuts");          // sanity check
    
    
    
    // TOP PAG sync selection v3
    
    set("Trigger", mbPar["trigger_cut"]);
    set("Primary vertex", mbPar["pv_cut"]);
    set("HBHE noise and scraping filter", mbPar["hbhe_cut"]);
    set("HBHE Iso noise filter", mbPar["hbheiso_cut"]); 
    set("CSC Tight Halo filter", mbPar["cscHalo_cut"]);
    set("EE Bad SC filter", mbPar["eesc_cut"]); 

    if (mbPar["muon_cuts"]){
        set("Min muon", miPar["min_muon"]);
        set("Max muon", miPar["max_muon"]);
    }
    else{
        set("Min muon", false);
        set("Max muon", false);
    }
    
    if (mbPar["electron_cuts"]){
        set("Min electron", miPar["min_electron"]);
        set("Max electron", miPar["max_electron"]);
    }
    else{
        set("Min electron", false);
        set("Max electron", false);
    }
    
    set("Min lepton", miPar["min_lepton"]);
    
    if (mbPar["jet_cuts"]){
        set("One jet or more", true);
        set("Two jets or more", true);
        set("Three jets or more", true);
        set("Min jet multiplicity", miPar["min_jet"]);
        set("Max jet multiplicity", miPar["max_jet"]);
    }
    else{
        set("One jet or more", false);
        set("Two jets or more", false);
        set("Three jets or more", false);
        set("Min jet multiplicity", false);
        set("Max jet multiplicity", false);
    }
    
    if (mbPar["met_cuts"]) set("Min MET", mdPar["min_met"]);
    
    set("All cuts", true);
    
} // initialize()




bool DileptonEventSelector::operator()( edm::EventBase const & event, pat::strbitset & ret){

  //load stuff for miniIso here
  //packed pf candidates and rho source needed miniIso
  edm::Handle<pat::PackedCandidateCollection> packedPFCands;
  edm::InputTag packedPFCandsLabel_("packedPFCandidates");
  event.getByLabel(packedPFCandsLabel_, packedPFCands);
  //rho isolation from susy recommendation
  edm::Handle<double> rhoJetsNC;
  event.getByLabel(edm::InputTag("fixedGridRhoFastjetCentralNeutral","") , rhoJetsNC);
  double myRhoJetsNC = *rhoJetsNC;
  
    
  pat::strbitset retJet       = jetSel_->getBitTemplate();
  
  while(1){ // standard infinite while loop trick to avoid nested ifs
        
    passCut(ret, "No selection");
    
    //
    //_____ Trigger cuts __________________________________
    //
    if ( considerCut("Trigger") ) {
      
      
      event.getByLabel( mtPar["trigger_collection"], mhEdmTriggerResults );
      //const edm::ParameterSetID ps = mhEdmTriggerResults->parameterSetID();
            const edm::TriggerNames trigNames = event.triggerNames(*mhEdmTriggerResults);
            
            unsigned int _tSize = mhEdmTriggerResults->size();
            
            // dump trigger names
            if (bFirstEntry && mbPar["dump_trigger"]){
                for (unsigned int i=0; i<_tSize; i++){
                    std::string trigName = trigNames.triggerName(i);
                    std::cout << i << "   " << trigName << std::endl;
                }
            }
            
            //Loop over each channel separately
            int passEE = 0;
            for (unsigned int ipath = 0; ipath < mvsPar["trigger_path_ee"].size(); ipath++){
                unsigned int _tIndex = trigNames.triggerIndex(mvsPar["trigger_path_ee"].at(ipath));
                if ( _tIndex<_tSize){
                    if (mhEdmTriggerResults->accept(_tIndex)){
                        passEE++;
                        break;
                    }
                }
            }
            
            int passEM = 0;
            for (unsigned int ipath = 0; ipath < mvsPar["trigger_path_em"].size(); ipath++){
                unsigned int _tIndex = trigNames.triggerIndex(mvsPar["trigger_path_em"].at(ipath));
                if ( _tIndex<_tSize){
                    if (mhEdmTriggerResults->accept(_tIndex)){
                        passEM++;
                        break;
                    }
                }
            }
            
            int passMM = 0;
            for (unsigned int ipath = 0; ipath < mvsPar["trigger_path_mm"].size(); ipath++){
                unsigned int _tIndex = trigNames.triggerIndex(mvsPar["trigger_path_mm"].at(ipath));
                if ( _tIndex<_tSize){
                    if (mhEdmTriggerResults->accept(_tIndex)){
                        passMM++;
                        break;
                    }
                }
            }
            
            mvSelTriggers.clear();
            mvSelTriggers.push_back(passEE);
            mvSelTriggers.push_back(passEM);
            mvSelTriggers.push_back(passMM);
            
            if ( ignoreCut("Trigger") || passEE + passEM + passMM > 0 ) passCut(ret, "Trigger");
            else break;
            
        } // end of trigger cuts
        
        
        
        //
        //_____ Primary vertex cuts __________________________________
        //
        if ( considerCut("Primary vertex") ) {
            
            if ( (*pvSel_)(event) ){
                passCut(ret, "Primary vertex"); // PV cuts total
            }
            
        } // end of PV cuts
        
        
        
        
        //
        //_____ HBHE noise and scraping filter________________________
        //
        if ( considerCut("HBHE noise and scraping filter") || considerCut("HBHE Iso noise filter") ) {
            
	  //set cut value considered
	  bool run1Cut=false;
	  bool run2LooseCut=false;
	  bool run2TightCut=false;
	  if(msPar["hbhe_cut_value"].find("Run1")!=std::string::npos){
	    run1Cut=true;
	  }
	  else if(msPar["hbhe_cut_value"].find("Run2Loose")!=std::string::npos){
	    run2LooseCut=true;
	  }
	  else if(msPar["hbhe_cut_value"].find("Run2Tight")!=std::string::npos){
	    run2TightCut=true;
	  }
	  else{
	    std::cout<<"HBHE cut not configured correctly!...exiting..."<<std::endl;
	    std::exit(-1);
	  }

	  if (mbPar["debug"]) std::cout<<"HBHE cuts..."<<std::endl;

	  //taken from http://cmslxr.fnal.gov/lxr/source/CommonTools/RecoAlgos/plugins/HBHENoiseFilterResultProducer.cc?v=CMSSW_7_4_1
	  //and http://cmslxr.fnal.gov/lxr/source/CommonTools/RecoAlgos/python/HBHENoiseFilterResultProducer_cfi.py?v=CMSSW_7_4_1

	  edm::InputTag noiselabel ("hcalnoise");
	  int minHPDHits = 17;
	  int minHPDNoOtherHits = 10;
	  //int minZeros = 10;
	  int minZeros = 999999.;
	  bool IgnoreTS4TS5ifJetInLowBVRegion = false;

	  //iso cuts
          int minNumIsolatedNoiseChannels = 10;
          double minIsolatedNoiseSumE = 50.0;
          double minIsolatedNoiseSumEt = 25.0;

	  // get the Noise summary object
	  edm::Handle<HcalNoiseSummary> summary_h;
	  event.getByLabel(noiselabel, summary_h);
	  if(!summary_h.isValid()) {
	    std::cout << "Could not find HcalNoiseSummary.\n";
	  }
	  const HcalNoiseSummary& summary(*summary_h);

	  bool goodJetFoundInLowBVRegion = false;
	  if (IgnoreTS4TS5ifJetInLowBVRegion) goodJetFoundInLowBVRegion = summary.goodJetFoundInLowBVRegion();

	  const bool failCommon = summary.maxHPDHits() >= minHPDHits || summary.maxHPDNoOtherHits() >= minHPDNoOtherHits || summary.maxZeros() >= minZeros;

	  //const bool failFull = failCommon || (summary.HasBadRBXTS4TS5() && !goodJetFoundInLowBVRegion);

	  const bool failRun1 = failCommon || (summary.HasBadRBXTS4TS5() && !goodJetFoundInLowBVRegion);

	  const bool failRun2Loose = failCommon || (summary.HasBadRBXRechitR45Loose() && !goodJetFoundInLowBVRegion);

	  const bool failRun2Tight = failCommon || (summary.HasBadRBXRechitR45Tight() && !goodJetFoundInLowBVRegion);

	  // Check isolation requirements
	  //const bool failIsolation = summary.numIsolatedNoiseChannels() >= minNumIsolatedNoiseChannels || summary.isolatedNoiseSumE() >= minIsolatedNoiseSumE || summary.isolatedNoiseSumEt() >= minIsolatedNoiseSumEt;
	  if(run1Cut){
	    if (!failRun1) passCut(ret, "HBHE noise and scraping filter"); // HBHE cuts total
	  }
	  else if(run2LooseCut){
	    if (!failRun2Loose) passCut(ret, "HBHE noise and scraping filter"); // HBHE cuts total
	  }
	  else if(run2TightCut){
	    if (!failRun2Tight) passCut(ret, "HBHE noise and scraping filter"); // HBHE cuts total
	  }
	  else {std::cout<<"No HBHE cut!"<<std::endl; passCut(ret, "HBHE noise and scraping filter");} // HBHE cuts total


	  // Check isolation requirements
	  const bool failIsolation = summary.numIsolatedNoiseChannels() >= minNumIsolatedNoiseChannels || summary.isolatedNoiseSumE() >= minIsolatedNoiseSumE || summary.isolatedNoiseSumEt() >= minIsolatedNoiseSumEt;
          if (considerCut("HBHE Iso noise filter")) {
	    if (!failIsolation) passCut(ret, "HBHE Iso noise filter"); // HBHE Iso cut
            else break;
          }


        } // end of hbhe cuts

	//_________CSC Halo Filter and EE Bad SuperCluster Filter_________

        if ( considerCut("CSC Tight Halo filter") || considerCut("EE Bad SC filter") ) {
	  edm::Handle<edm::TriggerResults > PatTriggerResults;
	  event.getByLabel( mtPar["flag_tag"], PatTriggerResults );
	  const edm::TriggerNames patTrigNames = event.triggerNames(*PatTriggerResults);
	  if ( considerCut("CSC Tight Halo filter") ) {
	        
	    bool cscpass = false;
    
	    for (unsigned int i=0; i<PatTriggerResults->size(); i++){
	      if (patTrigNames.triggerName(i) == "Flag_CSCTightHaloFilter") cscpass = PatTriggerResults->accept(patTrigNames.triggerIndex(patTrigNames.triggerName(i)));
	    }
    
	    if (cscpass) passCut(ret, "CSC Tight Halo filter"); // CSC cut
	    else break;
	  } // end of CSC cuts
    
	  if ( considerCut("EE Bad SC filter") ) {
	        
	    bool eescpass = false;
    
	    for (unsigned int i=0; i<PatTriggerResults->size(); i++){
	      if (patTrigNames.triggerName(i) == "Flag_eeBadScFilter") eescpass = PatTriggerResults->accept(patTrigNames.triggerIndex(patTrigNames.triggerName(i)));
	    }
    
	    if (eescpass) passCut(ret, "EE Bad SC filter"); // EE Bad SC cut
	    else break;
	  } // end of EE Bad SC cuts
        }

        //
        //_____ Muon cuts ________________________________
        //
        // loop over muons
        
        int _n_muons  = 0;
        int nSelMuons = 0;
        
        if ( mbPar["muon_cuts"] ) {
            
            //get muons
            event.getByLabel( mtPar["muon_collection"], mhMuons );
            
            mvSelMuons.clear();
            
            for (std::vector<pat::Muon>::const_iterator _imu = mhMuons->begin(); _imu != mhMuons->end(); _imu++){
                
                bool pass = false;
                
                //muon cuts
                while(1){
                    if (not _imu->globalTrack().isNonnull() or not _imu->globalTrack().isAvailable()) break;
                    if (not _imu->innerTrack().isNonnull()  or not _imu->innerTrack().isAvailable())  break;
                    
                    if (_imu->pt()>mdPar["muon_minpt"]){ }
                    else break;
                    
                    if ( fabs(_imu->eta())<mdPar["muon_maxeta"] ){ }
                    else break;
                    
                    pass = true; // success
                    break;
                }
                
                if ( pass ){
                    ++nSelMuons;
                    
                    // save every good muon
                    mvSelMuons.push_back( edm::Ptr<pat::Muon>( mhMuons, _n_muons) );
                }
                _n_muons++;
            } // end of the muon loop
            
            if( nSelMuons >= cut("Min muon", int()) || ignoreCut("Min muon") ) passCut(ret, "Min muon");
            else break;
            
            if( nSelMuons <= cut("Max muon", int()) || ignoreCut("Max muon") ) passCut(ret, "Max muon");
            else break;
            
        } // end of muon cuts
        
        //
        //_____ Electron cuts __________________________________
        //
        // loop over electrons
        
        int _n_electrons  = 0;
        int nSelElectrons = 0;
        //collection of electrons for cleaning
	std::vector<pat::Electron> electronsForCleaning;
	//read in pvs which we will need for dZ calculation:
	std::vector<reco::Vertex> goodPVs;
	edm::Handle<std::vector<reco::Vertex> > pvHandle;
	event.getByLabel(mtPar["pv_collection"], pvHandle);
	goodPVs = *(pvHandle.product());

        if ( mbPar["electron_cuts"] ) {

	  //get rho src
	  edm::InputTag rhoSrc_it("fixedGridRhoFastjetAll","");
	  edm::Handle<double> rhoHandle;
	  event.getByLabel(rhoSrc_it, rhoHandle);
	  double rhoIso = std::max(*(rhoHandle.product()), 0.0);

            
	  //get electrons
	  event.getByLabel( mtPar["electron_collection"], mhElectrons );
          
	  mvSelElectrons.clear();


          
	  for ( std::vector<pat::Electron>::const_iterator _iel = mhElectrons->begin(); _iel != mhElectrons->end(); _iel++){
	    
	    bool pass = false;
	    bool passLoose=false;
	    while(1){
	      if (not _iel->gsfTrack().isNonnull() or not _iel->gsfTrack().isAvailable()) break;
	      //skip if in barrel-endcap gap; doing it here means I never have to worry about it downstream since both electrons for analysis and those for cleaning are made here
	      if (_iel->isEBEEGap()) break;

	      //mva loose for cleaning
	      float mvaVal = mvaValue( *_iel,event);
	      pat::Electron* elptr = new pat::Electron(*_iel);
	      float miniIso = getPFMiniIsolation_EffectiveArea(packedPFCands, dynamic_cast<const reco::Candidate* > (elptr), 0.05, 0.2, 10., false, false,myRhoJetsNC);
	      

	      if(_iel->pt() < 10) passLoose=false;
	      else if(miniIso > 0.4) passLoose=false;
	      else{
		if(fabs(_iel->ecalDrivenMomentum().eta()) <0.8){
		  if(mvaVal>0.913286) passLoose = true;
		}
		else if(fabs(_iel->ecalDrivenMomentum().eta()) < 1.479){
		  if(mvaVal>0.805013) passLoose = true;
		}
		else if(fabs(_iel->ecalDrivenMomentum().eta())<2.4){
		  if(mvaVal > 0.358969) passLoose=true;
		}
	      }
	      /* NOT USING CUT BASED LOOSE ANYMORE
	      //get effective area to do pu correction for iso
	      double AEff;
	      if( fabs(_iel->ecalDrivenMomentum().eta())<1.0) AEff=0.1752;
	      else if(fabs(_iel->ecalDrivenMomentum().eta())<1.479) AEff=0.1862;
	      else if(fabs(_iel->ecalDrivenMomentum().eta())<2.0) AEff=0.1411;
	      else if(fabs(_iel->ecalDrivenMomentum().eta())<2.2) AEff=0.1534;
	      else if(fabs(_iel->ecalDrivenMomentum().eta())<2.3) AEff=0.1903;
	      else if(fabs(_iel->ecalDrivenMomentum().eta())<2.4) AEff=0.2243;
	      else if(fabs(_iel->ecalDrivenMomentum().eta())<2.5) AEff=0.2687;
	      
	      //calculate relIso
	      reco::GsfElectron::PflowIsolationVariables pfIso = _iel->pfIsolationVariables();
	      double chIso = pfIso.sumChargedHadronPt;
	      double nhIso = pfIso.sumNeutralHadronEt;
	      double phIso = pfIso.sumPhotonEt;
	      double PUIso = pfIso.sumPUPt;
	      double relIso = ( chIso + std::max(0.0, nhIso + phIso - rhoIso*AEff) ) / _iel->pt();
	      
	      //get d0 and dZ
	      double d0=_iel->dB();
	      double dZ;
	      if(goodPVs.size() > 0){
		dZ=_iel->gsfTrack()->dz(goodPVs.at(0).position());
	      } 
	      else {dZ=-999;}
	      //get 1/e -1/1p
	      float ooEmooP = 1.0/_iel->ecalEnergy() - _iel->eSuperClusterOverP()/_iel->ecalEnergy();



		if(fabs(_iel->ecalDrivenMomentum().eta()) <= 1.479){
		  if(_iel->full5x5_sigmaIetaIeta() >= 0.0103) {passLoose= false; }
		  else if(fabs(_iel->deltaEtaSuperClusterTrackAtVtx()) >= 0.0105)    {passLoose= false; }
		  else if(fabs(_iel->deltaPhiSuperClusterTrackAtVtx()) >= 0.115)    {passLoose= false; }
		  else if(_iel->hcalOverEcal() >= 0.104)         {passLoose= false; }
		  else if(relIso >= 0.0893)          {passLoose= false; }
		  else if(ooEmooP >= 0.102) {passLoose= false; }
		  else if(fabs(d0) >= 0.0261)      {passLoose= false; }
		  else if(fabs(dZ) >= 0.41)     {passLoose= false; }
		  else if(_iel->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 2)              {passLoose= false; }
		  else if(_iel->isGsfCtfScPixChargeConsistent() < 1)  {passLoose= false; }
		  else if(!_iel->passConversionVeto())        {passLoose= false; }
		  else passLoose=true;
		}
		
		//Endcap
		else{
		  if(_iel->full5x5_sigmaIetaIeta() >= 0.0301)  {passLoose= false; }
		  else if(fabs(_iel->deltaEtaSuperClusterTrackAtVtx()) >= 0.00814)    {passLoose= false; }
		  else if(fabs(_iel->deltaPhiSuperClusterTrackAtVtx()) >= 0.182)    {passLoose= false; }
		  else if(_iel->hcalOverEcal() >= 0.0897)        {passLoose= false; }
		  else if(relIso >= 0.121)        {passLoose= false; }
		  else if(ooEmooP >= 0.126) {passLoose= false; }
		  else if(fabs(d0) >= 0.118)       {passLoose= false; }
		  else if(fabs(dZ) >= 0.822)      {passLoose= false; }
		  else if(_iel->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1)              {passLoose= false; }
		  else if(_iel->isGsfCtfScPixChargeConsistent() < 1)  {passLoose= false; }
		  else if(!_iel->passConversionVeto())        {passLoose= false; }
		  else passLoose=true;
		}
		*/
	      // electron Et cut
	      if (_iel->pt()>mdPar["electron_minpt"]){ }
	      else break;
              
	      // electron eta cut
	      if ( fabs(_iel->eta())<mdPar["electron_maxeta"] ){ }
	      else break;
              
	      pass = true;
	      break;
	    }
	    //store loose electron for lepton-jet cleaning
	    if(passLoose){
	      electronsForCleaning.push_back(*_iel);
	    }
            
	    if ( pass ){
	      ++nSelElectrons;              
	      // save every good electron
	      mvSelElectrons.push_back( edm::Ptr<pat::Electron>( mhElectrons, _n_electrons) );
	      
	    }
	    _n_electrons++;
	  } // end of the electron loop
          
	  if( nSelElectrons >= cut("Min electron", int()) || ignoreCut("Min electron") ) passCut(ret, "Min electron");
	  else break;
          
	  if( nSelElectrons <= cut("Max electron", int()) || ignoreCut("Max electron") ) passCut(ret, "Max electron");
	  else break;
          
        } // end of electron cuts
        
        int nSelLeptons = nSelElectrons + nSelMuons;
        
        if( nSelLeptons >= cut("Min lepton", int()) || ignoreCut("Min lepton") ) passCut(ret, "Min lepton");
        else break;
        
        //     if (nSelLeptons < 2) std::cout<<"Too few leptons!"<<std::endl;
        
        //======================================================
        //
        // jet loop
        //
        //
        
        int _n_good_jets = 0;
	int _n_good_uncleaned_jets=0;
        int _n_jets = 0;
        //int njetsPF = 0;
        
        mvSelJets.clear();
	mvSelJetsCleaned.clear();
        mvAllJets.clear();
        
        event.getByLabel( mtPar["jet_collection"], mhJets );
        for (std::vector<pat::Jet>::const_iterator _ijet = mhJets->begin();
             _ijet != mhJets->end(); ++_ijet){
            
            retJet.set(false);
            // jet cuts for uncleaned jets
            if ( (*jetSel_)( *_ijet, retJet ) ){ 
                mvAllJets.push_back(edm::Ptr<pat::Jet>(mhJets, _n_jets)); 		
		//cut on corrected jet quantities
		TLorentzVector corJetP4 = correctJet(*_ijet,event);
                if (( corJetP4.Pt()>mdPar["jet_minpt"] ) && ( fabs(corJetP4.Eta())<mdPar["jet_maxeta"] )){ 
                    ++_n_good_uncleaned_jets;
                    mvSelJets.push_back(edm::Ptr<pat::Jet>(mhJets, _n_jets)); 
                }	       
            }
	    
	    //lepton jet cleaning
	    bool _cleaned = false;	    
	    pat::Jet cleanedJet = *_ijet;
	    TLorentzVector jetP4;
	    pat::Jet tmpJet = _ijet->correctedJet(0);

	    if ( mbPar["doLepJetCleaning"] ){

	      if (mbPar["debug"]) std::cout << "Checking Overlap" << std::endl;
	      //clean of muons
	      for(unsigned int ilep=0; ilep < mvSelMuons.size(); ilep++){
		bool looseMuon; //bool to loop over only loose muons can easily do here since loose muon id so easy to implement
		//get muon isolation
		//new definition of iso based on muon pog page
		const reco::MuonPFIsolation pfIsolationR04 = mvSelMuons[ilep]->pfIsolationR04();
		double chIso  = pfIsolationR04.sumChargedHadronPt;
		double nhIso  = pfIsolationR04.sumNeutralHadronEt;
		double gIso   = pfIsolationR04.sumPhotonEt;
		double puIso  = pfIsolationR04.sumPUPt;
		double relIso = (chIso + std::max(0.,nhIso + gIso - 0.5*puIso)) / mvSelMuons[ilep]->pt();

		//get miniIso
		pat::Muon* muptr = new pat::Muon(mvSelMuons[ilep]);
		float miniIso = getPFMiniIsolation_EffectiveArea(packedPFCands, dynamic_cast<const reco::Candidate* > (muptr), 0.05, 0.2, 10., false, false,myRhoJetsNC);

		if(miniIso > 0.4) looseMuon=false;
		else if(!(mvSelMuons[ilep]->isPFMuon())) looseMuon=false;
		else if(!( mvSelMuons[ilep]->isGlobalMuon() || mvSelMuons[ilep]->isTrackerMuon())) looseMuon=false;
		else looseMuon=true;
		if(!looseMuon) continue;
		if ( deltaR(mvSelMuons[ilep]->p4(),_ijet->p4()) < 0.4 ){
		  if (mbPar["debug"]) {
		    std::cout << "Jet Overlaps with the Muon... Cleaning jet..." << std::endl;
		    std::cout << "Lepton : pT = " << mvSelMuons[ilep]->pt() << " eta = " << mvSelMuons[ilep]->eta() << " phi = " << mvSelMuons[ilep]->phi() << std::endl;
		    std::cout << "      Raw Jet : pT = " << _ijet->pt() << " eta = " << _ijet->eta() << " phi = " << _ijet->phi() << std::endl;
		  }
		  //get the daughters of the muons
		  std::vector<reco::CandidatePtr> muDaughters;
		  for ( unsigned int isrc = 0; isrc < mvSelMuons[ilep]->numberOfSourceCandidatePtrs(); ++isrc ){
		    if (mvSelMuons[ilep]->sourceCandidatePtr(isrc).isAvailable()) muDaughters.push_back( mvSelMuons[ilep]->sourceCandidatePtr(isrc) );
		  }
		  const std::vector<edm::Ptr<reco::Candidate> > _ijet_consts = _ijet->daughterPtrVector();
		  for ( std::vector<edm::Ptr<reco::Candidate> >::const_iterator _i_const = _ijet_consts.begin(); _i_const != _ijet_consts.end(); ++_i_const){
		    /*if ( (*_i_const).key() == mvSelMuons[ilep]->originalObjectRef().key() ) {
		      cleanedJet.setP4( _ijet->p4() - mvSelMuons[ilep]->p4() );
		      //get the correction for the cleaned jet and apply it
		      jetP4 = correctJet(cleanedJet, event);
		      //annoying thing to convert our tlorentzstd::vector to root::math::lorentzstd::vector
		      ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double > > rlv;
		      rlv.SetXYZT(jetP4.X(),jetP4.Y(),jetP4.Z(),jetP4.T());
		      cleanedJet.setP4(rlv);
		      if (mbPar["debug"]) std::cout << "Corrected Jet : pT = " << jetP4.Pt() << " eta = " << jetP4.Eta() << " phi = " << jetP4.Phi() << std::endl;
		      _cleaned = true;
		      }*/ //switching to new lep-jet cleaning from dylan
		    for (unsigned int muI = 0; muI < muDaughters.size(); muI++) {
		      if ( (*_i_const).key() == muDaughters[muI].key() ) {
			

			tmpJet.setP4( tmpJet.p4() - muDaughters[muI]->p4() );
			
			if (mbPar["debug"]) std::cout << "Corrected Jet : pT = " << jetP4.Pt() << " eta = " << jetP4.Eta() << " phi = " << jetP4.Phi() << std::endl;
			_cleaned = true;
			muDaughters.erase( muDaughters.begin()+muI );
			break; //breaks out of mu daughters loop
		      }// end of if muon ref == jet ref
		    }//end loop over mu daughters
		   		 
		  }//end loop over jet constituents
		}//end check of muons being inside jet
	      }// end loop over muons for cleaning
	      	    
	      
	      //clean for electrons
	      for (unsigned int ilep=0; ilep <electronsForCleaning.size();ilep++){
		if ( deltaR(electronsForCleaning.at(ilep).p4(),_ijet->p4()) < 0.4 ){
		  //get the daughters of the electron
		  std::vector<reco::CandidatePtr> elDaughters;
		  for ( unsigned int isrc = 0; isrc < electronsForCleaning.at(ilep).numberOfSourceCandidatePtrs(); ++isrc ){
		    if (electronsForCleaning.at(ilep).sourceCandidatePtr(isrc).isAvailable()) elDaughters.push_back( electronsForCleaning.at(ilep).sourceCandidatePtr(isrc) );
		  }
		  if (mbPar["debug"]) {
		    std::cout << "Jet Overlaps with the Electron... Cleaning jet..." << std::endl;
		    std::cout << "Lepton : pT = " << electronsForCleaning.at(ilep).pt() << " eta = " << electronsForCleaning.at(ilep).eta() << " phi = " << electronsForCleaning.at(ilep).phi() << std::endl;
		    std::cout << "      Raw Jet : pT = " << _ijet->pt() << " eta = " << _ijet->eta() << " phi = " << _ijet->phi() << std::endl;
		  }
		  const std::vector<edm::Ptr<reco::Candidate> > _ijet_consts = _ijet->daughterPtrVector();
		  for ( std::vector<edm::Ptr<reco::Candidate> >::const_iterator _i_const = _ijet_consts.begin(); _i_const != _ijet_consts.end(); ++_i_const){
		    /*if ( (*_i_const).key() == electronsForCleaning.at(ilep).originalObjectRef().key() ) {
		      cleanedJet.setP4( _ijet->p4() - electronsForCleaning.at(ilep).p4() );
		      //get the correct 4std::vector
		      jetP4 = correctJet(cleanedJet, event);
		      //annoying thing to convert our tlorentzstd::vector to root::math::lorentzstd::vector
		      ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double > > rlv;
		      rlv.SetXYZT(jetP4.X(),jetP4.Y(),jetP4.Z(),jetP4.T());
		      cleanedJet.setP4(rlv);
		      if (mbPar["debug"]) std::cout << "Corrected Jet : pT = " << jetP4.Pt() << " eta = " << jetP4.Eta() << " phi = " << jetP4.Phi() << std::endl;
		      _cleaned = true;
		      }*/
		    for (unsigned int elI = 0; elI < elDaughters.size(); elI++) {
		      if ( (*_i_const).key() == elDaughters[elI].key() ) {
			tmpJet.setP4( tmpJet.p4() - elDaughters[elI]->p4() );

			if (mbPar["debug"]) std::cout << "Corrected Jet : pT = " << jetP4.Pt() << " eta = " << jetP4.Eta() << " phi = " << jetP4.Phi() << std::endl;
			_cleaned = true;
			elDaughters.erase( elDaughters.begin()+elI );
			break;
		      }
		    }
		  }
		}
	      }//end electron cleaning
	      


	      //if not cleaned just use first jet (remember if no cleaning then cleanedJet==*_ijet) to get corrected four std::vector and set the cleaned jet to have it
	      if (!_cleaned) {
		jetP4 = correctJet(cleanedJet, event);
		//annoying thing to convert our tlorentzstd::vector to root::math::lorentzstd::vector
		ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double > > rlv;
		rlv.SetXYZT(jetP4.X(),jetP4.Y(),jetP4.Z(),jetP4.T());
		cleanedJet.setP4(rlv);
	      }
	      else{
		//get the correct 4std::vector
		jetP4 = correctJet(tmpJet, event);
		//annoying thing to convert our tlorentzstd::vector to root::math::lorentzstd::vector
		ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double > > rlv;
		rlv.SetXYZT(jetP4.X(),jetP4.Y(),jetP4.Z(),jetP4.T());
		cleanedJet.setP4(rlv);     

	      }
	      
	      
	      //now do kinematic cuts and save them, cleaning and JEC could only mess up pfID so use original jet for it (though probably do nothing)
	      if ( (*jetSel_)( *_ijet, retJet ) ){ 
		
                if (( cleanedJet.pt()>mdPar["jet_minpt"] ) && ( fabs(cleanedJet.eta())<mdPar["jet_maxeta"] )){ 
		  ++_n_good_jets;
		  mvSelJetsCleaned.push_back(cleanedJet);
                }	       
	      }
	      
	    }//end lepton jet cleaning 
	    ++_n_jets;
	} // end of loop over jets
        
        if ( mbPar["jet_cuts"] ) {
            
            if ( ignoreCut("One jet or more") || _n_good_jets >= 1 ) passCut(ret, "One jet or more");
            else break; 
            
            if ( ignoreCut("Two jets or more") || _n_good_jets >= 2 ) passCut(ret, "Two jets or more");
            else break; 
            
            if ( ignoreCut("Three jets or more") || _n_good_jets >= 3 ) passCut(ret, "Three jets or more");
            else break; 
            
            if ( ignoreCut("Min jet multiplicity") || _n_good_jets >= cut("Min jet multiplicity",int()) ) passCut(ret, "Min jet multiplicity");
            else break; 
            
            if ( ignoreCut("Max jet multiplicity") || _n_good_jets <= cut("Max jet multiplicity",int()) ) passCut(ret, "Max jet multiplicity");
            else break; 
            
        } // end of jet cuts
        
        
        
        
        //
        //_____ MET cuts __________________________________
        //      
        event.getByLabel( mtPar["met_collection"], mhMet );
        mpMet = edm::Ptr<pat::MET>( mhMet, 0);
        
        if ( mbPar["met_cuts"] ) {
            
            
            // pfMet
            if ( mpMet.isNonnull() && mpMet.isAvailable() ) {
                
                pat::MET const & met = mhMet->at(0);
                if ( ignoreCut("Min MET") || met.et()>cut("Min MET", double()) ) passCut(ret, "Min MET");
            }
        } // end of MET cuts
        
        
        //
        //_____ Btagging selection _____________________
        //
        mvSelBtagJets.clear();
        bool _isTagged;
        
        for (std::vector<edm::Ptr<pat::Jet> >::const_iterator _ijet = mvSelJets.begin();
             _ijet != mvSelJets.end(); ++_ijet){
            
            _isTagged = isJetTagged(**_ijet, event);
            
            if (_isTagged) mvSelBtagJets.push_back(*_ijet); 
        }
        
        passCut(ret, "All cuts");
        break;
        
    } // end of while loop
    
    bFirstEntry = false;
    
    
    return (bool)ret;
    
    setIgnored(ret);
    
    return false;
}// end of operator()




void DileptonEventSelector::AnalyzeEvent( edm::EventBase const & event,
                                         LjmetEventContent & ec ){
    //
    // Compute analysis-specific quantities in the event,
    // return via event content
    //
    
    //   ec.SetValue("pi", 3.14);
    
    
    return;
}


#endif
