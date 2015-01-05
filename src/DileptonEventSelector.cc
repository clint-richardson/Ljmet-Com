// -*- C++ -*-
//
// FWLite PAT analyzer-selector for dilepton analyses
//
// Adapted from StopEventSelector
// Aram Avetisyan, September 2012
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
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFMuonSelector.h"
#include "PhysicsTools/SelectorUtils/interface/PVSelector.h"
#include "LJMet/Com/interface/BaseEventSelector.h"
#include "LJMet/Com/interface/LjmetFactory.h"



using trigger::TriggerObject;



class DileptonEventSelector : public BaseEventSelector {
    
public:
    
    
    DileptonEventSelector();
    ~DileptonEventSelector();
    
    
    // executes before loop over events
    virtual void BeginJob(std::map<std::string, edm::ParameterSet const > par);
    
    // main method where the cuts are applied
    virtual bool operator()( edm::Event const & event, pat::strbitset & ret);
    
    // executes after loop over events
    virtual void EndJob(){}
    
    
    virtual void AnalyzeEvent( edm::Event const & event, LjmetEventContent & ec );
    
    
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




bool DileptonEventSelector::operator()( edm::Event const & event, pat::strbitset & ret){
    
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
        if ( considerCut("HBHE noise and scraping filter") ) {
            
            passCut(ret, "HBHE noise and scraping filter"); // PV cuts total
            
        } // end of PV cuts
        
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
        
        if ( mbPar["electron_cuts"] ) {
            
            //get electrons
            event.getByLabel( mtPar["electron_collection"], mhElectrons );
            
            mvSelElectrons.clear();
            
            for ( std::vector<pat::Electron>::const_iterator _iel = mhElectrons->begin(); _iel != mhElectrons->end(); _iel++){
                
                bool pass = false;
                
                while(1){
                    if (not _iel->gsfTrack().isNonnull() or not _iel->gsfTrack().isAvailable()) break;
                    
                    // electron Et cut
                    if (_iel->pt()>mdPar["electron_minpt"]){ }
                    else break;
                    
                    // electron eta cut
                    if ( fabs(_iel->eta())<mdPar["electron_maxeta"] ){ }
                    else break;
                    
                    pass = true;
                    break;
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
        int _n_jets = 0;
        //int njetsPF = 0;
        
        mvSelJets.clear();
        mvAllJets.clear();
        
        event.getByLabel( mtPar["jet_collection"], mhJets );
        for (std::vector<pat::Jet>::const_iterator _ijet = mhJets->begin();
             _ijet != mhJets->end(); ++_ijet){
            
            retJet.set(false);
            // jet cuts
            if ( (*jetSel_)( *_ijet, retJet ) ){ 
                mvAllJets.push_back(edm::Ptr<pat::Jet>(mhJets, _n_jets)); 
                if (( _ijet->pt()>mdPar["jet_minpt"] ) && ( fabs(_ijet->eta())<mdPar["jet_maxeta"] )){ 
                    ++_n_good_jets;
                    mvSelJets.push_back(edm::Ptr<pat::Jet>(mhJets, _n_jets)); 
                }
            }
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




void DileptonEventSelector::AnalyzeEvent( edm::Event const & event,
                                         LjmetEventContent & ec ){
    //
    // Compute analysis-specific quantities in the event,
    // return via event content
    //
    
    //   ec.SetValue("pi", 3.14);
    
    
    return;
}


#endif
