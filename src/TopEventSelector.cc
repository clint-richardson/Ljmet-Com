// -*- C++ -*-
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

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObject.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/Provenance/interface/ParameterSetID.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "LJMet/Com/interface/BaseEventSelector.h"
#include "LJMet/Com/interface/ElectronSelector.h"
#include "LJMet/Com/interface/JetIDSelectionFunctor.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/MetSelectionFunctor.h"
#include "LJMet/Com/interface/MuonSelectionFunctor.h"
#include "LJMet/Com/interface/PvObjectSelector.h"
#include "PhysicsTools/SelectorUtils/interface/EventSelector.h"
#include "PhysicsTools/SelectorUtils/interface/PVSelector.h"
#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "PhysicsTools/SelectorUtils/interface/strbitset.h"
//#include "LJMet/Com/interface/PVSelector.h"



using namespace std;
using trigger::TriggerObject;



class TopEventSelector : public BaseEventSelector {
    
public:
    
    
    enum Selection_t { Top, TopBTag, EWK_Zmumu, DiLepton, TprimeInc };
    enum Version_t { Off, test, SelV1, SelV2, SelV3, Tight, Loose, DiLepton_MuMu, DiLepton_ElEl, DiLepton_MuEl, Wprime_Mu, Wprime_Mu_forPrompt, Wprime_El, Tprime_Inc};
    
    
    
    TopEventSelector();
    ~TopEventSelector();
    
    
    // executes before loop over events
    // we parse config parameters here
    virtual void BeginJob(std::map<std::string, edm::ParameterSet const > par);
    
    // main method where the cuts are applied
    bool operator()( edm::Event const & event, pat::strbitset & ret);
    
    // executes after loop over events
    virtual void EndJob(){}
    
    
protected:
    
    void initialize(Version_t version);
    
    std::string legend;
    
    boost::shared_ptr<MuonSelectionFunctor>  muonSel_;
    boost::shared_ptr<MuonSelectionFunctor>  looseMuonSel_;
    boost::shared_ptr<ElectronSelector>      electronSel_;
    boost::shared_ptr<ElectronSelector>      looseElectronSel_;
    boost::shared_ptr<JetIDSelectionFunctor> jetSel_;
    boost::shared_ptr<MetSelectionFunctor>   metSel_;
    boost::shared_ptr<PVSelector>            pvSel_;
    boost::shared_ptr<PvObjectSelector>      pvObjSel_;
    
    edm::InputTag pvSrc_;
    edm::InputTag muonSrc_;
    edm::InputTag electronSrc_;
    edm::InputTag jetSrc_;
    edm::InputTag jetSrcPF_;
    //edm::InputTag jetSrcJPT_;
    edm::InputTag metSrc_;
    edm::InputTag metTCSrc_;
    edm::InputTag metPFSrc_;
    edm::InputTag metPfTypeISrc_;
    edm::InputTag metMuonSrc_;
    edm::InputTag metTypeIISrc_;
    edm::InputTag rhoCorrectionSrc_;
    edm::InputTag triggerSrc_;
    edm::InputTag BeamspotSrc_;
    
    edm::Handle<std::vector<reco::Vertex> >     h_primVtx;
    edm::Handle<vector<pat::MET> >              h_met_;
    edm::Handle<vector<pat::MET> >              h_tcMet_;
    edm::Handle<vector<pat::MET> >              h_pfMet_;
    edm::Handle<vector<reco::MET> >             h_pfTypeIMet_;
    edm::Handle<vector<reco::CaloMET> >         h_caloMet_;
    edm::Handle<vector<reco::CaloMET> >         h_caloMet2_;
    edm::Handle<vector<pat::Jet> >              h_jets_;
    edm::Handle<vector<pat::Jet> >              h_jetsPF_;
    edm::Handle<vector<pat::Jet> >              h_jetsJPT_;
    edm::Handle<vector<pat::Muon> >             h_muons_;
    edm::Handle<vector<pat::Electron> >         h_electrons_;
    edm::Handle<double>                         h_rho;
    edm::Handle<std::vector<pat::TriggerPath> > hltresults_;
    edm::Handle<edm::TriggerResults>            hlt_;
    edm::Handle<pat::TriggerEvent>              triggerEvent;
    edm::Handle<pat::TriggerEvent >             h_triggerEvent_;
    edm::Handle<edm::TriggerResults >           mhEdmTriggerResults;
    edm::Handle<reco::BeamSpot >                Beamspot_;
    
    
    edm::Ptr<pat::MET>      met_;
    edm::Ptr<pat::MET>      tcMet_;
    edm::Ptr<pat::MET>      pfMet_;
    edm::Ptr<reco::MET>     pfTypeIMet_;
    edm::Ptr<reco::CaloMET> caloMet_;
    edm::Ptr<reco::CaloMET> caloMet2_;
    edm::Ptr<pat::Jet>      jet0_;
    edm::Ptr<pat::Jet>      jet1_;
    edm::Ptr<pat::Jet>      jet2_;
    edm::Ptr<pat::Jet>      jet3_;
    edm::Ptr<pat::Muon>     muon0_;
    edm::Ptr<pat::Muon>     muon1_;
    edm::Ptr<pat::Electron> electron0_;
    edm::Ptr<pat::Electron> electron1_;
    
    double                           _dimuon_mass;
    
    std::vector<edm::Ptr<reco::Vertex> > good_pvs_;
    std::map<std::string,bool>       _hlt;
    std::vector<edm::Ptr<pat::Jet> > good_jets_;
    std::vector<edm::Ptr<pat::Jet> > good_jetsPF_;
    std::vector<edm::Ptr<pat::Jet> > all_jetsPF_;
    std::vector<edm::Ptr<pat::Jet> > good_jetsJPT_;
    std::vector<edm::Ptr<pat::Muon> > good_muons_;
    std::vector<edm::Ptr<pat::Muon> > loose_muons_;
    std::vector<edm::Ptr<pat::Electron> > good_electrons_;
    std::vector<edm::Ptr<pat::Electron> > loose_electrons_;
    
    bool muonMatchedToHlt_firstpass;
    
    Selection_t selection;
    Version_t version;
};



// register plugin
//static int reg = LjmetFactory::GetInstance()->Register(new TopEventSelector(), "TopSelector");



TopEventSelector::TopEventSelector(){
}



TopEventSelector::~TopEventSelector(){
}



void TopEventSelector::BeginJob(std::map<std::string, edm::ParameterSet const> par){
    //
    // Actions needed before event loop happens
    // The map contains all parameters from the python config
    //
    
    
    std::string selectionStr = par["event_selector"].getParameter<std::string>("selection");
    std::string versionStr   = par["event_selector"].getParameter<std::string>("version");
    
    
    pvSrc_            = par["event_selector"].getParameter<edm::InputTag>("pvSrc");
    muonSrc_          = par["event_selector"].getParameter<edm::InputTag>("muonSrc");
    electronSrc_      = par["event_selector"].getParameter<edm::InputTag>("electronSrc");
    jetSrc_           = par["event_selector"].getParameter<edm::InputTag>("jetSrc");
    jetSrcPF_         = par["event_selector"].getParameter<edm::InputTag>("jetSrcPF");
    //jetSrcJPT_        = par["event_selector"].getParameter<edm::InputTag>("jetSrcJPT");
    metSrc_           = par["event_selector"].getParameter<edm::InputTag>("metSrc");
    metTCSrc_         = par["event_selector"].getParameter<edm::InputTag>("metTCSrc");
    metPFSrc_         = par["event_selector"].getParameter<edm::InputTag>("metPFSrc");
    metPfTypeISrc_    = par["event_selector"].getParameter<edm::InputTag>("metPfTypeISrc");
    metMuonSrc_       = par["event_selector"].getParameter<edm::InputTag>("metMuonSrc");
    metTypeIISrc_     = par["event_selector"].getParameter<edm::InputTag>("metTypeIISrc");
    rhoCorrectionSrc_ = par["event_selector"].getParameter<edm::InputTag>("rhoCorrection");
    triggerSrc_       = par["event_selector"].getParameter<edm::InputTag>("triggerSrc");
    BeamspotSrc_      = par["event_selector"].getParameter<edm::InputTag>("BeamspotSrc");
    
    
    // legend for identifying output messages
    legend = "[TopEventSelector]: ";
    
    
    // create object selectors
    muonSel_          = boost::shared_ptr<MuonSelectionFunctor> ( new MuonSelectionFunctor (par["muonSelector"]         ) );
    looseMuonSel_     = boost::shared_ptr<MuonSelectionFunctor> ( new MuonSelectionFunctor (par["looseMuonSelector"]    ) );
    electronSel_      = boost::shared_ptr<ElectronSelector>     ( new ElectronSelector     (par["electronSelector"]     ) );
    looseElectronSel_ = boost::shared_ptr<ElectronSelector>     ( new ElectronSelector     (par["looseElectronSelector"]) );
    jetSel_           = boost::shared_ptr<JetIDSelectionFunctor>( new JetIDSelectionFunctor(par["jetIDSelector"]        ) );
    metSel_           = boost::shared_ptr<MetSelectionFunctor>  ( new MetSelectionFunctor  (par["metSelector"]          ) );
    pvSel_            = boost::shared_ptr<PVSelector>           ( new PVSelector           (par["pvSelector"]           ) );
    pvObjSel_         = boost::shared_ptr<PvObjectSelector>     ( new PvObjectSelector     (par["pvObjectSelector"]     ) );
    
    
    
    
    // this event selector is for top
    selection = Top;
    
    if      ( versionStr == "Off" )   version = Off;
    else if ( versionStr == "test" ) version = test;
    else if ( versionStr == "SelV1" ) version = SelV1;
    else if ( versionStr == "SelV2" ) version = SelV2;
    //Final Selection for mu+jets TTbar events
    else if ( versionStr == "SelV3" ) version = SelV3;
    else if ( versionStr == "Tight" ) version = Tight;
    else if ( versionStr == "Wprime_Mu" ) version = Wprime_Mu;
    else if ( versionStr == "Wprime_Mu_forPrompt" ) version = Wprime_Mu_forPrompt;
    else if ( versionStr == "Wprime_El" ) version = Wprime_El;
    else {
        throw cms::Exception("InvalidInput") << "Expect version to be one of Off, test, SelV1, SelV2, SelV3" << std::endl;
    }
    
    std::cout << legend << "initializing top semileptonic selection" << std::endl;
    std::cout << legend << "selection version: " << versionStr << std::endl;
    
    initialize( version );
}



void TopEventSelector::initialize( Version_t version){
    
    push_back("No selection");
    set("No selection");
    
    push_back("HLT_Mu15");
    push_back("HLT_Ele");
    push_back("Trigger names"); // activate to get a list of available trigger names
    push_back("Prompt Reco cut"); // activate to cut out May10 Lumi Sections
    push_back("Trigger cuts");
    
    push_back("Min number of PV");
    push_back("Max number of PV");
    push_back("Primary vertex cuts");
    
    push_back("Min mu multiplicity");
    push_back("Max mu multiplicity");
    push_back("Loose muon veto");
    push_back("Muon veto");
    push_back("Muon cuts");
    
    push_back("Min ele multiplicity");
    push_back("Max ele multiplicity");
    push_back("Loose ele veto");
    push_back("Electron veto");
    push_back("Electron cuts");
    
    
    push_back("One jet or more");
    push_back("Two jets or more");
    push_back("Three jets or more");
    push_back("Min jet multiplicity");
    push_back("Max jet multiplicity");
    push_back("Leading jet pt");
    push_back("Jet cuts");
    
    push_back("MET pt");
    push_back("MET cuts");
    
    push_back("No btag jets");
    push_back("One btag jet or more");
    push_back("Two btag jets or more");
    push_back("Three btag jets or more");
    push_back("Btagging");
    
    
    set("HLT_Ele", false);
    set("Prompt Reco cut", false);
    set("Muon veto", false);
    
    set("Min ele multiplicity", false);
    set("Max ele multiplicity", false);
    set("Loose ele veto", false);
    set("Electron veto", true);
    set("Electron cuts", true);
    
    set("Leading jet pt", false);
    
    set("No btag jets", false);
    set("One btag jet or more", false);
    set("Two btag jets or more", false);
    set("Three btag jets or more", false);
    set("Btagging", false);
    
    
    // test cuts, for debugging
    if (version == test){
        
        set("HLT_Mu15", true);
        set("Trigger names", false);
        set("Trigger cuts", true);
        
        set("Min number of PV", 1);
        set("Max number of PV", false);
        set("Primary vertex cuts", true);
        
        set("Min mu multiplicity", 1);
        set("Max mu multiplicity", 1);
        set("Loose muon veto", true);
        set("Muon cuts", true);
        
        set("Electron veto", true);
        set("Electron cuts", true);
        
        set("One jet or more", true);
        set("Two jets or more", true);
        set("Three jets or more", true);
        set("Min jet multiplicity", 4);
        set("Max jet multiplicity", false);
        set("Jet cuts", true);
        
        set("MET pt", false);
        set("MET cuts", false);
        
    }
    
    // TOP PAG sync selection v1
    else if (version == SelV1){
        
        set("HLT_Mu15", true);
        set("Trigger names", false);
        set("Trigger cuts", true);
        
        set("Min number of PV", 1);
        set("Max number of PV", false);
        set("Primary vertex cuts", true);
        
        set("Min mu multiplicity", 1);
        set("Max mu multiplicity", 1);
        set("Loose muon veto", true);
        set("Muon cuts", true);
        
        set("Electron veto");
        set("Electron cuts");
        
        set("One jet or more");
        set("Two jets or more");
        set("Three jets or more");
        set("Min jet multiplicity", 4);
        set("Max jet multiplicity", std::numeric_limits<int>::max());
        set("Jet cuts");
        
        set("MET pt", false);
        set("MET cuts", true);
    }
    // TOP PAG sync selection v3
    else if (version == SelV3){
        
        //set("HLT_Mu15", true);
        set("HLT_Mu15", true);
        set("Trigger names", false);
        set("Trigger cuts", true);
        
        set("Min number of PV", 1);
        set("Max number of PV", false);
        set("Primary vertex cuts", true);
        
        set("Min mu multiplicity", 1);
        set("Max mu multiplicity", 1);
        set("Loose muon veto", true);
        set("Muon cuts");
        
        set("Min ele multiplicity", false);
        set("Max ele multiplicity", false);
        set("Loose ele veto", false);
        set("Electron veto", true);
        set("Electron cuts", true);
        
        set("One jet or more", true);
        set("Two jets or more", true);
        set("Three jets or more", true);
        set("Min jet multiplicity", 4);
        set("Max jet multiplicity", false);
        set("Leading jet pt", false);
        set("Jet cuts");
        
        set("MET pt", false);
        set("MET cuts", true); // keep true for access to met collections
        
        set("No btag jets", false);
        set("One btag jet or more", false);
        set("Two btag jets or more", false);
        set("Three btag jets or more", false);
        set("Btagging", false);
        
    }
    // TOP PAG sync selection v3
    else if (version == Tight){
        
        set("HLT_Mu15", true);
        set("Trigger names", false);
        set("Trigger cuts");
        
        set("Min number of PV", 1);
        set("Max number of PV", false);
        set("Primary vertex cuts", true);
        
        set("Min mu multiplicity", 1);
        set("Max mu multiplicity", 1);
        set("Loose muon veto");
        set("Muon cuts");
        
        set("Electron veto");
        set("Electron cuts");
        
        set("One jet or more", true);
        set("Two jets or more", true);
        set("Three jets or more", true);
        set("Min jet multiplicity", 1);
        set("Max jet multiplicity", 4);
        set("Jet cuts");
        
        set("MET pt", false);
        set("MET cuts", true); // keep true for access to met collections
        
        set("No btag jets", false);
        set("One btag jet or more", false);
        set("Two btag jets or more", false);
        set("Three btag jets or more", false);
        set("Btagging", false);
        
        
    }
    
    // Wprime Muon Selection
    else if (version == Wprime_Mu){
        
        set("HLT_Mu15", true);
        set("Trigger names", false);
        set("Trigger cuts", true);
        
        set("Min number of PV", 1);
        set("Max number of PV", false);
        set("Primary vertex cuts", true);
        
        set("Min mu multiplicity", 1);
        set("Max mu multiplicity", 1);
        set("Loose muon veto", true);
        set("Muon cuts");
        
        set("Min ele multiplicity", false);
        set("Max ele multiplicity", false);
        set("Loose ele veto", false);
        set("Electron veto", true);
        set("Electron cuts", true);
        
        set("One jet or more", true);
        set("Two jets or more", true);     //change back to true
        //set("One jet or more", false);
        //set("Two jets or more", false);
        set("Three jets or more", false);
        set("Min jet multiplicity", 2);     //change back to 2
        //set("Min jet multiplicity", 0);
        set("Max jet multiplicity", false);
        //set("Leading jet pt", 100.0);
        set("Leading jet pt", false);
        set("Jet cuts");
        
        set("MET pt", true);
        set("MET cuts", true); // keep true for access to met collections
        
        set("No btag jets", false);
        set("One btag jet or more", false);
        set("Two btag jets or more", false);
        set("Three btag jets or more", false);
        set("Btagging", false);
        
    }
    
    // Wprime Muon Selection
    else if (version ==  Wprime_Mu_forPrompt){
        
        set("HLT_Mu15", true);
        set("Trigger names", false);
        set("Prompt Reco cut", false);
        set("Trigger cuts", true);
        
        set("Min number of PV", 1);
        set("Max number of PV", false);
        set("Primary vertex cuts", true);
        
        set("Min mu multiplicity", 1);
        set("Max mu multiplicity", 1);
        set("Loose muon veto", true);
        set("Muon cuts");
        
        set("Min ele multiplicity", false);
        set("Max ele multiplicity", false);
        set("Loose ele veto", false);
        set("Electron veto", true);
        set("Electron cuts", true);
        
        set("One jet or more", true);
        set("Two jets or more", true);
        set("Three jets or more", false);
        set("Min jet multiplicity", 2);
        set("Max jet multiplicity", false);
        set("Jet cuts");
        
        set("MET pt", false);
        set("MET cuts", true); // keep true for access to met collections
        
        set("No btag jets", false);
        set("One btag jet or more", false);
        set("Two btag jets or more", false);
        set("Three btag jets or more", false);
        set("Btagging", false);
        
    }
    
    
    
    // Wprime Electron Selection
    else if (version == Wprime_El){
        
        set("HLT_Mu15", false);
        set("HLT_Ele", true);
        set("Trigger names", false);
        set("Trigger cuts", true);
        
        set("Min number of PV", 1);
        set("Max number of PV", false);
        set("Primary vertex cuts", true);
        
        set("Min mu multiplicity", false);
        set("Max mu multiplicity", false);
        set("Loose muon veto", false);
        set("Muon veto", true);
        set("Muon cuts", true);
        
        set("Min ele multiplicity", 1);
        set("Max ele multiplicity", 1);
        set("Loose ele veto", true);
        set("Electron veto", false);
        set("Electron cuts", true);
        
        set("One jet or more", true);
        set("Two jets or more", true);
        set("Three jets or more", false);
        set("Min jet multiplicity", 2);
        set("Max jet multiplicity", false);
        set("Leading jet pt", 80.0);
        set("Jet cuts");
        
        set("MET pt", true);
        set("MET cuts", true); // keep true for access to met collections
        
        set("No btag jets", false);
        set("One btag jet or more", false);
        set("Two btag jets or more", false);
        set("Three btag jets or more", false);
        set("Btagging", false);
        
    }
    
    
    else if (version == Off){
        std::cout << "Selection version: Off" << std::endl;
        
        set("HLT_Mu15", false);
        set("Trigger names", false);
        set("Trigger cuts", false);
        
        set("Min number of PV", false);
        set("Max number of PV", false);
        set("Primary vertex cuts", false);
        
        set("Min mu multiplicity", false);
        set("Max mu multiplicity", false);
        set("Loose muon veto", false);
        set("Muon cuts", false);
        
        set("Electron veto", false);
        set("Electron cuts", false);
        
        set("One jet or more", false);
        set("Two jets or more", false);
        set("Three jets or more", false);
        set("Min jet multiplicity", false);
        set("Jet cuts");
        
        set("MET pt", false);
        set("MET cuts", true);
    }
    
    
} // initialize()




bool TopEventSelector::operator()( edm::Event const & event, pat::strbitset & ret){
    
    pat::strbitset retMuon      = muonSel_->getBitTemplate();
    pat::strbitset retLooseMuon = looseMuonSel_->getBitTemplate();
    pat::strbitset retElectron  = electronSel_->getBitTemplate();
    pat::strbitset retLooseElectron = looseElectronSel_->getBitTemplate();
    pat::strbitset retJet       = jetSel_->getBitTemplate();
    pat::strbitset retMet       = metSel_->getBitTemplate();
    pat::strbitset retPv        = pvObjSel_->getBitTemplate();
    
    
    
    
    while(1){ // standard infinite while loop trick to avoid nested ifs
        
        passCut(ret, "No selection");
        
        //
        //_____ Trigger cuts __________________________________
        //
        if ( considerCut("Trigger cuts") ) {
            
            
            edm::InputTag _triggerEventSrc("TriggerResults::HLT");
            event.getByLabel( _triggerEventSrc, mhEdmTriggerResults );
            const edm::ParameterSetID ps = mhEdmTriggerResults->parameterSetID();
            const edm::TriggerNames trigNames = event.triggerNames(*mhEdmTriggerResults);
            
            bool passTrig = false;
            bool passTrigEle = false;
            
            
            //================================================
            /*
             TRIGGERS FOR MUONS
             */
            //=================================================
            
            std::string pathName;
            bool skip_event = false;
            int run_number = event.id().run();
            
            // For MC
            
            if( run_number < 159999)
                pathName = "HLT_IsoMu17_eta2p1_TriCentralPFJet30_v2";
            
            else if ( run_number >= 190456 && run_number <= 190738 )
                pathName = "HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2";
            
            else if ( run_number >= 190782 && run_number <= 191411 )
                pathName = "HLT_IsoMu20_eta2p1_TriCentralPFJet30_v3";
            
            else if ( run_number >= 191695 && run_number <= 193621 )
                pathName = "HLT_IsoMu20_eta2p1_TriCentralPFJet30_v4";
            
            else if ( run_number >= 193834 && run_number <= 194115 )
                pathName = "HLT_IsoMu17_eta2p1_TriCentralPFJet30_v5";
            
            else{
                //cout << "Unknown run for HLTpath selection: " << run_number << endl;
                skip_event = true;
            }
            
            if (skip_event) {
                passTrig = false;
            }
            else{
                //passTrig = mhEdmTriggerResults->accept(trigNames.triggerIndex(pathName));
                
                for (unsigned int i=0; i<mhEdmTriggerResults->size(); i++){
                    std::string trigName = trigNames.triggerName(i);
                    if (trigName == pathName){
                        passTrig = mhEdmTriggerResults->accept(trigNames.triggerIndex(trigName));
                    }
                }
            }
            
            // if we are only dumping trigger names, no need to go further
            if ( considerCut("Trigger names") ) {
                
                for (unsigned int i=0; i<mhEdmTriggerResults->size(); i++)
                {
                    std::string trigName = trigNames.triggerName(i);
                    if (trigName == pathName){
                        cout << "trigg Name = "<< trigName;
                        bool fired = mhEdmTriggerResults->accept(trigNames.triggerIndex(trigName));
                        cout <<", FIRED = "<<fired<<endl;
                    }
                }
                cout << endl;
                return 0;
            }
            
            
            if ( ignoreCut("HLT_Mu15") || passTrig ) passCut(ret, "HLT_Mu15");
            else break;
            
            if ( ignoreCut("HLT_Ele") || passTrigEle ) passCut(ret, "HLT_Ele");
            else break;
            
            passCut(ret, "Trigger cuts"); // trigger cuts total
            
        } // end of trigger cuts
        
        
        //
        //_____ Primary vertex cuts __________________________________
        //
        if ( considerCut("Primary vertex cuts") ) {
            event.getByLabel(pvSrc_, h_primVtx);
            
            int nVtx = 0;
            int nSelVtx = 0;
            good_pvs_.clear();
            for ( vector<reco::Vertex>::const_iterator v_ = h_primVtx->begin();
                 v_ != h_primVtx->end(); ++v_){
                
                retPv.set(false);
                bool pass = (*pvObjSel_)( *v_, retPv );
                
                if (pass){
                    ++nSelVtx;
                    
                    // return a vector of refs to good primary vertices
                    good_pvs_.push_back( edm::Ptr<reco::Vertex>(h_primVtx, nVtx) );
                }
                
                ++nVtx;
            }
            
            if( nSelVtx >= cut("Min number of PV", int()) || ignoreCut("Min number of PV") ){
                passCut(ret, "Min number of PV");
            }
            else break;
            
            
            if( nSelVtx <= cut("Max number of PV", int()) || ignoreCut("Max number of PV") ){
                passCut(ret, "Max number of PV");
            }
            else break;
            
            passCut(ret, "Primary vertex cuts"); // PV cuts total
            
        } // end of PV cuts
        
        //======================================================
        //
        // jet loop
        //
        // We have to do it here, before the muon cuts, because
        // the muon cuts require dR(mu, jet)
        //
        int nGoodJets = 0;
        //int njets = 0;
        int njetsPF = 0;
        double leadJetPt = 0;
        int cnt_leadJetPt = 0;
        
        if ( considerCut("Muon cuts") || considerCut("Jet cuts") ) {
            
            good_jets_.clear();
            good_jetsPF_.clear();
            all_jetsPF_.clear();
            
            ///////////////
            // For PF Jets
            //////////////
            
            event.getByLabel( jetSrcPF_, h_jetsPF_ );
            for (vector<pat::Jet>::const_iterator jetPF = h_jetsPF_->begin();
                 jetPF != h_jetsPF_->end(); ++jetPF){
                
                retJet.set(false);
                bool passPF = (*jetSel_)( *jetPF, retJet );
                
                if ( passPF ){
                    // save all the good jets
                    ++nGoodJets;
                    good_jetsPF_.push_back(edm::Ptr<pat::Jet>(h_jetsPF_, njetsPF));
                    
                    ++cnt_leadJetPt;
                    if ( cnt_leadJetPt == 1 ){leadJetPt = (*jetPF).pt();}
                    
                }
                
                //Fill all the jets regardless if they pass or not
                all_jetsPF_.push_back(edm::Ptr<pat::Jet>(h_jetsPF_, njetsPF));
                
                ++njetsPF;
                
            } // end of loop over PF jets
            
            //if (nGoodJets == 0) std::cout << "NO GOOD JETS IN THIS EVET" <<std::endl;
            
            
        } // end of if (consider mu or jet cuts)
        
        
        // get rho correction
        event.getByLabel( rhoCorrectionSrc_, h_rho);
        double rho = (*h_rho);
        //double rho = 0.0;
        
        //
        //_____ Muon cuts ________________________________
        //
        if ( considerCut("Muon cuts") ) {
            
            // now, finally, get muons
            event.getByLabel( muonSrc_, h_muons_ );
            
            // loop over muons
            int n_muons = 0;
            int nSelMuons = 0;
            int nLooseMuons = 0;
            
            good_muons_.clear();
            loose_muons_.clear();
            
            for ( vector<pat::Muon>::const_iterator mu = h_muons_->begin();
                 mu != h_muons_->end(); mu++){
                
                retMuon.set(false);
                bool pass = (*muonSel_)( *mu, retMuon, &good_jetsPF_, &good_pvs_, rho  );
                
                if ( pass ){
                    ++nSelMuons;
                    
                    // save the first good muon
                    if (nSelMuons == 1){
                        muon0_ = edm::Ptr<pat::Muon>( h_muons_, n_muons);
                        
                    }
                    
                    // save the second good muon
                    else if (nSelMuons == 2){
                        muon1_ = edm::Ptr<pat::Muon>( h_muons_, n_muons);
                    }
                    
                    // save every good muon
                    good_muons_.push_back( edm::Ptr<pat::Muon>( h_muons_, n_muons) );
                    
                }
                else{ // check for loose muon
                    
                    retLooseMuon.set(false);
                    bool pass_loose = (*looseMuonSel_)( *mu, retLooseMuon, &good_jetsPF_, &good_pvs_,rho );
                    
                    if ( pass_loose ){
                        loose_muons_.push_back( edm::Ptr<pat::Muon>( h_muons_, n_muons) );
                        ++nLooseMuons;
                    }
                    
                }
                
                
                ++n_muons;
            } // end of the muon loop
            
            if( nSelMuons >= cut("Min mu multiplicity", int()) || ignoreCut("Min mu multiplicity") ) passCut(ret, "Min mu multiplicity");
            else break;
            
            if( nSelMuons <= cut("Max mu multiplicity", int()) || ignoreCut("Max mu multiplicity") ) passCut(ret, "Max mu multiplicity");
            else break;
            
            if( (nLooseMuons == 0 && nSelMuons == 1)  || ignoreCut("Loose muon veto") ) passCut(ret, "Loose muon veto");
            else break;
            
            if( nSelMuons == 0 || ignoreCut("Muon veto")) passCut(ret, "Muon veto");
            else break;
            
            
            passCut(ret, "Muon cuts"); // muon cuts total
            
        } // end of muon cuts
        
        //
        //_____ Electron cuts __________________________________
        //
        
        int n_elec = 0;
        int nSelElectrons = 0;
        int nLooseElectrons = 0;
        
        good_electrons_.clear();
        loose_electrons_.clear();
        
        if ( considerCut("Electron cuts") ) {
            
            event.getByLabel( electronSrc_, h_electrons_ );
            
            // loop over electrons
            for ( vector<pat::Electron>::const_iterator el = h_electrons_->begin();
                 el != h_electrons_->end(); el++){
                
                retElectron.set(false);
                bool pass = (*electronSel_)( *el, retElectron, &good_jetsPF_, &good_pvs_, rho);
                
                if ( pass ){
                    ++nSelElectrons;
                    
                    // save the first good muon
                    if (nSelElectrons == 1){
                        electron0_ = edm::Ptr<pat::Electron>( h_electrons_, n_elec);	      
                    }
                    // save the second good muon
                    else if (nSelElectrons == 2){
                        electron1_ = edm::Ptr<pat::Electron>( h_electrons_, n_elec);	      
                    }
                    // save every good electorn
                    good_electrons_.push_back( edm::Ptr<pat::Electron>( h_electrons_, n_elec) );
                    
                }	  
                else{ // check for loose electron
                    
                    retLooseElectron.set(false);
                    bool pass_loose = (*looseElectronSel_)( *el, retLooseElectron, &good_jetsPF_, &good_pvs_, rho);
                    
                    if ( pass_loose ){
                        loose_electrons_.push_back( edm::Ptr<pat::Electron>( h_electrons_, n_elec) );
                        ++nLooseElectrons;
                    }
                    
                }
                
                ++n_elec;
            } // end of the electron loop
            
            if( nSelElectrons >= cut("Min ele multiplicity", int()) || ignoreCut("Min ele multiplicity") ) passCut(ret, "Min ele multiplicity");
            else break;
            
            if( nSelElectrons <= cut("Max ele multiplicity", int()) || ignoreCut("Max ele multiplicity") ) passCut(ret, "Max ele multiplicity");
            else break;
            
            if( (nLooseElectrons == 0 && nSelElectrons == 1) || ignoreCut("Loose ele veto")  ) passCut(ret, "Loose ele veto");
            else break;
            
            if( nSelElectrons == 0 || ignoreCut("Electron veto")) passCut(ret, "Electron veto");
            else break;
            
            passCut(ret, "Electron cuts"); // muon cuts total
            
        } // end of electron cuts
        
        
        
        
        //
        if ( considerCut("Jet cuts") ) {
            
            if ( ignoreCut("One jet or more") || nGoodJets >= 1 ) passCut(ret, "One jet or more");
            else break; 
            
            if ( ignoreCut("Two jets or more") || nGoodJets >= 2 ) passCut(ret, "Two jets or more");
            else break; 
            
            if ( ignoreCut("Three jets or more") || nGoodJets >= 3 ) passCut(ret, "Three jets or more");
            else break; 
            
            if ( ignoreCut("Min jet multiplicity") || nGoodJets >= cut("Min jet multiplicity",int()) ) passCut(ret, "Min jet multiplicity");
            else break; 
            
            if ( ignoreCut("Max jet multiplicity") || nGoodJets >= cut("Max jet multiplicity",int()) ) passCut(ret, "Max jet multiplicity");
            else break; 
            
            if ( ignoreCut("Leading jet pt") || leadJetPt >= cut("Leading jet pt",double()) ) passCut(ret, "Leading jet pt");
            else break;
            
            
            passCut(ret, "Jet cuts");
            
        } // end of jet cuts
        
        
        
        
        //
        //_____ MET cuts __________________________________
        //      
        if ( considerCut("MET cuts") ) {
            
            event.getByLabel( metPFSrc_, h_pfMet_ );
            
            // pfMet
            if ( signed (h_pfMet_->size()) >= 0 ) {
                
                pat::MET const & met = h_pfMet_->at(0);
                
                retMet.set(false);
                bool pass = (*metSel_)( met, retMet );
                
                if (pass)
                {
                    pfMet_ = edm::Ptr<pat::MET>( h_pfMet_, 0);
                }
                
                if ( ignoreCut("MET pt") || pass ){
                    passCut(ret, "MET pt");
                }
                else break;
                
            }
            
            passCut(ret, "MET cuts");
            
        } // end of MET cuts
        
        
        
        //
        //_____ Btagging cuts _____________________
        //
        if ( considerCut("Btagging") ) {
            
            
            int nBtagJets = 0;
            double bTagCut = 2.0;
            for (std::vector<edm::Ptr<pat::Jet> >::const_iterator jet = good_jetsPF_.begin(); jet != good_jetsPF_.end(); ++jet){
                
                if ( (*jet)->bDiscriminator( "simpleSecondaryVertexHighEffBJetTags" ) >= bTagCut ) ++nBtagJets;
                
            }
            
            
            if ( nBtagJets >= 0 || ignoreCut("No btag jets") )  passCut(ret, "No btag jets");
            else break; 
            if ( nBtagJets >= 1 || ignoreCut("One btag jet or more") )  passCut(ret, "One btag jet or more");
            else break; 
            if ( nBtagJets >= 2 || ignoreCut("Two btag jets or more") )  passCut(ret, "Two btag jets or more");
            else break; 
            if ( nBtagJets >= 3 || ignoreCut("Three btag jets or more") )  passCut(ret, "Three btag jets or more");
            else break; 
            
            
            passCut(ret, "Btagging");
        }
        
        
        
        break;
        
    } // end of while loop
    
    
    
    return (bool)ret;
    
    setIgnored(ret);
    
    return false;
}// end of operator()
