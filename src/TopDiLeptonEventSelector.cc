// -*- C++ -*-
/*   FWLite PAT analyzer-selector for DiLepton Top events
 */

#include <iostream>
#include <cmath>      // for fabs()

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
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/SelectorUtils/interface/EventSelector.h"
#include "PhysicsTools/SelectorUtils/interface/PVSelector.h"
#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "PhysicsTools/SelectorUtils/interface/strbitset.h"



using namespace std;



class TopDiLeptonEventSelector : public BaseEventSelector {
    
public:
    
    
    enum Selection_t { Top, TopBTag, EWK_Zmumu, DiLepton, TprimeInc };
    enum Version_t { Off, test, SelV1, SelV2, SelV3, Tight, Loose, DiLepton_MuMu, DiLepton_ElEl, DiLepton_MuEl, Wprime_Mu, Wprime_Mu_forPrompt, Wprime_El, Tprime_Inc};
    
    
    
    TopDiLeptonEventSelector();
    ~TopDiLeptonEventSelector();
    
    
    
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
//static int reg = LjmetFactory::GetInstance()->Register(new TopDiLeptonEventSelector(), "DiLeptonSelector");



TopDiLeptonEventSelector::TopDiLeptonEventSelector(){
}



TopDiLeptonEventSelector::~TopDiLeptonEventSelector(){
}




void TopDiLeptonEventSelector::BeginJob(std::map<std::string, edm::ParameterSet const> par){
    //
    // Actions needed before event loop happens
    // The map contains all parameters from the python config
    //
    
    std::string versionStr = par["event_selector"].getParameter<std::string>("version");
    
    
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
    legend = "[TopDiLeptonEventSelector]: ";
    
    
    // create object selectors
    muonSel_          = boost::shared_ptr<MuonSelectionFunctor> ( new MuonSelectionFunctor (par["muonSelector"]         ) );
    looseMuonSel_     = boost::shared_ptr<MuonSelectionFunctor> ( new MuonSelectionFunctor (par["looseMuonSelector"]    ) );
    electronSel_      = boost::shared_ptr<ElectronSelector>     ( new ElectronSelector     (par["electronSelector"]     ) );
    looseElectronSel_ = boost::shared_ptr<ElectronSelector>     ( new ElectronSelector     (par["looseElectronSelector"]) );
    jetSel_           = boost::shared_ptr<JetIDSelectionFunctor>( new JetIDSelectionFunctor(par["jetIDSelector"]        ) );
    metSel_           = boost::shared_ptr<MetSelectionFunctor>  ( new MetSelectionFunctor  (par["metSelector"]          ) );
    pvSel_            = boost::shared_ptr<PVSelector>           ( new PVSelector           (par["pvSelector"]           ) );
    pvObjSel_         = boost::shared_ptr<PvObjectSelector>     ( new PvObjectSelector     (par["pvObjectSelector"]     ) );
    
    
    // this event selector is for top dilepton
    selection = DiLepton;
    
    if ( versionStr == "DiLepton_MuMu" ) version = DiLepton_MuMu;
    else if ( versionStr == "DiLepton_ElEl" ) version = DiLepton_ElEl;
    else if ( versionStr == "DiLepton_MuEl" ) version = DiLepton_MuEl;
    else {
        throw cms::Exception("InvalidInput") << "Expect version to be one of DiLepton_MuMu, DiLepton_ElEl, DiLepton_MuEl" << std::endl;
    }
    
    std::cout << legend << "initializing top Dilepton selection" << std::endl;
    std::cout << legend << "selection version: " << versionStr << std::endl;
    
    initialize( version );
}



void TopDiLeptonEventSelector::initialize( Version_t version){
    
    push_back("No selection");
    set("No selection");
    
    push_back("HLT");
    push_back("Trigger names"); // activate to get a list of available trigger names
    push_back("Trigger cuts");
    
    push_back("Min number of PV");
    push_back("Max number of PV");
    push_back("Primary vertex cuts");
    
    push_back("Min mu multiplicity");
    push_back("Max mu multiplicity");
    push_back("Loose muon veto");
    push_back("Muon cuts");
    
    push_back("Min ele multiplicity");
    push_back("Max ele multiplicity");
    push_back("Loose ele veto");
    push_back("Electron cuts");
    
    push_back("Opposite charges");
    push_back("quarkonia");
    push_back("Dilepton mass");
    push_back("Z candidate cuts");
    
    push_back("One jet or more");
    push_back("Two jets or more");
    push_back("Three jets or more");
    push_back("Min jet multiplicity");
    push_back("Max jet multiplicity");
    push_back("Jet cuts");
    
    push_back("MET pt");
    push_back("MET cuts");
    
    push_back("No btag jets");
    push_back("One btag jet or more");
    push_back("Two btag jets or more");
    push_back("Btagging");
    
    
    set("HLT", false);
    set("Trigger names", false);
    set("Trigger cuts", false);
    
    set("Min number of PV", false);
    set("Max number of PV", false);
    set("Primary vertex cuts", false);
    
    set("Min mu multiplicity", false);
    set("Max mu multiplicity", false);
    set("Loose muon veto", false);
    set("Muon cuts", false);
    
    set("Min ele multiplicity", false);
    set("Max ele multiplicity", false);
    set("Loose ele veto", false);
    set("Electron cuts", false);
    
    
    set("Opposite charges", false);
    set("quarkonia", false);
    set("Dilepton mass", false);
    set("Z candidate cuts", false);
    
    set("One jet or more", false);
    set("Two jets or more", false);
    set("Three jets or more", false);
    set("Min jet multiplicity", false);
    set("Max jet multiplicity", false);
    set("Jet cuts", false);
    
    set("MET pt", false);
    set("MET cuts", false);
    
    set("No btag jets", false);
    set("One btag jet or more", false);
    set("Two btag jets or more", false);
    set("Btagging", false);
    
    if (version == DiLepton_MuMu){
        
        set("HLT", true);
        set("Trigger names", false);
        set("Trigger cuts", true);
        
        set("Min number of PV", 1);
        set("Max number of PV", false);
        set("Primary vertex cuts", true);
        
        set("Min mu multiplicity", 2);
        set("Max mu multiplicity", 2);
        set("Loose muon veto", false);
        set("Muon cuts", true);
        
        set("Min ele multiplicity", false);
        set("Max ele multiplicity", false);
        set("Loose ele veto", false);
        set("Electron cuts", false);
        
        set("Opposite charges", false);
        set("quarkonia", 12);
        //set("Dilepton mass", true);
        set("Dilepton mass", false);
        set("Z candidate cuts", true);
        
        set("One jet or more", false);
        set("Two jets or more", false);
        set("Three jets or more", false);
        set("Min jet multiplicity", false);
        set("Max jet multiplicity", false);
        set("Jet cuts");
        
        set("MET pt", false);
        set("MET cuts", true); // keep true for access to met collections
        
        set("No btag jets", false);
        set("One btag jet or more", false);
        set("Two btag jets or more", false);
        set("Btagging", false);
        
    }
    
    
    
    if (version == DiLepton_ElEl){
        
        set("HLT", true);
        set("Trigger names", false);
        set("Trigger cuts", true);
        
        set("Min number of PV", 1);
        set("Max number of PV", false);
        set("Primary vertex cuts", true);
        
        set("Min mu multiplicity", false);
        set("Max mu multiplicity", false);
        set("Loose muon veto", false);
        set("Muon cuts", false);
        
        set("Min ele multiplicity", 2);
        set("Max ele multiplicity", 2);
        set("Loose ele veto", false);
        set("Electron cuts", true);
        
        set("Opposite charges", false);
        set("quarkonia", 12);
        set("Dilepton mass", true);
        set("Z candidate cuts", true);
        
        set("One jet or more", false);
        set("Two jets or more", false);
        set("Three jets or more", false);
        set("Min jet multiplicity", false);
        set("Max jet multiplicity", false);
        set("Jet cuts");
        
        set("MET pt", false);
        set("MET cuts", true); // keep true for access to met collections
        
        set("No btag jets", false);
        set("One btag jet or more", false);
        set("Two btag jets or more", false);
        set("Btagging", false);
        
    }
    
    
    
    
    if (version == DiLepton_MuEl){
        
        set("HLT", true);
        set("Trigger names", false);
        set("Trigger cuts", true);
        
        set("Min number of PV", 1);
        set("Max number of PV", false);
        set("Primary vertex cuts", true);
        
        set("Min mu multiplicity", 1);
        set("Max mu multiplicity", 1);
        set("Loose muon veto", false);
        set("Muon cuts", true);
        
        set("Min ele multiplicity", 1);
        set("Max ele multiplicity", 1);
        set("Loose ele veto", false);
        set("Electron cuts", true);
        
        set("Opposite charges", false);
        set("quarkonia", 12);
        set("Dilepton mass", false);
        set("Z candidate cuts", true);
        
        set("One jet or more", false);
        set("Two jets or more", false);
        set("Three jets or more", false);
        set("Min jet multiplicity", false);
        set("Max jet multiplicity", false);
        set("Jet cuts");
        
        set("MET pt", false);
        set("MET cuts", true); // keep true for access to met collections
        
        set("No btag jets", false);
        set("One btag jet or more", false);
        set("Two btag jets or more", false);
        set("Btagging", false);
        
    }
    
    
} // initialize()




bool TopDiLeptonEventSelector::operator()( edm::Event const & event, pat::strbitset & ret){
    
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
            
            bool passTrig = true;
            
            edm::InputTag _triggerEventSrc("TriggerResults::HLT");
            event.getByLabel( _triggerEventSrc, mhEdmTriggerResults );
            const edm::ParameterSetID ps = mhEdmTriggerResults->parameterSetID();
            const edm::TriggerNames trigNames = event.triggerNames(*mhEdmTriggerResults);
            
            std::vector<std::string> pathName_mumu_data;
            pathName_mumu_data . push_back("HLT_DoubleMu7_v1");
            pathName_mumu_data . push_back("HLT_DoubleMu7_v2");
            
            pathName_mumu_data . push_back("HLT_Mu13_Mu8_v2");
            pathName_mumu_data . push_back("HLT_Mu13_Mu8_v3");
            pathName_mumu_data . push_back("HLT_Mu13_Mu8_v4");
            pathName_mumu_data . push_back("HLT_Mu13_Mu8_v6");
            pathName_mumu_data . push_back("HLT_Mu13_Mu8_v7");
            
            pathName_mumu_data . push_back("HLT_Mu17_Mu8_v10");
            pathName_mumu_data . push_back("HLT_Mu17_Mu8_v11");
            
            std::vector<std::string> pathName_elel_data;
            pathName_elel_data . push_back("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1");
            pathName_elel_data . push_back("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2");
            pathName_elel_data . push_back("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3");
            pathName_elel_data . push_back("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4");
            pathName_elel_data . push_back("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5");
            pathName_elel_data . push_back("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6");
            
            pathName_elel_data . push_back("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2");
            pathName_elel_data . push_back("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3");
            pathName_elel_data . push_back("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v4");
            pathName_elel_data . push_back("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v5");
            
            pathName_elel_data . push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5");
            pathName_elel_data . push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6");
            pathName_elel_data . push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7");
            pathName_elel_data . push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8");
            pathName_elel_data . push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9");
            pathName_elel_data . push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10");
            
            std::vector<std::string> pathName_muel_data;
            pathName_muel_data . push_back("HLT_Mu10_Ele10_CaloIdL_v2");
            pathName_muel_data . push_back("HLT_Mu10_Ele10_CaloIdL_v3");
            pathName_muel_data . push_back("HLT_Mu10_Ele10_CaloIdL_v4");
            
            pathName_muel_data . push_back("HLT_Mu17_Ele8_CaloIdL_v1");
            pathName_muel_data . push_back("HLT_Mu17_Ele8_CaloIdL_v2");
            pathName_muel_data . push_back("HLT_Mu17_Ele8_CaloIdL_v3");
            pathName_muel_data . push_back("HLT_Mu17_Ele8_CaloIdL_v4");
            pathName_muel_data . push_back("HLT_Mu17_Ele8_CaloIdL_v5");
            pathName_muel_data . push_back("HLT_Mu17_Ele8_CaloIdL_v6");
            pathName_muel_data . push_back("HLT_Mu17_Ele8_CaloIdL_v8");
            
            pathName_muel_data . push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4");
            pathName_muel_data . push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7");
            pathName_muel_data . push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8");
            
            pathName_muel_data . push_back("HLT_Mu8_Ele17_CaloIdL_v1");
            pathName_muel_data . push_back("HLT_Mu8_Ele17_CaloIdL_v2");
            pathName_muel_data . push_back("HLT_Mu8_Ele17_CaloIdL_v3");
            pathName_muel_data . push_back("HLT_Mu8_Ele17_CaloIdL_v4");
            pathName_muel_data . push_back("HLT_Mu8_Ele17_CaloIdL_v5");
            pathName_muel_data . push_back("HLT_Mu8_Ele17_CaloIdL_v6");
            
            pathName_muel_data . push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3");
            pathName_muel_data . push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4");
            pathName_muel_data . push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7");
            pathName_muel_data . push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8");
            
            //MC Selection is for Summer11 and will have to be updated when Fall11 becomes available
            std::vector<std::string> pathName_mumu_mc;
            pathName_mumu_mc . push_back("HLT_DoubleMu6_v1");
            pathName_mumu_mc . push_back("HLT_DoubleMu7_v1");
            
            std::vector<std::string> pathName_elel_mc;
            pathName_elel_mc . push_back("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2");
            pathName_elel_mc . push_back("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2");
            
            std::vector<std::string> pathName_muel_mc;
            pathName_muel_mc . push_back("HLT_Mu10_Ele10_CaloIdL_v3");
            pathName_muel_mc . push_back("HLT_Mu8_Ele17_CaloIdL_v2");
            pathName_muel_mc . push_back("HLT_Mu17_Ele8_CaloIdL_v2");
            
            //Generic to Loop over
            std::vector<std::string> pathName;
            
            bool data = true;
            if (event.id().run() < 160000) data = false;
            
            if (  version == DiLepton_MuMu ) {
                if (data) pathName = pathName_mumu_data;
                else  pathName = pathName_mumu_mc;
            }
            if (  version == DiLepton_ElEl ) {
                if (data) pathName = pathName_elel_data;
                else  pathName = pathName_elel_mc;
            }
            if (  version == DiLepton_MuEl ) {
                if (data) pathName = pathName_muel_data;
                else  pathName = pathName_muel_mc;
            }
            
            unsigned int ii = 0;
            while( ii < pathName.size()){
                
                if ( trigNames.size() != trigNames.triggerIndex(pathName[ii])  ) {
                    
                    bool passTrigTmp = mhEdmTriggerResults->accept(trigNames.triggerIndex(pathName[ii]));
                    
                    if (passTrigTmp == false) {
                        passTrig = false;
                        break;
                    }
                }
                ii++;
            }
            
            
            if ( ignoreCut("HLT") || passTrig ) passCut(ret, "HLT");
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
        int njetsPF = 0;
        
        if ( considerCut("Muon cuts") || considerCut("Jet cuts") ) {
            
            good_jets_.clear();
            good_jetsPF_.clear();
            
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
                }
                
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
        
        int n_muons = 0;
        int nSelMuons = 0;
        int nLooseMuons = 0;
        
        good_muons_.clear();
        loose_muons_.clear();
        
        if ( considerCut("Muon cuts") ) {
            
            // now, finally, get muons
            event.getByLabel( muonSrc_, h_muons_ );
            
            // loop over muons
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
                    bool pass_loose = (*looseMuonSel_)( *mu, retLooseMuon, &good_jetsPF_, &good_pvs_, rho );
                    
                    if ( pass_loose ){
                        loose_muons_.push_back( edm::Ptr<pat::Muon>( h_muons_, n_muons) );
                        ++nLooseMuons;
                    }
                    
                }
                
                ++n_muons;
            } // end of the muon loop
            
            if( nSelMuons >= cut("Min mu multiplicity", int()) || ignoreCut("Min mu multiplicity") ) passCut(ret, "Min mu multiplicity");
            else break;
            
            if( nSelMuons <= cut("Max mu multiplicity", int()) || ignoreCut("Max mu multiplicity")  ) passCut(ret, "Max mu multiplicity");
            else break;
            
            if( (nLooseMuons == 0 && nSelMuons == 1) || ignoreCut("Loose muon veto")  ) passCut(ret, "Loose muon veto");
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
                bool pass = (*electronSel_)( *el, retElectron, &good_jetsPF_, &good_pvs_, rho );
                
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
                    bool pass_loose = (*looseElectronSel_)( *el, retLooseElectron, &good_jetsPF_, &good_pvs_, rho );
                    
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
            
            passCut(ret, "Electron cuts"); // muon cuts total
            
        } // end of electron cuts
        
        
        
        //
        //_____ Z candidate cuts __________________________________
        //      
        if ( considerCut("Z candidate cuts") ) {
            
            // opposite lepton charges
            
            if (version == DiLepton_MuMu){
                if ( (n_muons + nLooseMuons > 1 && ( muon0_ -> charge() + muon1_ -> charge() ) == 0) || ignoreCut("Opposite charges") ) passCut(ret, "Opposite charges");
                else break;
            }
            if (version == DiLepton_ElEl){
                if ( (n_elec + nLooseElectrons > 1 && ( electron0_ -> charge() + electron1_ -> charge() ) == 0) || ignoreCut("Opposite charges") ) passCut(ret, "Opposite charges");
                else break;
            }
            if (version == DiLepton_MuEl){
                if ( (n_muons + n_elec > 1 && ( muon0_ -> charge() + electron0_ -> charge() ) == 0) || ignoreCut("Opposite charges") ) passCut(ret, "Opposite charges");
                else break;
            }
            
            
            // dimuon mass, reject events with Zmass +/- 15GeV 
            if ( considerCut("Dilepton mass") || considerCut("quarkonia")){
                double _m2 = -1;
                
                if (version == DiLepton_MuMu) _m2 = (muon0_->p4() + muon1_->p4()).M2();
                if (version == DiLepton_ElEl) _m2 = (electron0_->p4() + electron1_->p4()).M2();
                if (version == DiLepton_MuEl) _m2 = (muon0_->p4() + electron0_->p4()).M2();
                
                if ( (_m2 > 0 && sqrt(_m2) > cut("quarkonia", double())) || ignoreCut("quarkonia") ) passCut(ret, "quarkonia");
                else break;
                
                if ( ((_m2 > 0 && (sqrt(_m2) < 106 || sqrt(_m2) > 76))) || ignoreCut("Dilepton mass")  )  passCut(ret, "Dilepton mass");
                else break;
                
            }
            passCut(ret, "Z candidate cuts");
        }
        
        
        
        // Jet Cuts
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
            
            passCut(ret, "Jet cuts");
            
        } // end of jet cuts
        
        
        
        
        //
        //_____ MET cuts __________________________________
        //      
        if ( considerCut("MET cuts") ) {
            
            event.getByLabel( metSrc_, h_met_ );
            event.getByLabel( metTCSrc_, h_tcMet_ );
            event.getByLabel( metPFSrc_, h_pfMet_ );
            event.getByLabel( metPfTypeISrc_, h_pfTypeIMet_ );
            event.getByLabel( metMuonSrc_, h_caloMet_ );
            event.getByLabel( metTypeIISrc_, h_caloMet2_ );
            
            
            // pfMet
            
            if ( signed (h_pfMet_->size()) >= 0 ) {
                
                pat::MET const & met = h_pfMet_->at(0);
                
                retMet.set(false);
                bool pass = (*metSel_)( met, retMet );
                
                if (pass)
                {
                    pfMet_ = edm::Ptr<pat::MET>( h_pfMet_, 0);
                }
                
                if ( ignoreCut("MET pt") || retMet.test("pt") ){
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
            double bTagCut = 1.7;
            for (std::vector<edm::Ptr<pat::Jet> >::const_iterator jet = good_jetsPF_.begin(); jet != good_jetsPF_.end(); ++jet){
                
                if ( (*jet)->bDiscriminator( "trackCountingHighEffBJetTags" ) >= bTagCut ) ++nBtagJets;
                
            }
            
            
            if ( nBtagJets >= 0 || ignoreCut("No btag jets") )  passCut(ret, "No btag jets");
            else break; 
            if ( nBtagJets >= 1 || ignoreCut("One btag jet or more") )  passCut(ret, "One btag jet or more");
            else break; 
            if ( nBtagJets >= 2 || ignoreCut("Two btag jets or more") )  passCut(ret, "Two btag jets or more");
            else break; 
            
            passCut(ret, "Btagging");
        }
        
        
        
        
        break;
        
    } // end of while loop
    
    
    
    return (bool)ret;
    
    setIgnored(ret);
    
    return false;
}// end of operator()
