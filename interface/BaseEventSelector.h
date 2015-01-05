#ifndef LJMet_Com_interface_BaseEventSelector_h
#define LJMet_Com_interface_BaseEventSelector_h

/*
 Interface class for FWLite PAT analyzer-selectors
 Specific selectors must implement the () operator
 
 Author: Gena Kukartsev, 2010, 2012
         Orduna@2014
 */

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <utility> // std::pair

#include "FWCore/Framework/interface/Event.h"
#include "LJMet/Com/interface/LjmetEventContent.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "PhysicsTools/SelectorUtils/interface/EventSelector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "TLorentzVector.h"
#include "LJMet/Com/interface/BTagSFUtil.h"
#include "LJMet/Com/interface/BtagHardcodedConditions.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

//#include "TROOT.h"
//#include "TVector3.h"

class BaseEventSelector : public EventSelector {
    //
    // Base class for all event selector plugins
    //
    
    friend class LjmetFactory;
    
public:
    BaseEventSelector();
    virtual ~BaseEventSelector() { };
    virtual void BeginJob(std::map<std::string, edm::ParameterSet const > par);
    virtual bool operator()( edm::Event const & event, pat::strbitset & ret) = 0;
    virtual void EndJob() { }
    virtual void AnalyzeEvent( edm::Event const & event, LjmetEventContent & ec ) { }
    std::string GetName() { return mName; }
    /// Evaluates a signed perp components of v1 relative to v2. The sign is defined by Phi
    double GetPerp(TVector3 & v1, TVector3 & v2);
    bool AreMatched(const reco::Candidate & c1, const reco::Candidate & c2, double DR, double DPtRel) { return ((reco::deltaR(c1, c2) < DR) && (abs(c2.pt() - c1.pt())/ c2.pt() < DPtRel)); }
    
    std::vector<edm::Ptr<pat::Jet>> const & GetAllJets() const { return mvAllJets; }
    std::vector<edm::Ptr<pat::Jet>> const & GetSelectedJets() const { return mvSelJets; }
    std::vector<edm::Ptr<pat::Jet>> const & GetLooseJets() const { return mvSelJets; }
    std::vector<edm::Ptr<pat::Jet>> const & GetSelectedBtagJets() const { return mvSelBtagJets; }
    std::vector<std::pair<TLorentzVector, bool>> const & GetCorrJetsWithBTags() const { return mvCorrJetsWithBTags; }
    std::vector<edm::Ptr<pat::Muon>> const & GetAllMuons() const { return mvAllMuons; }
    std::vector<edm::Ptr<pat::Muon>> const & GetSelectedMuons() const { return mvSelMuons; }
    std::vector<edm::Ptr<pat::Muon>> const & GetLooseMuons() const { return mvLooseMuons; }
    std::vector<edm::Ptr<pat::Electron>> const & GetAllElectrons() const { return mvAllElectrons; }
    std::vector<edm::Ptr<pat::Electron>> const & GetSelectedElectrons() const { return mvSelElectrons; }
    std::vector<edm::Ptr<pat::Electron>> const & GetLooseElectrons() const { return mvLooseElectrons; }
    edm::Ptr<pat::MET> const & GetMet() const { return mpMet; }
    edm::Ptr<reco::PFMET> const & GetType1CorrMet() const { return mpType1CorrMet; }
    TLorentzVector const & GetCorrectedMet() const { return correctedMET_p4; }
    std::vector<unsigned int> const & GetSelectedTriggers() const { return mvSelTriggers; }
    std::vector<edm::Ptr<reco::Vertex>> const & GetSelectedPVs() const { return mvSelPVs; }
    double const & GetTestValue() const { return mTestValue; }
    void SetMc(bool isMc) { mbIsMc = isMc; }
    bool IsMc() { return mbIsMc; }
    
    // LJMET event content setters
    void Init( void );
    void SetEventContent(LjmetEventContent * pEc) { mpEc = pEc; }
    /// Declare a new histogram to be created for the module
    void SetHistogram(std::string name, int nbins, double low, double high) { mpEc->SetHistogram(mName, name, nbins, low, high); }
    void SetHistValue(std::string name, double value) { mpEc->SetHistValue(mName, name, value); }
    void SetTestValue(double & test) { mTestValue = test; }
    
    void SetCorrectedMet(TLorentzVector & met) { correctedMET_p4 = met; }
    void SetCorrJetsWithBTags(std::vector<std::pair<TLorentzVector, bool>> & jets) { mvCorrJetsWithBTags = jets; }
    
    bool isJetTagged(const pat::Jet &jet, edm::Event const & event, bool applySF = true);
    TLorentzVector correctJet(const pat::Jet & jet, edm::Event const & event);
    TLorentzVector correctMet(const pat::MET & met, edm::Event const & event);
    
protected:
    std::vector<edm::Ptr<pat::Jet>> mvAllJets;
    std::vector<edm::Ptr<pat::Jet>> mvSelJets;
    std::vector<edm::Ptr<pat::Jet>> mvLooseJets;
    std::vector<std::pair<TLorentzVector, bool>> mvCorrJetsWithBTags;
    std::vector<edm::Ptr<pat::Jet>> mvSelBtagJets;
    std::vector<edm::Ptr<pat::Muon>> mvAllMuons;
    std::vector<edm::Ptr<pat::Muon>> mvSelMuons;
    std::vector<edm::Ptr<pat::Muon>> mvLooseMuons;
    std::vector<edm::Ptr<pat::Electron>> mvAllElectrons;
    std::vector<edm::Ptr<pat::Electron>> mvSelElectrons;
    std::vector<edm::Ptr<pat::Electron>> mvLooseElectrons;
    edm::Ptr<pat::MET> mpMet;
    edm::Ptr<reco::PFMET> mpType1CorrMet;
    TLorentzVector correctedMET_p4;
    std::vector<unsigned int> mvSelTriggers;
    std::vector<edm::Ptr<reco::Vertex>> mvSelPVs;
    double mTestValue;
    
    // containers for config parameter values
    std::map<std::string, bool> mbPar;
    std::map<std::string, int> miPar;
    std::map<std::string, double> mdPar;
    std::map<std::string, std::string> msPar;
    std::map<std::string, edm::InputTag> mtPar;
    std::map<std::string, std::vector<std::string>> mvsPar;
    
    std::string mName;
    std::string mLegend;
    bool mbIsMc;
    
private:
    int mNCorrJets;
    int mNBtagSfCorrJets;
    double bTagCut;
    BTagSFUtil mBtagSfUtil;
    BtagHardcodedConditions mBtagCond;
    JetCorrectionUncertainty *jecUnc;
    LjmetEventContent * mpEc;
    
    /// Private init method to be called by LjmetFactory when registering the selector
    void init() { mLegend = "[" + mName + "]: "; std::cout << mLegend << "registering " << mName << std::endl; }
    void setName(std::string name) { mName = name; }
    /// Do what any event selector must do before event gets checked
    void BeginEvent(edm::Event const & event, LjmetEventContent & ec) { mNCorrJets = 0; mNBtagSfCorrJets = 0; }
    /// Do what any event selector must do after event processing is done, but before event content gets saved to file
    void EndEvent(edm::Event const & event, LjmetEventContent & ec) { SetHistValue("nBtagSfCorrections", mNBtagSfCorrJets); }
};

#endif
