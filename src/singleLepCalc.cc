/*
 Calculator for a generic single lepton analysis
 
 Author: Joshua Swanson, 2014
 */

#include <iostream>
#include <limits>   // std::numeric_limits

#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "TLorentzVector.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Math/interface/deltaR.h"

class LjmetFactory;

class singleLepCalc : public BaseCalc {
public:
    singleLepCalc();
    virtual ~singleLepCalc();
    virtual int BeginJob();
    virtual int AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector);
    virtual int EndJob();
    
private:
    edm::InputTag triggerSummary_;
    edm::InputTag triggerCollection_;
    edm::InputTag rhoSrc_;
    edm::InputTag muonSrc_;
    bool isWJets_;
    bool isTB_;
    bool isTT_;
};

//static int reg = LjmetFactory::GetInstance()->Register(new singleLepCalc(), "singleLepCalc");

singleLepCalc::singleLepCalc()
{
}

singleLepCalc::~singleLepCalc()
{
}

int singleLepCalc::BeginJob()
{
    if (mPset.exists("triggerSummary")) triggerSummary_ = mPset.getParameter<edm::InputTag>("triggerSummary");
    else                                triggerSummary_ = edm::InputTag("hltTriggerSummaryAOD");
    
    if (mPset.exists("triggerCollection")) triggerCollection_ = mPset.getParameter<edm::InputTag>("triggerCollection");
    else                                triggerCollection_ = edm::InputTag("TriggerResults::HLT");
    
    if (mPset.exists("rhoSrc")) rhoSrc_ = mPset.getParameter<edm::InputTag>("rhoSrc");
    else                        rhoSrc_ = edm::InputTag("fixedGridRhoAll");
    
    if (mPset.exists("isTB"))    isTB_ = mPset.getParameter<bool>("isTB");
    else                         isTB_ = false;
    
    if (mPset.exists("isTT"))    isTT_ = mPset.getParameter<bool>("isTT");
    else                         isTT_ = false;
    
    if (mPset.exists("isWJets")) isWJets_ = mPset.getParameter<bool>("isWJets");
    else                         isWJets_ = false;
    
    return 0;
}

int singleLepCalc::AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector)
{
    // ----- Get objects from the selector -----
    std::vector<edm::Ptr<pat::Jet> >            const & vSelJets = selector->GetSelectedJets();
    std::vector<edm::Ptr<pat::Jet> >            const & vSelBtagJets = selector->GetSelectedBtagJets();
    std::vector<edm::Ptr<pat::Jet> >            const & vAllJets = selector->GetAllJets();
    std::vector<std::pair<TLorentzVector,bool>> const & vCorrBtagJets = selector->GetCorrJetsWithBTags();
    std::vector<edm::Ptr<pat::Muon> >           const & vSelMuons = selector->GetSelectedMuons();
    std::vector<edm::Ptr<pat::Muon> >           const & vLooseMuons = selector->GetLooseMuons();
    std::vector<edm::Ptr<pat::Electron> >       const & vSelElectrons = selector->GetSelectedElectrons();
    std::vector<edm::Ptr<pat::Electron> >       const & vLooseElectrons = selector->GetLooseElectrons();
    edm::Ptr<pat::MET>                          const & pMet = selector->GetMet();
    edm::Ptr<reco::PFMET>                       const & pType1CorrMet = selector->GetType1CorrMet();
    TLorentzVector                              const & corrMET = selector->GetCorrectedMet();
    std::vector<edm::Ptr<reco::Vertex> >        const & vSelPVs = selector->GetSelectedPVs();
    
    // ----- Event kinematics -----
    int _nAllJets        = (int)vAllJets.size();
    int _nSelJets        = (int)vSelJets.size();
    int _nSelBtagJets    = (int)vSelBtagJets.size();
    int _nCorrBtagJets   = (int)vCorrBtagJets.size();
    int _nSelMuons       = (int)vSelMuons.size();
    int _nLooseMuons     = (int)vLooseMuons.size();
    int _nSelElectrons   = (int)vSelElectrons.size();
    int _nLooseElectrons = (int)vLooseElectrons.size();
    int _nSelPVs         = (int)vSelPVs.size();
    
    SetValue("nPV", _nSelPVs);
    
    // Trigger
    edm::Handle<edm::TriggerResults > mhEdmTriggerResults;
    event.getByLabel(triggerCollection_ , mhEdmTriggerResults);
    //const edm::ParameterSetID ps = mhEdmTriggerResults->parameterSetID();
    const edm::TriggerNames trigNames = event.triggerNames(*mhEdmTriggerResults);
    
    unsigned int _tSize = mhEdmTriggerResults->size();
    
    int passTrigEle27v10 = -1;
    int passTrigIsoMu24v13 = -1;
    
    unsigned int _tElMCIndex = trigNames.triggerIndex("HLT_Ele27_WP80_v10");
    if ( _tElMCIndex<_tSize){
        passTrigEle27v10 = mhEdmTriggerResults->accept(_tElMCIndex);
    }
    
    unsigned int _tMuMCIndex = trigNames.triggerIndex("HLT_IsoMu24_eta2p1_v13");
    if ( _tMuMCIndex<_tSize){
        passTrigIsoMu24v13 = mhEdmTriggerResults->accept(_tMuMCIndex);
    }
    
    // Muon
    double _muon_1_pt     = -std::numeric_limits<double>::max();
    double _muon_1_phi    = -std::numeric_limits<double>::max();
    double _muon_1_eta    = -std::numeric_limits<double>::max();
    double _muon_1_RelIso = -std::numeric_limits<double>::max();
    
    if (_nSelMuons>0) {
        std::cout << _muon_1_pt << std::endl;
        _muon_1_pt = vSelMuons[0]->pt();
        std::cout << _muon_1_pt << std::endl;
        int test ; std::cin >> test;
        _muon_1_phi = vSelMuons[0]->phi();
        _muon_1_eta = vSelMuons[0]->eta();
        
        double chIso = vSelMuons[0]->userIsolation(pat::PfChargedHadronIso);
        double nhIso = vSelMuons[0]->userIsolation(pat::PfNeutralHadronIso);
        double gIso  = vSelMuons[0]->userIsolation(pat::PfGammaIso);
        double puIso = vSelMuons[0]->userIsolation(pat::PfPUChargedHadronIso);
        _muon_1_RelIso = (chIso + std::max(0.,nhIso + gIso - 0.5*puIso))/_muon_1_pt;
        
    }
    
    SetValue("muon_1_pt",_muon_1_pt);
    SetValue("muon_1_phi",_muon_1_phi);
    SetValue("muon_1_eta",_muon_1_eta);
    SetValue("muon_1_RelIso",_muon_1_RelIso);
    
    // Electron
    
    edm::Handle<double> rhoHandle;
    event.getByLabel(rhoSrc_, rhoHandle);
    double rhoIso = std::max(*(rhoHandle.product()), 0.0);
    
    double _electron_1_pt     = -std::numeric_limits<double>::max();
    double _electron_1_phi    = -std::numeric_limits<double>::max();
    double _electron_1_eta    = -std::numeric_limits<double>::max();
    double _electron_1_RelIso = -std::numeric_limits<double>::max();
    
    if (_nSelElectrons>0) {
        _electron_1_pt = vSelElectrons[0]->ecalDrivenMomentum().pt();
        _electron_1_phi = vSelElectrons[0]->phi();
        _electron_1_eta = vSelElectrons[0]->eta();
        
        double AEff  = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, vSelElectrons[0]->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012);
        double chIso = vSelElectrons[0]->chargedHadronIso();
        double nhIso = vSelElectrons[0]->neutralHadronIso();
        double phIso = vSelElectrons[0]->photonIso();
        _electron_1_RelIso  = ( chIso + max(0.0, nhIso + phIso - rhoIso*AEff) )/ vSelElectrons[0]->ecalDrivenMomentum().pt();
        
    }
    
    SetValue("elec_1_pt", _electron_1_pt);
    SetValue("elec_1_phi", _electron_1_phi);
    SetValue("elec_1_eta", _electron_1_eta);
    SetValue("elec_1_RelIso", _electron_1_RelIso);
    
    double _electron_2_pt     = -std::numeric_limits<double>::max();
    double _electron_2_phi    = -std::numeric_limits<double>::max();
    double _electron_2_eta    = -std::numeric_limits<double>::max();
    double _electron_2_RelIso = -std::numeric_limits<double>::max();
    
    if (_nSelElectrons>1) {
        _electron_2_pt = vSelElectrons[1]->ecalDrivenMomentum().pt();
        _electron_2_phi = vSelElectrons[1]->phi();
        _electron_2_eta = vSelElectrons[1]->eta();
        
        double AEff  = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, vSelElectrons[1]->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012);
        double chIso = vSelElectrons[1]->chargedHadronIso();
        double nhIso = vSelElectrons[1]->neutralHadronIso();
        double phIso = vSelElectrons[1]->photonIso();
        _electron_2_RelIso  = ( chIso + max(0.0, nhIso + phIso - rhoIso*AEff) )/ vSelElectrons[1]->ecalDrivenMomentum().pt();
        
    }
    
    SetValue("elec_2_pt", _electron_2_pt);
    SetValue("elec_2_phi", _electron_2_phi);
    SetValue("elec_2_eta", _electron_2_eta);
    SetValue("elec_2_RelIso", _electron_2_RelIso);
    
    // Trigger Matching
    //    edm::Handle<trigger::TriggerEvent> mhEdmTriggerEvent;
    //    event.getByLabel(triggerSummary_,mhEdmTriggerEvent);
    //    trigger::TriggerObjectCollection allObjects = mhEdmTriggerEvent->getObjects();
    //
    //    int _electron_1_hltmatched =0;
    //    int _electron_2_hltmatched =0;
    //
    //    if (_nSelElectrons>0) {
    //      for(int i=0; i<mhEdmTriggerEvent->sizeFilters(); i++){
    //        if( mhEdmTriggerEvent->filterTag(i).label()!="hltEle27WP80TrackIsoFilter") continue;
    //	trigger::Keys keys = mhEdmTriggerEvent->filterKeys(i);
    //        double dR1 = 999.0;
    //        double dR2 = 999.0;
    //        for(size_t j=0; j<keys.size(); j++){
    //	  dR1=deltaR(allObjects[keys[j]].eta(),allObjects[keys[j]].phi(),_electron_1_eta,_electron_1_phi);
    //	  if (_nSelElectrons>1) dR2=deltaR(allObjects[keys[j]].eta(),allObjects[keys[j]].phi(),_electron_2_eta,_electron_2_phi);
    //	  if ( (dR1 < 0.5) || (dR2 < 0.5) ) {
    //	    if (dR1 < 0.5 && dR2 < 0.5 && dR1<dR2 ) _electron_1_hltmatched = 1;
    //	    else if (dR1 < 0.5 && dR2 < 0.5 && dR2<dR1 ) _electron_2_hltmatched = 1;
    //	    else if (dR1 < 0.5 && dR2 >= 0.5 ) _electron_1_hltmatched =1;
    //	    else if (dR2 < 0.5 && dR1 >= 0.5 ) _electron_2_hltmatched =1;
    //	  }
    //        }
    //      }
    //    }
    //
    //    SetValue("electron_1_hltmatched",_electron_1_hltmatched);
    //    SetValue("electron_2_hltmatched",_electron_2_hltmatched);
    //
    //    double _el1el2_mass = -999.0;
    //    if (_nSelElectrons>1) {
    //      math::XYZTLorentzVector lv_el1 = vSelElectrons[0]->p4();
    //      math::XYZTLorentzVector lv_el2 = vSelElectrons[1]->p4();
    //      math::XYZTLorentzVector lv_el1el2 = lv_el1+lv_el2;
    //      _el1el2_mass = lv_el1el2.M();
    //    }
    //    SetValue("el1el2_mass",_el1el2_mass);
    
    // MET
    double _met = -std::numeric_limits<double>::max();
    double _met_phi = -std::numeric_limits<double>::max();
    math::XYZTLorentzVector lv_met;
    if(pMet.isNonnull() && pMet.isAvailable()) {
        lv_met = pMet->p4();
        _met = lv_met.pt();
        _met_phi = lv_met.phi();
    }
    SetValue("met", _met);
    SetValue("met_phi", _met_phi);
    
    double _type1corrmet = -std::numeric_limits<double>::max();
    double _type1corrmet_phi = -std::numeric_limits<double>::max();
    math::XYZTLorentzVector lv_type1corrmet;
    if(pType1CorrMet.isNonnull() && pType1CorrMet.isAvailable()) {
        lv_type1corrmet = pType1CorrMet->p4();
        _type1corrmet = lv_type1corrmet.pt();
        _type1corrmet_phi = lv_type1corrmet.phi();
    }
    SetValue("type1corrmet", _type1corrmet);
    SetValue("type1corrmet_phi", _type1corrmet_phi);
    
    // Corrected MET
    Double_t _corr_met = -std::numeric_limits<double>::max();
    Double_t _corr_met_phi = -std::numeric_limits<double>::max();
    
    if(corrMET.Pt()>0) {
        _corr_met = corrMET.Pt();
        _corr_met_phi = corrMET.Phi();
    }
    
    SetValue("corr_met", _corr_met);
    SetValue("corr_met_phi", _corr_met_phi);
    
    // Store Jet Vectors for the selected jets and the bjets
    vector<Double_t> jet_pt;
    vector<Double_t> jet_eta;
    vector<Double_t> jet_phi;
    vector<int> jet_tag;
    vector<Double_t> bjet_pt;
    vector<Double_t> bjet_eta;
    vector<Double_t> bjet_phi;
    int nbtags = 0;
    
    
    for(int i=0; i<_nCorrBtagJets; i++){
        
        jet_pt.push_back(vCorrBtagJets[i].first.Pt());
        jet_eta.push_back(vCorrBtagJets[i].first.Eta());
        jet_phi.push_back(vCorrBtagJets[i].first.Phi());
        jet_tag.push_back(int(vCorrBtagJets[i].second));
        
        if( vCorrBtagJets[i].second ){
            nbtags++;
            bjet_pt.push_back(vCorrBtagJets[i].first.Pt());
            bjet_eta.push_back(vCorrBtagJets[i].first.Eta());
            bjet_phi.push_back(vCorrBtagJets[i].first.Phi());
        }
    }
    
    SetValue("jet_pt", jet_pt);
    SetValue("jet_eta", jet_eta);
    SetValue("jet_phi", jet_phi);
    SetValue("jet_tag", jet_tag);
    SetValue("bjet_pt", bjet_pt);
    SetValue("bjet_eta", bjet_eta);
    SetValue("bjet_phi", bjet_phi);
    SetValue("n_btags", nbtags);
    
    
    vector<int> jet_flavor;
    int nBjets = 0;
    int nCjets = 0;
    
    for(int i=0; i<_nSelJets; i++){
        jet_flavor.push_back( abs(vSelJets[i]->partonFlavour()) );
        if( abs(vSelJets[i]->partonFlavour())==5 ) nBjets++;
        if( abs(vSelJets[i]->partonFlavour())==4 ) nCjets++;
    }
    
    SetValue("n_Bjets", nBjets);
    SetValue("n_Cjets", nCjets);
    SetValue("jet_flavor", jet_flavor);
    
    //Event weights
    
    double _weight_WJets = 1.0;
    if (isWJets_) {
        if ( nBjets > 0)  _weight_WJets = 1.21*0.92;
        if ( nBjets == 0 and nCjets > 0) _weight_WJets = 1.66*0.92;
        if ( nBjets == 0 and nCjets == 0) _weight_WJets = 1.0*0.85;
    } else {
        _weight_WJets = 1.0;
    }
    SetValue("weight_WJets", _weight_WJets);
    
    
    double _weight_muon_eff = 1.0;
    //A+B+C Tight ID*ISO(dBeta<0.12)*TRIG(IsoMu24_eta2p1)
    //A+B = 5.32 C = 0.495
    if (_nSelMuons>0) {
        if (abs(_muon_1_eta) <= 0.9 ) _weight_muon_eff = 0.9925*1.0000*0.9815;
        if ( 0.9 < abs(_muon_1_eta) && abs(_muon_1_eta) <= 1.2 ) _weight_muon_eff = 0.9928*1.0006*0.9651;
        if ( 1.2 < abs(_muon_1_eta) && abs(_muon_1_eta) <= 2.1 ) _weight_muon_eff = 0.9960*1.0006*0.9968;
    }
    SetValue("weight_MuonEff", _weight_muon_eff);
    
    double _weight_electron_eff = 1.0;
    // Cut based Tight ID scale factors
    if (_nSelElectrons>0) {
        if ( abs(_electron_1_eta) <= 0.8 ) {
            if (_electron_1_pt > 30.0 && _electron_1_pt <= 40.0 ) _weight_electron_eff = 0.978;
            else if (_electron_1_pt > 40.0 && _electron_1_pt <= 50.0 ) _weight_electron_eff = 0.981;
            else if (_electron_1_pt > 50.0 ) _weight_electron_eff = 0.982;
        }
        else if (0.8 < abs(_electron_1_eta) && abs(_electron_1_eta) <= 1.442 ) {
            if (_electron_1_pt > 30.0 && _electron_1_pt <= 40.0 ) _weight_electron_eff = 0.958;
            else if (_electron_1_pt > 40.0 && _electron_1_pt <= 50.0 ) _weight_electron_eff = 0.969;
            else if (_electron_1_pt > 50.0 ) _weight_electron_eff = 0.969;
        }
        else if ( 1.442 < abs(_electron_1_eta) && abs(_electron_1_eta) <= 1.556 ) {
            if (_electron_1_pt > 30.0 && _electron_1_pt <= 40.0 ) _weight_electron_eff = 0.907;
            else if (_electron_1_pt > 40.0 && _electron_1_pt <= 50.0 ) _weight_electron_eff = 0.904;
            else if (_electron_1_pt > 50.0 ) _weight_electron_eff = 0.926;
        }
        else if ( 1.556 < abs(_electron_1_eta) && abs(_electron_1_eta) <=  2.0 ) {
            if (_electron_1_pt > 30.0 && _electron_1_pt <= 40.0 ) _weight_electron_eff = 0.909;
            else if (_electron_1_pt > 40.0 && _electron_1_pt <= 50.0 ) _weight_electron_eff = 0.942;
            else if (_electron_1_pt > 50.0 ) _weight_electron_eff = 0.957;
        }
        else if ( 2.0 < abs(_electron_1_eta) && abs(_electron_1_eta) < 2.5 ) {
            if (_electron_1_pt > 30.0 && _electron_1_pt <= 40.0 ) _weight_electron_eff = 0.987;
            else if (_electron_1_pt > 40.0 && _electron_1_pt <= 50.0 ) _weight_electron_eff = 0.991;
            else if (_electron_1_pt > 50.0 ) _weight_electron_eff = 0.999;
        }
    }
    SetValue("weight_ElectronEff", _weight_electron_eff);
    
    double _weight_electron_eff_53x = 1.0;
    // Cut based Tight ID scale factors for 53x
    if (_nSelElectrons>0) {
        if ( abs(_electron_1_eta) <= 0.8 ) {
            if (_electron_1_pt > 30.0 && _electron_1_pt <= 40.0 ) _weight_electron_eff_53x = 0.988;
            else if (_electron_1_pt > 40.0 && _electron_1_pt <= 50.0 ) _weight_electron_eff_53x = 0.989;
            else if (_electron_1_pt > 50.0 ) _weight_electron_eff_53x = 0.988;
        }
        else if (0.8 < abs(_electron_1_eta) && abs(_electron_1_eta) <= 1.442 ) {
            if (_electron_1_pt > 30.0 && _electron_1_pt <= 40.0 ) _weight_electron_eff_53x = 0.971;
            else if (_electron_1_pt > 40.0 && _electron_1_pt <= 50.0 ) _weight_electron_eff_53x = 0.977;
            else if (_electron_1_pt > 50.0 ) _weight_electron_eff_53x = 0.981;
        }
        else if ( 1.442 < abs(_electron_1_eta) && abs(_electron_1_eta) <= 1.556 ) {
            if (_electron_1_pt > 30.0 && _electron_1_pt <= 40.0 ) _weight_electron_eff_53x = 0.995;
            else if (_electron_1_pt > 40.0 && _electron_1_pt <= 50.0 ) _weight_electron_eff_53x = 0.963;
            else if (_electron_1_pt > 50.0 ) _weight_electron_eff_53x = 0.974;
        }
        else if ( 1.556 < abs(_electron_1_eta) && abs(_electron_1_eta) <=  2.0 ) {
            if (_electron_1_pt > 30.0 && _electron_1_pt <= 40.0 ) _weight_electron_eff_53x = 0.977;
            else if (_electron_1_pt > 40.0 && _electron_1_pt <= 50.0 ) _weight_electron_eff_53x = 0.994;
            else if (_electron_1_pt > 50.0 ) _weight_electron_eff_53x = 0.991;
        }
        else if ( 2.0 < abs(_electron_1_eta) && abs(_electron_1_eta) < 2.5 ) {
            if (_electron_1_pt > 30.0 && _electron_1_pt <= 40.0 ) _weight_electron_eff_53x = 1.027;
            else if (_electron_1_pt > 40.0 && _electron_1_pt <= 50.0 ) _weight_electron_eff_53x = 1.025;
            else if (_electron_1_pt > 50.0 ) _weight_electron_eff_53x = 1.015;
        }
    }
    SetValue("weight_ElectronEff_53x", _weight_electron_eff_53x);
    
    // Generator information
    double _genTopMass = -1.0;
    double _genTopPt = -1.0;
    double _genDrLeptonTopBjet = -std::numeric_limits<double>::max();
    double _genTopBjetPt = -1.0;
    double _genTopBjetEta = -1.0;
    double _genTopBjetPhi = -1.0;
    
    if (isTB_) {
        edm::Handle<reco::GenParticleCollection> genParticles;
        edm::InputTag genParticles_it = edm::InputTag("prunedGenParticles");
        event.getByLabel(genParticles_it, genParticles);
        
        int qLep = 0;
        math::XYZTLorentzVector lv_genLep;
        math::XYZTLorentzVector lv_genNu;
        math::XYZTLorentzVector lv_genB;
        math::XYZTLorentzVector lv_genBbar;
        
        for (size_t i = 0; i < genParticles->size(); i++) {
            const reco::GenParticle & p = (*genParticles).at(i);
            if (p.status() == 3) {
                if (fabs(p.pdgId())==11 or fabs(p.pdgId())==13) {
                    lv_genLep = p.p4();
                    if (p.pdgId()<0) qLep = 1;
                    else qLep = -1;
                }
                if (fabs(p.pdgId())==12 or fabs(p.pdgId())==14) lv_genNu = p.p4();
                if (p.pdgId()==5) lv_genB = p.p4();
                if (p.pdgId()==-5) lv_genBbar = p.p4();
            }
        }
        
        math::XYZTLorentzVector lv_genTop;
        math::XYZTLorentzVector lv_genOtherTop;
        
        if (qLep > 0) {
            lv_genTop = lv_genLep + lv_genNu + lv_genB;
            lv_genOtherTop = lv_genLep + lv_genNu + lv_genBbar;
        }
        else {
            lv_genTop = lv_genLep + lv_genNu + lv_genBbar;
            lv_genOtherTop = lv_genLep + lv_genNu + lv_genB;
        }
        
        
        
        _genTopMass = lv_genTop.M();
        _genTopPt = lv_genTop.Pt();
        _genTopBjetPt = (qLep<0?lv_genBbar.Pt():lv_genB.Pt());
        _genTopBjetEta = (qLep<0?lv_genBbar.Eta():lv_genB.Eta());
        _genTopBjetPhi = (qLep<0?lv_genBbar.Phi():lv_genB.Phi());
        
        double deta = -std::numeric_limits<double>::max();
        double dphi = -std::numeric_limits<double>::max();
        
        if (qLep<0) {
            deta = lv_genLep.Eta()-lv_genBbar.Eta();
            dphi = lv_genLep.Phi()-lv_genBbar.Phi();
            if ( dphi > TMath::Pi() ) dphi -= 2.*TMath::Pi();
            if ( dphi <= -TMath::Pi() ) dphi += 2.*TMath::Pi();
            _genDrLeptonTopBjet = TMath::Sqrt(deta*deta + dphi*dphi);
        } else {
            deta = lv_genLep.Eta()-lv_genB.Eta();
            dphi = lv_genLep.Phi()-lv_genB.Phi();
            if ( dphi > TMath::Pi() ) dphi -= 2.*TMath::Pi();
            if ( dphi <= -TMath::Pi() ) dphi += 2.*TMath::Pi();
            _genDrLeptonTopBjet = TMath::Sqrt(deta*deta + dphi*dphi);
        }
        
    }
    
    SetValue("genTopMass", _genTopMass);
    SetValue("genTopPt", _genTopPt);
    SetValue("genDrLeptonTopBjet", _genDrLeptonTopBjet);
    SetValue("genTopBjetPt", _genTopBjetPt);
    SetValue("genTopBjetEta", _genTopBjetEta);
    SetValue("genTopBjetPhi", _genTopBjetPhi);
    
    double _genTTMass = -1.0;
    
    if (isTT_) {
        edm::Handle<reco::GenParticleCollection> genParticles;
        edm::InputTag genParticles_it = edm::InputTag("prunedGenParticles");
        event.getByLabel(genParticles_it, genParticles);
        
        math::XYZTLorentzVector lv_genT;
        math::XYZTLorentzVector lv_genTbar;
        
        for (size_t i = 0; i < genParticles->size(); i++) {
            const reco::GenParticle & p = (*genParticles).at(i);
            if (p.pdgId()==6) lv_genT = p.p4();
            if (p.pdgId()==-6) lv_genTbar = p.p4();
        }
        math::XYZTLorentzVector lv_genTT =  lv_genT + lv_genTbar;
        _genTTMass = lv_genTT.M();
    }
    
    SetValue("genTTMass", _genTTMass);
    
    
    return 0;
}

int singleLepCalc::EndJob()
{
    return 0;
}
