/*
 // Calculator for charged Higgs to tb analysis
 //
 // Adapted by John Stupak from:
 //   WprimeEventSelector by David Sperka
 */


#include <iostream>
#include <math.h> /* M_PI and sqrt */

#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "TLorentzVector.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
//#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Math/interface/deltaR.h"

class LjmetFactory;


class ChargedHiggsCalc : public BaseCalc{
    
public:
    
    ChargedHiggsCalc();
    virtual ~ChargedHiggsCalc(){}
    
    virtual int BeginJob(){
        
        if (mPset.exists("triggerSummary")) triggerSummary_ = mPset.getParameter<edm::InputTag>("triggerSummary");
        else                                triggerSummary_ = edm::InputTag("hltTriggerSummaryAOD");
        
        if (mPset.exists("rhoSrc")) rhoSrc_ = mPset.getParameter<edm::InputTag>("rhoSrc");
        else                        rhoSrc_ = edm::InputTag("kt6PFJets", "rho");
        
        if (mPset.exists("isTB"))    isTB_ = mPset.getParameter<bool>("isTB");
        else                         isTB_ = false;
        
        if (mPset.exists("isTT"))    isTT_ = mPset.getParameter<bool>("isTT");
        else                         isTT_ = false;
        
        if (mPset.exists("isWJets")) isWJets_ = mPset.getParameter<bool>("isWJets");
        else                         isWJets_ = false;
        
        
        
        return 0;
    }
    
    virtual int AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector);
    virtual int EndJob(){return 0;}
    
    
private:
    
    edm::InputTag triggerSummary_;
    edm::InputTag rhoSrc_;
    bool isWJets_;
    bool isTB_;
    bool isTT_;
    
};



//static int reg = LjmetFactory::GetInstance()->Register(new ChargedHiggsCalc(), "ChargedHiggsCalc");



ChargedHiggsCalc::ChargedHiggsCalc(){
}

int ChargedHiggsCalc::AnalyzeEvent(edm::EventBase const & event,
                                   BaseEventSelector * selector){
    //
    // compute event variables here
    //
    
    //
    // _____ Get objects from the selector _____________________
    //
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
    
    
    //
    //_____Event kinematics _____________________
    //
    
    int _nAllJets      = (int)vAllJets.size();
    int _nSelJets      = (int)vSelJets.size();
    int _nSelBtagJets  = (int)vSelBtagJets.size();
    int _nCorrBtagJets  = (int)vCorrBtagJets.size();
    int _nSelMuons     = (int)vSelMuons.size();
    int _nLooseMuons   = (int)vLooseMuons.size();
    int _nSelElectrons = (int)vSelElectrons.size();
    int _nLooseElectrons = (int)vLooseElectrons.size();
    int _nSelPVs = (int)vSelPVs.size();
    
    SetValue("nPV",_nSelPVs);
    
    // Muon
    double _muon_1_pt = -9999.0;
    double _muon_1_phi = -9999.0;
    double _muon_1_eta = -9999.0;
    double _muon_1_RelIso = -9999.0;
    
    if (_nSelMuons>0) {
        _muon_1_pt = vSelMuons[0]->pt();
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
    
    double _electron_1_pt = -9999.0;
    double _electron_1_phi = -9999.0;
    double _electron_1_eta = -9999.0;
    double _electron_1_RelIso = -9999.0;
    
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
    
    double _electron_2_pt = -9999.0;
    double _electron_2_phi = -9999.0;
    double _electron_2_eta = -9999.0;
    double _electron_2_RelIso = -9999.0;
    
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
    edm::Handle<trigger::TriggerEvent> mhEdmTriggerEvent;
    event.getByLabel(triggerSummary_,mhEdmTriggerEvent);
    trigger::TriggerObjectCollection allObjects = mhEdmTriggerEvent->getObjects();
    
    int _electron_1_hltmatched =0;
    int _electron_2_hltmatched =0;
    
    if (_nSelElectrons>0) {
        for(int i=0; i<mhEdmTriggerEvent->sizeFilters(); i++){
            if( mhEdmTriggerEvent->filterTag(i).label()!="hltEle27WP80TrackIsoFilter") continue;
            trigger::Keys keys = mhEdmTriggerEvent->filterKeys(i);
            double dR1 = 999.0;
            double dR2 = 999.0;
            for(size_t j=0; j<keys.size(); j++){
                dR1=deltaR(allObjects[keys[j]].eta(),allObjects[keys[j]].phi(),_electron_1_eta,_electron_1_phi);
                if (_nSelElectrons>1) dR2=deltaR(allObjects[keys[j]].eta(),allObjects[keys[j]].phi(),_electron_2_eta,_electron_2_phi);
                if ( (dR1 < 0.5) || (dR2 < 0.5) ) {
                    if (dR1 < 0.5 && dR2 < 0.5 && dR1<dR2 ) _electron_1_hltmatched = 1;
                    else if (dR1 < 0.5 && dR2 < 0.5 && dR2<dR1 ) _electron_2_hltmatched = 1;
                    else if (dR1 < 0.5 && dR2 >= 0.5 ) _electron_1_hltmatched =1;
                    else if (dR2 < 0.5 && dR1 >= 0.5 ) _electron_2_hltmatched =1;
                }
            }
        }
    }
    
    SetValue("electron_1_hltmatched",_electron_1_hltmatched);
    SetValue("electron_2_hltmatched",_electron_2_hltmatched);
    
    double _el1el2_mass = -999.0;
    if (_nSelElectrons>1) {
        math::XYZTLorentzVector lv_el1 = vSelElectrons[0]->p4();
        math::XYZTLorentzVector lv_el2 = vSelElectrons[1]->p4();
        math::XYZTLorentzVector lv_el1el2 = lv_el1+lv_el2;
        _el1el2_mass = lv_el1el2.M();
    }
    SetValue("el1el2_mass",_el1el2_mass);
    
    // MET
    double _met = -9999.0;
    double _met_phi = -9999.0;
    math::XYZTLorentzVector lv_met;
    if(pMet.isNonnull() && pMet.isAvailable()) {
        lv_met = pMet->p4();
        _met = lv_met.pt();
        _met_phi = lv_met.phi();
    }
    SetValue("met", _met);
    SetValue("met_phi", _met_phi);
    
    double _type1corrmet = -9999.0;
    double _type1corrmet_phi = -9999.0;
    math::XYZTLorentzVector lv_type1corrmet;
    if(pType1CorrMet.isNonnull() && pType1CorrMet.isAvailable()) {
        lv_type1corrmet = pType1CorrMet->p4();
        _type1corrmet = lv_type1corrmet.pt();
        _type1corrmet_phi = lv_type1corrmet.phi();
    }
    SetValue("type1corrmet", _type1corrmet);
    SetValue("type1corrmet_phi", _type1corrmet_phi);
    
    // Corrected MET
    Double_t _corr_met = -9999.0;
    Double_t _corr_met_phi = -9999.0;
    
    if(corrMET.Pt()>0) {
        _corr_met = corrMET.Pt();
        _corr_met_phi = corrMET.Phi();
    }
    
    SetValue("corr_met", _corr_met);
    SetValue("corr_met_phi", _corr_met_phi);
    
    Double_t _jet_0_pt = -9999.0, _jet_0_eta = -9999.0, _jet_0_phi = -9999.0;
    Double_t _jet_1_pt = -9999.0, _jet_1_eta = -9999.0, _jet_1_phi = -9999.0;
    Double_t _jet_2_pt = -9999.0, _jet_2_eta = -9999.0, _jet_2_phi = -9999.0;
    Double_t _jet_3_pt = -9999.0, _jet_3_eta = -9999.0, _jet_3_phi = -9999.0;
    Double_t _jet_4_pt = -9999.0, _jet_4_eta = -9999.0, _jet_4_phi = -9999.0;
    Double_t _jet_5_pt = -9999.0, _jet_5_eta = -9999.0, _jet_5_phi = -9999.0;
    Double_t _jet_6_pt = -9999.0, _jet_6_eta = -9999.0, _jet_6_phi = -9999.0;
    Double_t _jet_7_pt = -9999.0, _jet_7_eta = -9999.0, _jet_7_phi = -9999.0;
    Double_t _jet_8_pt = -9999.0, _jet_8_eta = -9999.0, _jet_8_phi = -9999.0;
    Double_t _jet_9_pt = -9999.0, _jet_9_eta = -9999.0, _jet_9_phi = -9999.0;
    
    bool _jet_0_tag = false;
    bool _jet_1_tag = false;
    bool _jet_2_tag = false;
    bool _jet_3_tag = false;
    bool _jet_4_tag = false;
    bool _jet_5_tag = false;
    bool _jet_6_tag = false;
    bool _jet_7_tag = false;
    bool _jet_8_tag = false;
    bool _jet_9_tag = false;
    
    if (_nCorrBtagJets>0) {
        _jet_0_pt = vCorrBtagJets[0].first.Pt();
        _jet_0_eta = vCorrBtagJets[0].first.Eta();
        _jet_0_phi = vCorrBtagJets[0].first.Phi();
        _jet_0_tag = vCorrBtagJets[0].second;
    }
    if (_nCorrBtagJets>1) {
        _jet_1_pt = vCorrBtagJets[1].first.Pt();
        _jet_1_eta = vCorrBtagJets[1].first.Eta();
        _jet_1_phi = vCorrBtagJets[1].first.Phi();
        _jet_1_tag = vCorrBtagJets[1].second;
    }
    if (_nCorrBtagJets>2) {
        _jet_2_pt = vCorrBtagJets[2].first.Pt();
        _jet_2_eta = vCorrBtagJets[2].first.Eta();
        _jet_2_phi = vCorrBtagJets[2].first.Phi();
        _jet_2_tag = vCorrBtagJets[2].second;
    }
    if (_nCorrBtagJets>3) {
        _jet_3_pt = vCorrBtagJets[3].first.Pt();
        _jet_3_eta = vCorrBtagJets[3].first.Eta();
        _jet_3_phi = vCorrBtagJets[3].first.Phi();
        _jet_3_tag = vCorrBtagJets[3].second;
    }
    if (_nCorrBtagJets>4) {
        _jet_4_pt = vCorrBtagJets[4].first.Pt();
        _jet_4_eta = vCorrBtagJets[4].first.Eta();
        _jet_4_phi = vCorrBtagJets[4].first.Phi();
        _jet_4_tag = vCorrBtagJets[4].second;
    }
    if (_nCorrBtagJets>5) {
        _jet_5_pt = vCorrBtagJets[5].first.Pt();
        _jet_5_eta = vCorrBtagJets[5].first.Eta();
        _jet_5_phi = vCorrBtagJets[5].first.Phi();
        _jet_5_tag = vCorrBtagJets[5].second;
    }
    if (_nCorrBtagJets>6) {
        _jet_6_pt = vCorrBtagJets[6].first.Pt();
        _jet_6_eta = vCorrBtagJets[6].first.Eta();
        _jet_6_phi = vCorrBtagJets[6].first.Phi();
        _jet_6_tag = vCorrBtagJets[6].second;
    }
    if (_nCorrBtagJets>7) {
        _jet_7_pt = vCorrBtagJets[7].first.Pt();
        _jet_7_eta = vCorrBtagJets[7].first.Eta();
        _jet_7_phi = vCorrBtagJets[7].first.Phi();
        _jet_7_tag = vCorrBtagJets[7].second;
    }
    if (_nCorrBtagJets>8) {
        _jet_8_pt = vCorrBtagJets[8].first.Pt();
        _jet_8_eta = vCorrBtagJets[8].first.Eta();
        _jet_8_phi = vCorrBtagJets[8].first.Phi();
        _jet_8_tag = vCorrBtagJets[8].second;
    }
    if (_nCorrBtagJets>9) {
        _jet_9_pt = vCorrBtagJets[9].first.Pt();
        _jet_9_eta = vCorrBtagJets[9].first.Eta();
        _jet_9_phi = vCorrBtagJets[9].first.Phi();
        _jet_9_tag = vCorrBtagJets[9].second;
    }
    
    SetValue("jet_0_pt", _jet_0_pt);
    SetValue("jet_1_pt", _jet_1_pt);
    SetValue("jet_2_pt", _jet_2_pt);
    SetValue("jet_3_pt", _jet_3_pt);
    SetValue("jet_4_pt", _jet_4_pt);
    SetValue("jet_5_pt", _jet_5_pt);
    SetValue("jet_6_pt", _jet_6_pt);
    SetValue("jet_7_pt", _jet_7_pt);
    SetValue("jet_8_pt", _jet_8_pt);
    SetValue("jet_9_pt", _jet_9_pt);
    
    SetValue("jet_0_eta", _jet_0_eta);
    SetValue("jet_1_eta", _jet_1_eta);
    SetValue("jet_2_eta", _jet_2_eta);
    SetValue("jet_3_eta", _jet_3_eta);
    SetValue("jet_4_eta", _jet_4_eta);
    SetValue("jet_5_eta", _jet_5_eta);
    SetValue("jet_6_eta", _jet_6_eta);
    SetValue("jet_7_eta", _jet_7_eta);
    SetValue("jet_8_eta", _jet_8_eta);
    SetValue("jet_9_eta", _jet_9_eta);
    
    SetValue("jet_0_phi", _jet_0_phi);
    SetValue("jet_1_phi", _jet_1_phi);
    SetValue("jet_2_phi", _jet_2_phi);
    SetValue("jet_3_phi", _jet_3_phi);
    SetValue("jet_4_phi", _jet_4_phi);
    SetValue("jet_5_phi", _jet_5_phi);
    SetValue("jet_6_phi", _jet_6_phi);
    SetValue("jet_7_phi", _jet_7_phi);
    SetValue("jet_8_phi", _jet_8_phi);
    SetValue("jet_9_phi", _jet_9_phi);
    
    SetValue("jet_0_tag", _jet_0_tag);
    SetValue("jet_1_tag", _jet_1_tag);
    SetValue("jet_2_tag", _jet_2_tag);
    SetValue("jet_3_tag", _jet_3_tag);
    SetValue("jet_4_tag", _jet_4_tag);
    SetValue("jet_5_tag", _jet_5_tag);
    SetValue("jet_6_tag", _jet_6_tag);
    SetValue("jet_7_tag", _jet_7_tag);
    SetValue("jet_8_tag", _jet_8_tag);
    SetValue("jet_9_tag", _jet_9_tag);
    
    SetValue("n_btags", (_jet_0_tag+_jet_1_tag+_jet_2_tag+_jet_3_tag+_jet_4_tag+_jet_5_tag+_jet_6_tag+_jet_7_tag+_jet_8_tag+_jet_9_tag));
    
    int _jet_0_flavor = -999;
    int _jet_1_flavor = -999;
    int _jet_2_flavor = -999;
    int _jet_3_flavor = -999;
    int _jet_4_flavor = -999;
    int _jet_5_flavor = -999;
    int _jet_6_flavor = -999;
    int _jet_7_flavor = -999;
    int _jet_8_flavor = -999;
    int _jet_9_flavor = -999;
    
    if (_nSelJets>0) _jet_0_flavor = abs(vSelJets[0]->partonFlavour());
    if (_nSelJets>1) _jet_1_flavor = abs(vSelJets[1]->partonFlavour());
    if (_nSelJets>2) _jet_2_flavor = abs(vSelJets[2]->partonFlavour());
    if (_nSelJets>3) _jet_3_flavor = abs(vSelJets[3]->partonFlavour());
    if (_nSelJets>4) _jet_4_flavor = abs(vSelJets[4]->partonFlavour());
    if (_nSelJets>5) _jet_5_flavor = abs(vSelJets[5]->partonFlavour());
    if (_nSelJets>6) _jet_6_flavor = abs(vSelJets[6]->partonFlavour());
    if (_nSelJets>7) _jet_7_flavor = abs(vSelJets[7]->partonFlavour());
    if (_nSelJets>8) _jet_8_flavor = abs(vSelJets[8]->partonFlavour());
    if (_nSelJets>9) _jet_9_flavor = abs(vSelJets[9]->partonFlavour());
    
    SetValue("jet_0_flavor", _jet_0_flavor);
    SetValue("jet_1_flavor", _jet_1_flavor);
    SetValue("jet_2_flavor", _jet_2_flavor);
    SetValue("jet_3_flavor", _jet_3_flavor);
    SetValue("jet_4_flavor", _jet_4_flavor);
    SetValue("jet_5_flavor", _jet_5_flavor);
    SetValue("jet_6_flavor", _jet_6_flavor);
    SetValue("jet_7_flavor", _jet_7_flavor);
    SetValue("jet_8_flavor", _jet_8_flavor);
    SetValue("jet_9_flavor", _jet_9_flavor);
    
    int nBjets = 0;
    int nCjets = 0;
    
    if (_jet_0_flavor == 5) nBjets += 1;
    else if (_jet_0_flavor == 4) nCjets +=1;
    if (_jet_1_flavor == 5) nBjets += 1;
    else if (_jet_1_flavor == 4) nCjets +=1;
    if (_jet_2_flavor == 5) nBjets += 1;
    else if (_jet_2_flavor == 4) nCjets +=1;
    if (_jet_3_flavor == 5) nBjets += 1;
    else if (_jet_3_flavor == 4) nCjets +=1;
    if (_jet_4_flavor == 5) nBjets += 1;
    else if (_jet_4_flavor == 4) nCjets +=1;
    if (_jet_5_flavor == 5) nBjets += 1;
    else if (_jet_5_flavor == 4) nCjets +=1;
    if (_jet_6_flavor == 5) nBjets += 1;
    else if (_jet_6_flavor == 4) nCjets +=1;
    if (_jet_7_flavor == 5) nBjets += 1;
    else if (_jet_7_flavor == 4) nCjets +=1;
    if (_jet_8_flavor == 5) nBjets += 1;
    else if (_jet_8_flavor == 4) nCjets +=1;
    if (_jet_9_flavor == 5) nBjets += 1;
    else if (_jet_9_flavor == 4) nCjets +=1;
    
    SetValue("n_Bjets", nBjets);
    SetValue("n_Cjets", nCjets);
    
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
        if (abs(_muon_1_eta) <= 0.9 ) _weight_muon_eff = 0.994*0.993*0.976;
        if ( 0.9 < abs(_muon_1_eta) && abs(_muon_1_eta) <= 1.2 ) _weight_muon_eff = 0.992*0.998*0.961;
        if ( 1.2 < abs(_muon_1_eta) && abs(_muon_1_eta) <= 2.1 ) _weight_muon_eff = 0.998*1.002*0.983;
    }
    SetValue("weight_MuonEff", _weight_muon_eff);
    
    double _weight_electron_eff = 1.0;
    // Cut based Tight ID scale factors
    if (_nSelElectrons>0) {
        if ( abs(_electron_1_eta) <= 0.8 ) {
            if (_electron_1_pt > 30.0 && _electron_1_pt <= 40.0 ) _weight_electron_eff = 1.003;
            else if (_electron_1_pt > 40.0 && _electron_1_pt <= 50.0 ) _weight_electron_eff = 0.998;
            else if (_electron_1_pt > 50.0 ) _weight_electron_eff = 0.987;
        }
        else if (0.8 < abs(_electron_1_eta) && abs(_electron_1_eta) <= 1.442 ) {
            if (_electron_1_pt > 30.0 && _electron_1_pt <= 40.0 ) _weight_electron_eff = 0.986;
            else if (_electron_1_pt > 40.0 && _electron_1_pt <= 50.0 ) _weight_electron_eff = 0.990;
            else if (_electron_1_pt > 50.0 ) _weight_electron_eff = 0.979;
        }
        else if ( 1.442 < abs(_electron_1_eta) && abs(_electron_1_eta) <= 1.556 ) {
            if (_electron_1_pt > 30.0 && _electron_1_pt <= 40.0 ) _weight_electron_eff = 1.013;
            else if (_electron_1_pt > 40.0 && _electron_1_pt <= 50.0 ) _weight_electron_eff = 0.945;
            else if (_electron_1_pt > 50.0 ) _weight_electron_eff = 0.956;
        }
        else if ( 1.556 < abs(_electron_1_eta) && abs(_electron_1_eta) <=  2.0 ) {
            if (_electron_1_pt > 30.0 && _electron_1_pt <= 40.0 ) _weight_electron_eff = 1.000;
            else if (_electron_1_pt > 40.0 && _electron_1_pt <= 50.0 ) _weight_electron_eff = 1.005;
            else if (_electron_1_pt > 50.0 ) _weight_electron_eff = 0.999;
        }
        else if ( 2.0 < abs(_electron_1_eta) && abs(_electron_1_eta) < 2.5 ) {
            if (_electron_1_pt > 30.0 && _electron_1_pt <= 40.0 ) _weight_electron_eff = 1.091;
            else if (_electron_1_pt > 40.0 && _electron_1_pt <= 50.0 ) _weight_electron_eff = 1.055;
            else if (_electron_1_pt > 50.0 ) _weight_electron_eff = 1.035;
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
    
    double _genTopMass = -1.0;
    double _genTopPt = -1.0;
    double _genDrLeptonTopBjet = -9999.0;
    double _genChargedHiggsMass = -1.0;
    double _genChargedHiggsBjetPt = -1.0;
    double _genChargedHiggsBjetEta = -1.0;
    double _genChargedHiggsBjetPhi = -1.0;
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
        math::XYZTLorentzVector lv_genChargedHiggs;
        
        if (qLep > 0) {
            lv_genTop = lv_genLep + lv_genNu + lv_genB;
            lv_genOtherTop = lv_genLep + lv_genNu + lv_genBbar;
        }
        else {
            lv_genTop = lv_genLep + lv_genNu + lv_genBbar;
            lv_genOtherTop = lv_genLep + lv_genNu + lv_genB;
        }
        
        lv_genChargedHiggs = lv_genLep + lv_genNu + lv_genB + lv_genBbar;
        
        std::cout<<"genTopMass = "<<lv_genTop.M()<<" genOtherTopMass = " <<lv_genOtherTop.M()<<" genChargedHiggsMass = "<<lv_genChargedHiggs.M()<<std::endl;
        
        _genTopMass = lv_genTop.M();
        _genTopPt = lv_genTop.Pt();
        _genChargedHiggsMass = lv_genChargedHiggs.M();
        _genChargedHiggsBjetPt = (qLep>0?lv_genBbar.Pt():lv_genB.Pt());
        _genChargedHiggsBjetEta = (qLep>0?lv_genBbar.Eta():lv_genB.Eta());
        _genChargedHiggsBjetPhi = (qLep>0?lv_genBbar.Phi():lv_genB.Phi());
        _genTopBjetPt = (qLep<0?lv_genBbar.Pt():lv_genB.Pt());
        _genTopBjetEta = (qLep<0?lv_genBbar.Eta():lv_genB.Eta());
        _genTopBjetPhi = (qLep<0?lv_genBbar.Phi():lv_genB.Phi());
        
        double deta = -999.0;
        double dphi = -999.0;
        
        if (qLep<0) {
            deta = lv_genLep.Eta()-lv_genBbar.Eta();
            dphi = lv_genLep.Phi()-lv_genBbar.Phi();
            if ( dphi > M_PI ) dphi -= 2.*M_PI;
            if ( dphi <= -M_PI ) dphi += 2.*M_PI;
            _genDrLeptonTopBjet = sqrt(deta*deta + dphi*dphi);
        } else {
            deta = lv_genLep.Eta()-lv_genB.Eta();
            dphi = lv_genLep.Phi()-lv_genB.Phi();
            if ( dphi > M_PI ) dphi -= 2.*M_PI;
            if ( dphi <= -M_PI ) dphi += 2.*M_PI;
            _genDrLeptonTopBjet = sqrt(deta*deta + dphi*dphi);
        }
        
    }
    
    SetValue("genTopMass", _genTopMass);
    SetValue("genTopPt", _genTopPt);
    SetValue("genDrLeptonTopBjet", _genDrLeptonTopBjet);
    SetValue("genChargedHiggsMass", _genChargedHiggsMass);
    SetValue("genChargedHiggsBjetPt", _genChargedHiggsBjetPt);
    SetValue("genChargedHiggsBjetEta", _genChargedHiggsBjetEta);
    SetValue("genChargedHiggsBjetPhi", _genChargedHiggsBjetPhi);
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
