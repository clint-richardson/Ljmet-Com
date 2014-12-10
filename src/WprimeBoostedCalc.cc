/*
 Calculator for the Wprime to tb analysis
 
 Author: David Sperka, 2012
 */


#include <iostream>
#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "TLorentzVector.h"
//#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "AnalysisDataFormats/TopObjects/interface/CATopJetTagInfo.h"
#include "DataFormats/JetReco/interface/CATopJetTagInfo.h"

typedef math::XYZPoint Point;

using namespace std;

class LjmetFactory;


class WprimeBoostedCalc : public BaseCalc{
    
public:
    
    WprimeBoostedCalc();
    virtual ~WprimeBoostedCalc(){}
    
    virtual int BeginJob(){
        
        if (mPset.exists("triggerSummary")) triggerSummary_ = mPset.getParameter<edm::InputTag>("triggerSummary");
        else                                triggerSummary_ = edm::InputTag("selectedPatTrigger");
        
        if (mPset.exists("triggerCollection")) triggerCollection_ = mPset.getParameter<edm::InputTag>("triggerCollection");
        else                                triggerCollection_ = edm::InputTag("TriggerResults::HLT");
        
        
        if (mPset.exists("rhoSrc")) rhoSrc_ = mPset.getParameter<edm::InputTag>("rhoSrc");
        else                        rhoSrc_ = edm::InputTag("fixedGridRhoAll", "");
        
        if (mPset.exists("isTB"))    isTB_ = mPset.getParameter<bool>("isTB");
        else                         isTB_ = false;
        
        if (mPset.exists("isTT"))    isTT_ = mPset.getParameter<bool>("isTT");
        else                         isTT_ = false;
        
        if (mPset.exists("isWJets")) isWJets_ = mPset.getParameter<bool>("isWJets");
        else                         isWJets_ = false;
        
        if (mPset.exists("btagger")) btagger_ = mPset.getParameter<std::string>("btagger");
        else                         btagger_ = "slimmedSecondaryVertices";
        
        
        return 0;
    }
    
    virtual int AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector);
    virtual int EndJob(){return 0;}
    
    
private:
    
    edm::InputTag triggerSummary_;
    edm::InputTag triggerCollection_;
    edm::InputTag rhoSrc_;
    bool isWJets_;
    bool isTB_;
    bool isTT_;
    std::string btagger_;
    static bool my_compare(math::XYZTLorentzVector i, math::XYZTLorentzVector j){return i.Pt() > j.Pt();}
};



//static int reg = LjmetFactory::GetInstance()->Register(new WprimeBoostedCalc(), "WprimeBoostedCalc");



WprimeBoostedCalc::WprimeBoostedCalc(){
}

int WprimeBoostedCalc::AnalyzeEvent(edm::EventBase const & event,
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
    
    // Trigger
    edm::Handle<edm::TriggerResults > mhEdmTriggerResults;
    event.getByLabel( triggerCollection_ , mhEdmTriggerResults );
    //const edm::ParameterSetID ps = mhEdmTriggerResults->parameterSetID();
    const edm::TriggerNames trigNames = event.triggerNames(*mhEdmTriggerResults);
    
    unsigned int _tSize = mhEdmTriggerResults->size();
    
    int passTrigEle27v10 = -1;
    int passTrigIsoMu24v13 = -1;
    
    unsigned int _tElMCIndex = trigNames.triggerIndex("HLT_Ele30_CaloIdVT_TrkIdT_PFJet150_PFJet25_v9");
    if ( _tElMCIndex<_tSize){
        passTrigEle27v10 = mhEdmTriggerResults->accept(_tElMCIndex);
    }
    
    unsigned int _tMuMCIndex = trigNames.triggerIndex("HLT_Mu40_eta2p1_v12");
    if ( _tMuMCIndex<_tSize){
        passTrigIsoMu24v13 = mhEdmTriggerResults->accept(_tMuMCIndex);
    }
    
    //SetValue("passTrigEle27v10",passTrigEle27v10);
    //SetValue("passTrigIsoMu24v13",passTrigIsoMu24v13);
    
    // Muon
    double _muon_1_pt = -9999.0;
    double _muon_1_phi = -9999.0;
    double _muon_1_eta = -9999.0;
    double _muon_1_RelIso = -9999.0;
    
    double minDr_mu = 9999.0;
    int minDr_mu_index = -1;
    int jet_index_mu=0;
    double _ptrel_mu = -1.0;
    
    if (_nSelMuons>0) {
        _muon_1_pt = vSelMuons[0]->pt();
        _muon_1_phi = vSelMuons[0]->phi();
        _muon_1_eta = vSelMuons[0]->eta();
        
        TLorentzVector muP4;
        muP4.SetPtEtaPhiM(_muon_1_pt, _muon_1_eta, _muon_1_phi, 0.105658367);
        
        double chIso = vSelMuons[0]->userIsolation(pat::PfChargedHadronIso);
        double nhIso = vSelMuons[0]->userIsolation(pat::PfNeutralHadronIso);
        double gIso  = vSelMuons[0]->userIsolation(pat::PfGammaIso);
        double puIso = vSelMuons[0]->userIsolation(pat::PfPUChargedHadronIso);
        _muon_1_RelIso = (chIso + std::max(0.,nhIso + gIso - 0.5*puIso))/_muon_1_pt;
        
        for (unsigned int ijet = 0; ijet<vCorrBtagJets.size(); ijet++) {
            double deta = vCorrBtagJets[ijet].first.Eta()-vSelMuons[0]->eta();
            double dphi = vCorrBtagJets[ijet].first.Phi()-vSelMuons[0]->phi();
            if ( dphi > TMath::Pi() ) dphi -= 2.*TMath::Pi();
            if ( dphi <= -TMath::Pi() ) dphi += 2.*TMath::Pi();
            double dR = TMath::Sqrt(deta*deta + dphi*dphi);
            if (dR<minDr_mu) {
                minDr_mu = dR;
                minDr_mu_index=jet_index_mu;
            }
            jet_index_mu++;
        }
        
        //std::cout<<"muon px,py,pz,E = "<<muP4.Px()<<", "<<muP4.Py()<<", "<<muP4.Pz()<<", "<<muP4.E()<<", "<<std::endl;
        TVector3 p_mu(muP4.Px(),muP4.Py(),muP4.Pz());
        TVector3 p_jet(vCorrBtagJets[minDr_mu_index].first.Px(),vCorrBtagJets[minDr_mu_index].first.Py(),vCorrBtagJets[minDr_mu_index].first.Pz());
        
        double sin_alpha = (p_mu.Cross(p_jet)).Mag()/p_mu.Mag()/p_jet.Mag();
        _ptrel_mu = p_mu.Mag()*sin_alpha;
        
    }
    
    SetValue("muon_1_pt",_muon_1_pt);
    SetValue("muon_1_phi",_muon_1_phi);
    SetValue("muon_1_eta",_muon_1_eta);
    SetValue("muon_1_RelIso",_muon_1_RelIso);
    SetValue("muon_1_ptrel", _ptrel_mu);
    SetValue("muon_1_minDr", minDr_mu);
    
    // Electron
    edm::Handle<double> rhoHandle;
    event.getByLabel(rhoSrc_, rhoHandle);
    double rhoIso = std::max(*(rhoHandle.product()), 0.0);
    Point PVtx = vSelPVs[0]->position();
    
    double _electron_1_pt = -9999.0;
    double _electron_1_phi = -9999.0;
    double _electron_1_eta = -9999.0;
    double _electron_1_RelIso = -9999.0;
    
    double _electron_1_Deta  = -9999.0;
    double _electron_1_Dphi  = -9999.0;
    double _electron_1_sihih = -9999.0;
    double _electron_1_HoE   = -9999.0;
    double _electron_1_D0    = -9999.0;
    double _electron_1_DZ    = -9999.0;
    double _electron_1_Ooemoop = -9999.0;
    int _electron_1_mHits = -1;
    bool _electron_1_vtxFitConv = 0;
    
    SetValue("nPATElectrons",(int)vSelElectrons.size());
    
    double minDr_el = 9999.0;
    int minDr_el_index = -1;
    int jet_index_el=0;
    double _ptrel_el = -1.0;
    
    if (_nSelElectrons>0) {
        _electron_1_pt = vSelElectrons[0]->ecalDrivenMomentum().pt();
        _electron_1_phi = vSelElectrons[0]->phi();
        _electron_1_eta = vSelElectrons[0]->eta();
        
        TLorentzVector elP4;
        elP4.SetPtEtaPhiM(_electron_1_pt,_electron_1_eta,_electron_1_phi,0.00051099891);
        
        double AEff  = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, vSelElectrons[0]->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012);
        double chIso = vSelElectrons[0]->chargedHadronIso();
        double nhIso = vSelElectrons[0]->neutralHadronIso();
        double phIso = vSelElectrons[0]->photonIso();
        _electron_1_RelIso  = ( chIso + max(0.0, nhIso + phIso - rhoIso*AEff) )/ vSelElectrons[0]->ecalDrivenMomentum().pt();
        
        _electron_1_Deta  = vSelElectrons[0]->deltaEtaSuperClusterTrackAtVtx();
        _electron_1_Dphi  = vSelElectrons[0]->deltaPhiSuperClusterTrackAtVtx();
        _electron_1_sihih = vSelElectrons[0]->sigmaIetaIeta();
        _electron_1_HoE   = vSelElectrons[0]->hadronicOverEm();
        _electron_1_D0    = vSelElectrons[0]->dB();
        _electron_1_DZ    = vSelElectrons[0]->gsfTrack()->dz(PVtx);
        
        _electron_1_Ooemoop = (1.0/vSelElectrons[0]->ecalEnergy() - vSelElectrons[0]->eSuperClusterOverP()/vSelElectrons[0]->ecalEnergy());
        _electron_1_RelIso  = ( chIso + max(0.0, nhIso + phIso - rhoIso*AEff) )/ vSelElectrons[0]->ecalDrivenMomentum().pt();
        _electron_1_mHits   =  vSelElectrons[0]->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
        _electron_1_vtxFitConv = vSelElectrons[0]->passConversionVeto();
        
        for (unsigned int ijet = 0; ijet<vCorrBtagJets.size(); ijet++) {
            double deta = vCorrBtagJets[ijet].first.Eta()-vSelElectrons[0]->eta();
            double dphi = vCorrBtagJets[ijet].first.Phi()-vSelElectrons[0]->phi();
            if ( dphi > TMath::Pi() ) dphi -= 2.*TMath::Pi();
            if ( dphi <= -TMath::Pi() ) dphi += 2.*TMath::Pi();
            double dR = TMath::Sqrt(deta*deta + dphi*dphi);
            if (dR<minDr_el) {
                minDr_el = dR;
                minDr_el_index=jet_index_el;
            }
            jet_index_el++;
        }
        
        //std::cout<<"electron px,py,pz,E = "<<elP4.Px()<<", "<<elP4.Py()<<", "<<elP4.Pz()<<", "<<elP4.E()<<", "<<std::endl;
        TVector3 p_el(elP4.Px(),elP4.Py(),elP4.Pz());
        TVector3 p_jet(vCorrBtagJets[minDr_el_index].first.Px(),vCorrBtagJets[minDr_el_index].first.Py(),vCorrBtagJets[minDr_el_index].first.Pz());
        
        double sin_alpha = (p_el.Cross(p_jet)).Mag()/p_el.Mag()/p_jet.Mag();
        _ptrel_el = p_el.Mag()*sin_alpha;
        
    }
    
    SetValue("elec_1_pt", _electron_1_pt);
    SetValue("elec_1_phi", _electron_1_phi);
    SetValue("elec_1_eta", _electron_1_eta);
    SetValue("elec_1_RelIso", _electron_1_RelIso);
    SetValue("elec_1_ptrel", _ptrel_el);
    SetValue("elec_1_minDr", minDr_el);
    
    /*
     SetValue("elec_1_Deta", _electron_1_Deta);
     SetValue("elec_1_Dphi", _electron_1_Dphi);
     SetValue("elec_1_sihih", _electron_1_sihih);
     SetValue("elec_1_HoE", _electron_1_HoE);
     SetValue("elec_1_D0", _electron_1_D0);
     SetValue("elec_1_DZ", _electron_1_DZ);
     SetValue("elec_1_Ooemoop", _electron_1_Ooemoop);
     SetValue("elec_1_RelIso", _electron_1_RelIso);
     SetValue("elec_1_mHits", _electron_1_mHits);
     SetValue("elec_1_vtxFitConv", _electron_1_vtxFitConv);
     */
    
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
    edm::Handle<pat::TriggerObjectStandAloneCollection> mhEdmTriggerObjectColl;
    event.getByLabel(triggerSummary_,mhEdmTriggerObjectColl);
    
    const edm::TriggerNames &names = event.triggerNames(*mhEdmTriggerResults);
    
    int _electron_1_hltmatched =0;
    int _electron_2_hltmatched =0;
    
    if (_nSelElectrons>0) {
        for(pat::TriggerObjectStandAlone obj : *mhEdmTriggerObjectColl){
            obj.unpackPathNames(names);
            double dR1 = 999.0;
            double dR2 = 999.0;
            for(unsigned h = 0; h < obj.filterLabels().size(); ++h){
                if( obj.filterLabels()[h]!="hltEle30CaloIdVTTrkIdTDphiFilter") continue;
                dR1=deltaR(obj.eta(),obj.phi(),_electron_1_eta,_electron_1_phi);
                if (_nSelElectrons>1) dR2=deltaR(obj.eta(),obj.phi(),_electron_2_eta,_electron_2_phi);
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
    
    Double_t _jet_0_discrim = -999.;
    Double_t _jet_1_discrim = -999.;
    Double_t _jet_2_discrim = -999.;
    Double_t _jet_3_discrim = -999.;
    Double_t _jet_4_discrim = -999.;
    Double_t _jet_5_discrim = -999.;
    Double_t _jet_6_discrim = -999.;
    
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
    
    if (_nSelJets>0) _jet_0_discrim = vSelJets[0]->bDiscriminator(btagger_);
    if (_nSelJets>1) _jet_1_discrim = vSelJets[1]->bDiscriminator(btagger_);
    if (_nSelJets>2) _jet_2_discrim = vSelJets[2]->bDiscriminator(btagger_);
    if (_nSelJets>3) _jet_3_discrim = vSelJets[3]->bDiscriminator(btagger_);
    if (_nSelJets>4) _jet_4_discrim = vSelJets[4]->bDiscriminator(btagger_);
    if (_nSelJets>5) _jet_5_discrim = vSelJets[5]->bDiscriminator(btagger_);
    if (_nSelJets>6) _jet_6_discrim = vSelJets[6]->bDiscriminator(btagger_);
    
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
    
    SetValue("jet_0_discrim", _jet_0_discrim);
    SetValue("jet_1_discrim", _jet_1_discrim);
    SetValue("jet_2_discrim", _jet_2_discrim);
    SetValue("jet_3_discrim", _jet_3_discrim);
    SetValue("jet_4_discrim", _jet_4_discrim);
    SetValue("jet_5_discrim", _jet_5_discrim);
    SetValue("jet_6_discrim", _jet_6_discrim);
    
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
    
    /*  edm::InputTag topJetColl = edm::InputTag("goodPatJetsCATopTagPFPacked");
     edm::Handle<std::vector<pat::Jet> > topJets;
     event.getByLabel(topJetColl, topJets);
     
     //Four vector
     std::vector <double> CATopJetPt;
     std::vector <double> CATopJetEta;
     std::vector <double> CATopJetPhi;
     std::vector <double> CATopJetEnergy;
     
     //Check
     std::vector <double> CATopJetTag;
     
     //Top-like properties
     std::vector <double> CATopJetTopMass;
     std::vector <double> CATopJetMinPairMass;
     
     for (std::vector<pat::Jet>::const_iterator ijet = topJets->begin(); ijet != topJets->end(); ijet++){
     
     CATopJetPt     . push_back(ijet->pt());
     CATopJetEta    . push_back(ijet->eta());
     CATopJetPhi    . push_back(ijet->phi());
     CATopJetEnergy . push_back(ijet->energy());
     
     reco::CATopJetTagInfo* jetInfo = (reco::CATopJetTagInfo*) ijet->tagInfo("CATop");
     
     CATopJetTopMass     . push_back(jetInfo->properties().topMass);
     CATopJetMinPairMass . push_back(jetInfo->properties().minMass);
     
     
     //check if top jet should be tagged
     //nDaughters>=3, 250>topMass>140, minPairMass>50
     //should we make a pt cut? (tag efficiencies and tagger should only work for certain pt threshold/ranges)
     if (ijet->numberOfDaughters()<3 || jetInfo->properties().topMass<140 || jetInfo->properties().topMass>250 || jetInfo->properties().minMass <50) CATopJetTag      . push_back(false);
     else CATopJetTag . push_back(true);
     
     }
     */
    //Four vector
    //SetValue("CATopJetPt"     , CATopJetPt);
    SetValue("CATopJetPt"     , -1);
    //SetValue("CATopJetEta"    , CATopJetEta);
    SetValue("CATopJetEta"     , 9999);
    //SetValue("CATopJetPhi"    , CATopJetPhi);
    SetValue("CATopJetPhi"     , 9999);
    //SetValue("CATopJetEnergy" , CATopJetEnergy);
    SetValue("CATopJetEnergy"     , -1);
    
    //check
    //SetValue("CATopJetTag"    , CATopJetTag);
    SetValue("CATopJetTag"     , -1);
    
    //Properties
    //SetValue("CATopJetTopMass"     , CATopJetTopMass);
    SetValue("CATopJetMass"     , -1);
    //SetValue("CATopJetMinPairMass" , CATopJetMinPairMass);
    SetValue("CATopJetMinPairMass"     , -1);
    
    double _weight_WJets = 1.0;
    if (isWJets_) {
        if ( nBjets > 0)  _weight_WJets = 1.21*0.98;
        if ( nBjets == 0 and nCjets > 0) _weight_WJets = 1.66*0.98;
        if ( nBjets == 0 and nCjets == 0) _weight_WJets = 1.0*0.82;
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
    double _genDrWprimeBjetJet0 = -9999.0;
    double _genDrWprimeBjetJet1 = -9999.0;
    double _genDrTopBjetJet0 = -9999.0;
    double _genDrTopBjetJet1 = -9999.0;
    double _genWprimeMass = -1.0;
    double _genWprimeBjetPt = -1.0;
    double _genWprimeBjetEta = -1.0;
    double _genWprimeBjetPhi = -1.0;
    double _genTopBjetPt = -1.0;
    double _genTopBjetEta = -1.0;
    double _genTopBjetPhi = -1.0;
    double _genLeptonEta = -1.0;
    double _genLeptonPhi = -1.0;
    double _genLeptonPt = -1.0;
    double _genNuPz = -1.0;
    double _genJet0Pt = -1.0;
    double _genJet1Pt = -1.0;
    double _genJet2Pt = -1.0;
    double _genJet3Pt = -1.0;
    double _genrecoDRJet0 = 9999.0;
    double _genrecoDRJet1 = 9999.0;
    double _genrecoDRJet2 = 9999.0;
    double _genrecoDRJet3 = 9999.0;
    double _genrecoPTJet0 = 9999.0;
    double _genrecoPTJet1 = 9999.0;
    double _genrecoPTJet2 = 9999.0;
    double _genrecoPTJet3 = 9999.0;
    
    int _hasGenMu = -1;
    int _hasGenEl = -1;
    
    if (isTB_) {
        edm::Handle<reco::GenParticleCollection> genParticles;
        edm::InputTag genParticles_it = edm::InputTag("prunedGenParticles");
        event.getByLabel(genParticles_it, genParticles);
        
        edm::Handle<reco::GenJetCollection> genJets;
        edm::InputTag genJets_it = edm::InputTag("slimmedGenJets");
        event.getByLabel(genJets_it, genJets);
        
        int qLep = 0;
        math::XYZTLorentzVector lv_genLep;
        math::XYZTLorentzVector lv_genNu;
        math::XYZTLorentzVector lv_genB;
        math::XYZTLorentzVector lv_genBbar;
        std::vector<math::XYZTLorentzVector> lv_genJets;
        
        //std::cout << "-----------------Event start-------------------" << std::endl;
        
        for (size_t i = 0; i < genParticles->size(); i++) {
            const reco::GenParticle & p = (*genParticles).at(i);
            //std::cout << "Status = " << p.status() << "\tId = " << p.pdgId() << std::endl;
            if (p.status() == 23 || p.status() == 33) {
                if (fabs(p.pdgId())==11 or fabs(p.pdgId())==13) {
                    if (fabs(p.pdgId())==11) _hasGenEl = 1;
                    if (fabs(p.pdgId())==13) _hasGenMu = 1;
                    
                    lv_genLep = p.p4();
                    if (p.pdgId()<0) qLep = 1;
                    else qLep = -1;
                }
                if (fabs(p.pdgId())==12 or fabs(p.pdgId())==14) lv_genNu = p.p4();
                if (p.pdgId()==5) lv_genB = p.p4();
                if (p.pdgId()==-5) lv_genBbar = p.p4();
                //std::cout << "genParticle Id=" << p.pdgId() << " , pT=" << p.pt() << " , eta=" << p.eta() << " , phi=" << p.phi() << std::endl;
            }
            //      if (p.status() == 23 || p.status() == 33 || p.status() == 14 || p.status() == 15 || p.status() == 24 || p.status() == 34 || p.status() == 43 || p.status() == 51 || p.status() == 62 || p.status() == 63 || p.status() == 91) {
            //        if (fabs(p.pdgId())<=5 || p.pdgId()==9 || p.pdgId()==21) lv_genPartons.push_back(p.p4());
            //      }
        }
        for (size_t j = 0; j < genJets->size(); j++) {
            const reco::GenJet & je = (*genJets).at(j);
            //std::cout << "genJet Id=" << je.pdgId() << " , pT=" << je.pt() << " , eta=" << je.eta() << " , phi=" << je.phi() << std::endl;
            lv_genJets.push_back(je.p4());
        }
        /*if (_nSelJets>0) std::cout << "recoJet partonFlavour=" << _jet_0_flavor << " , pT=" << _jet_0_pt << " , eta=" << _jet_0_eta << " , phi=" << _jet_0_phi << std::endl;
         if (_nSelJets>1) std::cout << "recoJet partonFlavour=" << _jet_1_flavor << " , pT=" << _jet_1_pt << " , eta=" << _jet_1_eta << " , phi=" << _jet_1_phi << std::endl;
         if (_nSelJets>2) std::cout << "recoJet partonFlavour=" << _jet_2_flavor << " , pT=" << _jet_2_pt << " , eta=" << _jet_2_eta << " , phi=" << _jet_2_phi << std::endl;
         if (_nSelJets>3) std::cout << "recoJet partonFlavour=" << _jet_3_flavor << " , pT=" << _jet_3_pt << " , eta=" << _jet_3_eta << " , phi=" << _jet_3_phi << std::endl;
         if (_nSelJets>4) std::cout << "recoJet partonFlavour=" << _jet_4_flavor << " , pT=" << _jet_4_pt << " , eta=" << _jet_4_eta << " , phi=" << _jet_4_phi << std::endl;
         if (_nSelJets>5) std::cout << "recoJet partonFlavour=" << _jet_5_flavor << " , pT=" << _jet_5_pt << " , eta=" << _jet_5_eta << " , phi=" << _jet_5_phi << std::endl;
         if (_nSelJets>6) std::cout << "recoJet partonFlavour=" << _jet_6_flavor << " , pT=" << _jet_6_pt << " , eta=" << _jet_6_eta << " , phi=" << _jet_6_phi << std::endl;
         if (_nSelJets>7) std::cout << "recoJet partonFlavour=" << _jet_7_flavor << " , pT=" << _jet_7_pt << " , eta=" << _jet_7_eta << " , phi=" << _jet_7_phi << std::endl;
         if (_nSelJets>8) std::cout << "recoJet partonFlavour=" << _jet_8_flavor << " , pT=" << _jet_8_pt << " , eta=" << _jet_8_eta << " , phi=" << _jet_8_phi << std::endl;
         if (_nSelJets>9) std::cout << "recoJet partonFlavour=" << _jet_9_flavor << " , pT=" << _jet_9_pt << " , eta=" << _jet_9_eta << " , phi=" << _jet_9_phi << std::endl;*/
        
        
        double deta = -999.0;
        double dphi = -999.0;
        
        std::sort(lv_genJets.begin(), lv_genJets.end(), my_compare);
        
        if (_nSelJets>0 && lv_genJets.size()>0 ) {
            _genJet0Pt=lv_genJets[0].Pt();
            deta = lv_genJets[0].Eta()-_jet_0_eta;
            dphi = lv_genJets[0].Phi()-_jet_0_phi;
            if ( dphi > TMath::Pi() ) dphi -= 2.*TMath::Pi();
            if ( dphi <= -TMath::Pi() ) dphi += 2.*TMath::Pi();
            _genrecoDRJet0 = TMath::Sqrt(deta*deta + dphi*dphi);
            _genrecoPTJet0 = _genJet0Pt-_jet_0_pt;
        }
        if (_nSelJets>1 && lv_genJets.size()>1 ) {
            _genJet1Pt=lv_genJets[1].Pt();
            deta = lv_genJets[1].Eta()-_jet_1_eta;
            dphi = lv_genJets[1].Phi()-_jet_1_phi;
            if ( dphi > TMath::Pi() ) dphi -= 2.*TMath::Pi();
            if ( dphi <= -TMath::Pi() ) dphi += 2.*TMath::Pi();
            _genrecoDRJet1 = TMath::Sqrt(deta*deta + dphi*dphi);
            _genrecoPTJet1 = _genJet1Pt-_jet_1_pt;
        }
        if (_nSelJets>2 && lv_genJets.size()>2 ) {
            _genJet2Pt=lv_genJets[2].Pt();
            deta = lv_genJets[2].Eta()-_jet_0_eta;
            dphi = lv_genJets[2].Phi()-_jet_0_phi;
            if ( dphi > TMath::Pi() ) dphi -= 2.*TMath::Pi();
            if ( dphi <= -TMath::Pi() ) dphi += 2.*TMath::Pi();
            _genrecoDRJet2 = TMath::Sqrt(deta*deta + dphi*dphi);
            _genrecoPTJet2 = _genJet2Pt-_jet_2_pt;
        }
        if (_nSelJets>3 && lv_genJets.size()>3 ) {
            _genJet3Pt=lv_genJets[3].Pt();
            deta = lv_genJets[3].Eta()-_jet_0_eta;
            dphi = lv_genJets[3].Phi()-_jet_0_phi;
            if ( dphi > TMath::Pi() ) dphi -= 2.*TMath::Pi();
            if ( dphi <= -TMath::Pi() ) dphi += 2.*TMath::Pi();
            _genrecoDRJet3 = TMath::Sqrt(deta*deta + dphi*dphi);
            _genrecoPTJet3 = _genJet3Pt-_jet_3_pt;
        }
        
        _genNuPz = lv_genNu.Pz();
        
        math::XYZTLorentzVector lv_genTop;
        math::XYZTLorentzVector lv_genOtherTop;
        math::XYZTLorentzVector lv_genWprime;
        
        if (qLep > 0) {
            lv_genTop = lv_genLep + lv_genNu + lv_genB;
            lv_genOtherTop = lv_genLep + lv_genNu + lv_genBbar;
        }
        else {
            lv_genTop = lv_genLep + lv_genNu + lv_genBbar;
            lv_genOtherTop = lv_genLep + lv_genNu + lv_genB;
        }
        
        lv_genWprime = lv_genLep + lv_genNu + lv_genB + lv_genBbar;
        
        _genTopMass = lv_genTop.M();
        _genTopPt = lv_genTop.Pt();
        _genWprimeMass = lv_genWprime.M();
        _genWprimeBjetPt = (qLep>0?lv_genBbar.Pt():lv_genB.Pt());
        _genWprimeBjetEta = (qLep>0?lv_genBbar.Eta():lv_genB.Eta());
        _genWprimeBjetPhi = (qLep>0?lv_genBbar.Phi():lv_genB.Phi());
        _genTopBjetPt = (qLep<0?lv_genBbar.Pt():lv_genB.Pt());
        _genTopBjetEta = (qLep<0?lv_genBbar.Eta():lv_genB.Eta());
        _genTopBjetPhi = (qLep<0?lv_genBbar.Phi():lv_genB.Phi());
        _genLeptonEta = lv_genLep.Eta();
        _genLeptonPhi = lv_genLep.Phi();
        _genLeptonPt = lv_genLep.Pt();
        
        if (qLep<0) {
            deta = lv_genLep.Eta()-lv_genBbar.Eta();
            dphi = lv_genLep.Phi()-lv_genBbar.Phi();
            //std::cout<<"lepton pt,eta,phi,E = "<<lv_genLep.Pt()<<", "<<lv_genLep.Eta()<<", "<<lv_genLep.Phi()<<", "<<lv_genLep.E()<<", "<<std::endl;
            //std::cout<<"met pt,eta,phi,E = "<<lv_genNu.Pt()<<", "<<lv_genNu.Eta()<<", "<<lv_genNu.Phi()<<", "<<lv_genNu.E()<<", "<<std::endl;
            //std::cout<<"bjet pt,eta,phi,E = "<<lv_genBbar.Pt()<<", "<<lv_genBbar.Eta()<<", "<<lv_genBbar.Phi()<<", "<<lv_genBbar.E()<<", "<<std::endl;
            //std::cout<<"other bjet pt,eta,phi,E = "<<lv_genB.Pt()<<", "<<lv_genB.Eta()<<", "<<lv_genB.Phi()<<", "<<lv_genB.E()<<", "<<std::endl;
            if ( dphi > TMath::Pi() ) dphi -= 2.*TMath::Pi();
            if ( dphi <= -TMath::Pi() ) dphi += 2.*TMath::Pi();
            _genDrLeptonTopBjet = TMath::Sqrt(deta*deta + dphi*dphi);
            //std::cout<<"gen dR = "<<_genDrLeptonTopBjet<<std::endl;
        } else {
            deta = lv_genLep.Eta()-lv_genB.Eta();
            dphi = lv_genLep.Phi()-lv_genB.Phi();
            //std::cout<<"lepton pt,eta,phi,E = "<<lv_genLep.Pt()<<", "<<lv_genLep.Eta()<<", "<<lv_genLep.Phi()<<", "<<lv_genLep.E()<<", "<<std::endl;
            //std::cout<<"met pt,eta,phi,E = "<<lv_genNu.Pt()<<", "<<lv_genNu.Eta()<<", "<<lv_genNu.Phi()<<", "<<lv_genNu.E()<<", "<<std::endl;
            //std::cout<<"bjet pt,eta,phi,E = "<<lv_genB.Pt()<<", "<<lv_genB.Eta()<<", "<<lv_genB.Phi()<<", "<<lv_genB.E()<<", "<<std::endl;
            //std::cout<<"other bjet pt,eta,phi,E = "<<lv_genBbar.Pt()<<", "<<lv_genBbar.Eta()<<", "<<lv_genBbar.Phi()<<", "<<lv_genBbar.E()<<", "<<std::endl;
            if ( dphi > TMath::Pi() ) dphi -= 2.*TMath::Pi();
            if ( dphi <= -TMath::Pi() ) dphi += 2.*TMath::Pi();
            _genDrLeptonTopBjet = TMath::Sqrt(deta*deta + dphi*dphi);
            //std::cout<<"gen dR = "<<_genDrLeptonTopBjet<<std::endl;
        }
        if (_nCorrBtagJets>0) {
            deta = _jet_0_eta-_genTopBjetEta;
            dphi = _jet_0_phi-_genTopBjetPhi;
            if ( dphi > TMath::Pi() ) dphi -= 2.*TMath::Pi();
            if ( dphi <= -TMath::Pi() ) dphi += 2.*TMath::Pi();
            _genDrTopBjetJet0 = TMath::Sqrt(deta*deta + dphi*dphi);
            deta = _jet_0_eta-_genWprimeBjetEta;
            dphi = _jet_0_phi-_genWprimeBjetPhi;
            if ( dphi > TMath::Pi() ) dphi -= 2.*TMath::Pi();
            if ( dphi <= -TMath::Pi() ) dphi += 2.*TMath::Pi();
            _genDrWprimeBjetJet0 = TMath::Sqrt(deta*deta + dphi*dphi);
        }
        if (_nCorrBtagJets>1) {
            deta = _jet_1_eta-_genTopBjetEta;
            dphi = _jet_1_phi-_genTopBjetPhi;
            if ( dphi > TMath::Pi() ) dphi -= 2.*TMath::Pi();
            if ( dphi <= -TMath::Pi() ) dphi += 2.*TMath::Pi();
            _genDrTopBjetJet1 = TMath::Sqrt(deta*deta + dphi*dphi);    
            deta = _jet_1_eta-_genWprimeBjetEta;
            dphi = _jet_1_phi-_genWprimeBjetPhi;
            if ( dphi > TMath::Pi() ) dphi -= 2.*TMath::Pi();
            if ( dphi <= -TMath::Pi() ) dphi += 2.*TMath::Pi();
            _genDrWprimeBjetJet1 = TMath::Sqrt(deta*deta + dphi*dphi);    
        }
        
    }
    
    SetValue("hasGenEl", _hasGenEl);
    SetValue("hasGenMu", _hasGenMu);
    SetValue("genTopMass", _genTopMass);
    SetValue("genNuPz", _genNuPz);
    
    SetValue("genTopPt", _genTopPt);
    SetValue("genDrLeptonTopBjet", _genDrLeptonTopBjet);
    SetValue("genWprimeMass", _genWprimeMass);
    SetValue("genWprimeBjetPt", _genWprimeBjetPt);
    SetValue("genWprimeBjetEta", _genWprimeBjetEta);
    SetValue("genWprimeBjetPhi", _genWprimeBjetPhi);
    SetValue("genTopBjetPt", _genTopBjetPt);
    SetValue("genTopBjetEta", _genTopBjetEta);
    SetValue("genTopBjetPhi", _genTopBjetPhi);
    SetValue("genLeptonPt", _genLeptonPt);
    SetValue("genLeptonEta", _genLeptonEta);
    SetValue("genLeptonPhi", _genLeptonPhi);
    SetValue("genDrWprimeBjetJet0", _genDrWprimeBjetJet0);
    SetValue("genDrWprimeBjetJet1", _genDrWprimeBjetJet1);
    SetValue("genDrTopBjetJet0", _genDrTopBjetJet0);
    SetValue("genDrTopBjetJet1", _genDrTopBjetJet1);
    SetValue("genJet0Pt", _genJet0Pt);
    SetValue("genJet1Pt", _genJet1Pt);
    SetValue("genJet2Pt", _genJet2Pt);
    SetValue("genJet3Pt", _genJet3Pt);
    SetValue("genrecoDRJet0", _genrecoDRJet0);
    SetValue("genrecoDRJet1", _genrecoDRJet1);
    SetValue("genrecoDRJet2", _genrecoDRJet2);
    SetValue("genrecoDRJet3", _genrecoDRJet3);
    SetValue("genrecoPTJet0", _genrecoPTJet0);
    SetValue("genrecoPTJet1", _genrecoPTJet1);
    SetValue("genrecoPTJet2", _genrecoPTJet2);
    SetValue("genrecoPTJet3", _genrecoPTJet3);
    
    double _genTTMass = -1.0;
    double _genTPt = -1.0;
    double _genTbarPt = -1.0;
    double _wTPt = -1.0;
    double _wTbarPt = -1.0;
    double _sfTopPt = 1.0;
    
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
        _genTPt = lv_genT.Pt();
        _genTbarPt = lv_genTbar.Pt();
        
        math::XYZTLorentzVector lv_genTT =  lv_genT + lv_genTbar;
        _genTTMass = lv_genTT.M();
        
        // TOP-12-030 / AN-12-475
        _wTPt = TMath::Exp(0.156-0.00137*_genTPt);
        _wTbarPt = TMath::Exp(0.156-0.00137*_genTbarPt);
        _sfTopPt = TMath::Sqrt(_wTPt*_wTbarPt);
        
    }
    
    SetValue("genTTMass", _genTTMass);
    SetValue("genTPt", _genTPt);
    SetValue("genTbarPt", _genTbarPt);
    SetValue("wTPt", _wTPt);
    SetValue("wTbarPt", _wTbarPt);
    SetValue("weight_TopPt", _sfTopPt);  
    
    return 0;
}
