/*
  Calculator for the Inclusive Tprime analysis

  Author: Michael Luk, 2012
*/


#include <iostream>
#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "TLorentzVector.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/JetReco/interface/CATopJetTagInfo.h"
class LjmetFactory;


class TprimeCalc : public BaseCalc{
  
public:
  
    TprimeCalc();
    virtual ~TprimeCalc(){}

    virtual int BeginJob(){

        if (mPset.exists("triggerSummary")) triggerSummary_ = mPset.getParameter<edm::InputTag>("triggerSummary");
        else                                triggerSummary_ = edm::InputTag("hltTriggerSummaryAOD");

        if (mPset.exists("rhoSrc")) rhoSrc_ = mPset.getParameter<edm::InputTag>("rhoSrc");
        else                        rhoSrc_ = edm::InputTag("kt6PFJets", "rho");


	
	if (mPset.exists("isTB")) isTB_ = mPset.getParameter<bool>("isTB");
	else                         isTB_ = false;

	if (mPset.exists("genParticles")) genParticles_it = mPset.getParameter<edm::InputTag>("genParticles");
	else                              genParticles_it = edm::InputTag("prunedGenParticles");

        if (mPset.exists("isWJets")) isWJets_ = mPset.getParameter<bool>("isWJets");
        else                         isWJets_ = false;

	
	if (mPset.exists("isHiggs")) isHiggs_ = mPset.getParameter<bool>("isHiggs");
        else                         isHiggs_ = false;



        return 0;
    }

    virtual int AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector);
    virtual int EndJob(){return 0;}

  
private:
  
    edm::InputTag triggerSummary_; 
    edm::InputTag rhoSrc_;
    edm::InputTag             genParticles_it;
    bool isWJets_;
    bool isHiggs_; 
    bool isTB_;
};



//static int reg = LjmetFactory::GetInstance()->Register(new TprimeCalc(), "TprimeCalc");



TprimeCalc::TprimeCalc(){
}

int TprimeCalc::AnalyzeEvent(edm::EventBase const & event,
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
    double _electron_1_SCeta = -9999.0;
    double _electron_1_RelIso = -9999.0;
    if (_nSelElectrons>0) {
        _electron_1_pt = vSelElectrons[0]->ecalDrivenMomentum().pt();
        _electron_1_phi = vSelElectrons[0]->phi();
        _electron_1_eta = vSelElectrons[0]->eta();
	_electron_1_SCeta = vSelElectrons[0]->superCluster()->eta();
	
        double AEff  = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, vSelElectrons[0]->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012);
        double chIso = vSelElectrons[0]->chargedHadronIso();
        double nhIso = vSelElectrons[0]->neutralHadronIso();
        double phIso = vSelElectrons[0]->photonIso();
        _electron_1_RelIso  = ( chIso + max(0.0, nhIso + phIso - rhoIso*AEff) )/ vSelElectrons[0]->ecalDrivenMomentum().pt();

    }

    SetValue("elec_1_pt", _electron_1_pt);
    SetValue("elec_1_phi", _electron_1_phi);
    SetValue("elec_1_eta", _electron_1_eta);
    SetValue("elec_1_SCeta",_electron_1_SCeta);
    SetValue("elec_1_RelIso", _electron_1_RelIso);

    double _electron_2_pt = -9999.0;
    double _electron_2_phi = -9999.0;
    double _electron_2_eta = -9999.0;
    double _electron_2_SCeta = -9999.0;
    double _electron_2_RelIso = -9999.0;

    if (_nSelElectrons>1) {
        _electron_2_pt = vSelElectrons[1]->ecalDrivenMomentum().pt();
        _electron_2_phi = vSelElectrons[1]->phi();
        _electron_2_eta = vSelElectrons[1]->eta();
        _electron_2_SCeta = vSelElectrons[1]->superCluster()->eta();

        double AEff  = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, vSelElectrons[1]->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012);
        double chIso = vSelElectrons[1]->chargedHadronIso();
        double nhIso = vSelElectrons[1]->neutralHadronIso();
        double phIso = vSelElectrons[1]->photonIso();
        _electron_2_RelIso  = ( chIso + max(0.0, nhIso + phIso - rhoIso*AEff) )/ vSelElectrons[1]->ecalDrivenMomentum().pt();

    }

    SetValue("elec_2_pt", _electron_2_pt);
    SetValue("elec_2_phi", _electron_2_phi);
    SetValue("elec_2_eta", _electron_2_eta);
    SetValue("elec_2_SCeta", _electron_2_SCeta);
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
        if ( 0.9 < abs(_muon_1_eta) && abs(_muon_1_eta) <= 1.2 ) _weight_muon_eff = 0.992*0.998*0.964;
        if ( 1.2 < abs(_muon_1_eta) && abs(_muon_1_eta) <= 2.1 ) _weight_muon_eff = 0.998*1.002*0.994;
    }
    SetValue("weight_MuonEff", _weight_muon_eff);


    double _weight_muon_eff_ABCD = 1.0;
    //A+B+C+D Tight ID*ISO(dBeta<0.12)*TRIG(IsoMu24_eta2p1)
    //
    // trigger SF have been updated
    // --> A(4.59%) + B(23.62) + C(35.18) + D(36.61)
    if (_nSelMuons>0) {
      if (abs(_muon_1_eta) <= 0.9 ) _weight_muon_eff_ABCD = 0.9939*1.0004*0.981;
        if ( 0.9 < abs(_muon_1_eta) && abs(_muon_1_eta) <= 1.2 ) _weight_muon_eff_ABCD = 0.9902*1.0031*0.961;
        if ( 1.2 < abs(_muon_1_eta) && abs(_muon_1_eta) <= 2.1 ) _weight_muon_eff_ABCD = 0.9970*1.0050*0.983;
    }
    SetValue("weight_MuonEff_ABCD", _weight_muon_eff_ABCD);


    



  
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
	if (_electron_1_pt > 30.0 && _electron_1_pt <= 40.0 ) _weight_electron_eff_53x = 0.979;
	else if (_electron_1_pt > 40.0 && _electron_1_pt <= 50.0 ) _weight_electron_eff_53x = 0.984;
	else if (_electron_1_pt > 50.0 ) _weight_electron_eff_53x = 0.983;
      } 
      else if (0.8 < abs(_electron_1_eta) && abs(_electron_1_eta) <= 1.442 ) {
	if (_electron_1_pt > 30.0 && _electron_1_pt <= 40.0 ) _weight_electron_eff_53x = 0.961;
	else if (_electron_1_pt > 40.0 && _electron_1_pt <= 50.0 ) _weight_electron_eff_53x = 0.972;
	else if (_electron_1_pt > 50.0 ) _weight_electron_eff_53x = 0.977;
      }
      else if ( 1.442 < abs(_electron_1_eta) && abs(_electron_1_eta) <= 1.556 ) {
	if (_electron_1_pt > 30.0 && _electron_1_pt <= 40.0 ) _weight_electron_eff_53x = 0.983;
	else if (_electron_1_pt > 40.0 && _electron_1_pt <= 50.0 ) _weight_electron_eff_53x = 0.957;
	else if (_electron_1_pt > 50.0 ) _weight_electron_eff_53x = 0.978;
      }
      else if ( 1.556 < abs(_electron_1_eta) && abs(_electron_1_eta) <=  2.0 ) {
	if (_electron_1_pt > 30.0 && _electron_1_pt <= 40.0 ) _weight_electron_eff_53x = 0.962;
	else if (_electron_1_pt > 40.0 && _electron_1_pt <= 50.0 ) _weight_electron_eff_53x = 0.985;
	else if (_electron_1_pt > 50.0 ) _weight_electron_eff_53x = 0.986;
      }
      else if ( 2.0 < abs(_electron_1_eta) && abs(_electron_1_eta) < 2.5 ) {
	if (_electron_1_pt > 30.0 && _electron_1_pt <= 40.0 ) _weight_electron_eff_53x = 1.002;
	else if (_electron_1_pt > 40.0 && _electron_1_pt <= 50.0 ) _weight_electron_eff_53x =0.999;
	else if (_electron_1_pt > 50.0 ) _weight_electron_eff_53x = 0.995;
      }
    }
    SetValue("weight_ElectronEff_53x", _weight_electron_eff_53x);



    

    // 
    //
    //_____ Jets ______________________________
    //
    
    //Get Top-like jets
    edm::InputTag topJetColl = edm::InputTag("goodPatJetsCATopTagPF");
    edm::Handle<std::vector<pat::Jet> > topJets;
    event.getByLabel(topJetColl, topJets);

    //Four vector
    std::vector <double> CATopJetPt;
    std::vector <double> CATopJetEta;
    std::vector <double> CATopJetPhi;
    std::vector <double> CATopJetEnergy;

    std::vector <double> CATopJetCSV;
    //   std::vector <double> CATopJetRCN;

    //Identity
    std::vector <int> CATopJetIndex;
    std::vector <int> CATopJetnDaughters;

    //Top-like properties
    std::vector <double> CATopJetTopMass;
    std::vector <double> CATopJetMinPairMass;

    //Daughter four vector and index
    std::vector <double> CATopDaughterPt;
    std::vector <double> CATopDaughterEta;
    std::vector <double> CATopDaughterPhi;
    std::vector <double> CATopDaughterEnergy;

    std::vector <int> CATopDaughterMotherIndex;

    for (std::vector<pat::Jet>::const_iterator ijet = topJets->begin(); ijet != topJets->end(); ijet++){

      int index = (int)(ijet-topJets->begin());

      CATopJetPt     . push_back(ijet->pt());
      CATopJetEta    . push_back(ijet->eta());
      CATopJetPhi    . push_back(ijet->phi());
      CATopJetEnergy . push_back(ijet->energy());
    
      CATopJetCSV    . push_back(ijet->bDiscriminator( "combinedSecondaryVertexBJetTags"));
      //     CATopJetRCN    . push_back((ijet->chargedEmEnergy()+ijet->chargedHadronEnergy()) / (ijet->neutralEmEnergy()+ijet->neutralHadronEnergy()));    
             
      CATopJetIndex      . push_back(index);
      CATopJetnDaughters . push_back((int)ijet->numberOfDaughters());

      reco::CATopJetTagInfo* jetInfo = (reco::CATopJetTagInfo*) ijet->tagInfo("CATop");

      CATopJetTopMass     . push_back(jetInfo->properties().topMass);
      CATopJetMinPairMass . push_back(jetInfo->properties().minMass);

      for (size_t ui = 0; ui < ijet->numberOfDaughters(); ui++){
	CATopDaughterPt     . push_back(ijet->daughter(ui)->pt());
	CATopDaughterEta    . push_back(ijet->daughter(ui)->eta());
	CATopDaughterPhi    . push_back(ijet->daughter(ui)->phi());
	CATopDaughterEnergy . push_back(ijet->daughter(ui)->energy());        
      
	CATopDaughterMotherIndex . push_back(index);      
      }
    }
    

    
    //Four vector
    SetValue("CATopJetPt"     , CATopJetPt);
    SetValue("CATopJetEta"    , CATopJetEta);
    SetValue("CATopJetPhi"    , CATopJetPhi);
    SetValue("CATopJetEnergy" , CATopJetEnergy);

    SetValue("CATopJetCSV"    , CATopJetCSV);
    //   SetValue("CATopJetRCN"    , CATopJetRCN);

    //Identity
    SetValue("CATopJetIndex"      , CATopJetIndex);
    SetValue("CATopJetnDaughters" , CATopJetnDaughters);

    //Properties
    SetValue("CATopJetTopMass"     , CATopJetTopMass);
    SetValue("CATopJetMinPairMass" , CATopJetMinPairMass);

    //Daughter four vector and index
    SetValue("CATopDaughterPt"     , CATopDaughterPt);
    SetValue("CATopDaughterEta"    , CATopDaughterEta);
    SetValue("CATopDaughterPhi"    , CATopDaughterPhi);
    SetValue("CATopDaughterEnergy" , CATopDaughterEnergy);
    
    SetValue("CATopDaughterMotherIndex"      , CATopDaughterMotherIndex);

    //Get CA8 jets for W's
    edm::InputTag CAWJetColl = edm::InputTag("goodPatJetsCA8PrunedPF");
    edm::Handle<std::vector<pat::Jet> > CAWJets;
    event.getByLabel(CAWJetColl, CAWJets);

    //Four vector
    std::vector <double> CAWJetPt;
    std::vector <double> CAWJetEta;
    std::vector <double> CAWJetPhi;
    std::vector <double> CAWJetEnergy;

    std::vector <double> CAWJetCSV;
    //   std::vector <double> CAWJetRCN;

    //Identity
    std::vector <int> CAWJetIndex;
    std::vector <int> CAWJetnDaughters;

    //Mass
    std::vector <double> CAWJetMass;
  
    //Daughter four vector and index
    std::vector <double> CAWDaughterPt;
    std::vector <double> CAWDaughterMass;
    
    std::vector <double> CAWDaughterEta;
    std::vector <double> CAWDaughterPhi;
    std::vector <double> CAWDaughterEnergy;
    //std::vector <double> CAWDaughterMass;

    std::vector <int> CAWDaughterMotherIndex;

    //
    // counter in 60 130, mass window
    //
    // counter in 60 130 && 
    // highest pt W
    // highest pt top
    std::vector<double> SelWJetsIndex;
    std::vector<double> SelHighPtWJetsIndex;
    int wjetIndex_ = -1;
    std::vector<double> nSelTopJets;
    
    
    for (std::vector<pat::Jet>::const_iterator ijet = CAWJets->begin(); ijet != CAWJets->end(); ijet++){
      wjetIndex_+=1;
      int index = (int)(ijet-CAWJets->begin());

      //Four vector
      CAWJetPt     . push_back(ijet->pt());
      CAWJetEta    . push_back(ijet->eta());
      CAWJetPhi    . push_back(ijet->phi());
      CAWJetEnergy . push_back(ijet->energy());        

      CAWJetCSV    . push_back(ijet->bDiscriminator( "combinedSecondaryVertexBJetTags"));
      //     CAWJetRCN    . push_back((ijet->chargedEmEnergy()+ijet->chargedHadronEnergy()) / (ijet->neutralEmEnergy()+ijet->neutralHadronEnergy()));    

      //Identity
      CAWJetIndex      . push_back(index);
      CAWJetnDaughters . push_back((int)ijet->numberOfDaughters());

      //Mass
      CAWJetMass . push_back(ijet->mass());
    
      
      if( (ijet->mass() > 60.) && (ijet->mass() < 130.)){
	SelWJetsIndex.push_back(wjetIndex_);
	if(ijet->pt() > 200.){
	  SelHighPtWJetsIndex.push_back(wjetIndex_);
	}
      }
      

      for (size_t ui = 0; ui < ijet->numberOfDaughters(); ui++){
	CAWDaughterPt     . push_back(ijet->daughter(ui)->pt());
	CAWDaughterEta    . push_back(ijet->daughter(ui)->eta());
	CAWDaughterPhi    . push_back(ijet->daughter(ui)->phi());
	CAWDaughterEnergy . push_back(ijet->daughter(ui)->energy());        
	
	CAWDaughterMotherIndex . push_back(index);      
	CAWDaughterMass     . push_back(ijet->daughter(ui)->mass());
      }
    }

    //Four vector
    SetValue("CAWJetPt"     , CAWJetPt);
    SetValue("CAWJetMass"     , CAWJetMass);
    SetValue("CAWJetEta"    , CAWJetEta);
    SetValue("CAWJetPhi"    , CAWJetPhi);
    SetValue("CAWJetEnergy" , CAWJetEnergy);

    SetValue("CAWJetCSV"    , CAWJetCSV);
    //   SetValue("CAWJetRCN"    , CAWJetRCN);

    //Identity
    SetValue("CAWJetIndex"      , CAWJetIndex);
    SetValue("CAWJetnDaughters" , CAWJetnDaughters);
  
    //Mass
    SetValue("CAWJetMass"     , CAWJetMass);

    //Daughter four vector and index
    
    SetValue("CAWDaughterPt"     , CAWDaughterPt);
    SetValue("CAWDaughterMass"     , CAWDaughterMass);
    SetValue("CAWDaughterEta"    , CAWDaughterEta);
    SetValue("CAWDaughterPhi"    , CAWDaughterPhi);
    SetValue("CAWDaughterEnergy" , CAWDaughterEnergy);

    SetValue("CAWDaughterMotherIndex" , CAWDaughterMotherIndex);
    
    
    SetValue("SelWJetsIndex",SelWJetsIndex);
    SetValue("SelHighPtWJetsIndex",SelHighPtWJetsIndex);
    
    int _nSelWJets = (int)SelWJetsIndex.size();
    SetValue("nSelWJets",_nSelWJets);
        
    int _wjet_0_pt = -999;
    int _wjet_1_pt = -999;
    int _wjet_2_pt = -999;
    int _wjet_3_pt = -999;
    int _wjet_4_pt = -999;

    if (_nSelWJets>0) _wjet_0_pt = CAWJetPt[SelWJetsIndex[0]];
    if (_nSelWJets>1) _wjet_1_pt = CAWJetPt[SelWJetsIndex[1]];
    if (_nSelWJets>2) _wjet_2_pt = CAWJetPt[SelWJetsIndex[2]];
    if (_nSelWJets>3) _wjet_3_pt = CAWJetPt[SelWJetsIndex[3]];
    if (_nSelWJets>4) _wjet_4_pt = CAWJetPt[SelWJetsIndex[4]];
    
    SetValue("wjet_0_pt", _wjet_0_pt);
    SetValue("wjet_1_pt", _wjet_1_pt);
    SetValue("wjet_2_pt", _wjet_2_pt);
    SetValue("wjet_3_pt", _wjet_3_pt);
    SetValue("wjet_4_pt", _wjet_4_pt);


    int _nSelHighPtWJets = (int)SelHighPtWJetsIndex.size();
    SetValue("nSelHighPtWJets",_nSelHighPtWJets);
    
    int _wjet_highpt_0_pt = -999;
    int _wjet_highpt_1_pt = -999;
    int _wjet_highpt_2_pt = -999;
    int _wjet_highpt_3_pt = -999;
    int _wjet_highpt_4_pt = -999;

    if (_nSelHighPtWJets>0) _wjet_highpt_0_pt = CAWJetPt[SelHighPtWJetsIndex[0]];
    if (_nSelHighPtWJets>1) _wjet_highpt_1_pt = CAWJetPt[SelHighPtWJetsIndex[1]];
    if (_nSelHighPtWJets>2) _wjet_highpt_2_pt = CAWJetPt[SelHighPtWJetsIndex[2]];
    if (_nSelHighPtWJets>3) _wjet_highpt_3_pt = CAWJetPt[SelHighPtWJetsIndex[3]];
    if (_nSelHighPtWJets>4) _wjet_highpt_4_pt = CAWJetPt[SelHighPtWJetsIndex[4]];
    
    SetValue("wjet_highpt_0_pt", _wjet_highpt_0_pt);
    SetValue("wjet_highpt_1_pt", _wjet_highpt_1_pt);
    SetValue("wjet_highpt_2_pt", _wjet_highpt_2_pt);
    SetValue("wjet_highpt_3_pt", _wjet_highpt_3_pt);
    SetValue("wjet_highpt_4_pt", _wjet_highpt_4_pt);


    //Get all CA8 jets (not just for W and Top)
    edm::InputTag CA8JetColl = edm::InputTag("goodPatJetsCA8PF");
    edm::Handle<std::vector<pat::Jet> > CA8Jets;
    event.getByLabel(CA8JetColl, CA8Jets);

    //Four vector
    std::vector <double> CA8JetPt;
    std::vector <double> CA8JetEta;
    std::vector <double> CA8JetPhi;
    std::vector <double> CA8JetEnergy;

    std::vector <double> CA8JetCSV;
    //   std::vector <double> CA8JetRCN;

    for (std::vector<pat::Jet>::const_iterator ijet = CA8Jets->begin(); ijet != CA8Jets->end(); ijet++){
    
      //Four vector
      CA8JetPt     . push_back(ijet->pt());
      CA8JetEta    . push_back(ijet->eta());
      CA8JetPhi    . push_back(ijet->phi());
      CA8JetEnergy . push_back(ijet->energy());

      CA8JetCSV    . push_back(ijet->bDiscriminator( "combinedSecondaryVertexBJetTags"));
      //     CA8JetRCN    . push_back((ijet->chargedEmEnergy()+ijet->chargedHadronEnergy()) / (ijet->neutralEmEnergy()+ijet->neutralHadronEnergy()));    
    }

    //Four vector
    SetValue("CA8JetPt"     , CA8JetPt);
    SetValue("CA8JetEta"    , CA8JetEta);
    SetValue("CA8JetPhi"    , CA8JetPhi);
    SetValue("CA8JetEnergy" , CA8JetEnergy);

    SetValue("CA8JetCSV"    , CA8JetCSV);
    //   SetValue("CA8JetRCN"    , CA8JetRCN);
    



    // gen info to reweight Higgs decays - in tprime MC only. 
    double higgsWeight = 1.;
    double higgsWeight_test = 0.;
    double higgsMassGenerator = 0.;
    
    if (isHiggs_){
      // scale factors used to scale BR of 120 GeV higgs -> 125 GeV higgs (i.e. BR(H125->XX)/BR(H120->XX))
      double higgsBBSf     = 1.47415;
      double higgsTauTauSf = 1.32511;
      double higgsMuMuSf   = 1.30178;
      double higgsCCSf     = 1.35842;
      
      double higgsGGSf         = 0.25024;
      double higgsGammaGammaSf = 5.55457;
      double higgsZGammaSf     = 1.61765;
      double higgsWWSf = 1.22012;
      double higgsZZSf = 1.38307;
      
      edm::Handle<reco::GenParticleCollection> genParticles;
      event.getByLabel(genParticles_it, genParticles);
      // loop over all gen particles in event
      for(size_t i = 0; i < genParticles->size(); i++){
	const reco::GenParticle & p = (*genParticles).at(i);
	int id = p.pdgId();
	
	
	// find higgs (25)
	if(abs(id) == 25 && p.status() == 3){
	  higgsMassGenerator = p.mass();
	  std::cout<<higgsMassGenerator<<std::endl;
	  int absDauIds = 0;
	  size_t nDaughters = p.numberOfDaughters();
	  // get all daughters
	  //std::cout<<"number of daughters "<<nDaughters<<std::endl;
	  for(size_t j = 0; j < nDaughters; ++ j) {
	    int dauId = (p.daughter(j))->pdgId();
	    absDauIds += abs(dauId);
	    ///std::cout<<" daughid "<<dauId<<std::endl;
	  }// daughters
	  //std::cout<<" total value "<<absDauIds<<std::endl;
	  // for each higgs, find the decay products and weight accordingly
	  if(absDauIds==10){
	    higgsWeight *= higgsBBSf;
	    higgsWeight_test = 1;
	  }
	    
	  else if(absDauIds==30){
	    higgsWeight *= higgsTauTauSf;
	    higgsWeight_test = 2;
	  }
	  else if(absDauIds==26){
	    higgsWeight *= higgsMuMuSf;
	    higgsWeight_test = 3;
	  }
	  else if(absDauIds==8){
	    higgsWeight *= higgsCCSf;
	    higgsWeight_test = 4;
	  }
	  else if(absDauIds==42){
	    higgsWeight *= higgsGGSf;
	    higgsWeight_test = 5;
	  }
	  else if(absDauIds==44){
	    higgsWeight *= higgsGammaGammaSf;
	    higgsWeight_test = 6;
	  }
	  else if(absDauIds==45){
	    higgsWeight *= higgsZGammaSf;
	    higgsWeight_test = 7;
	  }
	  else if(absDauIds==48){
	    higgsWeight *= higgsWWSf;
	    higgsWeight_test = 8;
	  }
	  else if(absDauIds==46){
	    higgsWeight *= higgsZZSf;
	    higgsWeight_test = 9;
	  }
	  else{
	    higgsWeight_test = - absDauIds;
	  }
	} // if higgs
      }  // end gen particles loop
    }// ends if Higgs
    SetValue("weight_Higgs",higgsWeight);
    SetValue("weight_Higgs_test",higgsWeight_test);


    
    double ttbarWeight = 1.;
    bool isTTbar_ = false;
    if (isTTbar_){
      // scale factors used to scale BR of 120 GeV higgs -> 125 GeV higgs (i.e. BR(H125->XX)/BR(H120->XX))
      edm::Handle<reco::GenParticleCollection> genParticles;
      event.getByLabel(genParticles_it, genParticles);
      // loop over all gen particles in event
      for(size_t i = 0; i < genParticles->size(); i++){
	const reco::GenParticle & p = (*genParticles).at(i);
	reco::Candidate* mother = (reco::Candidate*) p.mother();
	int id = p.pdgId();
	
	// find bottom (5)
	if(abs(id) == 5 && p.status() == 3){
	  int motherId = mother->pdgId();//(p.mother).pdgId();
	  if(abs(motherId) != 6){
	    ttbarWeight += 1;
	  }
	} 
      }  
    }// ends if Higgs
    SetValue("weight_TTbar",ttbarWeight);


    return 0;
}
