/*
 Calculator for the SLiTT analysis
 
 Author: Gena Kukartsev, 2012
 */



#include <iostream>
#include <limits>
#include <stdio.h>
#include "TFile.h"
#include "TMath.h"
#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "DataFormats/FWLite/interface/Record.h"
#include "DataFormats/FWLite/interface/EventSetup.h"
#include "DataFormats/FWLite/interface/ESHandle.h"



class LjmetFactory;



class StopCalc : public BaseCalc{
    
public:
    
    StopCalc();
    virtual ~StopCalc();
    
    virtual int BeginJob();
    virtual int AnalyzeEvent(edm::Event const & event, BaseEventSelector * selector);
    virtual int EndJob();
    
    
private:
    
    void SetBestCandidateVars(math::XYZTLorentzVector vLep,
                              math::XYZTLorentzVector vMet,
                              std::vector<math::XYZTLorentzVector> vLJets,
                              std::vector<math::XYZTLorentzVector> vBJets,
                              std::string suffix = "");
    math::XYZTLorentzVector TlvToXyzt(TLorentzVector tlv);
    
    double GetLikelihood(double mWlep, double mWhad, double mTlep, double mThad);
    
    std::vector<double> GetNeutrinoPz(math::XYZTLorentzVector lv_mu,
                                      math::XYZTLorentzVector lv_met,
                                      int & success);
    
    
    TFile * f_BTag;
    fwlite::EventSetup * es;
    fwlite::RecordID testRecID;
    
    double mSigWhad;
    double mSigWlep;
    double mSigThad;
    double mSigTlep;
    double mMw, mMtop;
    
};



//static int reg = LjmetFactory::GetInstance()->Register(new StopCalc(), "StopCalc");



StopCalc::StopCalc(){
}



StopCalc::~StopCalc(){
}



double StopCalc::GetLikelihood(double mWlep, double mWhad,
                               double mTlep, double mThad){
    
    double _likelihood =
    //TMath::Gaus(mWlep, mMw, mSigWlep)*
    TMath::Gaus(mWhad, mMw, mSigWhad)*
    TMath::Gaus(mTlep, mMtop, mSigTlep)*
    TMath::Gaus(mThad, mMtop, mSigThad);
    
    return _likelihood;
}



std::vector<double> StopCalc::GetNeutrinoPz(math::XYZTLorentzVector lv_mu,
                                            math::XYZTLorentzVector lv_met,
                                            int & success){
    //
    // reconstruct the neutrino from W
    // returns both solutions or pz=0 if fail
    // failure will set success=0
    //
    
    
    std::vector<double> pz(2, 0.0);
    
    // lepton
    double plx = lv_mu.px();
    double ply = lv_mu.py();
    double plz = lv_mu.pz();
    double el  = lv_mu.E();
    
    // neutrino
    double pnx = lv_met.px();
    double pny = lv_met.py();
    
    double a = mMw*mMw/2.0+plx*pnx+ply*pny;
    double b = el*el*(pnx*pnx+pny*pny);
    double c = el*el-plz*plz;
    double d = 1.0+c*(a*a-b)/a/a/plz/plz;
    
    if (d>=0){
        pz[0]  = a*plz/c*(1+sqrt(d));
        pz[1] = a*plz/c*(1-sqrt(d));
        
        success = 1;
    }
    else{
        // debug
        //std::cout << mLegend << "no solution for neutrino pz" << std::endl;
        
        success = 0;
    }
    
    return pz;
}



int StopCalc::BeginJob(){
    
    //f_BTag = new TFile("cond/btag_performance_2012.root","READ");
    //es = new fwlite::EventSetup(f_BTag);
    //if ( !es->exists("BTagPerformanceRecord") ) {
    //std::cout << "Can't find BTagPerformanceRecord in EventSetup" << std::endl;
    //}
    //testRecID = es->recordID("BTagPerformanceRecord");
    
    // SMP-12-015
    // AN-12/224, sec 14.1 fit
    mSigWhad =  15.3; // gev
    
    mSigWlep =  15.0;
    mSigThad = 100.0;
    mSigTlep =  50.0;
    mMw      =  80.4;
    mMtop    = 172.5;
    
    return 0;
}



int StopCalc::EndJob(){
    
    //f_BTag->Close();
    //delete es;
    //delete f_BTag;
    
    return 0;
}



int StopCalc::AnalyzeEvent(edm::Event const & event,
                           BaseEventSelector * selector){
    //
    // compute event variables here
    //
    
    //
    // _____ Get objects from the selector _________________________
    //
    TLorentzVector                              const & corrMET = selector->GetCorrectedMet();
    edm::Ptr<pat::MET>                          const & pMet = selector->GetMet();
    edm::Ptr<reco::PFMET>                       const & pType1CorrMet = selector->GetType1CorrMet();
    std::vector<edm::Ptr<pat::Electron> >       const & vLooseElectrons = selector->GetLooseElectrons();
    std::vector<edm::Ptr<pat::Electron> >       const & vSelElectrons = selector->GetSelectedElectrons();
    std::vector<edm::Ptr<pat::Jet> >            const & vAllJets = selector->GetAllJets();
    std::vector<edm::Ptr<pat::Jet> >            const & vSelBtagJets = selector->GetSelectedBtagJets();
    std::vector<edm::Ptr<pat::Jet> >            const & vSelJets = selector->GetSelectedJets();
    std::vector<edm::Ptr<pat::Muon> >           const & vLooseMuons = selector->GetLooseMuons();
    std::vector<edm::Ptr<pat::Muon> >           const & vSelMuons = selector->GetSelectedMuons();
    std::vector<edm::Ptr<reco::Vertex> >        const & vSelPVs = selector->GetSelectedPVs();
    std::vector<std::pair<TLorentzVector,bool>> const & vCorrBtagJets = selector->GetCorrJetsWithBTags();
    
    
    double mu_pt  = -10.0;
    double mu_eta = -10.0;
    double mu_phi = -10.0;
    
    // jets as come from TLBSM
    double vJetPt[10];
    double vJetEta[10];
    double vJetPhi[10];
    
    // jets with extra corrections (usually better)
    double vCorrJetPt[10];
    double vCorrJetEta[10];
    double vCorrJetPhi[10];
    bool   vCorrJetTag[10];
    
    
    int _nAllJets      = (int)vAllJets.size();
    int _nSelJets      = (int)vSelJets.size();
    int _nSelBtagJets  = (int)vSelBtagJets.size();
    int _nSelMuons     = (int)vSelMuons.size();
    int _nCorrBtagJets  = (int)vCorrBtagJets.size();
    int _nLooseMuons   = (int)vLooseMuons.size();
    int _nSelElectrons = (int)vSelElectrons.size();
    int _nLooseElectrons = (int)vLooseElectrons.size();
    int _nSelPVs = (int)vSelPVs.size();
    
    
    
    //_____ Muons __________________________________________________
    
    if ( vSelMuons.size()>0 ){
        
        mu_pt  = vSelMuons[0]->pt();
        mu_eta = vSelMuons[0]->eta();
        mu_phi = vSelMuons[0]->phi();
        //add tracker muon flag
        
    }
    SetValue("muon_0_pt", vSelMuons[0]->pt());
    SetValue("muon_0_eta", vSelMuons[0]->eta());
    SetValue("muon_0_phi", vSelMuons[0]->phi());
    
    
    //_____ Jets __________________________________________________
    std::vector<int> vCorrLightJetIndex;
    std::vector<int> vCorrTaggedJetIndex;
    for (unsigned int i=0; i<10; ++i){
        
        if ( vSelJets.size()>i ){
            
            vJetPt[i]  = vSelJets[i]->pt();
            vJetEta[i] = vSelJets[i]->eta();
            vJetPhi[i] = vSelJets[i]->phi();
            
        }
        else{
            
            vJetPt[i]  =  -1.0;
            vJetEta[i] = -10.0;
            vJetPhi[i] = -10.0;
            
        }
        
        char buf[128];
        sprintf(buf, "jet_%d_pt", i);
        SetValue(buf, vJetPt[i]);
        sprintf(buf, "jet_%d_eta", i);
        SetValue(buf, vJetEta[i]);
        sprintf(buf, "jet_%d_phi", i);
        SetValue(buf, vJetPhi[i]);
        
        
        // extra corrected jets
        if ( _nCorrBtagJets>(int)i ){
            
            vCorrJetPt[i]  = vCorrBtagJets[i].first.Pt();
            vCorrJetEta[i] = vCorrBtagJets[i].first.Eta();
            vCorrJetPhi[i] = vCorrBtagJets[i].first.Phi();
            vCorrJetTag[i] = vCorrBtagJets[i].second;
            
            
        }
        else{
            
            vCorrJetPt[i]  =  -1.0;
            vCorrJetEta[i] = -10.0;
            vCorrJetPhi[i] = -10.0;
            vCorrJetTag[i] = false;
            
        }
        
        sprintf(buf, "corr_jet_%d_pt", i);
        SetValue(buf, vCorrJetPt[i]);
        sprintf(buf, "corr_jet_%d_eta", i);
        SetValue(buf, vCorrJetEta[i]);
        sprintf(buf, "corr_jet_%d_phi", i);
        SetValue(buf, vCorrJetPhi[i]);
        sprintf(buf, "corr_jet_%d_tag", i);
        SetValue(buf, (int)(vCorrJetTag[i]));
    }
    
    for (unsigned int i=0; i<vCorrBtagJets.size(); ++i){
        // cache indices of tagged and light jets
        if (vCorrBtagJets[i].second) vCorrTaggedJetIndex.push_back(i);
        else vCorrLightJetIndex.push_back(i);
    }
    
    
    //
    //_____ Event kinematics _______________________________________
    //
    
    SetValue("nPV",_nSelPVs);
    
    
    // m(l,b)
    math::XYZTLorentzVector lv_mu;
    math::XYZTLorentzVector lv_b1;
    math::XYZTLorentzVector lv_b2;
    math::XYZTLorentzVector lv_lb1;
    math::XYZTLorentzVector lv_lb2;
    
    // corrected
    math::XYZTLorentzVector lv_b1_corr;
    math::XYZTLorentzVector lv_b2_corr;
    math::XYZTLorentzVector lv_lb1_corr;
    math::XYZTLorentzVector lv_lb2_corr;
    
    if (vCorrTaggedJetIndex.size()>0){
        lv_b1_corr.SetPxPyPzE(
                              vCorrBtagJets[vCorrTaggedJetIndex[0]].first.Px(),
                              vCorrBtagJets[vCorrTaggedJetIndex[0]].first.Py(),
                              vCorrBtagJets[vCorrTaggedJetIndex[0]].first.Pz(),
                              vCorrBtagJets[vCorrTaggedJetIndex[0]].first.Energy()
                              );
    }
    if (vCorrTaggedJetIndex.size()>1){
        lv_b2_corr.SetPxPyPzE(
                              vCorrBtagJets[vCorrTaggedJetIndex[1]].first.Px(),
                              vCorrBtagJets[vCorrTaggedJetIndex[1]].first.Py(),
                              vCorrBtagJets[vCorrTaggedJetIndex[1]].first.Pz(),
                              vCorrBtagJets[vCorrTaggedJetIndex[1]].first.Energy()
                              );
    }
    
    if (_nSelMuons>0)    lv_mu = vSelMuons[0]->p4();
    if (_nSelBtagJets>0) lv_b1 = vSelBtagJets[0]->p4();
    if (_nSelBtagJets>1) lv_b2 = vSelBtagJets[1]->p4();
    if (_nSelMuons>0 && _nSelBtagJets>0) lv_lb1 = lv_mu+lv_b1;
    if (_nSelMuons>0 && _nSelBtagJets>1) lv_lb2 = lv_mu+lv_b2;
    
    SetValue("mlb1", sqrt(lv_lb1.M2()));
    SetValue("mlb2", sqrt(lv_lb2.M2()));
    
    // b-jet pT
    SetValue("bjet_1_pt", lv_b1.pt());
    SetValue("bjet_1_eta", lv_b1.eta());
    SetValue("bjet_1_phi", lv_b1.phi());
    SetValue("bjet_2_pt", lv_b2.pt());
    SetValue("bjet_2_eta", lv_b2.eta());
    SetValue("bjet_2_phi", lv_b2.phi());
    
    // MET
    double _met = -1.0;
    double _met_phi = -10.0;
    math::XYZTLorentzVector lv_met;
    if(pMet.isNonnull() && pMet.isAvailable()) lv_met = pMet->p4();
    _met = lv_met.pt();
    _met_phi = lv_met.phi();
    SetValue("met", _met);
    SetValue("met_phi", _met_phi);
    
    
    // MET Type 1 Corrected
    double _type1corrmet = -1.0;
    double _type1corrmet_phi = -10.0;
    math::XYZTLorentzVector lv_type1corrmet;
    if(pType1CorrMet.isNonnull() && pType1CorrMet.isAvailable()){
        lv_type1corrmet = pType1CorrMet->p4();
        _type1corrmet = lv_type1corrmet.pt();
        _type1corrmet_phi = lv_type1corrmet.phi();
    }
    SetValue("met_type1Corr", _type1corrmet);
    SetValue("met_phi_type1Corr", _type1corrmet_phi);
    
    // Corrected MET
    Double_t _corr_met = -1.0;
    Double_t _corr_met_phi = -10.0;
    
    if(corrMET.Pt()>0.0) {
        _corr_met = corrMET.Pt();
        _corr_met_phi = corrMET.Phi();
    }
    
    SetValue("met_corr", _corr_met);
    SetValue("met_phi_corr", _corr_met_phi);
    
    
    // mT
    double _mt = -1.0;
    if (_met>0 && _nSelMuons>0){
        _mt = sqrt( 2.0*mu_pt*_met*(1.0-cos(reco::deltaPhi(mu_phi,_met_phi))) );
    }
    SetValue("mt", _mt);
    
    // mT corr
    double _mt_corr = -1.0;
    if (_corr_met>0 && _nSelMuons>0){
        _mt_corr = sqrt( 2.0*mu_pt*_corr_met*(1.0-cos(reco::deltaPhi(mu_phi,_corr_met_phi))) );
    }
    SetValue("mt_corr", _mt_corr);
    
    // HT
    double _ht = _met;
    for (std::vector<edm::Ptr<pat::Jet> >::const_iterator i=vSelJets.begin();
         i!=vSelJets.end(); ++i){
        _ht += (*i)->pt();
    }
    for (std::vector<edm::Ptr<pat::Muon> >::const_iterator i=vSelMuons.begin();
         i!=vSelMuons.end(); ++i){
        _ht += (*i)->pt();
    }
    SetValue("ht", _ht);
    
    // HT corr
    double _ht_corr = _corr_met;
    for (std::vector<std::pair<TLorentzVector,bool> >::const_iterator i=vCorrBtagJets.begin();
         i!=vCorrBtagJets.end(); ++i){
        _ht_corr += i->first.Pt();
    }
    for (std::vector<edm::Ptr<pat::Muon> >::const_iterator i=vSelMuons.begin();
         i!=vSelMuons.end(); ++i){
        _ht_corr += (*i)->pt();
    }
    SetValue("ht_corr", _ht_corr);
    
    
    
    //
    //_____ candidate reconstruction _______________________________
    //
    
    //_____ identify untagged jets only
    
    std::vector<edm::Ptr<pat::Jet> > vLightJets;
    
    // debug
    //std::cout << mLegend << "jets: " << std::endl;
    
    for (std::vector<edm::Ptr<pat::Jet> >::const_iterator i=vSelJets.begin();
         i!=vSelJets.end(); ++i){
        
        bool _isLight = true;
        
        // only consider first two b jets
        //for (std::vector<edm::Ptr<pat::Jet> >::const_iterator j=vSelBtagJets.begin();
        // j!=vSelBtagJets.begin(); ++j){
        for (unsigned int j=0; j<2; ++j){
            
            if ( (*i)==vSelBtagJets[j] ) _isLight=false;
            //if ( (*i)==(*j) ) _isLight=false;
            
        }
        
        if (_isLight) vLightJets.push_back(*i);
        
        // debug
        //std::cout << mLegend << "is light jet: " << _isLight << std::endl;
    }
    
    
    //
    // _____ best candidate with extra corrected objects ___________
    //
    math::XYZTLorentzVector lvCorrMet = TlvToXyzt(corrMET);
    std::vector<math::XYZTLorentzVector> vCorrLJets;
    for (std::vector<int>::const_iterator i=vCorrLightJetIndex.begin();
         i!=vCorrLightJetIndex.end();++i){
        vCorrLJets.push_back(TlvToXyzt(vCorrBtagJets[*i].first));
    }
    std::vector<math::XYZTLorentzVector> vCorrBJets;
    if (vCorrTaggedJetIndex.size()>1){
        vCorrBJets.push_back(TlvToXyzt(vCorrBtagJets[vCorrTaggedJetIndex[0]].first));
        vCorrBJets.push_back(TlvToXyzt(vCorrBtagJets[vCorrTaggedJetIndex[1]].first));
    }
    
    if (vCorrLJets.size()<2 || vCorrBJets.size()<2){
        std::cout << "DEBUG1: light jets " << vCorrLJets.size() << std::endl;
        std::cout << "DEBUG1: b jets " << vCorrBJets.size() << std::endl;
    }
    
    SetBestCandidateVars(lv_mu, lvCorrMet, vCorrLJets, vCorrBJets, "_corr");
    
    
    //
    // _____ best candidate with standard objects __________________
    //
    std::vector<math::XYZTLorentzVector> vDefLJets;
    for (unsigned int i=0; i<vLightJets.size() ;++i){
        vDefLJets.push_back(vLightJets[i]->p4());
    }
    std::vector<math::XYZTLorentzVector> vDefBJets;
    if (vSelBtagJets.size()>1){
        vDefBJets.push_back(vSelBtagJets[0]->p4());
        vDefBJets.push_back(vSelBtagJets[1]->p4());
    }
    
    if (vDefLJets.size()<2 || vDefBJets.size()<2){
        std::cout << "DEBUG2: light jets " << vDefLJets.size() << std::endl;
        std::cout << "DEBUG2: b jets " << vDefBJets.size() << std::endl;
    }
    
    SetBestCandidateVars(lv_mu, lv_met, vDefLJets, vDefBJets, "_def");
    
    
    //_____ find best candidate
    
    double _bestL = -1.0;
    
    double _bestMlb    = -1.0;
    double _bestMlbhad = -1.0;
    
    double _bestSubSmin = -1.0;
    
    double _bestNuPz   = -100000.0;
    double _bestNuPt   = -100.0;
    double _bestNuE    = -100.0;
    double _bestNuM2   = -1.0;
    double _bestNuEta  = -100.0;
    double _bestNuPhi  = -10.0;
    
    double _bestWlepPz = -100000.0;
    double _bestWlepPt = -100.0;
    double _bestWlepE  = -100.0;
    double _bestWlepM  = -1.0;
    double _bestWlepEta= -100.0;
    double _bestWlepPhi= -10.0;
    
    double _bestWhadPz = -100000.0;
    double _bestWhadPt = -100.0;
    double _bestWhadE  = -100.0;
    double _bestWhadM  = -1.0;
    double _bestWhadEta= -100.0;
    double _bestWhadPhi= -10.0;
    
    double _bestTlepPz = -100000.0;
    double _bestTlepPt = -100.0;
    double _bestTlepE  = -100.0;
    double _bestTlepM  = -1.0;
    double _bestTlepEta= -100.0;
    double _bestTlepPhi= -10.0;
    
    double _bestThadPz = -100000.0;
    double _bestThadPt = -100.0;
    double _bestThadE  = -100.0;
    double _bestThadM  = -1.0;
    double _bestThadEta= -100.0;
    double _bestThadPhi= -10.0;
    
    double _bestBlepPz = -100000.0;
    double _bestBlepPt = -100.0;
    double _bestBlepE  = -100.0;
    double _bestBlepEta= -100.0;
    double _bestBlepPhi= -10.0;
    
    double _bestBhadPz = -100000.0;
    double _bestBhadPt = -100.0;
    double _bestBhadE  = -100.0;
    double _bestBhadEta= -100.0;
    double _bestBhadPhi= -10.0;
    
    double _bestJet1Pz = -100000.0;
    double _bestJet1Pt = -100.0;
    double _bestJet1E  = -100.0;
    double _bestJet1Eta= -100.0;
    double _bestJet1Phi= -10.0;
    
    double _bestJet2Pz = -100000.0;
    double _bestJet2Pt = -100.0;
    double _bestJet2E  = -100.0;
    double _bestJet2Eta= -100.0;
    double _bestJet2Phi= -10.0;
    
    //_____ neutrino solution
    int _neuSuccess = 0;
    std::vector<double> pPz = GetNeutrinoPz(lv_mu, lv_met, _neuSuccess);
    //std::cout << mLegend << "first  solution for neutrino pz: " << pPz[0] << std::endl;
    //std::cout << mLegend << "second solution for neutrino pz: " << pPz[1] << std::endl;
    
    SetValue("neutrinoSuccess", _neuSuccess); // were able to reconstruct neutrino
    
    
    //_____ loop over two neutrino solutions
    for (unsigned int n=0; n<2; ++n){
        
        math::XYZTLorentzVector lv_neu(lv_met);
        lv_neu.SetPz(pPz[n]);
        lv_neu.SetE(sqrt(pPz[n]*pPz[n]+lv_met.pt()*lv_met.pt()));
        //std::cout << mLegend << "neutrino mass squared: " << lv_neu.M2() << std::endl;
        
        
        //_____ leptonic W
        math::XYZTLorentzVector lv_Wlep = lv_mu+lv_neu;
        double _mWlep = sqrt(lv_Wlep.M2());
        
        
        //_____ loop over all two-light-jet combinations
        for (unsigned int i=0; i<vLightJets.size(); ++i){
            for (unsigned int j=i+1; j<vLightJets.size(); ++j){
                // i,j make up the hadronic W
                
                //std::cout << mLegend << "total/first/second: " << vLightJets.size()
                //	  << " / " << i << " / " << j << std::endl;
                
                math::XYZTLorentzVector lv_jet1 = vLightJets[i]->p4();
                math::XYZTLorentzVector lv_jet2 = vLightJets[j]->p4();
                
                
                //_____ hadronic W
                math::XYZTLorentzVector lv_Whad = lv_jet1+lv_jet2;
                double _mWhad = sqrt(lv_Whad.M2());
                
                
                //_____ loop over two b jets combinations
                for (unsigned int k=0; k<2; ++k){
                    // loop over the two b jets
                    // k - lepton-side b jet
                    
                    math::XYZTLorentzVector lv_blep = vSelBtagJets[k]->p4();
                    math::XYZTLorentzVector lv_bhad = vSelBtagJets[1-k]->p4();
                    
                    
                    //_____ hadronic top
                    math::XYZTLorentzVector lv_Thad = lv_Whad+lv_bhad;
                    double _mThad = sqrt(lv_Thad.M2());
                    
                    
                    //_____ leptonic top
                    math::XYZTLorentzVector lv_Tlep = lv_Wlep+lv_blep;
                    double _mTlep = sqrt(lv_Tlep.M2());
                    
                    
                    //_____ candidate likelihood
                    double _candL = GetLikelihood(_mWlep, _mWhad, _mTlep, _mThad);
                    
                    //std::cout << mLegend << "lep W mass: " << _mWlep << std::endl;
                    //std::cout << mLegend << "candidate likelihood: " << _candL << std::endl;
                    
                    
                    //_____ update values if better candidate found
                    if (_candL>_bestL){
                        
                        _bestL = _candL;
                        
                        
                        // lepton+b(lep)
                        math::XYZTLorentzVector lv_lb = lv_mu+lv_blep;
                        _bestMlb = sqrt(lv_lb.M2());
                        
                        
                        // lepton+b(had)
                        math::XYZTLorentzVector lv_lbhad = lv_mu+lv_bhad;
                        _bestMlbhad = sqrt(lv_lbhad.M2());
                        
                        
                        // "subsmin" variable:
                        // JHEP 1106 (2011) 041
                        math::XYZTLorentzVector lv_sub = lv_mu +lv_blep+lv_bhad+lv_jet1+lv_jet2;
                        math::XYZTLorentzVector lv_tot = lv_sub+lv_neu;
                        _bestSubSmin = sqrt(
                                            ( sqrt(lv_sub.M2()+lv_sub.pt()*lv_sub.pt())
                                             +
                                             sqrt(lv_neu.M2()+lv_neu.pt()*lv_neu.pt()) )
                                            *
                                            ( sqrt(lv_sub.M2()+lv_sub.pt()*lv_sub.pt())
                                             +
                                             sqrt(lv_neu.M2()+lv_neu.pt()*lv_neu.pt()) )
                                            -
                                            lv_tot.pt()*lv_tot.pt()
                                            );
                        
                        
                        _bestNuPz   = lv_neu.pz();
                        _bestNuPt   = lv_neu.pt();
                        _bestNuE    = lv_neu.E();
                        _bestNuM2   = lv_neu.M2();
                        _bestNuEta  = lv_neu.Eta();
                        _bestNuPhi  = lv_neu.Phi();
                        
                        _bestWlepPz = lv_Wlep.pz();
                        _bestWlepPt = lv_Wlep.pt();
                        _bestWlepE  = lv_Wlep.E();
                        _bestWlepM  = _mWlep;
                        _bestWlepEta= lv_Wlep.Eta();
                        _bestWlepPhi= lv_Wlep.Phi();
                        
                        _bestWhadPz = lv_Whad.pz();
                        _bestWhadPt = lv_Whad.pt();
                        _bestWhadE  = lv_Whad.E();
                        _bestWhadM  = _mWhad;
                        _bestWhadEta= lv_Whad.Eta();
                        _bestWhadPhi= lv_Whad.Phi();
                        
                        _bestTlepPz = lv_Tlep.pz();
                        _bestTlepPt = lv_Tlep.pt();
                        _bestTlepE  = lv_Tlep.E();
                        _bestTlepM  = _mTlep;
                        _bestTlepEta= lv_Tlep.Eta();
                        _bestTlepPhi= lv_Tlep.Phi();
                        
                        _bestThadPz = lv_Thad.pz();
                        _bestThadPt = lv_Thad.pt();
                        _bestThadE  = lv_Thad.E();
                        _bestThadM  = _mThad;
                        _bestThadEta= lv_Thad.Eta();
                        _bestThadPhi= lv_Thad.Phi();
                        
                        _bestBlepPz = lv_blep.pz();
                        _bestBlepPt = lv_blep.pt();
                        _bestBlepE  = lv_blep.E();
                        _bestBlepEta= lv_blep.Eta();
                        _bestBlepPhi= lv_blep.Phi();
                        
                        _bestBhadPz = lv_bhad.pz();
                        _bestBhadPt = lv_bhad.pt();
                        _bestBhadE  = lv_bhad.E();
                        _bestBhadEta= lv_bhad.Eta();
                        _bestBhadPhi= lv_bhad.Phi();
                        
                        _bestJet1Pz = lv_jet1.pz();
                        _bestJet1Pt = lv_jet1.pt();
                        _bestJet1E  = lv_jet1.E();
                        _bestJet1Eta= lv_jet1.Eta();
                        _bestJet1Phi= lv_jet1.Phi();
                        
                        _bestJet2Pz = lv_jet2.pz();
                        _bestJet2Pt = lv_jet2.pt();
                        _bestJet2E  = lv_jet2.E();
                        _bestJet2Eta= lv_jet2.Eta();
                        _bestJet2Phi= lv_jet2.Phi();
                        
                    }
                    else{ // should never happen
                        if (_candL<0) std::cout << mLegend << "cand likelihood: " << _candL << std::endl;
                    }
                    
                    
                } // end of b jet loop
                
            } // inner light jet loop
        } // outer light jet loop
        
    } // neutrino loop
    
    SetValue("bestL", _bestL);
    
    SetValue("bestMlb", _bestMlb);
    SetValue("bestMlbhad", _bestMlbhad);
    
    SetValue("bestSubSmin", _bestSubSmin);
    
    SetValue("bestNuPz", _bestNuPz);
    SetValue("bestNuPt", _bestNuPt);
    SetValue("bestNuE", _bestNuE);
    SetValue("bestNuM2", _bestNuM2);
    SetValue("bestNuEta", _bestNuEta);
    SetValue("bestNuPhi", _bestNuPhi);
    
    SetValue("bestWlepPz", _bestWlepPz);
    SetValue("bestWlepPt", _bestWlepPt);
    SetValue("bestWlepE", _bestWlepE);
    SetValue("bestWlepM", _bestWlepM);
    SetValue("bestWlepEta", _bestWlepEta);
    SetValue("bestWlepPhi", _bestWlepPhi);
    
    SetValue("bestWhadPz", _bestWhadPz);
    SetValue("bestWhadPt", _bestWhadPt);
    SetValue("bestWhadE", _bestWhadE);
    SetValue("bestWhadM", _bestWhadM);
    SetValue("bestWhadEta", _bestWhadEta);
    SetValue("bestWhadPhi", _bestWhadPhi);
    
    SetValue("bestTlepPz", _bestTlepPz);
    SetValue("bestTlepPt", _bestTlepPt);
    SetValue("bestTlepE", _bestTlepE);
    SetValue("bestTlepM", _bestTlepM);
    SetValue("bestTlepEta", _bestTlepEta);
    SetValue("bestTlepPhi", _bestTlepPhi);
    
    SetValue("bestThadPz", _bestThadPz);
    SetValue("bestThadPt", _bestThadPt);
    SetValue("bestThadE", _bestThadE);
    SetValue("bestThadM", _bestThadM);
    SetValue("bestThadEta", _bestThadEta);
    SetValue("bestThadPhi", _bestThadPhi);
    
    SetValue("bestBlepPz", _bestBlepPz);
    SetValue("bestBlepPt", _bestBlepPt);
    SetValue("bestBlepE", _bestBlepE);
    SetValue("bestBlepEta", _bestBlepEta);
    SetValue("bestBlepPhi", _bestBlepPhi);
    
    SetValue("bestBhadPz", _bestBhadPz);
    SetValue("bestBhadPt", _bestBhadPt);
    SetValue("bestBhadE", _bestBhadE);
    SetValue("bestBhadEta", _bestBhadEta);
    SetValue("bestBhadPhi", _bestBhadPhi);
    
    SetValue("bestJet1Pz", _bestJet1Pz);
    SetValue("bestJet1Pt", _bestJet1Pt);
    SetValue("bestJet1E", _bestJet1E);
    SetValue("bestJet1Eta", _bestJet1Eta);
    SetValue("bestJet1Phi", _bestJet1Phi);
    
    SetValue("bestJet2Pz", _bestJet2Pz);
    SetValue("bestJet2Pt", _bestJet2Pt);
    SetValue("bestJet2E", _bestJet2E);
    SetValue("bestJet2Eta", _bestJet2Eta);
    SetValue("bestJet2Phi", _bestJet2Phi);
    
    return 0;
}


math::XYZTLorentzVector StopCalc::TlvToXyzt(TLorentzVector tlv){
    math::XYZTLorentzVector vec;
    vec.SetPxPyPzE(tlv.Px(),tlv.Py(),tlv.Pz(),tlv.Energy());
    return vec;
}



void StopCalc::SetBestCandidateVars(math::XYZTLorentzVector vLep,
                                    math::XYZTLorentzVector vMet,
                                    std::vector<math::XYZTLorentzVector> vLJets,
                                    std::vector<math::XYZTLorentzVector> vBJets,
                                    std::string suffix){
    //
    // Reconstructs the best mu+jets ttbar candidate
    // and saves corresponding quantities to file
    //
    
    double _bestL = -1.0;
    
    double _bestMlb    = -1.0;
    double _bestMlbhad = -1.0;
    
    double _bestSubSmin = -1.0;
    
    double _bestNuPz   = -100000.0;
    double _bestNuPt   = -100.0;
    double _bestNuE    = -100.0;
    double _bestNuM2   = -1.0;
    double _bestNuEta  = -100.0;
    double _bestNuPhi  = -10.0;
    
    double _bestWlepPz = -100000.0;
    double _bestWlepPt = -100.0;
    double _bestWlepE  = -100.0;
    double _bestWlepM  = -1.0;
    double _bestWlepEta= -100.0;
    double _bestWlepPhi= -10.0;
    
    double _bestWhadPz = -100000.0;
    double _bestWhadPt = -100.0;
    double _bestWhadE  = -100.0;
    double _bestWhadM  = -1.0;
    double _bestWhadEta= -100.0;
    double _bestWhadPhi= -10.0;
    
    double _bestTlepPz = -100000.0;
    double _bestTlepPt = -100.0;
    double _bestTlepE  = -100.0;
    double _bestTlepM  = -1.0;
    double _bestTlepEta= -100.0;
    double _bestTlepPhi= -10.0;
    
    double _bestThadPz = -100000.0;
    double _bestThadPt = -100.0;
    double _bestThadE  = -100.0;
    double _bestThadM  = -1.0;
    double _bestThadEta= -100.0;
    double _bestThadPhi= -10.0;
    
    double _bestBlepPz = -100000.0;
    double _bestBlepPt = -100.0;
    double _bestBlepE  = -100.0;
    double _bestBlepEta= -100.0;
    double _bestBlepPhi= -10.0;
    
    double _bestBhadPz = -100000.0;
    double _bestBhadPt = -100.0;
    double _bestBhadE  = -100.0;
    double _bestBhadEta= -100.0;
    double _bestBhadPhi= -10.0;
    
    double _bestJet1Pz = -100000.0;
    double _bestJet1Pt = -100.0;
    double _bestJet1E  = -100.0;
    double _bestJet1Eta= -100.0;
    double _bestJet1Phi= -10.0;
    
    double _bestJet2Pz = -100000.0;
    double _bestJet2Pt = -100.0;
    double _bestJet2E  = -100.0;
    double _bestJet2Eta= -100.0;
    double _bestJet2Phi= -10.0;
    
    //_____ neutrino solution
    int _neuSuccess = 0;
    
    math::XYZTLorentzVector lv_mu = vLep;
    math::XYZTLorentzVector lv_met = vMet;
    
    std::vector<double> pPz = GetNeutrinoPz(lv_mu, lv_met, _neuSuccess);
    //std::cout << mLegend << "first  solution for neutrino pz: " << pPz[0] << std::endl;
    //std::cout << mLegend << "second solution for neutrino pz: " << pPz[1] << std::endl;
    
    SetValue("neutrinoSuccess", _neuSuccess); // were able to reconstruct neutrino
    
    
    //_____ loop over two neutrino solutions
    for (unsigned int n=0; n<2; ++n){
        
        math::XYZTLorentzVector lv_neu(lv_met);
        lv_neu.SetPz(pPz[n]);
        lv_neu.SetE(sqrt(pPz[n]*pPz[n]+lv_met.pt()*lv_met.pt()));
        //std::cout << mLegend << "neutrino mass squared: " << lv_neu.M2() << std::endl;
        
        
        //_____ leptonic W
        math::XYZTLorentzVector lv_Wlep = lv_mu+lv_neu;
        double _mWlep = sqrt(lv_Wlep.M2());
        
        
        //_____ loop over all two-light-jet combinations
        for (unsigned int i=0; i<vLJets.size(); ++i){
            for (unsigned int j=i+1; j<vLJets.size(); ++j){
                // i,j make up the hadronic W
                
                //std::cout << mLegend << "total/first/second: " << vLightJets.size()
                //	  << " / " << i << " / " << j << std::endl;
                
                math::XYZTLorentzVector lv_jet1 = vLJets[i];
                math::XYZTLorentzVector lv_jet2 = vLJets[j];
                
                
                //_____ hadronic W
                math::XYZTLorentzVector lv_Whad = lv_jet1+lv_jet2;
                double _mWhad = sqrt(lv_Whad.M2());
                
                
                //_____ loop over two b jets combinations
                for (unsigned int k=0; k<2; ++k){
                    // loop over the two b jets
                    // k - lepton-side b jet
                    
                    math::XYZTLorentzVector lv_blep = vBJets[k];
                    math::XYZTLorentzVector lv_bhad = vBJets[1-k];
                    
                    
                    //_____ hadronic top
                    math::XYZTLorentzVector lv_Thad = lv_Whad+lv_bhad;
                    double _mThad = sqrt(lv_Thad.M2());
                    
                    
                    //_____ leptonic top
                    math::XYZTLorentzVector lv_Tlep = lv_Wlep+lv_blep;
                    double _mTlep = sqrt(lv_Tlep.M2());
                    
                    
                    //_____ candidate likelihood
                    double _candL = GetLikelihood(_mWlep, _mWhad, _mTlep, _mThad);
                    
                    //std::cout << mLegend << "lep W mass: " << _mWlep << std::endl;
                    //std::cout << mLegend << "candidate likelihood: " << _candL << std::endl;
                    
                    
                    //_____ update values if better candidate found
                    if (_candL>_bestL){
                        
                        _bestL = _candL;
                        
                        
                        // lepton+b(lep)
                        math::XYZTLorentzVector lv_lb = lv_mu+lv_blep;
                        _bestMlb = sqrt(lv_lb.M2());
                        
                        
                        // lepton+b(had)
                        math::XYZTLorentzVector lv_lbhad = lv_mu+lv_bhad;
                        _bestMlbhad = sqrt(lv_lbhad.M2());
                        
                        
                        // "subsmin" variable:
                        // JHEP 1106 (2011) 041
                        math::XYZTLorentzVector lv_sub = lv_mu +lv_blep+lv_bhad+lv_jet1+lv_jet2;
                        math::XYZTLorentzVector lv_tot = lv_sub+lv_neu;
                        _bestSubSmin = sqrt(
                                            ( sqrt(lv_sub.M2()+lv_sub.pt()*lv_sub.pt())
                                             + 
                                             sqrt(lv_neu.M2()+lv_neu.pt()*lv_neu.pt()) )
                                            *
                                            ( sqrt(lv_sub.M2()+lv_sub.pt()*lv_sub.pt()) 
                                             + 
                                             sqrt(lv_neu.M2()+lv_neu.pt()*lv_neu.pt()) )
                                            -
                                            lv_tot.pt()*lv_tot.pt() 
                                            );
                        
                        
                        _bestNuPz   = lv_neu.pz();
                        _bestNuPt   = lv_neu.pt();
                        _bestNuE    = lv_neu.E();
                        _bestNuM2   = lv_neu.M2();
                        _bestNuEta  = lv_neu.Eta();
                        _bestNuPhi  = lv_neu.Phi();
                        
                        _bestWlepPz = lv_Wlep.pz();
                        _bestWlepPt = lv_Wlep.pt();
                        _bestWlepE  = lv_Wlep.E();
                        _bestWlepM  = _mWlep;
                        _bestWlepEta= lv_Wlep.Eta();
                        _bestWlepPhi= lv_Wlep.Phi();
                        
                        _bestWhadPz = lv_Whad.pz();
                        _bestWhadPt = lv_Whad.pt();
                        _bestWhadE  = lv_Whad.E();
                        _bestWhadM  = _mWhad;
                        _bestWhadEta= lv_Whad.Eta();
                        _bestWhadPhi= lv_Whad.Phi();
                        
                        _bestTlepPz = lv_Tlep.pz();	
                        _bestTlepPt = lv_Tlep.pt();	
                        _bestTlepE  = lv_Tlep.E();	
                        _bestTlepM  = _mTlep;	
                        _bestTlepEta= lv_Tlep.Eta();
                        _bestTlepPhi= lv_Tlep.Phi();
                        
                        _bestThadPz = lv_Thad.pz();	
                        _bestThadPt = lv_Thad.pt();	
                        _bestThadE  = lv_Thad.E();	
                        _bestThadM  = _mThad;	
                        _bestThadEta= lv_Thad.Eta();
                        _bestThadPhi= lv_Thad.Phi();
                        
                        _bestBlepPz = lv_blep.pz();	
                        _bestBlepPt = lv_blep.pt();	
                        _bestBlepE  = lv_blep.E();	
                        _bestBlepEta= lv_blep.Eta();
                        _bestBlepPhi= lv_blep.Phi();
                        
                        _bestBhadPz = lv_bhad.pz();	
                        _bestBhadPt = lv_bhad.pt();	
                        _bestBhadE  = lv_bhad.E();	
                        _bestBhadEta= lv_bhad.Eta();
                        _bestBhadPhi= lv_bhad.Phi();
                        
                        _bestJet1Pz = lv_jet1.pz();	
                        _bestJet1Pt = lv_jet1.pt();	
                        _bestJet1E  = lv_jet1.E();	
                        _bestJet1Eta= lv_jet1.Eta();
                        _bestJet1Phi= lv_jet1.Phi();
                        
                        _bestJet2Pz = lv_jet2.pz();	
                        _bestJet2Pt = lv_jet2.pt();	
                        _bestJet2E  = lv_jet2.E();	
                        _bestJet2Eta= lv_jet2.Eta();
                        _bestJet2Phi= lv_jet2.Phi();
                        
                    }
                    else{ // should never happen
                        if (_candL<0) std::cout << mLegend << "cand likelihood: " << _candL << std::endl;
                    }
                    
                    
                } // end of b jet loop
                
            } // inner light jet loop
        } // outer light jet loop
        
    } // neutrino loop
    
    char buf[128];
    sprintf(buf, "bestL%s", suffix.c_str());
    SetValue(buf, _bestL);
    //  SetValue("bestL", _bestL);
    
    sprintf(buf, "bestMlb%s", suffix.c_str());
    SetValue(buf, _bestMlb);
    sprintf(buf, "bestMlbhad%s", suffix.c_str());
    SetValue(buf, _bestMlbhad);
    
    sprintf(buf, "bestSubSmin%s", suffix.c_str());
    SetValue(buf, _bestSubSmin);
    
    sprintf(buf, "bestNuPz%s", suffix.c_str());
    SetValue(buf, _bestNuPz);
    sprintf(buf, "bestNuPt%s", suffix.c_str());
    SetValue(buf, _bestNuPt);
    sprintf(buf, "bestNuE%s", suffix.c_str());
    SetValue(buf, _bestNuE);
    sprintf(buf, "bestNuM2%s", suffix.c_str());
    SetValue(buf, _bestNuM2);
    sprintf(buf, "bestNuEta%s", suffix.c_str());
    SetValue(buf, _bestNuEta);
    sprintf(buf, "bestNuPhi%s", suffix.c_str());
    SetValue(buf, _bestNuPhi);
    
    sprintf(buf, "bestWlepPz%s", suffix.c_str());
    SetValue(buf, _bestWlepPz);
    sprintf(buf, "bestWlepPt%s", suffix.c_str());
    SetValue(buf, _bestWlepPt);
    sprintf(buf, "bestWlepE%s", suffix.c_str());
    SetValue(buf, _bestWlepE);
    sprintf(buf, "bestWlepM%s", suffix.c_str());
    SetValue(buf, _bestWlepM);
    sprintf(buf, "bestWlepEta%s", suffix.c_str());
    SetValue(buf, _bestWlepEta);
    sprintf(buf, "bestWlepPhi%s", suffix.c_str());
    SetValue(buf, _bestWlepPhi);
    
    sprintf(buf, "bestWhadPz%s", suffix.c_str());
    SetValue(buf, _bestWhadPz);
    sprintf(buf, "bestWhadPt%s", suffix.c_str());
    SetValue(buf, _bestWhadPt);
    sprintf(buf, "bestWhadE%s", suffix.c_str());
    SetValue(buf, _bestWhadE);
    sprintf(buf, "bestWhadM%s", suffix.c_str());
    SetValue(buf, _bestWhadM);
    sprintf(buf, "bestWhadEta%s", suffix.c_str());
    SetValue(buf, _bestWhadEta);
    sprintf(buf, "bestWhadPhi%s", suffix.c_str());
    SetValue(buf, _bestWhadPhi);
    
    sprintf(buf, "bestTlepPz%s", suffix.c_str());
    SetValue(buf, _bestTlepPz);
    sprintf(buf, "bestTlepPt%s", suffix.c_str());
    SetValue(buf, _bestTlepPt);
    sprintf(buf, "bestTlepE%s", suffix.c_str());
    SetValue(buf, _bestTlepE);
    sprintf(buf, "bestTlepM%s", suffix.c_str());
    SetValue(buf, _bestTlepM);
    sprintf(buf, "bestTlepEta%s", suffix.c_str());
    SetValue(buf, _bestTlepEta);
    sprintf(buf, "bestTlepPhi%s", suffix.c_str());
    SetValue(buf, _bestTlepPhi);
    
    sprintf(buf, "bestThadPz%s", suffix.c_str());
    SetValue(buf, _bestThadPz);
    sprintf(buf, "bestThadPt%s", suffix.c_str());
    SetValue(buf, _bestThadPt);
    sprintf(buf, "bestThadE%s", suffix.c_str());
    SetValue(buf, _bestThadE);
    sprintf(buf, "bestThadM%s", suffix.c_str());
    SetValue(buf, _bestThadM);
    sprintf(buf, "bestThadEta%s", suffix.c_str());
    SetValue(buf, _bestThadEta);
    sprintf(buf, "bestThadPhi%s", suffix.c_str());
    SetValue(buf, _bestThadPhi);
    
    sprintf(buf, "bestBlepPz%s", suffix.c_str());
    SetValue(buf, _bestBlepPz);
    sprintf(buf, "bestBlepPt%s", suffix.c_str());
    SetValue(buf, _bestBlepPt);
    sprintf(buf, "bestBlepE%s", suffix.c_str());
    SetValue(buf, _bestBlepE);
    sprintf(buf, "bestBlepEta%s", suffix.c_str());
    SetValue(buf, _bestBlepEta);
    sprintf(buf, "bestBlepPhi%s", suffix.c_str());
    SetValue(buf, _bestBlepPhi);
    
    sprintf(buf, "bestBhadPz%s", suffix.c_str());
    SetValue(buf, _bestBhadPz);
    sprintf(buf, "bestBhadPt%s", suffix.c_str());
    SetValue(buf, _bestBhadPt);
    sprintf(buf, "bestBhadE%s", suffix.c_str());
    SetValue(buf, _bestBhadE);
    sprintf(buf, "bestBhadEta%s", suffix.c_str());
    SetValue(buf, _bestBhadEta);
    sprintf(buf, "bestBhadPhi%s", suffix.c_str());
    SetValue(buf, _bestBhadPhi);
    
    sprintf(buf, "bestJet1Pz%s", suffix.c_str());
    SetValue(buf, _bestJet1Pz);
    sprintf(buf, "bestJet1Pt%s", suffix.c_str());
    SetValue(buf, _bestJet1Pt);
    sprintf(buf, "bestJet1E%s", suffix.c_str());
    SetValue(buf, _bestJet1E);
    sprintf(buf, "bestJet1Eta%s", suffix.c_str());
    SetValue(buf, _bestJet1Eta);
    sprintf(buf, "bestJet1Phi%s", suffix.c_str());
    SetValue(buf, _bestJet1Phi);
    
    sprintf(buf, "bestJet2Pz%s", suffix.c_str());
    SetValue(buf, _bestJet2Pz);
    sprintf(buf, "bestJet2Pt%s", suffix.c_str());
    SetValue(buf, _bestJet2Pt);
    sprintf(buf, "bestJet2E%s", suffix.c_str());
    SetValue(buf, _bestJet2E);
    sprintf(buf, "bestJet2Eta%s", suffix.c_str());
    SetValue(buf, _bestJet2Eta);
    sprintf(buf, "bestJet2Phi%s", suffix.c_str());
    SetValue(buf, _bestJet2Phi);
    
    return;
}
