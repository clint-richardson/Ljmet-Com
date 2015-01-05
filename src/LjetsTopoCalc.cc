/*
 Calculator for the old D0-style LJets topo variables
 
 Author: Gena Kukartsev, 2012
 */



#include <iostream>

#include "DataFormats/Math/interface/LorentzVector.h"
#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LJetsTopoVars.h"
//#include "LJMet/Com/interface/LJetsTopoVarsNew.h" // needs work on compatibility and refactoring
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"



class LjmetFactory;



class LjetsTopoCalc : public BaseCalc{
    
public:
    
    LjetsTopoCalc();
    virtual ~LjetsTopoCalc(){}
    
    virtual int BeginJob(){return 0;}
    virtual int AnalyzeEvent(edm::Event const & event, BaseEventSelector * selector);
    virtual int EndJob(){return 0;}
    
    
private:
    
    int FillLjetsBranches( std::vector<edm::Ptr<pat::Muon> > const & vTightMuons,
                          std::vector<edm::Ptr<pat::Jet> >  const & vSelJets,
                          edm::Ptr<pat::MET>                const & pMet,
                          bool isMuon );
    
};



//static int reg = LjmetFactory::GetInstance()->Register(new LjetsTopoCalc(), "LjetsTopoCalc");



LjetsTopoCalc::LjetsTopoCalc(){
}



int LjetsTopoCalc::AnalyzeEvent(edm::Event const & event,
                                BaseEventSelector * selector){
    //
    // compute event variables here
    //
    
    
    
    //
    // _____ Get objects from the selector _____________________
    //
    std::vector<edm::Ptr<pat::Jet> >      const & vSelJets = selector->GetSelectedJets();
    std::vector<edm::Ptr<pat::Jet> >      const & vSelBtagJets = selector->GetSelectedBtagJets();
    std::vector<edm::Ptr<pat::Jet> >      const & vAllJets = selector->GetAllJets();
    std::vector<edm::Ptr<pat::Muon> >     const & vSelMuons = selector->GetSelectedMuons();
    std::vector<edm::Ptr<pat::Muon> >     const & vLooseMuons = selector->GetLooseMuons();
    std::vector<edm::Ptr<pat::Electron> > const & vSelElectrons = selector->GetSelectedElectrons();
    edm::Ptr<pat::MET>                    const & pMet = selector->GetMet();
    
    
    
    
    //
    //_____ Old D0-style LJetsTopoVars: mu ______
    //
    // (safety checks inside)
    FillLjetsBranches(vSelMuons,
                      vSelJets,
                      pMet,
                      true);
    
    
    
    return 0;
}



int LjetsTopoCalc::FillLjetsBranches( std::vector<edm::Ptr<pat::Muon> > const & vSelMuons,
                                     std::vector<edm::Ptr<pat::Jet> >  const & vSelJets,
                                     edm::Ptr<pat::MET>                const & pMet,
                                     bool isMuon ){
    //
    // Fills old D0-style LJetsTopoVars
    // Returns 0 if success,
    //        -1 if failed (no needed objects in the event)
    //
    
    int result = -1;
    
    while(1){
        
        // check if there is a muon
        if ( vSelMuons.size() == 0 ) break;
        TLorentzVector tlv_muon( vSelMuons[0]->px(),
                                vSelMuons[0]->py(),
                                vSelMuons[0]->pz(),
                                vSelMuons[0]->energy() );
        
        
        // make a vector of jets
        std::vector<TLorentzVector> tlv_jets;
        for ( std::vector<edm::Ptr<pat::Jet> >::const_iterator iJet = vSelJets.begin();
             iJet != vSelJets.end();
             ++iJet ){
            TLorentzVector tlv_jet( (*iJet)->px(),
                                   (*iJet)->py(),
                                   (*iJet)->pz(),
                                   (*iJet)->energy() );
            tlv_jets.push_back(tlv_jet);
        }
        if (tlv_jets.size() < 4) break;
        
        // check if MET exists
        if ( pMet.isNonnull() && pMet.isAvailable() ){ }
        else break;
        TLorentzVector tlv_met( pMet->px(),
                               pMet->py(),
                               pMet->pz(),
                               pMet->energy() );
        
        
        // topovars calculator
        LJetsTopoVars topovars(tlv_jets, tlv_muon, tlv_met, isMuon);
        
        
        // compute branches
        SetValue("aplanarity", topovars.aplanarity());
        SetValue("centrality", topovars.centrality());
        SetValue("sphericity", topovars.sphericity());
        SetValue("ht", topovars.ht());
        SetValue("htpluslepton", topovars.htpluslepton());
        SetValue("methtpluslepton", topovars.methtpluslepton());
        SetValue("h", topovars.h());
        SetValue("ktMinPrime", topovars.ktMinPrime());
        SetValue("dphiLepMet", topovars.dphiLepMet());
        SetValue("minDijetMass", topovars.minDijetMass());
        SetValue("maxJetEta", topovars.maxJetEta());
        SetValue("Et3", topovars.Et3());
        SetValue("minDijetDeltaR", topovars.minDijetDeltaR());
        SetValue("LeptonJet_DeltaR", topovars.LeptonJet_DeltaR());
        SetValue("Jet1Jet2_DeltaR", topovars.Jet1Jet2_DeltaR());
        SetValue("Jet1Jet2_DeltaPhi", topovars.Jet1Jet2_DeltaPhi());
        SetValue("Jet1Jet2_M", topovars.Jet1Jet2_M());
        SetValue("Jet1Jet2_Pt", topovars.Jet1Jet2_Pt());
        SetValue("Jet1Jet2W_M", topovars.Jet1Jet2W_M());
        SetValue("Jet1Jet2W_Pt", topovars.Jet1Jet2W_Pt());
        SetValue("Hz", topovars.Hz());
        SetValue("HT2", topovars.HT2());
        SetValue("HT2prime", topovars.HT2prime());
        SetValue("W_MT", topovars.W_MT());
        SetValue("W_M", topovars.W_M());
        SetValue("W_Pt", topovars.W_Pt());
        SetValue("DphiJMET", topovars.DphiJMET());
        SetValue("Muon_DeltaR", topovars.Muon_DeltaR());
        SetValue("getHt", topovars.getHt());
        SetValue("getHtp", topovars.getHtp());
        SetValue("getHtpp", topovars.getHtpp());
        SetValue("getHt2", topovars.getHt2());
        SetValue("getHt2p", topovars.getHt2p());
        SetValue("getHt2pp", topovars.getHt2pp());
        SetValue("getHt3", topovars.getHt3());
        SetValue("getHt3p", topovars.getHt3p());
        SetValue("getHt3pp", topovars.getHt3pp());
        SetValue("getCen", topovars.getCen());
        SetValue("getNJW", topovars.getNJW());
        SetValue("getJetEtaMax", topovars.getJetEtaMax());
        SetValue("getMdijetMin", topovars.getMdijetMin());
        SetValue("getMtjets", topovars.getMtjets());
        SetValue("getSqrtsT", topovars.getSqrtsT());
        SetValue("getMtAurelio", topovars.getMtAurelio());
        SetValue("getPzOverHT", topovars.getPzOverHT());
        SetValue("getMevent", topovars.getMevent());
        SetValue("getM123inv", topovars.getM123inv());
        SetValue("getEta2Sum", topovars.getEta2Sum());
        SetValue("getMwRec", topovars.getMwRec());
        SetValue("getH", topovars.getH());
        SetValue("getSph", topovars.getSph());
        SetValue("getApl", topovars.getApl());
        SetValue("getAplMu", topovars.getAplMu());
        SetValue("getKtminp", topovars.getKtminp());
        SetValue("getKtminpReduced", topovars.getKtminpReduced());
        SetValue("getDrMinJetJet", topovars.getDrMinJetJet());
        SetValue("getDphiMuMet", topovars.getDphiMuMet());
        SetValue("getMt", topovars.getMt());
        
        
        result = 0;
        break;
        
    } // end of while(1)
    
    return result;
}
