/*
 Calculator for the old D0-style LJets topo variables
 
 Author: Gena Kukartsev, 2012
 */



#include <iostream>

#include "DataFormats/Math/interface/LorentzVector.h"
#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LJetsTopoVarsNew.h" // needs work on compatibility and refactoring
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

class LjmetFactory;

class LjetsTopoCalcNew : public BaseCalc{
    
public:
    
    LjetsTopoCalcNew();
    virtual ~LjetsTopoCalcNew(){}
    
    virtual int BeginJob(){
        if (mPset.exists("useBestTop")) bestTop_ = mPset.getParameter<bool>("useBestTop");
        else                            bestTop_ = false;
        
        if (mPset.exists("debug"))      debug_ = mPset.getParameter<bool>("debug");
        else                            debug_ = false;
        
        return 0;
    }
    virtual int AnalyzeEvent(edm::Event const & event, BaseEventSelector * selector);
    virtual int EndJob(){return 0;}
    
    
private:
    
    bool bestTop_;
    bool debug_;
    
    int FillLjetsBranches( std::vector<edm::Ptr<pat::Muon> > const & vTightMuons,
                          std::vector<edm::Ptr<pat::Electron> > const & vTightElectrons,
                          std::vector<std::pair<TLorentzVector,bool> >  const & vCorrBtagJets,
                          //edm::Ptr<pat::MET> const & pMet,
                          TLorentzVector const & corrMET,
                          bool isMuon,
                          bool bestTop );
    
};



//static int reg = LjmetFactory::GetInstance()->Register(new LjetsTopoCalcNew(), "LjetsTopoCalcNew");



LjetsTopoCalcNew::LjetsTopoCalcNew(){
    mLegend = "[LjetsTopoCalcNew]: ";
}


int LjetsTopoCalcNew::AnalyzeEvent(edm::Event const & event,
                                   BaseEventSelector * selector){
    //
    // compute event variables here
    //
    
    //
    // _____ Get objects from the selector _____________________
    //
    std::vector<edm::Ptr<pat::Jet> >      const & vSelJets = selector->GetSelectedJets();
    std::vector<edm::Ptr<pat::Jet> >      const & vSelBtagJets = selector->GetSelectedBtagJets();
    std::vector<std::pair<TLorentzVector,bool>> const & vCorrBtagJets = selector->GetCorrJetsWithBTags();
    std::vector<edm::Ptr<pat::Jet> >      const & vAllJets = selector->GetAllJets();
    std::vector<edm::Ptr<pat::Muon> >     const & vSelMuons = selector->GetSelectedMuons();
    std::vector<edm::Ptr<pat::Electron> > const & vSelElectrons = selector->GetSelectedElectrons();
    edm::Ptr<pat::MET>                    const & pMet = selector->GetMet();
    TLorentzVector                        const & corrMET = selector->GetCorrectedMet();
    
    //
    //_____ Old D0-style LJetsTopoVarsNew: mu ______
    //
    // (safety checks inside)
    
    bool muonchannel = false;
    bool electronchannel = false;
    if (vSelMuons.size()>0) muonchannel = true;
    if (vSelElectrons.size()>0) electronchannel = true;
    if (muonchannel && electronchannel && debug_)
        std::cout << mLegend <<"WARNING: Two Leptons in the event"<<std::endl;
    
    if (vSelMuons.size()==0 and vSelElectrons.size()==0) {
        if (debug_)
            std::cout << mLegend << "No Lepton in event! "<<std::endl;
        return 0;
    }
    
    if (vCorrBtagJets.size()==0) {
        if (debug_) std::cout << mLegend <<"No Corrected jets in event! "<<std::endl;
        return 0;
    }
    
    if (corrMET.Pt()==0) {
        if (debug_) std::cout << mLegend <<"No Corrected MET in event! "<<std::endl;
        return 0;
    }
    
    
    FillLjetsBranches(vSelMuons,
                      vSelElectrons,
                      vCorrBtagJets,
                      corrMET,
                      muonchannel, // isMuon
                      bestTop_); //bestTop for neutrino pz
    
    
    return 0;
}



int LjetsTopoCalcNew::FillLjetsBranches( std::vector<edm::Ptr<pat::Muon> > const & vSelMuons,
                                        std::vector<edm::Ptr<pat::Electron> > const & vSelElectrons,
                                        std::vector<std::pair<TLorentzVector,bool>> const & vCorrBtagJets,
                                        //edm::Ptr<pat::MET> const & pMet,
                                        TLorentzVector const & corrMET,
                                        bool isMuon,
                                        bool bestTop ){
    //
    // Fills old D0-style LJetsTopoVarsNew
    // Returns 0 if success,
    //        -1 if failed (no needed objects in the event)
    //
    
    int result = -1;
    
    while(1){
        
        TLorentzVector tlv_lepton;
        
        if ( vSelMuons.size() == 0 && vSelElectrons.size() == 0) break;
        if ( vSelMuons.size() > 0 && vSelElectrons.size() > 0) {
            std::cout<<mLegend<<"Two Leptons, using the Muon...."<<std::endl;
            tlv_lepton.SetPxPyPzE( vSelMuons[0]->px(),
                                  vSelMuons[0]->py(),
                                  vSelMuons[0]->pz(),
                                  vSelMuons[0]->energy() );
        }
        if ( vSelMuons.size() > 0 ) {
            tlv_lepton.SetPxPyPzE( vSelMuons[0]->px(),
                                  vSelMuons[0]->py(),
                                  vSelMuons[0]->pz(),
                                  vSelMuons[0]->energy() );
        }
        if ( vSelElectrons.size() > 0 ) {
            tlv_lepton.SetPxPyPzE( vSelElectrons[0]->px(),
                                  vSelElectrons[0]->py(),
                                  vSelElectrons[0]->pz(),
                                  vSelElectrons[0]->energy() );
        }
        
        // make a vector of jets
        /*
         std::vector<TLorentzVector> tlv_jets;
         for ( std::vector<std::pair<TLorentzVector,bool>>::const_iterator iJet = vCorrBtagJets.begin();
         iJet != vCorrBtagJets.end();
         ++iJet ){
         //
         
         //
         TLorentzVector tlv_jet( (*iJet)->first.px(),
         (*iJet)->first.py(),
         (*iJet)->first.pz(),
         (*iJet)->first.energy() );
         
         tlv_jets.push_back(tlv_jet);
         }
         */
        
        // check if MET exists
        //if ( pMet.isNonnull() && pMet.isAvailable() ){ }
        //else break;
        
        // check if corrected MET exists
        if ( corrMET.Pt() > 0 ){ }
        else break;
        
        TLorentzVector tlv_met( corrMET.Px(),
                               corrMET.Py(),
                               corrMET.Pz(),
                               corrMET.Energy() );
        
        // topovars calculator
        //LJetsTopoVarsNew topovars(tlv_jets, tlv_muon, tlv_met, isMuon, bestTop);
        LJetsTopoVarsNew topovars(vCorrBtagJets, tlv_lepton, tlv_met, isMuon,  bestTop);
        
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
        SetValue("Ht", topovars.getHt());
        SetValue("Htp", topovars.getHtp());
        SetValue("Htpp", topovars.getHtpp());
        SetValue("Ht2", topovars.getHt2());
        SetValue("Ht2p", topovars.getHt2p());
        SetValue("Ht2pp", topovars.getHt2pp());
        SetValue("Ht3", topovars.getHt3());
        SetValue("Ht3p", topovars.getHt3p());
        SetValue("Ht3pp", topovars.getHt3pp());
        SetValue("Centrality", topovars.getCen());
        SetValue("NJW", topovars.getNJW());
        SetValue("JetEtaMax", topovars.getJetEtaMax());
        SetValue("MdijetMin", topovars.getMdijetMin());
        SetValue("Mtjets", topovars.getMtjets());
        SetValue("SqrtsT", topovars.getSqrtsT());
        SetValue("MtAurelio", topovars.getMtAurelio());
        SetValue("PzOverHT", topovars.getPzOverHT());
        SetValue("Mevent", topovars.getMevent());
        SetValue("M123inv", topovars.getM123inv());
        SetValue("Eta2Sum", topovars.getEta2Sum());
        SetValue("MwRec", topovars.getMwRec());
        SetValue("H", topovars.getH());
        SetValue("Sphericity", topovars.getSph());
        SetValue("Aplanarity", topovars.getApl());
        SetValue("AplanarityMu", topovars.getAplMu());
        SetValue("Ktminp", topovars.getKtminp());
        SetValue("KtminpReduced", topovars.getKtminpReduced());
        SetValue("DrMinJetJet", topovars.getDrMinJetJet());
        SetValue("DphiMuMet", topovars.getDphiMuMet());
        SetValue("Mt", topovars.getMt());
        SetValue("dphiLepJ1", topovars.dphiLepJ1());
        SetValue("dphiLepJ2", topovars.dphiLepJ2());
        SetValue("dphiLepJ3", topovars.dphiLepJ3());
        SetValue("dphiLepJ4", topovars.dphiLepJ4());
        SetValue("dphiLepLeadBTagJet", topovars.dphiLepLeadBTagJet());
        SetValue("dphiLepSecLeadBTagJet", topovars.dphiLepSecLeadBTagJet());
        SetValue("dphiLepLightJet", topovars.dphiLepLightJet());
        //SetValue("LeadBTagJet_DiscVal", topovars.LeadBTagJet_DiscVal());
        //SetValue("SecLeadBTagJet_DiscVal", topovars.SecLeadBTagJet_DiscVal());
        //SetValue("Lead12BTagJet_DiscVal", topovars.Lead12BTagJet_DiscVal());
        //SetValue("LightJet_DiscVal", topovars.LightJet_DiscVal());
        SetValue("BestTop", topovars.BestTop());
        SetValue("BestJet_Pt", topovars.BestJet_Pt());
        SetValue("BestJet_Eta", topovars.BestJet_Eta());
        SetValue("BestJet_Phi", topovars.BestJet_Phi());
        SetValue("BestTop_Pt", topovars.BestTop_Pt());
        SetValue("BestTopBJet_Eta", topovars.BestTopBJet_Eta());
        SetValue("BestTopBJet_Phi", topovars.BestTopBJet_Phi());
        SetValue("BestTopBJet_Pt", topovars.BestTopBJet_Pt());
        SetValue("SecBestTop", topovars.SecBestTop());
        SetValue("SecBestBTagTop", topovars.SecBestBTagTop());
        SetValue("BTagTopMass", topovars.BTagTopMass());
        SetValue("BTagTop_Pt", topovars.BTagTop_Pt());
        SetValue("SecBTagTopMass", topovars.SecBTagTopMass());
        SetValue("SecBTagTop_Pt", topovars.SecBTagTop_Pt());
        SetValue("BestJetJet2W_M", topovars.BestJetJet2W_M());
        SetValue("Jet1TagJet2TagW_M", topovars.Jet1TagJet2TagW_M());
        SetValue("Cos_BestJetLepton_BestTop", topovars.Cos_BestJetLepton_BestTop());
        SetValue("Cos_LightjetJetLepton_BestTop", topovars.Cos_LightjetJetLepton_BestTop());
        SetValue("Cos_LightjetJetLepton_BTagTop", topovars.Cos_LightjetJetLepton_BTagTop());
        SetValue("HT_AllJets_MinusBestJet", topovars.HT_AllJets_MinusBestJet());
        SetValue("H_AllJets_MinusBestJet", topovars.H_AllJets_MinusBestJet());
        SetValue("J1_NotBestJet_Eta", topovars.J1_NotBestJet_Eta());
        SetValue("J2_NotBestJet_Eta", topovars.J2_NotBestJet_Eta());
        SetValue("AllJets_MinusBestJet_Pt", topovars.AllJets_MinusBestJet_Pt());
        SetValue("AllJets_M", topovars.AllJets_M());
        SetValue("AllJetsW_M", topovars.AllJetsW_M());//sqrt_shat
        SetValue("J1_NotBestJet_Pt", topovars.J1_NotBestJet_Pt());
        SetValue("J2_NotBestJet_Pt", topovars.J2_NotBestJet_Pt());
        SetValue("J1_NotBestJet_Eta", topovars.J1_NotBestJet_Eta());
        SetValue("J1_NotBestJet_Phi", topovars.J1_NotBestJet_Phi());
        SetValue("Muon_DeltaR", topovars.Muon_DeltaR());
        SetValue("PzNu", topovars.pznu());
        SetValue("PzOtherNu", topovars.pzothernu());
        SetValue("LepTopBJet_DeltaR", topovars.LepTopBJet_DeltaR());
        SetValue("LepTopBJet_Dphi", topovars.LepTopBJet_Dphi());
        SetValue("LepTopBJet_Deta", topovars.LepTopBJet_Deta());
        
        result = 0;
        break;
        
    } // end of while(1)
    
    return result;
}
