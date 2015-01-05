/*
  Calculator for the SLiTT analysis
  
  Author: Gena Kukartsev, 2012
*/



#include <iostream>
#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"

#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "LJMet/Com/interface/TopElectronSelector.h"

#include "DataFormats/JetReco/interface/CATopJetTagInfo.h"
#include "PhysicsTools/SelectorUtils/interface/PFElectronSelector.h"
#include "PhysicsTools/SelectorUtils/interface/PFMuonSelector.h"
#include "TFile.h"
#include "TH1F.h"
#include "LJMet/Com/interface/FileExists.h"

#define IS_BHADRON_PDGID(id) ( ((abs(id)/100)%10 == 5) || (abs(id) >= 5000 && abs(id) <= 5999) )
#define IS_NEUTRINO_PDGID(id) ( (abs(id) == 12) || (abs(id) == 14) || (abs(id) == 16) )

using std::cout;
using std::endl;

class LjmetFactory;
namespace TopInitID{
  static const int status = 3;
  static const int tID    = 6; 
}
namespace TopDecayID{
  /// identification of top decays; used for following
  /// the decay chain in TopDecaySubset
  static const int stable = 2;
  static const int unfrag = 3;
  static const int tID    = 6;
  static const int bID    = 5;
  static const int glueID = 21;
  static const int photID = 22;
  static const int ZID    = 23;
  static const int WID    = 24;
  static const int elecID = 11;
  static const int muonID = 13;
  static const int tauID  = 15;
}

namespace WDecay{
  /// classification of leptons in the decay channel 
  /// of the W boson used in several places throughout 
  /// the package
  enum LeptonType {kNone, kElec, kMuon, kTau};
}
  /// supported modes to fill the new vectors 
  /// of gen particles
  enum  FillMode {kStable, kME};
  /// classification of potential shower types
  enum ShowerModel{kStart=-1, kNone, kPythia, kHerwig};

class TopEventReweightCalc : public BaseCalc{
  
public:
  
  TopEventReweightCalc();
  virtual ~TopEventReweightCalc(){}
  
  virtual int BeginJob();
  virtual int AnalyzeEvent(edm::Event const & event, BaseEventSelector * selector);
  virtual int EndJob(){return 0;}

private:
  
  bool                      isMc;
  edm::InputTag             genParticles_it;
  edm::InputTag genJetsIT_;
  bool ptReweight, reweightBSemiLeptDecyas, reweightBfragmentation;
  double nuDecayFractionSource_;
  double nuDecayFractionTarget_;

  string fragSourceFile_;
  string fragTargetFile_;

  TFile* sourceFile;
  TFile* targetFile;

  TH1F* sourceHist;
  TH1F* targetHist, *weightHist;

  double eventWeightBJES(edm::Handle<reco::GenParticleCollection> &genParticles,
	edm::Handle<std::vector< reco::GenJet > > & genJets);
  void printDaughters(const reco::Candidate * p)const;

};


static int reg = LjmetFactory::GetInstance()->Register(new TopEventReweightCalc(), "TopEventReweightCalc");


TopEventReweightCalc::TopEventReweightCalc(){
}

int TopEventReweightCalc::BeginJob(){
  cout << "TopEventReweightCalc\n";
  if (mPset.exists("isMc"))         isMc = mPset.getParameter<bool>("isMc");
  else                              isMc = false;
  
  if (mPset.exists("genParticles")) genParticles_it = mPset.getParameter<edm::InputTag>("genParticles");
  else                              genParticles_it = edm::InputTag("genParticles");

  if (isMc) {
    if (mPset.exists("ptReweight"))
	 ptReweight = mPset.getParameter<bool>("ptReweight");
    else ptReweight = false;

    if (mPset.exists("reweightBSemiLeptDecyas"))
	 reweightBSemiLeptDecyas = mPset.getParameter<bool>("reweightBSemiLeptDecyas");
    else reweightBSemiLeptDecyas = false;

    if (mPset.exists("reweightBfragmentation"))
	 reweightBfragmentation = mPset.getParameter<bool>("reweightBfragmentation");
    else reweightBfragmentation = false;
    
    if (reweightBSemiLeptDecyas) {
      nuDecayFractionSource_ = mPset.getParameter<double>("nuDecayFractionSource");
      nuDecayFractionTarget_ = mPset.getParameter<double>("nuDecayFractionTarget");
    }
    if (reweightBfragmentation){
      fragSourceFile_ = mPset.getParameter<string>("fragSourceFile");
      fragTargetFile_ = mPset.getParameter<string>("fragTargetFile");
      fexists(fragSourceFile_, true);
      fexists(fragTargetFile_, true);
      genJetsIT_ = edm::InputTag("selectedPatJetsPFlow","genJets");
    TH1::AddDirectory(kFALSE);
      TFile sourceFile(fragSourceFile_.c_str(),"READ");
      if (!sourceFile.IsOpen()) {
	cout << "PUWeighting::fragSourceFile "<< fragSourceFile_ << "could not be openend\n";
	assert(false);
      }
      sourceHist = (TH1F*) sourceFile.Get("EventWeightBJES/genBHadronPtFraction")->Clone();
      cout << fragSourceFile_<<endl;
      sourceHist->Print();
//       sourceFile.Close();
      TFile targetFile(fragTargetFile_.c_str(),"READ");
      if (!targetFile.IsOpen()) {
	cout << "PUWeighting::fragTargetFile "<< fragSourceFile_ << "could not be openend\n";
	assert(false);
      }
      targetHist = (TH1F*) targetFile.Get("EventWeightBJES/genBHadronPtFraction")->Clone();
      cout << fragTargetFile_<<endl;
      targetHist->Print();
//       targetFile.Close();

      if (sourceHist->GetNbinsX() != targetHist->GetNbinsX()) 
        std::cout << "Incompatible b-fragmentation histograms: Number of bins not equal" << std::endl;

      sourceHist->Scale(1./sourceHist->Integral());
      targetHist->Scale(1./targetHist->Integral());

      weightHist = (TH1F*) targetHist->Clone();
      weightHist->Divide(sourceHist);

      std::cout << "Weights for b-fragmentation" << std::endl;
      for (int i = 0; i < sourceHist->GetNbinsX(); ++i) {
      std::cout << 
      weightHist->GetXaxis()->GetBinLowEdge(i)
      <<std::setiosflags(std::ios::left)
      << " [" << i << "] "
      << std::setw(10) << weightHist->GetBinContent(i)<<endl;
      }
      std::cout << std::endl;


    }

  }
  


//   hists["eventWeight"] = fs->make<TH1F>("eventWeight", "eventWeight", 1000, 0, 10 );
//   hists["genBHadronNuDecay"] = fs->make<TH1F>("genBHadronNuDecay", "genBHadronNuDecay", 2, 0, 2);
//   hists["genBHadronPtFraction"] = fs->make<TH1F>("genBHadronPtFraction", "genBHadronPtFraction", 100, 0, 2);
// 
  

  return 0;
}

int TopEventReweightCalc::AnalyzeEvent(edm::Event const & event,
			       BaseEventSelector * selector){

  //
  //_____ Gen Info ______________________________
  //

  //Four vector
  std::vector <double> genPt;
  std::vector <double> genEta;
  std::vector <double> genPhi;
  std::vector <double> genEnergy;

  //Identity
  std::vector <int> genID;
  std::vector <int> genIndex;
  std::vector <int> genStatus;
  std::vector <int> genMotherID;
  std::vector <int> genMotherIndex;

  double eventWeight = 1.0;

  if (isMc){
    
////////////////////////////////////////////////////

    edm::Handle<reco::GenParticleCollection> genParticles;
    event.getByLabel(genParticles_it, genParticles);

    if ( reweightBSemiLeptDecyas || reweightBfragmentation) {
      edm::Handle<std::vector< reco::GenJet > > genJets;
      if (reweightBfragmentation) event.getByLabel(genJetsIT_, genJets);
      eventWeight  = eventWeightBJES(genParticles, genJets);
    }

    for(reco::GenParticleCollection::const_iterator t=genParticles->begin(); t!=genParticles->end(); ++t){
      if( std::abs(t->pdgId())==TopDecayID::tID && t->status()==TopDecayID::unfrag ){
//     cout << " id "<< t->pdgId()<<" "<<t->status();
//     cout << " ok "<< t->pt() << " "<<t->eta()<< " "<<t->phi();
//     cout << " " << exp(0.148-0.00129*t->pt())<<endl;
      if (ptReweight) eventWeight *= exp(0.148-0.00129*t->pt());
	//Four vector
      genPt     . push_back(t->pt());
      genEta    . push_back(t->eta());
      genPhi    . push_back(t->phi());
      genEnergy . push_back(t->energy());

      //Identity
      genID            . push_back(t->pdgId());
// 	genIndex         . push_back((int) i);
      genStatus        . push_back(t->status());
      genMotherID      . push_back(t->mother()->pdgId());
//     cout << endl;  
      }
    }
// 	      cout << "eventWeight: "<<eventWeight<<endl;

  }  //End MC-only if

  //Four vector
  SetValue("genPt"     , genPt);
  SetValue("genEta"    , genEta);
  SetValue("genPhi"    , genPhi);
  SetValue("genEnergy" , genEnergy);

  //Identity
  SetValue("genID"            , genID);
  SetValue("genIndex"         , genIndex);
  SetValue("genStatus"        , genStatus);
  SetValue("genMotherID"      , genMotherID);
  SetValue("genMotherIndex"   , genMotherIndex);
  SetValue("eventWeight",eventWeight);

  return 0;
}

void TopEventReweightCalc::printDaughters(const reco::Candidate * p) const
{
    int n = p->numberOfDaughters();
    for(int j = 0; j < n; ++j) {
      const reco::Candidate * d = p->daughter( j );
      int dauId = d->pdgId();
      cout << dauId<<" ";
      if (d->numberOfDaughters()) {
        cout <<" ( ";
	printDaughters(d);
        cout <<" ) ";
      };
    }
}

double TopEventReweightCalc::eventWeightBJES(edm::Handle<reco::GenParticleCollection> &genParticles,
	edm::Handle<std::vector< reco::GenJet > > &genJets)
{  


  double weight = 1.;

  //////////////////////////////////////////////////////////////////////////
  // GENPARTICLES
  ////////////////////////////////////////////////////////////////////////
    
  for(size_t i = 0; i < genParticles->size(); ++ i) {
    const reco::GenParticle & p = (*genParticles)[i];
    if (p.pt() == 0) continue;
    
    int id = p.pdgId();
    int n = p.numberOfDaughters();
//     if ((abs(id)==5)||(abs(id)==6)) {
//     cout << "Found B hadron " << id <<" "<<n<<endl;
//     printDaughters(&p);
//     cout <<endl<<endl;}
    if (!IS_BHADRON_PDGID(id)) continue;
    
    bool hasBDaughter = false;
    bool hasNuDaughter = false;
//     double genBHadronNuDecay = 0.;
    
    for(int j = 0; j < n; ++j) {
      const reco::Candidate * d = p.daughter( j );
      int dauId = d->pdgId();
//       cout << dauId<<" "<<d->numberOfDaughters()<<" ; ";
      if (IS_BHADRON_PDGID(dauId)) {
        hasBDaughter = true;
        break;
      }
      if (IS_NEUTRINO_PDGID(dauId)) hasNuDaughter = true;
    }
//     cout << endl;
    // Weakly decaying B hadron
    if (!hasBDaughter) {
    
      // Neutrino fraction weight
      if (reweightBSemiLeptDecyas) {
	if (hasNuDaughter) {
          // genBHadronNuDecay = 1.;
          weight *= nuDecayFractionTarget_/nuDecayFractionSource_;
//           cout << "hasNuDaughter " << weight<<endl;
	}
	else {
          weight *= (1.-nuDecayFractionTarget_)/(1.-nuDecayFractionSource_);
//           cout << "not hasNuDaughter " << weight<<endl;
	}
//       hists.find("genBHadronNuDecay")->second->Fill( genBHadronNuDecay );
      }
      
      // Fragmentation weight
      if (reweightBfragmentation) {
	for (std::vector< reco::GenJet >::const_iterator ijet = genJets->begin(); ijet != genJets->end(); ++ijet) {
          if (p.pt() == 0 || ijet->pt() == 0) continue;

          double deta = p.eta() - ijet->eta();
          double dphi = TVector2::Phi_mpi_pi(p.phi() - ijet->phi());
          double dr   = sqrt( deta*deta + dphi*dphi );

          // Simple dR match of hadron and GenJet
          if (dr < 0.5) {
            double xb = p.pt()/ijet->pt();
  //           hists.find("genBHadronPtFraction")->second->Fill(xb);
            if (xb < 2.) {
              weight *= weightHist->GetBinContent(weightHist->GetXaxis()->FindFixBin(xb));
            }
            break;
          }
	}
      }
    }
  }
  return weight;
}

