/*
  FWLite port of ElectroWeakAnalysis/Utilities/src/PdfWeightProducer

   Author: Gena Kukartsev, 2013
*/


#include "LJMet/Com/interface/LjmetPdfWeightProducer.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"



namespace LHAPDF {
      void initPDFSet(int nset, const std::string& filename, int member=0);
      int numberPDF(int nset);
      void usePDFMember(int nset, int member);
      double xfx(int nset, double x, double Q, int fl);
      double getXmin(int nset, int member);
      double getXmax(int nset, int member);
      double getQ2min(int nset, int member);
      double getQ2max(int nset, int member);
      void extrapolate(bool extrapolate=true);
}



/////////////////////////////////////////////////////////////////////////////////////

LjmetPdfWeightProducer::LjmetPdfWeightProducer(const edm::ParameterSet& pset) :
 fixPOWHEG_(pset.getUntrackedParameter<std::string> ("FixPOWHEG", "")),
 //genTag_(pset.getUntrackedParameter<edm::InputTag> ("GenTag", edm::InputTag("genParticles"))),
 genTag_(pset.getUntrackedParameter<edm::InputTag> ("GenTag", edm::InputTag("prunedGenParticles"))),
 pdfInfoTag_(pset.getUntrackedParameter<edm::InputTag> ("PdfInfoTag", edm::InputTag("generator"))),
 pdfSetNames_(pset.getUntrackedParameter<std::vector<std::string> > ("PdfSetNames"))
{
      if (fixPOWHEG_ != "") pdfSetNames_.insert(pdfSetNames_.begin(),fixPOWHEG_);

      if (pdfSetNames_.size()>3) {
            edm::LogWarning("") << pdfSetNames_.size() << " PDF sets requested on input. Using only the first 3 sets and ignoring the rest!!";
            pdfSetNames_.erase(pdfSetNames_.begin()+3,pdfSetNames_.end());
      }

      for (unsigned int k=0; k<pdfSetNames_.size(); k++) {
            size_t dot = pdfSetNames_[k].find_first_of('.');
            size_t underscore = pdfSetNames_[k].find_first_of('_');
            if (underscore<dot) {
                  pdfShortNames_.push_back(pdfSetNames_[k].substr(0,underscore));
            } else {
                  pdfShortNames_.push_back(pdfSetNames_[k].substr(0,dot));
            }
            //produces<std::vector<double> >(pdfShortNames_[k].data());
      }
} 



/////////////////////////////////////////////////////////////////////////////////////

LjmetPdfWeightProducer::~LjmetPdfWeightProducer(){}



/////////////////////////////////////////////////////////////////////////////////////

void LjmetPdfWeightProducer::beginJob() {
      for (unsigned int k=1; k<=pdfSetNames_.size(); k++) {
            LHAPDF::initPDFSet(k,pdfSetNames_[k-1]);
      }
}



/////////////////////////////////////////////////////////////////////////////////////
std::map<std::string,std::vector<double> >
LjmetPdfWeightProducer::produce(edm::Event const & event) {

  // FWLite
  std::map<std::string,std::vector<double> > result;

      if (event.isRealData()) return result;

      edm::Handle<GenEventInfoProduct> pdfstuff;
      if (!event.getByLabel(pdfInfoTag_, pdfstuff)) {
            edm::LogError("LjmetPdfWeightProducer") << ">>> PdfInfo not found: " << pdfInfoTag_.encode() << " !!!";
            return result;
      }

      float Q = pdfstuff->pdf()->scalePDF;

      int id1 = pdfstuff->pdf()->id.first;
      double x1 = pdfstuff->pdf()->x.first;
      double pdf1 = pdfstuff->pdf()->xPDF.first;

      int id2 = pdfstuff->pdf()->id.second;
      double x2 = pdfstuff->pdf()->x.second;
      double pdf2 = pdfstuff->pdf()->xPDF.second; 

      // debug
      //std::cout << "GenEventInfoProduct:    Q: " << Q << std::endl;
      //std::cout << "GenEventInfoProduct:  id1: " << id1 << std::endl;
      //std::cout << "GenEventInfoProduct:   x1: " << x1 << std::endl;
      //std::cout << "GenEventInfoProduct: pdf1: " << pdf1 << std::endl;
      //std::cout << "GenEventInfoProduct:  id2: " << id2 << std::endl;
      //std::cout << "GenEventInfoProduct:   x2: " << x2 << std::endl;
      //std::cout << "GenEventInfoProduct: pdf2: " << pdf2 << std::endl;

      // Ad-hoc fix for POWHEG
      if (fixPOWHEG_!="") {
            edm::Handle<reco::GenParticleCollection> genParticles;
            if (!event.getByLabel(genTag_, genParticles)) {
                  edm::LogError("LjmetPdfWeightProducer") << ">>> genParticles  not found: " << genTag_.encode() << " !!!";
                  return result;
            }
            unsigned int gensize = genParticles->size();
            double mboson = 0.;
            for(unsigned int i = 0; i<gensize; ++i) {
                  const reco::GenParticle& part = (*genParticles)[i];
                  int status = part.status();
                  if (status!=3) continue;
                  int id = part.pdgId();
                  if (id!=23 && abs(id)!=24) continue;
                  mboson = part.mass();
                  break;
            }
            Q = sqrt(mboson*mboson+Q*Q);
            LHAPDF::usePDFMember(1,0);
            pdf1 = LHAPDF::xfx(1, x1, Q, id1)/x1;
            pdf2 = LHAPDF::xfx(1, x2, Q, id2)/x2;
      }
      
      // GENA: nominal PDFs in GenInfo seem missing, so redo them
      //LHAPDF::usePDFMember(1,0);
      //pdf1 = LHAPDF::xfx(1, x1, Q, id1)/x1;
      //pdf2 = LHAPDF::xfx(1, x2, Q, id2)/x2;


      // Put PDF weights in the event
      for (unsigned int k=1; k<=pdfSetNames_.size(); ++k) {
	
	// GENA: nominal PDFs in GenInfo seem missing, so redo them
	LHAPDF::usePDFMember(k,0);
	pdf1 = LHAPDF::xfx(k, x1, Q, id1)/x1;
	pdf2 = LHAPDF::xfx(k, x2, Q, id2)/x2;

            std::auto_ptr<std::vector<double> > weights (new std::vector<double>);
            unsigned int nweights = 1;
            if (LHAPDF::numberPDF(k)>1) nweights += LHAPDF::numberPDF(k);
            weights->reserve(nweights);
        
	    // debug
	    //std::cout << "PDF short name: " << pdfShortNames_[k-1] << std::endl;

            for (unsigned int i=0; i<nweights; ++i) {
                  LHAPDF::usePDFMember(k,i);
                  double newpdf1 = LHAPDF::xfx(k, x1, Q, id1)/x1;
                  double newpdf2 = LHAPDF::xfx(k, x2, Q, id2)/x2;
                  weights->push_back(newpdf1/pdf1*newpdf2/pdf2);
		  // debug
		  //std::cout << "PDF weight: " << newpdf1/pdf1*newpdf2/pdf2 << std::endl;
            }

	    // FWLite
            //event.put(weights,pdfShortNames_[k-1]);
	    result.insert(std::pair<std::string,std::vector<double> >(pdfShortNames_[k-1], *weights));
      }

      
      return result;

}
