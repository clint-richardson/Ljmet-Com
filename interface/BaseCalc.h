#ifndef LJMet_Com_interface_BaseCalc_h
#define LJMet_Com_interface_BaseCalc_h

/*
 Base class for all calculators
 
 Author: Gena Kukartsev, 2012
 */

#include <iostream>
#include <vector>

#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
//#include "LJMet/Com/interface/VVString.h"

class BaseEventSelector;
class LjmetEventContent;

namespace edm {
    class EventBase;
}

class BaseCalc {
    //
    // Base class for all calculators
    //
    
    friend class LjmetFactory;
    
public:
    BaseCalc();
    virtual ~BaseCalc() { }
    BaseCalc(const BaseCalc &); // stop default
    std::string GetName() { return mName; }
    virtual int BeginJob() = 0;
    virtual int ProduceEvent(edm::EventBase const & event, BaseEventSelector * selector) { return 0; }
    virtual int AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector) { return 0; }
    virtual int EndJob() { return 0; }
    
    std::string mName;
    std::string mLegend;
    
    // LJMET event content setters
    /// Declare a new histogram to be created for the module
    void SetHistogram(std::string name, int nbins, double low, double high);
    void SetHistValue(std::string name, double value);
    void SetValue(std::string name, bool value);
    void SetValue(std::string name, int value);
    void SetValue(std::string name, long long value);
    void SetValue(std::string name, double value);
    void SetValue(std::string name, std::vector<bool> value);
    void SetValue(std::string name, std::vector<int> value);
    void SetValue(std::string name, std::vector<double> value);
    void SetValue(std::string name, std::vector<std::string> value);

protected:
    edm::ParameterSet mPset;
    
private:
    /// Private init method to be called by LjmetFactory when registering the calculator
    virtual void init();
    void setName(std::string name) { mName = name; }
    void SetEventContent(LjmetEventContent * pEc) { mpEc = pEc; }
    void SetPSet(edm::ParameterSet pset) { mPset = pset; }
    LjmetEventContent * mpEc;
};

#endif
