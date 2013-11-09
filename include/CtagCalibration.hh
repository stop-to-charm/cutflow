#ifndef CTAG_CALIBRATION_HH
#define CTAG_CALIBRATION_HH

#include <boost/noncopyable.hpp>
#include <string> 
#include <map>
#include "ctag_defs.hh"

class BaselineJet; 

namespace Analysis { 
  class CalibrationDataInterfaceROOT; 
  class CalibrationDataVariables; 
}

// input structure for both cuts and SF
struct JetTagFactorInputs { 
  double pt; 			// in MeV
  double eta; 			
  double anti_b; 		// log(pc/pb)
  double anti_u; 		// log(pc/pu)
  // these are the flavor_truth_label values, translated to enums
  ctag::Flavor flavor; 		// B, C, U, T, or DATA
};

// translator from the int value in D3PDs to enums used here
ctag::Flavor get_flavor(int flavor_truth_label); 

// the output structure. 'up' and 'down' variations may be symmetric in 
// many cases, but we keep both for flexibility. 
struct JetTagSF { 
  JetTagSF(const std::pair<double, double>&); 
  JetTagSF(); 
  double nominal; 
  double up; 
  double down; 
}; 

class CtagCalibration : boost::noncopyable 
{
public: 
  CtagCalibration(std::string calibration_file, std::string file_path); 
  ~CtagCalibration(); 
  // We should apply the product of the scale factors for every jet in 
  // the event. This is essentially the 'continuous' approach, where 
  // the scale factors depend on the tagger
  JetTagSF scale_factor(const JetTagFactorInputs& jet_tf_inputs) const; 
  bool pass(const JetTagFactorInputs& jet_tf_inputs, 
	    ctag::OperatingPoint) const; 
private: 
  typedef std::pair<double, double> CalResult; 
  typedef std::map<ctag::OperatingPoint, std::string> OPStrings; 
  typedef std::map<ctag::OperatingPoint, double> CutValue; 
  typedef std::pair<ctag::Flavor,ctag::OperatingPoint> FOPIndex; 
  void check_cdi() const; 
  Analysis::CalibrationDataVariables get_vars(double pt, double eta) const; 
  std::string get_label(ctag::Flavor) const; 
  void set_indices(ctag::Flavor, ctag::OperatingPoint); 
  unsigned get_index(const std::map<FOPIndex, unsigned>&,  
		     ctag::Flavor, ctag::OperatingPoint) const; 
  Analysis::CalibrationDataInterfaceROOT* m_cdi; 
  OPStrings m_op_strings; 
  std::string m_jet_author; 
  CutValue m_anti_u_cuts; 
  CutValue m_anti_b_cuts; 
  std::map<FOPIndex, unsigned> m_flav_op_eff_index; 
  std::map<FOPIndex, unsigned> m_flav_op_sf_index; 
}; 


#endif 
