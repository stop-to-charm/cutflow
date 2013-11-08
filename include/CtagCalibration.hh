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

struct JetTagFactorInputs { 
  double pt; 
  double eta; 
  double anti_b; 
  double anti_u; 
  ctag::Flavor flavor; 
};

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
  JetTagSF scale_factor(const JetTagFactorInputs& jet_tf_inputs) const; 
  bool pass(const JetTagFactorInputs& jet_tf_inputs, 
	    ctag::OperatingPoint) const; 
private: 
  typedef std::pair<double, double> CalResult; 
  typedef Analysis::CalibrationDataInterfaceROOT CDI; 
  typedef Analysis::CalibrationDataVariables CalVars; 
  typedef std::map<ctag::OperatingPoint, std::string> OPStrings; 
  typedef std::map<ctag::OperatingPoint, double> CutValue; 
  typedef std::pair<ctag::Flavor,ctag::OperatingPoint> FOPIndex; 
  void check_cdi() const; 
  std::string get_op(ctag::OperatingPoint) const; 
  CalVars get_vars(double pt, double eta) const; 
  std::string get_label(ctag::Flavor) const; 
  void set_indices(ctag::Flavor, ctag::OperatingPoint); 
  unsigned get_index(const std::map<FOPIndex, unsigned>&,  
		     ctag::Flavor, ctag::OperatingPoint) const; 
  CDI* m_cdi; 
  OPStrings m_op_strings; 
  std::string m_jet_author; 
  CutValue m_anti_u_cuts; 
  CutValue m_anti_b_cuts; 
  std::map<FOPIndex, unsigned> m_flav_op_eff_index; 
  std::map<FOPIndex, unsigned> m_flav_op_sf_index; 
}; 


#endif 
