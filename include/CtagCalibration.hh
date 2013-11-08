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

class CtagCalibration : boost::noncopyable 
{
public: 
  typedef std::pair<double, double> CalResult; 
  CtagCalibration(std::string calibration_file, std::string file_path); 
  ~CtagCalibration(); 
  CalResult applied_factor(const JetTagFactorInputs& jet_tf_inputs,  
			   ctag::OperatingPoint tagger, 
			   ctag::Uncertainty = ctag::Total) const; 
  bool pass_anti_u(const JetTagFactorInputs& jet_tf_inputs, 
		   ctag::OperatingPoint) const; 
  bool pass_anti_b(const JetTagFactorInputs& jet_tf_inputs, 
		   ctag::OperatingPoint) const; 
private: 
  CalResult pass_factor(double pt, double eta, 
			ctag::Flavor flavor, 
			ctag::OperatingPoint tagger, 
			ctag::Uncertainty = ctag::Total) const;
  CalResult fail_factor(double pt, double eta, 
			ctag::Flavor flavor,
			ctag::OperatingPoint tagger, 
			ctag::Uncertainty = ctag::Total) const;
  typedef Analysis::CalibrationDataInterfaceROOT CDI; 
  typedef Analysis::CalibrationDataVariables CalVars; 
  typedef std::map<ctag::OperatingPoint, CDI*> Interfaces;
  typedef std::map<ctag::OperatingPoint, std::string> Names; 
  typedef std::map<ctag::OperatingPoint, double> CutValue; 
  void check_cdi() const; 
  std::string get_op(ctag::OperatingPoint) const; 
  CDI* get_cdi(ctag::OperatingPoint) const; 
  CalVars get_vars(double pt, double eta) const; 
  std::string get_label(ctag::Flavor) const; 
  void set_indices(ctag::Flavor, ctag::OperatingPoint); 
  CDI* m_cnn; 
  Interfaces m_interfaces; 
  Names m_ops; 
  std::string m_jet_author; 
  CutValue m_anti_u_cuts; 
  CutValue m_anti_b_cuts; 
  typedef std::pair<ctag::Flavor,ctag::OperatingPoint> FOPIndex; 
  std::map<FOPIndex, unsigned> m_flav_op_eff_index; 
  std::map<FOPIndex, unsigned> m_flav_op_sf_index; 
}; 


#endif 
