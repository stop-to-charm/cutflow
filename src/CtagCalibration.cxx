#include "CtagCalibration.hh"
#include "CalibrationDataInterface/CalibrationDataInterfaceROOT.h"
#include <stdexcept> 
#include <cassert> 

Analysis::Uncertainty get_unct(btag::Uncertainty); 

CtagCalibration::CtagCalibration(std::string clb_file, 
				 std::string file_path): 
  m_cnn(0), 
  m_jet_author("AntiKt4TopoLCJVF")
{
  using namespace ctag; 
  if (clb_file.size() && file_path.size()) { 
    m_cnn = new CDI("JetFitterCOMBCharm", clb_file, file_path); 
  }
  m_ops[CNN_LOOSE] = "-1_0_0_0"; 
  m_ops[CNN_MEDIUM] = "-1_0_-0_82"; 
  m_ops[CNN_TIGHT] = "-1_0_1_0"; 
  m_interfaces[CNN_LOOSE] = m_cnn; 
  m_interfaces[CNN_MEDIUM] = m_cnn; 
  m_interfaces[CNN_TIGHT] = m_cnn; 

  // ACHTUNG: these are hacks until we get a better CDI
  m_ops[JFC_LOOSE] = m_ops[CNN_LOOSE]; 
  m_ops[JFC_MEDIUM] = m_ops[CNN_MEDIUM]; 
  m_ops[JFC_TIGHT] = m_ops[CNN_TIGHT]; 
  m_interfaces[JFC_LOOSE] = m_cnn; 
  m_interfaces[JFC_MEDIUM] = m_cnn; 
  m_interfaces[JFC_TIGHT] = m_cnn; 

  if (m_cnn) { 
    check_cdi(); 
  }
  m_anti_u_cuts[CNN_LOOSE] = -999; 
  m_anti_u_cuts[CNN_MEDIUM] = -0.82; 
  m_anti_u_cuts[CNN_TIGHT] = 1.0; 
  m_anti_b_cuts[CNN_LOOSE]  = -1.0; 
  m_anti_b_cuts[CNN_MEDIUM] = -1.0;
  m_anti_b_cuts[CNN_TIGHT]  = -1.0;

  m_anti_u_cuts[JFC_LOOSE] = -999; 
  m_anti_u_cuts[JFC_MEDIUM] = 0.95; 
  m_anti_u_cuts[JFC_TIGHT] = 2.5; 
  m_anti_b_cuts[JFC_LOOSE]  = -0.9; 
  m_anti_b_cuts[JFC_MEDIUM] = -0.9;
  m_anti_b_cuts[JFC_TIGHT]  = -0.9;

  for (auto flav: {B, C, U, T} ) { 
    for (auto op: m_ops) { 
      set_indices(flav, op.first); 
    }
  }

}

CtagCalibration::~CtagCalibration() { 
  delete m_cnn; 
  m_cnn = 0; 
}

CtagCalibration::CalResult 
CtagCalibration::pass_factor(double pt, double eta, 
			      ctag::Flavor flavor, 
			      ctag::OperatingPoint tagger, 
			      ctag::Uncertainty uncert) const { 
  CalVars vars = get_vars(pt, eta); 
  // std::string op = get_op(tagger); 
  // std::string label = get_label(flavor); 
  unsigned sf_index = m_flav_op_sf_index.at({flavor, tagger}); 
  Analysis::Uncertainty unct = get_unct(uncert); 
  return get_cdi(tagger)->getScaleFactor(vars, sf_index, unct);
}
CtagCalibration::CalResult 
CtagCalibration::fail_factor(double pt, double eta, 
			     ctag::Flavor flavor, 
			     ctag::OperatingPoint tagger, 
			     ctag::Uncertainty uncert) const { 
  CalVars vars = get_vars(pt, eta); 
  // std::string op = get_op(tagger); 
  // std::string label = get_label(flavor); 
  unsigned sf_index = m_flav_op_sf_index.at({flavor, tagger}); 
  unsigned eff_index = m_flav_op_eff_index.at({flavor, tagger}); 
  Analysis::Uncertainty unct = get_unct(uncert); 
  return get_cdi(tagger)->getInefficiencyScaleFactor(
    vars, sf_index, eff_index, unct);
}

CtagCalibration::CalResult
CtagCalibration::applied_factor(const JetTagFactorInputs& tf_inputs, 
				ctag::OperatingPoint tagger, 
				ctag::Uncertainty uncertainty) const { 

  if (pass_anti_b(tf_inputs, tagger) && pass_anti_u(tf_inputs, tagger)) { 
    return pass_factor(tf_inputs.pt, tf_inputs.eta, tf_inputs.flavor, 
		       tagger, uncertainty); 
  }
  else { 
    return fail_factor(tf_inputs.pt, tf_inputs.eta, tf_inputs.flavor, 
		       tagger, uncertainty); 
  }
}

bool CtagCalibration::pass_anti_u(const JetTagFactorInputs& tf_inputs, 
				  ctag::OperatingPoint tagger) 
  const { 
  CutValue::const_iterator ucut = m_anti_u_cuts.find(tagger); 
  if (ucut == m_anti_u_cuts.end() ) { 
    throw std::logic_error("asked for undefined tagger in " __FILE__); 
  }
  return tf_inputs.anti_u > ucut->second; 
}
bool CtagCalibration::pass_anti_b(const JetTagFactorInputs& tf_inputs, 
				  ctag::OperatingPoint tagger) 
  const { 
  CutValue::const_iterator bcut = m_anti_b_cuts.find(tagger); 
  if (bcut == m_anti_b_cuts.end() ) { 
    throw std::logic_error("asked for undefined tagger in " __FILE__); 
  }
  return tf_inputs.anti_b > bcut->second; 
}


// private stuff
void CtagCalibration::check_cdi() const { 
  for (Names::const_iterator itr = m_ops.begin(); itr != m_ops.end(); itr++){
    if (! m_cnn->getBinnedScaleFactors(m_jet_author, 
				       get_label(ctag::B), 
				       itr->second)) { 
      throw std::runtime_error("ctag calibration information not found"); 
    }
  }
}

CtagCalibration::CDI* CtagCalibration::get_cdi(ctag::OperatingPoint t)
  const { 
  Interfaces::const_iterator cdi = m_interfaces.find(t); 
  if (cdi == m_interfaces.end()) { 
    throw std::logic_error("didn't find cdi in " __FILE__); 
  }
  return cdi->second; 
}

std::string CtagCalibration::get_op(ctag::OperatingPoint t) const { 
  Names::const_iterator oper = m_ops.find(t); 
  if (oper == m_ops.end()) { 
    throw std::logic_error("didn't find OP in " __FILE__); 
  }
  return oper->second; 
}

Analysis::CalibrationDataVariables CtagCalibration::get_vars(double pt, 
							     double eta) 
  const { 
  CalVars vars; 
  vars.jetAuthor = m_jet_author; 
  vars.jetPt = pt; 
  vars.jetEta = eta; 
  return vars; 
}

std::string CtagCalibration::get_label(ctag::Flavor flavor) const { 
  switch (flavor) { 
  case ctag::B: return std::string("B"); 
  case ctag::C: return std::string("C"); 
  case ctag::U: return std::string("Light"); 
  case ctag::T: return std::string("T"); 
  default: 
    assert(false); 
  }
}

void CtagCalibration::set_indices(ctag::Flavor flav, ctag::OperatingPoint op) 
{ 
  const std::string& label = get_label(flav); 
  FOPIndex ind_key(flav, op); 
  bool ok_sf = m_interfaces.at(op)->retrieveCalibrationIndex(
    label, m_ops.at(op), m_jet_author, true, m_flav_op_sf_index[ind_key]); 
  bool ok_eff = m_interfaces.at(op)->retrieveCalibrationIndex(
    label, m_ops.at(op), m_jet_author, false, m_flav_op_eff_index[ind_key]); 
  if (!ok_eff || !ok_sf) { 
    std::string problem = "problem setting op " + m_ops.at(op) + " " + 
      get_label(flav); 
    throw std::runtime_error(problem); 
  }
}



