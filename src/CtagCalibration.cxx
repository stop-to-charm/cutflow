#include "CtagCalibration.hh"
#include "CalibrationDataInterface/CalibrationDataInterfaceROOT.h"
#include <stdexcept> 
#include <cassert> 

// Analysis::Uncertainty get_unct(btag::Uncertainty); 

// fill from the CalibrationDataInterfaceROOT output
JetTagSF::JetTagSF(const std::pair<double, double>& cal_result): 
  nominal(cal_result.first), 
  up(cal_result.first + cal_result.second), 
  down(cal_result.first - cal_result.second)
{}

// fill with NaN
JetTagSF::JetTagSF(): 
  nominal(0.0/0.0), 
  up(0.0/0.0), 
  down(0.0/0.0)
{}

// shorthand for variations
namespace { 
  double up(const std::pair<double, double>& res) { 
    return res.first + res.second;
  } 
  double down(const std::pair<double, double>& res) { 
    return res.first - res.second;
  } 
}

CtagCalibration::CtagCalibration(std::string clb_file, 
				 std::string file_path): 
  m_cdi(new CDI("JetFitterCOMBCharm", clb_file, file_path)),
  m_jet_author("AntiKt4TopoLCJVF")
{
  using namespace ctag; 
  // WARNING: these are hacks until we get a better CDI
  m_op_strings[LOOSE] = "-1_0_0_0"; 
  m_op_strings[MEDIUM] = "-1_0_-0_82"; 

  if (m_cdi) { 
    check_cdi(); 
  }

  m_anti_u_cuts[LOOSE] = -999; 
  m_anti_u_cuts[MEDIUM] = 0.95; 
  m_anti_b_cuts[LOOSE]  = -0.9; 
  m_anti_b_cuts[MEDIUM] = -0.9;

  // flavors happen to start with B and end with DATA
  for (int flavor = ctag::B; flavor < ctag::DATA; flavor++) { 
    for (OPStrings::const_iterator itr = m_op_strings.begin(); 
	 itr != m_op_strings.end(); itr++) { 
      set_indices(static_cast<ctag::Flavor>(flavor), itr->first); 
    }
  }

}

CtagCalibration::~CtagCalibration() { 
  delete m_cdi; 
  m_cdi = 0; 
}

JetTagSF CtagCalibration::scale_factor(
  const JetTagFactorInputs& tf_inputs) const { 

  CalVars vars = get_vars(tf_inputs.pt, tf_inputs.eta); 
  Analysis::Uncertainty unct = Analysis::Total; 

  // hack in some trickery here... 
  // if we fail loose, return the inefficiency SF for loose
  if (!pass(tf_inputs, ctag::LOOSE) ) { 
    unsigned sf_index = get_index(
      m_flav_op_sf_index, tf_inputs.flavor, ctag::LOOSE); 
    unsigned eff_index = get_index(
      m_flav_op_eff_index, tf_inputs.flavor, ctag::LOOSE); 
    
    return JetTagSF(m_cdi->getInefficiencyScaleFactor(
		      vars, sf_index, eff_index, unct));
  }
  // if we pass medium return that SF
  unsigned med_sf_index = get_index(
    m_flav_op_sf_index, tf_inputs.flavor, ctag::MEDIUM); 
  if (pass(tf_inputs, ctag::MEDIUM)) { 
    return JetTagSF(m_cdi->getScaleFactor(vars, med_sf_index, unct));
  }
  // if we pass loose but fail medium we need to be clever
  if (pass(tf_inputs, ctag::LOOSE) && !pass(tf_inputs, ctag::LOOSE)) { 
    unsigned loose_eff_index = get_index(
      m_flav_op_eff_index, tf_inputs.flavor, ctag::LOOSE); 
    unsigned med_eff_index = get_index(
      m_flav_op_eff_index, tf_inputs.flavor, ctag::MEDIUM); 
    unsigned loose_sf_index = get_index(
      m_flav_op_sf_index, tf_inputs.flavor, ctag::LOOSE); 
    unsigned med_sf_index = get_index(
      m_flav_op_sf_index, tf_inputs.flavor, ctag::MEDIUM); 

    CalResult data_med_eff = m_cdi->getEfficiency(
      vars, med_sf_index, med_eff_index, unct); 
    CalResult data_loose_eff = m_cdi->getEfficiency(
      vars, loose_sf_index, loose_eff_index, unct); 
    CalResult mc_med_eff = m_cdi->getMCEfficiency(
      vars, med_eff_index, unct); 
    CalResult mc_loose_eff = m_cdi->getMCEfficiency(
      vars, loose_eff_index, unct); 
    
    JetTagSF sf; 
    sf.nominal = ((data_loose_eff.first - data_med_eff.first) / 
		  (mc_loose_eff.first - mc_med_eff.first)); 

    // not sure of the best way to do the variation...
    sf.up = ( (up(data_loose_eff) - up(data_med_eff) ) / 
	      (up(mc_loose_eff) - up(mc_med_eff))); 
    sf.down = ( (down(data_loose_eff) - down(data_med_eff) ) / 
		(down(mc_loose_eff) - down(mc_med_eff))); 
    return sf; 
  }
  
  // if we're here it's an error
  throw std::logic_error("something went horribly wrong in " __FILE__); 
}

bool CtagCalibration::pass(const JetTagFactorInputs& tf_inputs, 
			   ctag::OperatingPoint op) 
  const { 
  CutValue::const_iterator ucut = m_anti_u_cuts.find(op); 
  if (ucut == m_anti_u_cuts.end() ) { 
    throw std::logic_error("asked for undefined op in " __FILE__); 
  }
  CutValue::const_iterator bcut = m_anti_b_cuts.find(op); 
  if (bcut == m_anti_b_cuts.end() ) { 
    throw std::logic_error("asked for undefined op in " __FILE__); 
  }
  bool pass_u = tf_inputs.anti_u > ucut->second; 
  bool pass_b = tf_inputs.anti_b > bcut->second; 
  return pass_b && pass_u; 
}


// ----- private stuff ---------
void CtagCalibration::check_cdi() const { 
  for (OPStrings::const_iterator itr = m_op_strings.begin(); itr != m_op_strings.end(); itr++){
    if (! m_cdi->getBinnedScaleFactors(m_jet_author, 
				       get_label(ctag::B), 
				       itr->second)) { 
      throw std::runtime_error("ctag calibration information not found"); 
    }
  }
}

std::string CtagCalibration::get_op(ctag::OperatingPoint t) const { 
  OPStrings::const_iterator oper = m_op_strings.find(t); 
  if (oper == m_op_strings.end()) { 
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
  bool ok_sf = m_cdi->retrieveCalibrationIndex(
    label, m_op_strings.at(op), m_jet_author, true, m_flav_op_sf_index[ind_key]); 
  bool ok_eff = m_cdi->retrieveCalibrationIndex(
    label, m_op_strings.at(op), m_jet_author, false, m_flav_op_eff_index[ind_key]); 
  if (!ok_eff || !ok_sf) { 
    std::string problem = "problem setting op " + m_op_strings.at(op) + " " + 
      get_label(flav); 
    throw std::runtime_error(problem); 
  }
}


unsigned CtagCalibration::get_index(
  const std::map<FOPIndex, unsigned>& vec, 
  ctag::Flavor flav, ctag::OperatingPoint op) const { 
  FOPIndex index = std::make_pair(flav, op); 
  std::map<FOPIndex, unsigned>::const_iterator pos = vec.find(index);
  if (pos == vec.end()) throw std::logic_error(
    "asked for unset operating point index in " __FILE__); 
  return pos->second; 
}



