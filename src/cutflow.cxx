#include <cstdio>
#include "SmartChain.hh"
#include "SusyBuffer.h"
#include "SUSYTools/SUSYObjDef.h"
#include "TLorentzVector.h"
#include "CutCounter.hh"

#include <cstdlib>
#include <vector> 
#include <algorithm>

// minimal class to keep track of particles
class IdLorentzVector : public TLorentzVector
{
public: 
  int index; 
  bool pass; 
};

// these functions check to see if the object passed the SUSYObjDef cuts
std::vector<IdLorentzVector> filter_pass(const std::vector<IdLorentzVector>&); 
std::vector<IdLorentzVector> filter_fail(const std::vector<IdLorentzVector>&); 
bool has_higher_pt(const TLorentzVector& v1, const TLorentzVector& v2); 

template<typename M, typename A>
A remove_overlaping(const M& mask, A altered, const float delta_r); 

int main (int narg, const char* argv[]) { 
  SmartChain* chain = new SmartChain("susy"); 
  for (int iii = 1; iii < narg; iii++) { 
    printf("file: %s, %i of %i\n", argv[iii], iii, narg - 1); 
    chain->add(argv[iii]); 
  }
  SusyBuffer buffer(chain); 

  // ------ initialize susytools here -----------------
  SUSYObjDef* def = new SUSYObjDef; 
  def->initialize(false, true); // not data, atlfast
  printf("initalized\n"); 

  CutCounter counter; 

  const unsigned n_entries = std::min(chain->GetEntries(),100000000LL); 
  printf("looping over %i entries\n", n_entries); 
  for (unsigned nnn = 0; nnn < n_entries; nnn++) {
    if (nnn % 1000 == 0) { 
      printf("%ik entries processed\n", nnn/1000); 
    }
    def->Reset(); 
    chain->GetEntry(nnn); 

    // ---- start filling objects here -----
    std::vector<IdLorentzVector> all_jets; 
    for (int jeti = 0; jeti < buffer.jet_n; jeti++) { 
      def->FillJet(
	jeti, 
	buffer.jet_pt                 ->at(jeti), 
	buffer.jet_eta                ->at(jeti), 
	buffer.jet_phi                ->at(jeti),
	buffer.jet_E                  ->at(jeti), 
	buffer.jet_constscale_eta        ->at(jeti), 
	buffer.jet_constscale_phi        ->at(jeti), 
	buffer.jet_constscale_E        ->at(jeti), 
	buffer.jet_constscale_m        ->at(jeti),
	buffer.jet_ActiveAreaPx->at(jeti), 
	buffer.jet_ActiveAreaPy->at(jeti), 
	buffer.jet_ActiveAreaPz->at(jeti), 
	buffer.jet_ActiveAreaE->at(jeti), 
	buffer.Eventshape_rhoKt4LC, 
	buffer.averageIntPerXing,
	buffer.vx_nTracks); 

      def->ApplyJetSystematics(
	jeti, 
	buffer.jet_constscale_eta        ->at(jeti), 
	buffer.jet_flavor_truth_label ->at(jeti), 
	buffer.averageIntPerXing,
	buffer.vx_nTracks, 
	SystErr::NONE); 

      bool good_jet = def->IsGoodJet(
	jeti, 
	buffer.jet_constscale_eta        ->at(jeti), 
	buffer.jet_emfrac             ->at(jeti), 
	buffer.jet_hecf               ->at(jeti),
	buffer.jet_LArQuality         ->at(jeti), 
	buffer.jet_HECQuality         ->at(jeti), 
	buffer.jet_AverageLArQF       ->at(jeti), 
	buffer.jet_Timing             ->at(jeti), 
	buffer.jet_sumPtTrk           ->at(jeti), 
	buffer.jet_fracSamplingMax    ->at(jeti),
	buffer.jet_SamplingMax        ->at(jeti), 
	buffer.jet_NegativeE          ->at(jeti), 
	buffer.RunNumber, 
	20e3, 	
	10,	
	JetID::VeryLooseBad);
      TLorentzVector jet_tlv = def->GetJetTLV(); 
      IdLorentzVector jet; 
      jet.SetPxPyPzE(jet_tlv.Px(), jet_tlv.Py(), jet_tlv.Pz(), jet_tlv.E()); 
      jet.index = jeti; 
      jet.pass = good_jet; 
      all_jets.push_back(jet); 
    } // end jet filling loop
    std::sort(all_jets.begin(),all_jets.end(),has_higher_pt); 

    std::vector<IdLorentzVector> all_electrons; 
    for (int eli = 0; eli < buffer.el_n; eli++) { 
      bool good_el = def->FillElectron(
	eli,
	buffer.el_eta                   ->at(eli), 
	buffer.el_phi                   ->at(eli), 
	buffer.el_cl_eta                ->at(eli),
	buffer.el_cl_phi                ->at(eli),
	buffer.el_cl_E                  ->at(eli),
	buffer.el_tracketa              ->at(eli),
	buffer.el_trackphi              ->at(eli),
	buffer.el_author                ->at(eli),
	buffer.el_mediumPP              ->at(eli),
	buffer.el_OQ                    ->at(eli),
	buffer.el_nPixHits              ->at(eli),
	buffer.el_nSCTHits              ->at(eli),
	buffer.el_MET_Egamma10NoTau_wet->at(eli).at(0), 
	10e3,			// et cut
	2.47);
      TLorentzVector el_tlv = def->GetElecTLV(); 
      IdLorentzVector electron; 
      electron.SetPxPyPzE(el_tlv.Px(), el_tlv.Py(), el_tlv.Pz(), el_tlv.E()); 
      electron.index = eli; 
      electron.pass = good_el; 
      all_electrons.push_back(electron); 
    } // end electron filling loop

    std::vector<IdLorentzVector> all_muons; 
    for (int mui = 0; mui < buffer.mu_staco_n; mui++) { 
      bool good_muon = def->FillMuon
	(mui,
	 buffer.mu_staco_pt                           ->at(mui),
	 buffer.mu_staco_eta                          ->at(mui),
	 buffer.mu_staco_phi                          ->at(mui),
	 buffer.mu_staco_me_qoverp_exPV               ->at(mui),
	 buffer.mu_staco_id_qoverp_exPV               ->at(mui),
	 buffer.mu_staco_me_theta_exPV                ->at(mui),
	 buffer.mu_staco_id_theta_exPV                ->at(mui),
	 buffer.mu_staco_id_theta                     ->at(mui),
	 buffer.mu_staco_charge                       ->at(mui), 
	 buffer.mu_staco_isCombinedMuon               ->at(mui),
	 buffer.mu_staco_isSegmentTaggedMuon          ->at(mui),
	 buffer.mu_staco_loose                        ->at(mui),
	 buffer.mu_staco_nPixHits                     ->at(mui),
	 buffer.mu_staco_nPixelDeadSensors            ->at(mui),
	 buffer.mu_staco_nPixHoles                    ->at(mui),
	 buffer.mu_staco_nSCTHits                     ->at(mui),
	 buffer.mu_staco_nSCTDeadSensors              ->at(mui),
	 buffer.mu_staco_nSCTHoles                    ->at(mui),
	 buffer.mu_staco_nTRTHits                     ->at(mui),
	 buffer.mu_staco_nTRTOutliers                 ->at(mui),
	 10e3, 
	 2.4);
      TLorentzVector muon_tlv = def->GetMuonTLV(mui); 
      IdLorentzVector muon; 
      muon.SetPxPyPzE(
	muon_tlv.Px(), muon_tlv.Py(), muon_tlv.Pz(), muon_tlv.E()); 
      muon.index = mui; 
      muon.pass = good_muon; 
      all_muons.push_back(muon); 
    } // end muon filling loop

    //  ----- object preselection ------
    std::vector<IdLorentzVector> preselected_el = filter_pass(all_electrons); 
    std::vector<IdLorentzVector> preselected_mu = filter_pass(all_muons); 
    std::vector<IdLorentzVector> preselected_jets; 
    for (std::vector<IdLorentzVector>::const_iterator jitr = all_jets.begin(); 
	 jitr != all_jets.end(); jitr++) { 
      bool is_good_pt = jitr->Pt() > 20e3; 
      bool is_good_eta = std::abs(jitr->Eta()) < 2.8; 
      if (is_good_eta && is_good_pt) { 
	preselected_jets.push_back(*jitr); 
      }
    }
    counter["preselected_el"] += preselected_el.size(); 
    counter["preselected_mu"] += preselected_mu.size(); 
    counter["preselected_jets"] += preselected_jets.size(); 

    // ---- overlap removal ------
    std::vector<IdLorentzVector> after_overlap_jets = remove_overlaping(
      preselected_el, preselected_jets, 0.2); 
    std::vector<IdLorentzVector> after_overlap_el = remove_overlaping(
      after_overlap_jets, preselected_el, 0.4); 
    std::vector<IdLorentzVector> after_overlap_mu = remove_overlaping(
      after_overlap_jets, preselected_mu, 0.4); 
    counter["after_overlap_el"]   += after_overlap_el.size(); 
    counter["after_overlap_mu"]   += after_overlap_mu.size(); 
    counter["after_overlap_jets"] += after_overlap_jets.size(); 

    // ---- veto object selection -----
    std::vector<IdLorentzVector> veto_jets = filter_fail(after_overlap_jets); 
    std::vector<IdLorentzVector> veto_electrons; 
    for (std::vector<IdLorentzVector>::const_iterator 
	   itr = after_overlap_el.begin(); 
	 itr != after_overlap_el.end(); itr++) { 
      float rel_isolation = (
	buffer.el_ptcone20->at(itr->index) / itr->Pt()); 
      bool isolated_el = rel_isolation < 0.1; 
      if (isolated_el) {
	veto_electrons.push_back(*itr); 
      }
    }
    std::vector<IdLorentzVector> veto_muons; 
    for (std::vector<IdLorentzVector>::const_iterator
	   itr = after_overlap_mu.begin(); 
	 itr != after_overlap_mu.end(); itr++) { 
      float isolation = buffer.mu_staco_ptcone20->at(itr->index); 
      bool isolated_mu = isolation < 1.8e3; 
      if (isolated_mu) { 
	veto_muons.push_back(*itr); 
      }
    }
    counter["veto_jets"] += veto_jets.size(); 
    counter["veto_electrons"] += veto_electrons.size(); 
    counter["veto_muons"] += veto_muons.size(); 

    // ---- signal object selection -----
    std::vector<IdLorentzVector> good_jets = filter_pass(after_overlap_jets);
    std::vector<IdLorentzVector> signal_jets; 
    for (std::vector<IdLorentzVector>::const_iterator
	   itr = good_jets.begin(); 
	 itr != good_jets.end(); itr++) { 
      bool signal_pt = itr->Pt() > 30e3; 
      bool tag_eta = std::abs(itr->Eta()) < 2.5; 
      float jet_jvf = buffer.jet_jvtxf->at(itr->index); 
      bool ok_jvf = (jet_jvf > 0.5) || (itr->Pt() > 50e3); 
      if (signal_pt && tag_eta && ok_jvf) { 
	signal_jets.push_back(*itr); 
      }
    }
    std::vector<IdLorentzVector> control_electrons; 
    for (std::vector<IdLorentzVector>::const_iterator
	   itr = veto_electrons.begin(); 
	 itr != veto_electrons.end(); itr++) { 
      bool control_pt = itr->Pt() > 20e3; 
      if (control_pt) { 
	control_electrons.push_back(*itr); 
      }
    }
    std::vector<IdLorentzVector> control_muons; 
    for (std::vector<IdLorentzVector>::const_iterator
	   itr = veto_muons.begin(); 
	 itr != veto_muons.end(); itr++) { 
      bool control_pt = itr->Pt() > 20e3; 
      if (control_pt) { 
	control_muons.push_back(*itr); 
      }
    }
    counter["good_jets"] += good_jets.size(); 
    counter["signal_jets"] += signal_jets.size(); 
    counter["control_electrons"] += control_electrons.size(); 
    counter["control_muons"] += control_muons.size(); 
  } // end of event loop


  // ------ dump results ------
  typedef std::vector<std::pair<std::string, int> > OrdCuts; 
  OrdCuts ordered_cuts = counter.get_ordered_cuts(); 
  for (OrdCuts::const_iterator itr = ordered_cuts.begin(); 
       itr != ordered_cuts.end(); itr++) { 
    printf("%s: %i\n", itr->first.c_str(), itr->second); 
  }

}


std::vector<IdLorentzVector> filter_pass(
  const std::vector<IdLorentzVector>& in) { 
  std::vector<IdLorentzVector> out; 
  for (std::vector<IdLorentzVector>::const_iterator itr = in.begin(); 
       itr != in.end(); itr++) { 
    if (itr->pass) { 
      out.push_back(*itr); 
    }
  }
  return out; 
}

std::vector<IdLorentzVector> filter_fail(
  const std::vector<IdLorentzVector>& in) { 
  std::vector<IdLorentzVector> out; 
  for (std::vector<IdLorentzVector>::const_iterator itr = in.begin(); 
       itr != in.end(); itr++) { 
    if (!itr->pass) { 
      out.push_back(*itr); 
    }
  }
  return out; 
}

template<typename M, typename A>
A remove_overlaping(const M& mask, A altered, const float delta_r) { 
  for (typename M::const_iterator mitr = mask.begin(); 
       mitr != mask.end(); mitr++) { 
    A new_container; 
    for (typename M::const_iterator vic = altered.begin(); 
	 vic != altered.end(); vic++) { 
      assert(mitr->Pt() > 0); 
      double delr = mitr->DeltaR(*vic); 
      if (delr > delta_r) { 
	new_container.push_back(*vic); 
      }
    }
    altered = new_container; 
  }
  return altered; 
} 

bool has_higher_pt(const TLorentzVector& v1, const TLorentzVector& v2) { 
  return v1.Pt() > v2.Pt(); 
}
