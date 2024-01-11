//FJSTARTHEADER
//FJENDHEADER
#include "fastjet/DynamicRJetPlugin.hh"

// fastjet stuff
#include "fastjet/ClusterSequence.hh"
#include "fastjet/NNH.hh"
#include "fastjet/PseudoJet.hh"
//#include <iostream>

// other stuff
#include <list>
#include <memory>
#include <cmath>
#include <vector>
#include <sstream>

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
//namespace contrib{

class DeltaRJetBJ {
public:
  void init(const PseudoJet & jet, const ExtraInfo * extraInfo){
    _p = extraInfo->p;
    _r = extraInfo->radius;
    _pti = jet.pt();
    _pti2 = _pti*_pti;
    _jeti = jet;
    _avg_R = jet.user_info<MoreInfo>().mean_R();
    _sd_R  = jet.user_info<MoreInfo>().rms_R() - _avg_R * _avg_R;
    _sd_R  = sqrt(_sd_R);
  }

  double distance(const DeltaRJetBJ * jet) const {
    double distance2 = jet->_jeti.plain_distance(_jeti);
    double ptmax = jet->_pti2;
    if (ptmax < _pti2) ptmax = _pti2;
    double dij = distance2/ptmax;

    return dij;
  }

  double beam_distance() const {
    
    double num = (_r + _sd_R);
    double diB = num * num/_pti2; 
    return diB; 
  }

private:
  int _p ;
  double _r, _pti, _pti2 ;
  PseudoJet _jeti;
  double _avg_R, _sd_R;

};

//------------------------------------------------------------------
// implementation of the DeltaRJet plugin
//------------------------------------------------------------------

bool DeltaRJetPlugin::_first_time = true;

string DeltaRJetPlugin::description () const {
  ostringstream desc;
  desc << "DeltaRJet algorithm with R = " << R();
  return desc.str();
}

void DeltaRJetPlugin::run_clustering(ClusterSequence & clust_seq) const {

  int njets = clust_seq.jets().size();
  int insize = njets;
  //std::cout << njets << "run_clustering" << std::endl;
  int algo_pow  = algorithm();
  double radius = R();
  PseudoJet newjet, jet_i, jet_j;
  double mean_DR, rms_DR, wt;
  double mean_DR1, rms_DR1, wt1;
  double jet_i_pt, jet_j_pt;
  double DRij;
  double pti, ptj;

  //MoreInfo minfo(0.4, 0.0, 1.0);

  //for (int ii =0; ii<njets; ii++){
    //clust_seq.jets()[ii].set_user_info(new MoreInfo(0.4, 0.0, 1.0));
    //PseudoJet particle = clust_seq.jets()[ii];
    //particle.set_user_info(new MoreInfo(0.4,0.0,1.0));
    //std::cout << "test test test" << std::endl;
    //std::cout << particle.user_info<MoreInfo>().mean_R() << std::endl;
  //}

  ExtraInfo extraInfo(algo_pow, radius);
  //DeltaRJetPlugin::_actual_run_clustering<NNH<DeltaRJetBJ, ExtraInfo> >(clust_seq, &extraInfo);

  NNH<DeltaRJetBJ, ExtraInfo> nn(clust_seq.jets(), &extraInfo);

  while (njets > 0) {
    int i, j, k;
    double dij = nn.dij_min(i, j);
    if (j >= 0) {
      jet_i  = clust_seq.jets()[i];
      jet_j  = clust_seq.jets()[j];

      jet_i_pt = jet_i.pt();
      jet_j_pt = jet_j.pt();
      
      //newjet_constituents.resize(0);
      //vector<PseudoJet> newjet_constituents = clust_seq.jets()[i].constituents() 
      //newjet_constituents.insert(newjet_constituents.end(), clust_seq.jets()[j].constituents());
      if (i < insize && j < insize) {
        mean_DR = jet_i.delta_R(jet_j);
	rms_DR = mean_DR*mean_DR;
	wt = jet_i_pt*jet_j_pt;
      }
      else if ( i < insize && j >= insize) {
	mean_DR = 0;
	rms_DR  = 0;
	wt      = 0;
	for (unsigned int jt=0; jt < jet_j.constituents().size(); jt++){
          DRij = jet_i.delta_R(jet_j.constituents()[jt]);
	  ptj  = jet_j.constituents()[jt].pt();
          mean_DR += DRij*ptj;
	  rms_DR  += DRij*DRij*ptj;
	  wt      += ptj;
        }
	mean_DR = mean_DR*jet_i_pt;
	rms_DR  = rms_DR *jet_i_pt;
	wt      = wt     *jet_i_pt;

	mean_DR = (mean_DR + jet_j.user_info<MoreInfo>().mean_R()*jet_j.user_info<MoreInfo>().wt())/(wt+jet_j.user_info<MoreInfo>().wt());
	rms_DR  = (rms_DR + jet_j.user_info<MoreInfo>().rms_R()*jet_j.user_info<MoreInfo>().wt())/(wt+jet_j.user_info<MoreInfo>().wt());
	wt      = wt + jet_j.user_info<MoreInfo>().wt();
      }
      else if ( i >= insize && j < insize) {
	mean_DR = 0;
	rms_DR  = 0;
	wt      = 0;
	for (unsigned int it=0; it < jet_i.constituents().size(); it++){
          DRij = jet_j.delta_R(jet_i.constituents()[it]);
	  pti  = jet_i.constituents()[it].pt();
          mean_DR += DRij*pti;
	  rms_DR  += DRij*DRij*pti;
	  wt      += pti;
        }
	mean_DR = mean_DR*jet_j_pt;
	rms_DR  = rms_DR *jet_j_pt;
	wt      = wt     *jet_j_pt;

	mean_DR = (mean_DR + jet_i.user_info<MoreInfo>().mean_R()*jet_i.user_info<MoreInfo>().wt())/(wt+jet_i.user_info<MoreInfo>().wt());
	rms_DR  = (rms_DR + jet_i.user_info<MoreInfo>().rms_R()*jet_i.user_info<MoreInfo>().wt())/(wt+jet_i.user_info<MoreInfo>().wt());
	wt      = wt + jet_i.user_info<MoreInfo>().wt();
      }
      else  {
	mean_DR = 0;
	rms_DR  = 0;
	wt      = 0;
	for (unsigned int it=0; it < jet_i.constituents().size(); it++){
	  pti  = jet_i.constituents()[it].pt();
	  mean_DR1 = 0;
	  rms_DR1  = 0;
	  wt1      = 0;
	  for (unsigned int jt=0; jt < jet_j.constituents().size(); jt++) {
            DRij = jet_i.constituents()[it].delta_R(jet_j.constituents()[jt]);
	    ptj  = jet_j.constituents()[jt].pt();
            mean_DR1 += DRij*ptj;
	    rms_DR1  += DRij*DRij*ptj;
	    wt1      += ptj;
	  }
	  mean_DR = mean_DR + mean_DR1*pti;
	  rms_DR  = rms_DR  + rms_DR1 *pti;
	  wt      = wt      + wt1     *pti;
        }

	mean_DR = (mean_DR + jet_i.user_info<MoreInfo>().mean_R()*jet_i.user_info<MoreInfo>().wt() + jet_j.user_info<MoreInfo>().mean_R()*jet_j.user_info<MoreInfo>().wt())/(wt + jet_i.user_info<MoreInfo>().wt() + jet_j.user_info<MoreInfo>().wt());
	rms_DR  = (rms_DR + jet_i.user_info<MoreInfo>().rms_R()*jet_i.user_info<MoreInfo>().wt() + jet_j.user_info<MoreInfo>().rms_R()*jet_j.user_info<MoreInfo>().wt())/(wt + jet_i.user_info<MoreInfo>().wt() + jet_j.user_info<MoreInfo>().wt());
	wt      = wt + jet_i.user_info<MoreInfo>().wt() + jet_j.user_info<MoreInfo>().wt();
      }
      newjet = jet_i + jet_j;
      newjet.set_user_info(new MoreInfo(mean_DR, rms_DR, wt));
      clust_seq.plugin_record_ij_recombination(i, j, dij, newjet, k);
      nn.merge_jets(i, j, clust_seq.jets()[k], k);
    }
    else {
      clust_seq.plugin_record_iB_recombination(i, dij);
      nn.remove_jet(i);
    }
    njets--;
  }
    
}

// print a banner for reference to the 3rd-party code
void DeltaRJetPlugin::_print_banner(ostream *ostr) const{
  if (! _first_time) return;
  _first_time=false;

  // make sure the user has not set the banner stream to NULL
  if (!ostr) return;  

  (*ostr) << "#-------------------------------------------------------------------------" << endl;
  (*ostr) << "# You are running the TrackJet plugin for FastJet. It is based on         " << endl;
  (*ostr) << "# the implementation by Andy Buckley and Manuel Bahr that is to be        " << endl;
  (*ostr) << "# found in Rivet 1.1.2. See http://www.hepforge.org/downloads/rivet.      " << endl;
  (*ostr) << "#-------------------------------------------------------------------------" << endl;

  // make sure we really have the output done.
  ostr->flush();
}

//} // namespace contrib
FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
