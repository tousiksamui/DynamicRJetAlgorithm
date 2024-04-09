// fastjet stuff
#include "fastjet/ClusterSequence.hh"
#include "fastjet/NNH.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/DynamicRJetPlugin.hh"

// other stuff
#include <list>
#include <memory>
#include <cmath>
#include <vector>
#include <sstream>

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
//namespace contrib{

class DRAKBriefJet {
public:
  void init(const PseudoJet & jet, const ExtraInfo * extraInfo){
    //_p = extraInfo->p;
    //_pti = jet.pt();
    _jeti = jet;
    _pti2 = jet.pt2();
    _avg_R = jet.user_info<MoreInfo>().mean_R();
    _sd_R  = jet.user_info<MoreInfo>().rms_R() - _avg_R * _avg_R;
    _sd_R  = (_sd_R > 0.0) ? sqrt(_sd_R) : -sqrt(-_sd_R);
    _r = extraInfo->R0();
  }

  double distance(const DRAKBriefJet * jet) const {
    double distance2 = jet->_jeti.plain_distance(_jeti);
    double pt2max = jet->_pti2;
    if (pt2max < _pti2) pt2max = _pti2;
    return distance2/pt2max;
    //return dij;
  }

  double beam_distance() const {
    double num = (_r + _sd_R);
    double diB = num * num/_pti2; 
    return diB; 
  }

private:
  PseudoJet _jeti;
  double _pti2, _avg_R, _sd_R, _r;

};

class DRCABriefJet {
public:
  void init(const PseudoJet & jet, const ExtraInfo * extraInfo){
    _jeti = jet;
    _avg_R = jet.user_info<MoreInfo>().mean_R();
    _sd_R  = jet.user_info<MoreInfo>().rms_R() - _avg_R * _avg_R;
    _sd_R  = (_sd_R > 0.0) ? sqrt(_sd_R) : -sqrt(-_sd_R);
    _r = extraInfo->R0();
  }

  double distance(const DRCABriefJet * jet) const {
    double distance2 = jet->_jeti.plain_distance(_jeti);
    return distance2;
  }

  double beam_distance() const {
    double num = (_r + _sd_R);
    double diB = num * num; 
    return diB; 
  }

private:
  PseudoJet _jeti;
  double _avg_R, _sd_R, _r;

};

class DRKTBriefJet {
public:
  void init(const PseudoJet & jet, const ExtraInfo * extraInfo){
    _jeti = jet;
    _pti2 = jet.pt2();
    _avg_R = jet.user_info<MoreInfo>().mean_R();
    _sd_R  = jet.user_info<MoreInfo>().rms_R() - _avg_R * _avg_R;
    _sd_R  = (_sd_R > 0.0) ? sqrt(_sd_R) : -sqrt(-_sd_R);
    _r = extraInfo->R0();
  }

  double distance(const DRKTBriefJet * jet) const {
    double distance2 = jet->_jeti.plain_distance(_jeti);
    double pt2min = jet->_pti2;
    if (pt2min > _pti2) pt2min = _pti2;
    return distance2*pt2min;
  }

  double beam_distance() const {
    double num = (_r + _sd_R);
    double diB = num * num * _pti2; 
    return diB; 
  }

private:
  PseudoJet _jeti;
  double _pti2, _avg_R, _sd_R, _r;

};


//------------------------------------------------------------------
// implementation of the DynamicRJet plugin
//------------------------------------------------------------------

bool DynamicRJetPlugin::_first_time = true;

string DynamicRJetPlugin::description () const {
  ostringstream desc;
  switch(_jet_algorithm) {
    case DRKT_algorithm:
      desc << "Dynamic Radius (KT) Jet algorithm with R0 = " << R();
      break;
    case DRCA_algorithm:
      desc << "Dynamic Radius (CA) Jet algorithm with R0 = " << R();
      break;
    case DRAK_algorithm:
      desc << "Dynamic Radius (AK) Jet algorithm with R0 = " << R();
      break;
    default:
      throw Error("unrecognized jet_algorithm");
  }

  return desc.str();
}



void DynamicRJetPlugin::run_clustering(ClusterSequence & clust_seq) const {

  switch(_jet_algorithm) {
    case DRKT_algorithm:
      actual_run<DRKTBriefJet>(clust_seq);
      break;
    case DRCA_algorithm:
      actual_run<DRCABriefJet>(clust_seq);
      break;
    case DRAK_algorithm:
      actual_run<DRAKBriefJet>(clust_seq);
      break;
    default:
      throw Error("unrecognized jet_algorithm");
  }
}

template<class BJ> void DynamicRJetPlugin::actual_run(ClusterSequence & clust_seq) const {

  int njets = clust_seq.jets().size();
  int insize = njets;
  PseudoJet newjet, jet_i, jet_j;
  double mean_DR, rms_DR, wt;
  double mean_DR1, rms_DR1, wt1;
  double jet_i_pt, jet_j_pt;
  double DRij, DR2ij;
  double pti, ptj;

  double radius = R();
  ExtraInfo extraInfo(radius);

  NNH<BJ, ExtraInfo> nn(clust_seq.jets(), &extraInfo);
  //NNBase<ExtraInfo> *nn = nullptr ;
  //NNHInfo<ExtraInfo> *nn = nullptr ;
  //NNH<BJ, ExtraInfo> *nn = nullptr;


  while (njets > 0) {
    int i, j, k;
    double dij = nn.dij_min(i, j);
    if (j >= 0) {
      jet_i  = clust_seq.jets()[i];
      jet_j  = clust_seq.jets()[j];
      
      if (i < insize && j < insize) {
        jet_i_pt = jet_i.pt();
        jet_j_pt = jet_j.pt();
        rms_DR = jet_i.plain_distance(jet_j);
        mean_DR = sqrt(rms_DR);
        wt = jet_i_pt*jet_j_pt;
      }
      else if ( i < insize && j >= insize) {
        jet_i_pt = jet_i.pt();
        mean_DR = 0;
        rms_DR  = 0;
        wt      = 0;
        vector<PseudoJet> jet_j_constituents = jet_j.constituents();
        for (unsigned int jt=0; jt < jet_j_constituents.size(); jt++){
          DR2ij = jet_i.plain_distance(jet_j_constituents[jt]);
          ptj  = jet_j_constituents[jt].pt();
          rms_DR  += DR2ij*ptj;
          mean_DR += sqrt(DR2ij)*ptj;
          wt      += ptj;
        }
        mean_DR = mean_DR*jet_i_pt;
        rms_DR  = rms_DR *jet_i_pt;
        wt      = wt     *jet_i_pt;

        wt      = wt + jet_j.user_info<MoreInfo>().wt();
        mean_DR = (mean_DR + jet_j.user_info<MoreInfo>().mean_R()*jet_j.user_info<MoreInfo>().wt())/wt;
        rms_DR  = (rms_DR + jet_j.user_info<MoreInfo>().rms_R()*jet_j.user_info<MoreInfo>().wt())/wt;
      }
      else if ( i >= insize && j < insize) {
        jet_j_pt = jet_j.pt();
        mean_DR = 0;
        rms_DR  = 0;
        wt      = 0;
        vector<PseudoJet> jet_i_constituents = jet_i.constituents();
        for (unsigned int it=0; it < jet_i_constituents.size(); it++){
          DR2ij = jet_j.plain_distance(jet_i_constituents[it]);
          pti  = jet_i_constituents[it].pt();
          rms_DR  += DR2ij*pti;
          mean_DR += sqrt(DR2ij)*pti;
          wt      += pti;
        }
        mean_DR = mean_DR*jet_j_pt;
        rms_DR  = rms_DR *jet_j_pt;
        wt      = wt     *jet_j_pt;

        wt      = wt + jet_i.user_info<MoreInfo>().wt();
        mean_DR = (mean_DR + jet_i.user_info<MoreInfo>().mean_R()*jet_i.user_info<MoreInfo>().wt())/wt;
        rms_DR  = (rms_DR + jet_i.user_info<MoreInfo>().rms_R()*jet_i.user_info<MoreInfo>().wt())/wt;
      }
      else  {
        mean_DR = 0;
        rms_DR  = 0;
        wt      = 0;
        vector<PseudoJet> jet_i_constituents = jet_i.constituents();
        vector<PseudoJet> jet_j_constituents = jet_j.constituents();
        vector<double> ptjs, ptis;
        ptjs.clear();
        ptis.clear();

        double wti = 0;
        double wtj = 0;

        for (unsigned int it=0; it < jet_i_constituents.size(); it++){
          pti  = jet_i_constituents[it].pt();
          wti = wti + pti;
          ptis.push_back(pti);
        }

        for (unsigned int jt=0; jt < jet_j_constituents.size(); jt++){
          ptj  = jet_j_constituents[jt].pt();
          wtj = wtj + ptj;
          ptjs.push_back(ptj);
        }

        for (unsigned int it=0; it < jet_i_constituents.size(); it++){
          mean_DR1 = 0;
          rms_DR1  = 0;
          for (unsigned int jt=0; jt < jet_j_constituents.size(); jt++) {
            DR2ij = jet_i_constituents[it].plain_distance(jet_j_constituents[jt]);
            rms_DR1  += DR2ij*ptjs[jt];
            mean_DR1 += sqrt(DR2ij)*ptjs[jt];
          }
          mean_DR = mean_DR + mean_DR1*ptis[it];
          rms_DR  = rms_DR  + rms_DR1 *ptis[it];
        }

        wt = wti*wtj + jet_i.user_info<MoreInfo>().wt() + jet_j.user_info<MoreInfo>().wt();
        mean_DR = (mean_DR + jet_i.user_info<MoreInfo>().mean_R()*jet_i.user_info<MoreInfo>().wt() + jet_j.user_info<MoreInfo>().mean_R()*jet_j.user_info<MoreInfo>().wt())/wt;
        rms_DR  = (rms_DR + jet_i.user_info<MoreInfo>().rms_R()*jet_i.user_info<MoreInfo>().wt() + jet_j.user_info<MoreInfo>().rms_R()*jet_j.user_info<MoreInfo>().wt())/wt;
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
void DynamicRJetPlugin::_print_banner(ostream *ostr) const{
  if (! _first_time) return;
  _first_time=false;

  // make sure the user has not set the banner stream to NULL
  if (!ostr) return;  

  (*ostr) << "#-------------------------------------------------------------------------" << endl;
  (*ostr) << "# You are running the Dyanamic Radius Jet Clustering Algorithm. It is     " << endl;
  (*ostr) << "# impleted within the framework of FastJet3. Please cite the FastJet3     " << endl;
  (*ostr) << "# references. It is based on arXiv:2301.13074. This implementation is     " << endl;
  (*ostr) << "# done by Tousik Samui.                                                   " << endl;
  (*ostr) << "#-------------------------------------------------------------------------" << endl;

  // make sure we really have the output done.
  ostr->flush();
}

//} // namespace contrib
FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
