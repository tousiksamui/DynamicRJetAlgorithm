//FJSTARTHEADER
//----------------------------------------------------------------------
//FJENDHEADER

#ifndef __DYNAMICRJETPLUGIN_HH__
#define __DYNAMICRJETPLUGIN_HH__

//#include "fastjet/internal/base.hh"
//#include "fastjet/config.h"

#include "fastjet/JetDefinition.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//namespace contrib{

// another forward declaration to reduce includes
class PseudoJet;

//----------------------------------------------------------------------

class ExtraInfo {
  public:
  ExtraInfo (double radin) {
    _r0   = radin;
  }
  double R0() const {return _r0;}

  private:
  double _r0;
};


class MoreInfo: public PseudoJet::UserInfoBase {
  public:
  MoreInfo(double mean_radin, double rms_radin, double wtin) { 
  _mean_rad = mean_radin;
  _rms_rad = rms_radin;
  _wt = wtin;
  }

  double mean_R() const {return _mean_rad;}
  double rms_R() const {return _rms_rad;}
  double wt() const {return _wt;}

  private:
  double _mean_rad, _rms_rad, _wt;
};


enum DynamicRJetAlgorithm {
  DRKT_algorithm = 0,
  KTlike=DRKT_algorithm,

  DRCA_algorithm = 1,
  CAlike=DRCA_algorithm,

  DRAK_algorithm = 2,
  AKlike=DRAK_algorithm
};

//----------------------------------------------------------------------

/// @ingroup plugins
/// \class DynamicRJetPlugin
/// Implementation of the DynamicRJet algorithm (plugin for fastjet v2.4 upwards)
//
class DynamicRJetPlugin : public JetDefinition::Plugin {
public:
  /// Main constructor for the DynamicRJet Plugin class.  
  ///
  DynamicRJetPlugin (double radius, DynamicRJetAlgorithm jet_algorithm_in=DRAK_algorithm,
		  RecombinationScheme jet_recombination_scheme=E_scheme){
    _radius  = radius;
    _radius2 = radius*radius;
    _jet_algorithm = jet_algorithm_in;
    _jet_recombiner = JetDefinition::DefaultRecombiner(jet_recombination_scheme);
  }

  /// copy constructor
  DynamicRJetPlugin (const DynamicRJetPlugin & plugin) {
    *this = plugin;
  }

  // the things that are required by base class
  virtual std::string description () const;
  virtual void run_clustering(ClusterSequence &) const;

  /// the plugin mechanism's standard way of accessing the jet radius
  /// here we return the R of the last alg in the list
  virtual double R() const {return _radius;}
  virtual DynamicRJetAlgorithm  algorithm() const {return _jet_algorithm;}

private:
  double _radius, _radius2;

  DynamicRJetAlgorithm _jet_algorithm;

  JetDefinition::DefaultRecombiner _jet_recombiner;

  static bool _first_time;

  template<class BJ> void actual_run(ClusterSequence &) const;

  /// print a banner for reference to the 3rd-party code
  void _print_banner(std::ostream *ostr) const;
};

class DRAK : public DynamicRJetPlugin {
  public:
    DRAK(double radius, RecombinationScheme recombination_scheme=E_scheme) : DynamicRJetPlugin(radius, DRAK_algorithm, recombination_scheme) {}
};

class DRCA : public DynamicRJetPlugin {
  public:
    DRCA(double radius, RecombinationScheme recombination_scheme=E_scheme) : DynamicRJetPlugin(radius, DRCA_algorithm, recombination_scheme) {}
};

class DRKT : public DynamicRJetPlugin {
  public:
    DRKT(double radius, RecombinationScheme recombination_scheme=E_scheme) : DynamicRJetPlugin(radius, DRKT_algorithm, recombination_scheme) {}
};
//} // namespace contrib

FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh

#endif // __TRACKJETPLUGIN_HH__

