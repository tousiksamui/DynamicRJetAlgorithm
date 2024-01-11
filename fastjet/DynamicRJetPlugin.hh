//FJSTARTHEADER
//----------------------------------------------------------------------
//FJENDHEADER

#ifndef __DELTARJETPLUGIN_HH__
#define __DELTARJETPLUGIN_HH__

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
  ExtraInfo (int pin, double radin) {
    p        = pin;
    radius   = radin;
  }
  int p;
  double radius;
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
  double _mean_rad, _rms_rad, _wt;
};

//----------------------------------------------------------------------


/// @ingroup plugins
/// \class DeltaRJetPlugin
/// Implementation of the DeltaRJet algorithm (plugin for fastjet v2.4 upwards)
//
class DeltaRJetPlugin : public JetDefinition::Plugin {
public:
  /// Main constructor for the DeltaRJet Plugin class.  
  ///
  DeltaRJetPlugin (int p, double radius, 
		  RecombinationScheme jet_recombination_scheme=E_scheme){
    _p       = p;
    _radius  = radius;
    _radius2 = radius*radius;
    _jet_recombiner = JetDefinition::DefaultRecombiner(jet_recombination_scheme);
  }

  /// copy constructor
  DeltaRJetPlugin (const DeltaRJetPlugin & plugin) {
    *this = plugin;
  }

  // the things that are required by base class
  virtual std::string description () const;
  virtual void run_clustering(ClusterSequence &) const;

  /// the plugin mechanism's standard way of accessing the jet radius
  /// here we return the R of the last alg in the list
  virtual double R() const {return _radius;}
  virtual int    algorithm() const {return _p;}

private:
  int _p;
  double _radius, _radius2;

  JetDefinition::DefaultRecombiner _jet_recombiner;

  static bool _first_time;

  //template<class N> void _actual_run_clustering(ClusterSequence &, ExtraInfo *) const;

  /// print a banner for reference to the 3rd-party code
  void _print_banner(std::ostream *ostr) const;
};
//} // namespace contrib

FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh

#endif // __TRACKJETPLUGIN_HH__

