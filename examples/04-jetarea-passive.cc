//---------------------------------------------------------------------
// file: 04-jetarea-passive.cc
// run:  ./run04.sh
// run it with sample.dat
//---------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/DynamicRJetPlugin.hh"

// Header and Namespace.
using namespace std;
using namespace fastjet;

int main() {
  // input file
  ifstream infile;
  string filepath = "sample.dat";

  infile.open(filepath);

  // No. of events.
  int nEvent = 10 ;

  // FastJet parameters
  double radius = 0.5;
  double pTjetMin = 10.0;
  double etaMax = 5.0;

  vector<PseudoJet> particles;

  // some useful variables
  double px, py, pz, e;

  // Loop over events 
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    particles.resize(0);

    // Loop over event record
    while(true) {
      infile >> px >> py >> pz >> e; 
      if (e < 0.0) break;
      // construct the particle PseudoJet
      PseudoJet particle(px, py, pz, e);

      particles.push_back(particle);
    }

    // continue to next event if no final state particle.
    if (int(particles.size())==0) continue; 

    AreaDefinition area_def(passive_area);

    // DR-AK algorithm
    DRAK DRAKJP(radius);
    JetDefinition jet_def_DRAK(&DRAKJP);

    ClusterSequenceArea cs_DRAK(particles, jet_def_DRAK, area_def);

    vector<PseudoJet> DRAKjets = sorted_by_pt(cs_DRAK.inclusive_jets(pTjetMin));

    cout << "############################################################" << endl;
    cout << "Event # " << iEvent+1 << "." << endl;
    cout << "------------------------------------------------------------" << endl;

    cout << "Output of " << jet_def_DRAK.description() << "." << endl;
    cout << setw(5) << right << "jet #" << setw(10) << right << "rapidity" << setw(10) << right << "phi" << setw(15) << right << "pt (GeV)" << setw(22) << right << "jet area (Passive)" << endl; 
    // Rd can't be obtained CluasterSequenceArea in case of passive_area.
    // Please use ClusterSequence to get Rd (see 01-basic.cc). 

    for (int iD = 0; iD < DRAKjets.size(); iD++) {

      double jet_area = DRAKjets[iD].area();

      cout << setw(5) << right << iD+1;
      cout << setw(10) << right << fixed << setprecision(2) << DRAKjets[iD].rap();
      cout << setw(10) << right << fixed << setprecision(2) << DRAKjets[iD].phi();
      cout << setw(15) << right << fixed << setprecision(2) << DRAKjets[iD].pt();
      cout << setw(16) << right << fixed << setprecision(2) << jet_area << endl;
    }
    cout << "------------------------------------------------------------" << endl;
    cout << endl;

    // DR-CA algorithm
    DRCA DRCAJP(radius);
    JetDefinition jet_def_DRCA(&DRCAJP);

    ClusterSequenceArea cs_DRCA(particles, jet_def_DRCA, area_def);
    vector<PseudoJet> DRCAjets = sorted_by_pt(cs_DRCA.inclusive_jets(pTjetMin));

    cout << "Output of " << jet_def_DRCA.description() << "." << endl;
    cout << setw(5) << right << "jet #" << setw(10) << right << "rapidity" << setw(10) << right << "phi" << setw(15) << right << "pt (GeV)" << setw(22) << right << "jet area (Passive)" << endl; 
    // Rd can't be obtained CluasterSequenceArea in case of passive_area.
    // Please use ClusterSequence to get Rd (see 01-basic.cc). 

    for (int iD = 0; iD < DRCAjets.size(); iD++) {

     double jet_area = DRCAjets[iD].area();

      cout << setw(5) << right << iD+1;
      cout << setw(10) << right << fixed << setprecision(2) << DRCAjets[iD].rap();
      cout << setw(10) << right << fixed << setprecision(2) << DRCAjets[iD].phi();
      cout << setw(15) << right << fixed << setprecision(2) << DRCAjets[iD].pt();
      cout << setw(16) << right << fixed << setprecision(2) << jet_area << endl;
    }
    cout << "------------------------------------------------------------" << endl;
    cout << endl;

    // DR-KT algorithm
    DRKT DRKTJP(radius);
    JetDefinition jet_def_DRKT(&DRKTJP);

    ClusterSequenceArea cs_DRKT(particles, jet_def_DRKT, area_def);
    vector<PseudoJet> DRKTjets = sorted_by_pt(cs_DRKT.inclusive_jets(pTjetMin));

    cout << "Output of " << jet_def_DRKT.description() << "." << endl;
    cout << setw(5) << right << "jet #" << setw(10) << right << "rapidity" << setw(10) << right << "phi" << setw(15) << right << "pt (GeV)" << setw(22) << right << "jet area (Passive)" << endl; 
    // Rd can't be obtained CluasterSequenceArea in case of passive_area.
    // Please use ClusterSequence to get Rd (see 01-basic.cc). 

    for (int iD = 0; iD < DRKTjets.size(); iD++) {

      double jet_area = DRKTjets[iD].area();

      cout << setw(5) << right << iD+1;
      cout << setw(10) << right << fixed << setprecision(2) << DRKTjets[iD].rap();
      cout << setw(10) << right << fixed << setprecision(2) << DRKTjets[iD].phi();
      cout << setw(15) << right << fixed << setprecision(2) << DRKTjets[iD].pt();
      cout << setw(16) << right << fixed << setprecision(2) << jet_area << endl;
    }
    cout << "------------------------------------------------------------" << endl;
    cout << endl;

  }

  return 0;
}
