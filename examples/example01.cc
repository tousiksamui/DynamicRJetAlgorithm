//---------------------------------------------------------------------
// file: example01.cc
// run it with sample.dat
//---------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/DynamicRJetPlugin.hh"

// Header and Namespace.
using namespace std;
using namespace fastjet;

int main() {
  // input file
  ifstream infile;
  string filepath = "sample.dat";

  infile.open(filepath);

  // No. of events to generate.
  int nEvent = 10 ;

  // FastJet parameters
  double radius = 0.5;
  double pTjetMin = 10.0;
  double etaMax = 5.0;

  vector<PseudoJet> particles;

  // some useful variables
  double part_pt2;
  double Rd;
  double px, py, pz, e;

  // Loop over events 
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    cout << "############################################################" << endl;
    cout << "Event # " << iEvent+1 << "." << endl;
    cout << "------------------------------------------------------------" << endl;
    particles.resize(0);

    // Loop over event record
    while(true) {
      infile >> px >> py >> pz >> e; 
      if (e < 0.0) break;
      // construct the particle PseudoJet
      PseudoJet particle(px, py, pz, e);
      part_pt2 = px*px + py*py;
      particle.set_user_info(new MoreInfo(0.0, 0.0, part_pt2));
      particles.push_back(particle);
    }

    // continue to next event if no final state particle.
    if (int(particles.size())==0) continue; 

    // DR-AK algorithm
    DRAK DRAKJP(radius);
    JetDefinition jet_def_DRAK(&DRAKJP);
    ClusterSequence cs_DRAK(particles, jet_def_DRAK);
    vector<PseudoJet> DRAKjets = sorted_by_pt(cs_DRAK.inclusive_jets(pTjetMin));

    cout << "Output of " << jet_def_DRAK.description() << "." << endl;
    //printf("%5s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt");
    cout << setw(5) << right << "jet #" << setw(10) << right << "rapidity" << setw(10) << right << "phi" << setw(15) << right << "pt (GeV)" << setw(10) << right << "Rd" << endl; 

    for (int iD = 0; iD < DRAKjets.size(); iD++) {
      Rd = DRAKjets[iD].user_info<MoreInfo>().mean_R();
      Rd = DRAKjets[iD].user_info<MoreInfo>().rms_R() - Rd*Rd;
      Rd = (Rd > 0.0) ? radius+sqrt(Rd) : radius-sqrt(-Rd);

      cout << setw(5) << right << iD+1;
      cout << setw(10) << right << fixed << setprecision(2) << DRAKjets[iD].rap();
      cout << setw(10) << right << fixed << setprecision(2) << DRAKjets[iD].phi();
      cout << setw(15) << right << fixed << setprecision(2) << DRAKjets[iD].pt();
      cout << setw(10) << right << fixed << setprecision(2) << Rd << endl;
    }
    cout << "------------------------------------------------------------" << endl;
    cout << endl;

    // DR-CA algorithm
    DRCA DRCAJP(radius);
    JetDefinition jet_def_DRCA(&DRCAJP);
    ClusterSequence cs_DRCA(particles, jet_def_DRCA);
    vector<PseudoJet> DRCAjets = sorted_by_pt(cs_DRCA.inclusive_jets(pTjetMin));

    cout << "Output of " << jet_def_DRCA.description() << "." << endl;
    //printf("%5s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt");
    cout << setw(5) << right << "jet #" << setw(10) << right << "rapidity" << setw(10) << right << "phi" << setw(15) << right << "pt (GeV)" << setw(10) << right << "Rd" << endl; 

    for (int iD = 0; iD < DRCAjets.size(); iD++) {
      Rd = DRCAjets[iD].user_info<MoreInfo>().mean_R();
      Rd = DRCAjets[iD].user_info<MoreInfo>().rms_R() - Rd*Rd;
      Rd = (Rd > 0.0) ? radius+sqrt(Rd) : radius-sqrt(-Rd);

      cout << setw(5) << right << iD+1;
      cout << setw(10) << right << fixed << setprecision(2) << DRCAjets[iD].rap();
      cout << setw(10) << right << fixed << setprecision(2) << DRCAjets[iD].phi();
      cout << setw(15) << right << fixed << setprecision(2) << DRCAjets[iD].pt();
      cout << setw(10) << right << fixed << setprecision(2) << Rd << endl;
    }
    cout << "------------------------------------------------------------" << endl;
    cout << endl;

    // DR-KT algorithm
    DRKT DRKTJP(radius);
    JetDefinition jet_def_DRKT(&DRKTJP);
    ClusterSequence cs_DRKT(particles, jet_def_DRKT);
    vector<PseudoJet> DRKTjets = sorted_by_pt(cs_DRKT.inclusive_jets(pTjetMin));

    cout << "Output of " << jet_def_DRKT.description() << "." << endl;
    //printf("%5s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt");
    cout << setw(5) << right << "jet #" << setw(10) << right << "rapidity" << setw(10) << right << "phi" << setw(15) << right << "pt (GeV)" << setw(10) << right << "Rd" << endl; 

    for (int iD = 0; iD < DRKTjets.size(); iD++) {
      Rd = DRKTjets[iD].user_info<MoreInfo>().mean_R();
      Rd = DRKTjets[iD].user_info<MoreInfo>().rms_R() - Rd*Rd;
      Rd = (Rd > 0.0) ? radius+sqrt(Rd) : radius-sqrt(-Rd);

      cout << setw(5) << right << iD+1;
      cout << setw(10) << right << fixed << setprecision(2) << DRKTjets[iD].rap();
      cout << setw(10) << right << fixed << setprecision(2) << DRKTjets[iD].phi();
      cout << setw(15) << right << fixed << setprecision(2) << DRKTjets[iD].pt();
      cout << setw(10) << right << fixed << setprecision(2) << Rd << endl;
    }
    cout << "------------------------------------------------------------" << endl;
    cout << endl;

  }

  return 0;
}
