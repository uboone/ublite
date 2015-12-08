////////////////////////////////////////////////////////////////////////
// Class:       FEMemulator
// Module Type: FEMemulator
// File:        FEMemulator_module.cc
//
// Generated at Wed Dec  2 14:11:27 2015 by Taritree Wongjirad using artmod
// from cetpkgsupport v1_10_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ libraries
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

// ROOT
#include "TTree.h"
#include "TH1D.h"

// LArSoft 
#include "SimpleTypesAndConstants/RawTypes.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "Utilities/TimeService.h"

//Optical Channel Maps
#include "uboone/Geometry/UBOpChannelTypes.h"
#include "uboone/Geometry/UBOpReadoutMap.h"

//RawDigits
#include "RawData/raw.h"
#include "RawData/RawDigit.h"
#include "RawData/OpDetWaveform.h"

// Pulse finding
#include "uboone/OpticalDetectorAna/OpticalSubEvents/cfdiscriminator_algo/cfdiscriminator.hh"


class FEMemulator;

class FEMemulator : public art::EDAnalyzer {
public:
  explicit FEMemulator(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FEMemulator(FEMemulator const &) = delete;
  FEMemulator(FEMemulator &&) = delete;
  FEMemulator & operator = (FEMemulator const &) = delete;
  FEMemulator & operator = (FEMemulator &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;


private:

  // basic trigger
  void basicTrigger( int BeamWinSize, int NChannels, const std::vector< std::vector<int> >& chwfms, std::vector<int>& vmaxdiff, std::vector<int>& vmulti );

  // Declare member data here.

  // Configuration Parameters
  std::string fOpDataModule;
  int fFEMslot;
  int fNChannels;
  int fBeamWinSize;
  int fFrontBuffer;
  int fMinReadoutTicks;
  int fDiscr0precount;

  int fDiscr0threshold;
  int fDiscr0width;
  int fDiscr0delay;
  int fDiscr0deadtime;

  int fDiscr3threshold;
  int fDiscr3width;
  int fDiscr3delay;
  int fDiscr3deadtime;

  // Output tree and variables

  // Beam Window Tree
  TTree* fTwindow;
  int run;
  int subrun;
  int event;
  int winid;
  int maxdiff;
  int multiplicity;
  
  // records configuration
  TTree* fTconfig;

};


FEMemulator::FEMemulator(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  
  // Declare handle to analyzer file
  art::ServiceHandle<art::TFileService> out_file;

  // Load Parameters
  fOpDataModule    = p.get<std::string>( "OpDataModule",  "pmtreadout" );
  fFEMslot         = p.get<int>( "FEMslot",       5 );
  fNChannels       = p.get<int>( "NChannels",    32 );
  fBeamWinSize     = p.get<int>( "BeamWinSize", 103 );
  fMinReadoutTicks = p.get<int>( "MinReadoutTicks", 500 );
  fFrontBuffer     = p.get<int>( "FrontPorchBuffer", 0 );
  fDiscr0precount  = p.get<int>( "Discr0precount", 2 );

  fDiscr0threshold = p.get<int>( "Discr0threshold", 5 );
  fDiscr0width     = p.get<int>( "Discr0width", 5 );
  fDiscr0delay     = p.get<int>( "Discr0delay", 5 );
  fDiscr0deadtime  = p.get<int>( "Discr0deadtime", 5 );

  fDiscr3threshold = p.get<int>( "Discr3threshold", 10 );
  fDiscr3width     = p.get<int>( "Discr3width", 4 );
  fDiscr3delay     = p.get<int>( "Discr3delay", 4 );
  fDiscr3deadtime  = p.get<int>( "Discr3deadtime", 4 );

  // Setup Output Trees

  // save configuration
  fTconfig = out_file->make<TTree>( "config", "FEM emulator config. parameters" );
  fTconfig->Branch( "femslot", &fFEMslot, "femslot/I" );
  fTconfig->Branch( "beamwinsize", &fBeamWinSize, "beamwinsize/I" );
  fTconfig->Branch( "frontbuffer", &fFrontBuffer, "frontbuffer/fTconfig" );
  fTconfig->Branch( "discr0precount", &fDiscr0precount, "discr0precount/I" );
  fTconfig->Branch( "discr0threshold", &fDiscr0threshold, "discr0threshold/I" );
  fTconfig->Branch( "discr0width", &fDiscr0width, "discr0width/I" );
  fTconfig->Branch( "discr0delay", &fDiscr0delay, "discr0delay/I" );
  fTconfig->Branch( "discr0deadtime", &fDiscr0deadtime, "discr0deadtime/I" );
  fTconfig->Branch( "discr3threshold", &fDiscr3threshold, "discr3threshold/I" );
  fTconfig->Branch( "discr3width", &fDiscr3width, "discr3width/I" );
  fTconfig->Branch( "discr3delay", &fDiscr3delay, "discr3delay/I" );
  fTconfig->Branch( "discr3deadtime", &fDiscr3deadtime, "discr3deadtime/I" );
  fTconfig->Fill();

  // output
  fTwindow = out_file->make<TTree>( "windowtree", "Trigger Variables per window" );
  fTwindow->Branch( "run", &run, "run/I" );
  fTwindow->Branch( "subrun", &subrun, "subrun/I" );
  fTwindow->Branch( "event", &event, "event/I" );
  fTwindow->Branch( "winid", &winid, "winid/I" );
  fTwindow->Branch( "maxdiff", &maxdiff, "maxdiff/I" );
  fTwindow->Branch( "multiplicity", &multiplicity, "multiplicity/I" );

}

void FEMemulator::analyze(art::Event const & evt)
{
  // Implementation of required member function here.
  run    = (int)evt.run();
  subrun = (int)evt.subRun();
  event    = (int)evt.event();

  // initialize data handles and services
  art::ServiceHandle<geo::UBOpReadoutMap> ub_PMT_channel_map;
  art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
  //art::Handle< std::vector< raw::OpDetWaveform > > LogicHandle;
  ub_PMT_channel_map->SetOpMapRun( evt.run() );

  evt.getByLabel( fOpDataModule, "OpdetBeamHighGain", wfHandle);
  //evt.getByLabel( fOpDataModule, "UnspecifiedLogic" , LogicHandle);
  std::vector<raw::OpDetWaveform> const& opwfms(*wfHandle);
  //std::vector<raw::OpDetWaveform> const& logwfms(*LogicHandle);
  
  // first accumulate waveforms
  std::vector< std::vector<int> > wfms;
  wfms.resize( fNChannels );

  for(auto &wfm : opwfms)  {

    unsigned int readout_ch = wfm.ChannelNumber();
    unsigned int c,s,f;
    ub_PMT_channel_map->GetCrateSlotFEMChFromReadoutChannel(readout_ch, c, s, f);
    int slot = (int)s;
    int ch = (int)f%100;

    if ( slot!=fFEMslot )
      continue;
    if ( ch%100>=fNChannels )
      continue;
    if ( (int)wfm.size()<fMinReadoutTicks )
      continue;

    wfms[ch].resize( wfm.size(), 0 );
    for (int i=0; i<(int)wfm.size(); i++)
      wfms[ch][i] = (int)wfm[i];
  }
  
  // below is for the basic trigger
  std::vector<int> vmaxdiff;
  std::vector<int> vmulti;
  basicTrigger( fBeamWinSize, fNChannels, wfms, vmaxdiff, vmulti );
  for (int i=0; i<(int)vmaxdiff.size(); i++) {
    winid = i;
    maxdiff = vmaxdiff.at(i);
    multiplicity = vmulti.at(i);
    fTwindow->Fill();
  }

  // NN trigger


  // ZO trigger

}

// Maybe these various triggers should go elsewhere
void FEMemulator::basicTrigger( int BeamWinSize, int NChannels, const std::vector< std::vector<int> >& chwfms, std::vector<int>& vmaxdiff, std::vector<int>& vmulti ) {

  //std::cout << __PRETTY_FUNCTION__ << std::endl;

  // first calculate accumulators for each waveform
  std::vector<int> chdiff[NChannels];
  std::vector<int> chhit[NChannels];
  for (int ch=0; ch<NChannels; ch++) {

    const std::vector<int>& wfm = chwfms[ch];

    chdiff[ch].resize( wfm.size(), 0 );
    chhit[ch].resize( wfm.size(), 0 );

    // memory for diff vectors
    std::vector<int> diff0( (int)wfm.size(), 0 );
    std::vector<int> diff3( (int)wfm.size(), 0 );
    for (int tick=fDiscr0delay; tick<(int)wfm.size(); tick++)
      diff0.at(tick) = wfm.at(tick)-wfm.at(tick-fDiscr0delay);
    for (int tick=fDiscr3delay; tick<(int)wfm.size(); tick++)
      diff3.at(tick) = wfm.at(tick)-wfm.at(tick-fDiscr3delay);


    // determine triggers and fill accumulators
    std::vector<int> ttrig0;
    std::vector<int> ttrig3;
    
    for (int tick=0; tick<(int)wfm.size(); tick++) {
      // discr0 must fire first: acts as pre-trigger. won't fire again until all discs are no longer active
      if ( diff0.at(tick)>=fDiscr0threshold ) {
	if ( (ttrig0.size()==0 || ttrig0.at( ttrig0.size()-1 )+fDiscr0precount<tick )
	     && ( ttrig3.size()==0 || ttrig3.at( ttrig3.size()-1 )+fDiscr3deadtime < tick ) ) {
	  // form discr0 trigger
	  ttrig0.push_back( tick );
	}
      } // end of if discr0 fires

      // discr3 fire
      if ( diff3.at(tick)>=fDiscr3threshold ) {
	// must be within discr0 prewindow and outside of past discr3 deadtime
	if ( ( ttrig0.size()>0 && tick-ttrig0.at( ttrig0.size()-1 ) < fDiscr0deadtime )
	     && ( ttrig3.size()==0 || ttrig3.at( ttrig3.size()-1 ) + fDiscr3deadtime < tick ) ) {
	  ttrig3.push_back( tick );
	  // find maxdiff
	  int tmaxdiff = diff3.at(tick);
	  for (int t=tick; t<std::min( tick+fDiscr3width, (int)diff3.size() ); t++) {
	    if ( tmaxdiff<diff3.at(t) )
	      tmaxdiff = diff3.at(t);
	  }
	  // fill the accumulators
	  int tend = std::min( tick+fDiscr3deadtime, (int)diff3.size() );
	  for (int t=tick; t<tend; t++) {
	    chdiff[ch].at( t ) = tmaxdiff;
	    chhit[ch].at( t ) = 1;
	  }
	}
      }
    }//end of wfm loop for trigger and accumulators
  }//end of channel loop

  // break up waveform into windows and calculate trigger vars for each window
  int wfmsize = (int)chwfms.at(0).size();
  if ( wfmsize < fMinReadoutTicks ) {
    std::cout << "Beam readout window size is too small! (" << wfmsize << " < " << fMinReadoutTicks << ")" << std::endl;
    return;
  }

  int nwindows = (wfmsize-1-fFrontBuffer)/fBeamWinSize;
  vmaxdiff.clear();
  vmulti.clear();

  for (int iwin=0; iwin<nwindows; iwin++) {
    int winstart = fFrontBuffer + iwin*fBeamWinSize;
    int winend   = fFrontBuffer + (iwin+1)*fBeamWinSize;
    winid = iwin;
    int winmaxmulti = 0;
    int winmaxdiff = 0;
    for (int tick=winstart; tick<winend; tick++) {
      int maxdiff_ = 0;
      int nhit_ = 0;
      for (int ch=0; ch<NChannels; ch++) {
	maxdiff_ += chdiff[ch].at(tick);
	nhit_    += chhit[ch].at(tick);
      }
      if ( winmaxdiff < maxdiff_ )
	winmaxdiff = maxdiff_;
      if ( winmaxmulti < nhit_ )
	winmaxmulti = nhit_;
    }

    // store for the window
    vmaxdiff.push_back( winmaxdiff );
    vmulti.push_back( winmaxmulti );
    
  }
  
}

DEFINE_ART_MODULE(FEMemulator)
