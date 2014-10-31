#ifndef SCANNERALGO_TEMPLATE_H
#define SCANNERALGO_TEMPLATE_H

/*
  This file defines certain specilization of templated functions.
  In particular it implements:

  ScannerAlgo::ScanData
  ScannerAlgo::GetPtrMap
  ScannerAlgo::LiteDataType

  One has to implement the relevant specilization for introducing a new data product!

 */

namespace larlite {

  //
  // ScanData specialization
  //
  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::simb::MCTruth> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;
    auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_mctruth*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      const art::Ptr<::simb::MCTruth> mct_ptr(dh,i);
      
      larlite::mctruth lite_mct;
      
      // Generator (Origin) ID
      lite_mct.SetOrigin((larlite::simb::Origin_t) mct_ptr->Origin() );

      // Particle Information
      for(size_t j=0; j<(size_t)(mct_ptr->NParticles()); ++j) {
	
	const ::simb::MCParticle lar_part(mct_ptr->GetParticle(j));
	
	larlite::mcpart lite_part(lar_part.TrackId(),
				  lar_part.PdgCode(),
				  lar_part.Process(), 
				  lar_part.Mother(),
				  lar_part.Mass(),
				  lar_part.StatusCode());
	
	lite_part.SetPolarization(lar_part.Polarization());
	lite_part.SetRescatter(lar_part.Rescatter());
	lite_part.SetWeight(lar_part.Weight());
	
	for(size_t k=0; k<(size_t)(lar_part.NumberDaughters()); ++k)
	  
	  lite_part.AddDaughter(lar_part.Daughter(k));
	
	// Add trajectory points
	larlite::mctrajectory lite_track;
	lite_track.reserve(lar_part.NumberTrajectoryPoints());
	
	bool   loopInFV=false;
	size_t start_FV=0;
	for(size_t l=0; l<(size_t)lar_part.NumberTrajectoryPoints(); ++l) {
	  
	  lite_track.push_back(lar_part.Position(l),lar_part.Momentum(l));
	  
	  // Record fiducial volume tracking
	  bool inFV = false;//IsFV(lar_part.Vx(l),lar_part.Vy(l),lar_part.Vz(l),lar_part.T(l));
	  
	  if(!loopInFV) {
	    
	    if(inFV) { loopInFV=true; start_FV=l; }
	    
	  }else if(!inFV) {
	    lite_part.AddFiducialTrack(start_FV,l-1);
	    loopInFV=false;
	  }
	}
	if(loopInFV){ lite_part.AddFiducialTrack(start_FV,lar_part.NumberTrajectoryPoints()-1); }
	
	lite_part.SetTrajectory(lite_track);
	
	lite_mct.Add(lite_part);
      }
      
      // Neutrino Information
      const ::simb::MCNeutrino lar_nu(mct_ptr->GetNeutrino());
      
      lite_mct.SetNeutrino( lar_nu.CCNC(),
			    lar_nu.Mode(),
			    lar_nu.InteractionType(),
			    lar_nu.Target(),
			    lar_nu.HitNuc(),
			    lar_nu.HitQuark(),
			    lar_nu.W(),
			    lar_nu.X(),
			    lar_nu.Y(),
			    lar_nu.QSqr() );

      // Store address map for downstream association
      fPtrIndex_mctruth[mct_ptr] = std::make_pair(lite_data->size(),name_index);

      // Save
      lite_data->push_back(lite_mct);
    }
  }
  
  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::simb::GTruth> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_gtruth*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      const art::Ptr<::simb::GTruth> gtruth_ptr(dh,i);
      
      larlite::gtruth lite_gtruth;
      
      lite_gtruth.fGint        = gtruth_ptr->fGint;
      lite_gtruth.fGscatter    = gtruth_ptr->fGscatter;
      lite_gtruth.fweight      = gtruth_ptr->fweight;
      lite_gtruth.fprobability = gtruth_ptr->fprobability;
      lite_gtruth.fXsec        = gtruth_ptr->fXsec;
      lite_gtruth.fDiffXsec    = gtruth_ptr->fDiffXsec;
      lite_gtruth.fNumPiPlus   = gtruth_ptr->fNumPiPlus;
      lite_gtruth.fNumPiMinus  = gtruth_ptr->fNumPiMinus;
      lite_gtruth.fNumPi0      = gtruth_ptr->fNumPi0;
      lite_gtruth.fNumProton   = gtruth_ptr->fNumProton;
      lite_gtruth.fNumNeutron  = gtruth_ptr->fNumNeutron;
      lite_gtruth.fIsCharm     = gtruth_ptr->fIsCharm;
      lite_gtruth.fResNum      = gtruth_ptr->fResNum;
      lite_gtruth.fgQ2         = gtruth_ptr->fgQ2;
      lite_gtruth.fgq2         = gtruth_ptr->fgq2;
      lite_gtruth.fgW          = gtruth_ptr->fgW;
      lite_gtruth.fgT          = gtruth_ptr->fgT;
      lite_gtruth.fgX          = gtruth_ptr->fgX;
      lite_gtruth.fgY          = gtruth_ptr->fgY;
      lite_gtruth.fFShadSystP4 = gtruth_ptr->fFShadSystP4;
      lite_gtruth.fIsSeaQuark  = gtruth_ptr->fIsSeaQuark;
      lite_gtruth.fHitNucP4    = gtruth_ptr->fHitNucP4;
      lite_gtruth.ftgtZ        = gtruth_ptr->ftgtZ;
      lite_gtruth.ftgtA        = gtruth_ptr->ftgtA;
      lite_gtruth.ftgtPDG      = gtruth_ptr->ftgtPDG;

      fPtrIndex_gtruth[gtruth_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_gtruth);
      
    }
    
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::simb::MCParticle> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_mcpart*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      art::Ptr<::simb::MCParticle> mcparticle_ptr(dh,i);
      
      larlite::mcpart mcpart_lite(mcparticle_ptr->TrackId(),
				  mcparticle_ptr->PdgCode(),
				  mcparticle_ptr->Process(), 
				  mcparticle_ptr->Mother(),
				  mcparticle_ptr->Mass(),
				  mcparticle_ptr->StatusCode());
      
      mcpart_lite.SetPolarization(mcparticle_ptr->Polarization());
      mcpart_lite.SetRescatter(mcparticle_ptr->Rescatter());
      mcpart_lite.SetWeight(mcparticle_ptr->Weight());
      
      for(size_t k=0; k<(size_t)(mcparticle_ptr->NumberDaughters()); ++k)
	mcpart_lite.AddDaughter(mcparticle_ptr->Daughter(k));
      
      // Add trajectory points
      larlite::mctrajectory lite_track;
      lite_track.reserve(mcparticle_ptr->NumberTrajectoryPoints());
      
      bool   loopInFV=false;
      size_t start_FV=0;
      for(size_t l=0; l<(size_t)mcparticle_ptr->NumberTrajectoryPoints(); ++l) {
	
	lite_track.push_back(mcparticle_ptr->Position(l),
			     mcparticle_ptr->Momentum(l));
	
	// Record fiducial volume tracking
	bool inFV = false;//IsFV(mcparticle_ptr->Vx(l), mcparticle_ptr->Vy(l),mcparticle_ptr->Vz(l),mcparticle_ptr->T(l));
	
	if(!loopInFV) {
	  
	  if(inFV) { loopInFV=true; start_FV=l; }
	  
	}else if(!inFV) {
	  mcpart_lite.AddFiducialTrack(start_FV,l-1);
	  loopInFV=false;
	}
      }
      if(loopInFV){ mcpart_lite.AddFiducialTrack(start_FV,mcparticle_ptr->NumberTrajectoryPoints()-1); }
      
      mcpart_lite.SetTrajectory(lite_track);
      
      // Store address map for downstream association
      fPtrIndex_mcpart[mcparticle_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(mcpart_lite);
    }
  }
  
  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::simb::MCFlux> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_mcflux*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      const art::Ptr<::simb::MCFlux> mcflux_ptr(dh,i);
      
      larlite::mcflux lite_mcflux;

      lite_mcflux.frun = mcflux_ptr->frun;
      lite_mcflux.fevtno = mcflux_ptr->fevtno;
      lite_mcflux.fndxdz = mcflux_ptr->fndxdz;
      lite_mcflux.fndydz = mcflux_ptr->fndydz;
      lite_mcflux.fnpz = mcflux_ptr->fnpz;
      lite_mcflux.fnenergy = mcflux_ptr->fnenergy;
      lite_mcflux.fndxdznea = mcflux_ptr->fndxdznea;
      lite_mcflux.fndydznea = mcflux_ptr->fndydznea;
      lite_mcflux.fnenergyn = mcflux_ptr->fnenergyn;
      lite_mcflux.fnwtnear = mcflux_ptr->fnwtnear;
      lite_mcflux.fndxdzfar = mcflux_ptr->fndxdzfar;
      lite_mcflux.fndydzfar = mcflux_ptr->fndydzfar;
      lite_mcflux.fnenergyf = mcflux_ptr->fnenergyf;
      lite_mcflux.fnwtfar = mcflux_ptr->fnwtfar;
      lite_mcflux.fnorig = mcflux_ptr->fnorig;
      lite_mcflux.fndecay = mcflux_ptr->fndecay;
      lite_mcflux.fntype = mcflux_ptr->fntype;
      lite_mcflux.fvx = mcflux_ptr->fvx;
      lite_mcflux.fvy = mcflux_ptr->fvy;
      lite_mcflux.fvz = mcflux_ptr->fvz;
      lite_mcflux.fpdpx = mcflux_ptr->fpdpx;
      lite_mcflux.fpdpy = mcflux_ptr->fpdpy;
      lite_mcflux.fpdpz = mcflux_ptr->fpdpz;
      lite_mcflux.fppdxdz = mcflux_ptr->fppdxdz;
      lite_mcflux.fppdydz = mcflux_ptr->fppdydz;
      lite_mcflux.fpppz = mcflux_ptr->fpppz;
      lite_mcflux.fppenergy = mcflux_ptr->fppenergy;
      lite_mcflux.fppmedium = mcflux_ptr->fppmedium;
      lite_mcflux.fptype = mcflux_ptr->fptype;
      lite_mcflux.fppvx = mcflux_ptr->fppvx;
      lite_mcflux.fppvy = mcflux_ptr->fppvy;
      lite_mcflux.fppvz = mcflux_ptr->fppvz;
      lite_mcflux.fmuparpx = mcflux_ptr->fmuparpx;
      lite_mcflux.fmuparpy = mcflux_ptr->fmuparpy;
      lite_mcflux.fmuparpz = mcflux_ptr->fmuparpz;
      lite_mcflux.fmupare = mcflux_ptr->fmupare;
      lite_mcflux.fnecm = mcflux_ptr->fnecm;
      lite_mcflux.fnimpwt = mcflux_ptr->fnimpwt;
      lite_mcflux.fxpoint = mcflux_ptr->fxpoint;
      lite_mcflux.fypoint = mcflux_ptr->fypoint;
      lite_mcflux.fzpoint = mcflux_ptr->fzpoint;
      lite_mcflux.ftvx = mcflux_ptr->ftvx;
      lite_mcflux.ftvy = mcflux_ptr->ftvy;
      lite_mcflux.ftvz = mcflux_ptr->ftvz;
      lite_mcflux.ftpx = mcflux_ptr->ftpx;
      lite_mcflux.ftpy = mcflux_ptr->ftpy;
      lite_mcflux.ftpz = mcflux_ptr->ftpz;
      lite_mcflux.ftptype = mcflux_ptr->ftptype;

      lite_mcflux.ftgen = mcflux_ptr->ftgen;
      lite_mcflux.ftgptype = mcflux_ptr->ftgptype;

      lite_mcflux.ftgppx = mcflux_ptr->ftgppx;
      lite_mcflux.ftgppy = mcflux_ptr->ftgppy;
      lite_mcflux.ftgppz = mcflux_ptr->ftgppz;
      lite_mcflux.ftprivx = mcflux_ptr->ftprivx;
      lite_mcflux.ftprivy = mcflux_ptr->ftprivy;
      lite_mcflux.ftprivz = mcflux_ptr->ftprivz;
      lite_mcflux.fbeamx = mcflux_ptr->fbeamx;
      lite_mcflux.fbeamy = mcflux_ptr->fbeamy;
      lite_mcflux.fbeamz = mcflux_ptr->fbeamz;
      lite_mcflux.fbeampx = mcflux_ptr->fbeampx;
      lite_mcflux.fbeampy = mcflux_ptr->fbeampy;
      lite_mcflux.fbeampz = mcflux_ptr->fbeampz;

      lite_mcflux.fFluxType = (::larlite::simb::flux_code_)(mcflux_ptr->fFluxType);

      lite_mcflux.fgenx = mcflux_ptr->fgenx;

      lite_mcflux.fgeny = mcflux_ptr->fgeny;
      lite_mcflux.fgenz = mcflux_ptr->fgenz;
      lite_mcflux.fdk2gen = mcflux_ptr->fdk2gen;

      lite_mcflux.fgen2vtx = mcflux_ptr->fgen2vtx;

      lite_mcflux.SetFluxGen( mcflux_ptr->Flux( 12,::simb::kGenerator), 
			      mcflux_ptr->Flux(-12,::simb::kGenerator),
			      mcflux_ptr->Flux( 14,::simb::kGenerator),
			      mcflux_ptr->Flux(-14,::simb::kGenerator),
			      mcflux_ptr->Flux( 16,::simb::kGenerator),
			      mcflux_ptr->Flux(-16,::simb::kGenerator) );

      lite_mcflux.SetFluxPos( mcflux_ptr->Flux( 12,::simb::kHistPlusFocus), 
			      mcflux_ptr->Flux(-12,::simb::kHistPlusFocus),
			      mcflux_ptr->Flux( 14,::simb::kHistPlusFocus),
			      mcflux_ptr->Flux(-14,::simb::kHistPlusFocus),
			      mcflux_ptr->Flux( 16,::simb::kHistPlusFocus),
			      mcflux_ptr->Flux(-16,::simb::kHistPlusFocus) );

      lite_mcflux.SetFluxPos( mcflux_ptr->Flux( 12,::simb::kHistMinusFocus), 
			      mcflux_ptr->Flux(-12,::simb::kHistMinusFocus),
			      mcflux_ptr->Flux( 14,::simb::kHistMinusFocus),
			      mcflux_ptr->Flux(-14,::simb::kHistMinusFocus),
			      mcflux_ptr->Flux( 16,::simb::kHistMinusFocus),
			      mcflux_ptr->Flux(-16,::simb::kHistMinusFocus) );

      // Register pointer to association look-up map
      fPtrIndex_mcflux[mcflux_ptr] = std::make_pair(lite_data->size(),name_index);

      lite_data->push_back(lite_mcflux);
    }
  }
  
  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::sim::MCShower> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_mcshower*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {

      art::Ptr<sim::MCShower> mcs_ptr(dh,i);
      
      larlite::mcshower lite_mcs;
      lite_mcs.SetMotherID(mcs_ptr->MotherPDGID(), mcs_ptr->MotherTrackID());
      
      lite_mcs.SetMotherPoint(mcs_ptr->MotherPosition());
      
      lite_mcs.SetMotherProcess(mcs_ptr->MotherCreationProcess());
      
      lite_mcs.SetMotherMomentum(mcs_ptr->MotherMomentum());
      
      lite_mcs.SetDaughterTrackList(mcs_ptr->DaughterTrackID());
      
      lite_mcs.SetDaughterMomentum(mcs_ptr->DaughterMomentum());
      lite_mcs.SetDaughterPosition(mcs_ptr->DaughterPosition());
      
      lite_mcs.SetPlaneCharge(mcs_ptr->Charge());
      
      fPtrIndex_mcshower[mcs_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_mcs);
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::sim::SimChannel> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_simch*)lite_dh;

    for(size_t i=0; i<dh->size(); ++i ) {
      
      const art::Ptr<::sim::SimChannel> sch_ptr(dh,i);
      
      std::map<unsigned short,std::vector<sim::IDE> > sch_map(sch_ptr->TDCIDEMap());
      
      larlite::simch lite_sch;
      lite_sch.set_channel(sch_ptr->Channel());
      
      for(auto sch_iter = sch_map.begin(); sch_iter!=sch_map.end(); ++sch_iter) {
	
	unsigned short tdc = (*sch_iter).first;
	
	for(auto const this_ide : (*sch_iter).second) {
	  
	  larlite::ide lite_ide;
	  lite_ide.trackID = this_ide.trackID;
	  lite_ide.numElectrons = this_ide.numElectrons;
	  lite_ide.energy = this_ide.energy;
	  lite_ide.x = this_ide.x;
	  lite_ide.y = this_ide.y;
	  lite_ide.z = this_ide.z;
	  
	  lite_sch.add_ide(tdc,lite_ide);
	}
      }

      fPtrIndex_simch[sch_ptr] = std::make_pair(lite_data->size(),name_index);

      lite_data->push_back(lite_sch);
    }

  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::Wire> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_wire*)lite_dh;

    for(size_t i=0; i<dh->size(); i++){

      const art::Ptr<::recob::Wire> wire_ptr(dh,i);
    
      larlite::wire wire_lite(wire_ptr->Signal(),
			      wire_ptr->Channel(),
			      (larlite::geo::View_t)(wire_ptr->View()),
			      (larlite::geo::SigType_t)(wire_ptr->SignalType()));

      fPtrIndex_wire[wire_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(wire_lite);
    }  
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::Hit> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;
    auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());    
    auto lite_data = (::larlite::event_hit*)lite_dh;
    art::ServiceHandle<::geo::Geometry> geo;

    for(size_t i=0; i<dh->size(); i++){
      
      art::Ptr<::recob::Hit> hit_ptr(dh,i);
      
      larlite::hit lite_hit;
      
      lite_hit.set_waveform(hit_ptr->fHitSignal);
      lite_hit.set_times(hit_ptr->StartTime(),
		       hit_ptr->PeakTime(),
			 hit_ptr->EndTime());
      lite_hit.set_times_err(hit_ptr->SigmaStartTime(),
			     hit_ptr->SigmaPeakTime(),
			     hit_ptr->SigmaEndTime());
      lite_hit.set_charge(hit_ptr->Charge(),hit_ptr->Charge(true));
      lite_hit.set_charge_err(hit_ptr->SigmaCharge(),hit_ptr->SigmaCharge(true));
      lite_hit.set_multiplicity(hit_ptr->Multiplicity());
      lite_hit.set_channel(geo->PlaneWireToChannel(hit_ptr->WireID()));
      //lite_hit.set_channel(hit_ptr->Channel());
      lite_hit.set_wire(hit_ptr->WireID().Wire);
      lite_hit.set_fit_goodness(hit_ptr->GoodnessOfFit());
      lite_hit.set_view((larlite::geo::View_t)(hit_ptr->View()));
      lite_hit.set_sigtype((larlite::geo::SigType_t)(hit_ptr->SignalType()));
      
      // Store address map for downstream association
      fPtrIndex_hit[hit_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_hit);
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::OpHit> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_ophit*)lite_dh;

    for(size_t i=0; i<dh->size(); ++i) {
      
      art::Ptr<::recob::OpHit> hit_ptr(dh,i);
      
      ::larlite::ophit lite_hit( hit_ptr->OpChannel(),
				 hit_ptr->PeakTime(),
				 hit_ptr->PeakTimeAbs(),
				 hit_ptr->Frame(),
				 hit_ptr->Width(),
				 hit_ptr->Area(),
				 hit_ptr->Amplitude(),
				 hit_ptr->PE(),
				 hit_ptr->FastToTotal() );

      fPtrIndex_ophit[hit_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_hit);
    }
    
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::OpFlash> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_opflash*)lite_dh;
    art::ServiceHandle<::geo::Geometry> geo;

    for(size_t i=0; i<dh->size(); ++i) {

      art::Ptr<::recob::OpFlash> flash_ptr(dh,i);
      
      std::vector<double> pe_per_opdet;
      pe_per_opdet.reserve(geo->NOpChannels());
      for(size_t j=0; j<geo->NOpChannels(); ++j)
	pe_per_opdet.push_back(flash_ptr->PE(j));
      
      ::larlite::opflash lite_flash( flash_ptr->Time(),
				     flash_ptr->TimeWidth(),
				     flash_ptr->AbsTime(),
				     flash_ptr->Frame(),
				     pe_per_opdet,
				     flash_ptr->InBeamFrame(),
				     flash_ptr->OnBeamTime(),
				     flash_ptr->FastToTotal(),
				     flash_ptr->YCenter(),
				     flash_ptr->YWidth(),
				     flash_ptr->ZCenter(),
				     flash_ptr->ZWidth(),
				     flash_ptr->WireCenters(),
				     flash_ptr->WireWidths());
      
      fPtrIndex_opflash[flash_ptr] = std::make_pair(lite_data->size(),name_index);

      lite_data->push_back(lite_flash);
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::anab::CosmicTag> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_cosmictag*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      art::Ptr<::anab::CosmicTag> tag_ptr(dh,i);
      
      larlite::cosmictag lite_tag(tag_ptr->endPt1,
				  tag_ptr->endPt2,
				  tag_ptr->fCosmicScore,
				  (::larlite::anab::CosmicTagID_t)(tag_ptr->fCosmicType));
      
      // store product ptr for association
      fPtrIndex_cosmictag[tag_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_tag);
      
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::Cluster> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_cluster*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      art::Ptr<::recob::Cluster> cluster_ptr(dh,i);
      
      larlite::cluster lite_cluster;
      lite_cluster.set_charge(cluster_ptr->Charge());
      lite_cluster.set_dtdw(cluster_ptr->dTdW());
      lite_cluster.set_dqdw(cluster_ptr->dQdW());
      lite_cluster.set_dtdw_err(cluster_ptr->SigmadTdW());
      lite_cluster.set_dqdw_err(cluster_ptr->SigmadQdW());
      lite_cluster.set_id(cluster_ptr->ID());
      lite_cluster.set_view((larlite::geo::View_t)(cluster_ptr->View()));
      lite_cluster.set_start_vtx(cluster_ptr->StartPos());
      lite_cluster.set_end_vtx(cluster_ptr->EndPos());
      lite_cluster.set_start_vtx_err(cluster_ptr->SigmaStartPos());
      lite_cluster.set_end_vtx_err(cluster_ptr->SigmaEndPos());
      
      // Store address map for downstream association
      fPtrIndex_cluster[cluster_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_cluster);
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::Seed> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_seed*)lite_dh;

    double fSeedPoint[3];
    double fSeedDirection[3];
    double fSeedPointError[3];
    double fSeedDirectionError[3];

    for(size_t i=0; i<dh->size(); ++i) {
      
      art::Ptr<::recob::Seed> seed_ptr(dh,i);
      larlite::seed lite_seed;

      seed_ptr->GetPoint(fSeedPoint,fSeedPointError);
      seed_ptr->GetDirection(fSeedDirection,fSeedDirectionError);

      lite_seed.SetPoint(fSeedPoint,fSeedPointError);
      lite_seed.SetDirection(fSeedDirection,fSeedDirectionError);
      lite_seed.SetValidity(seed_ptr->IsValid());
      
      // Store address map for downstream association
      fPtrIndex_seed[seed_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_seed);
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::EndPoint2D> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_endpoint2d*)lite_dh;
    for(size_t i=0; i < dh->size(); ++i) {
      
      const art::Ptr<::recob::EndPoint2D> end2d_ptr(dh,i);
      
      // get vertex info
      larlite::endpoint2d lite_end2d(end2d_ptr->DriftTime(),
				     end2d_ptr->WireID().Wire,
				     end2d_ptr->Strength(),
				     end2d_ptr->ID(),
				     (larlite::geo::View_t)(end2d_ptr->View()),
				     end2d_ptr->Charge());
      
      fPtrIndex_end2d[end2d_ptr] = std::make_pair(lite_data->size(),name_index);
      
      // Store data
      lite_data->push_back(lite_end2d);
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::SpacePoint> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_spacepoint*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      art::Ptr<::recob::SpacePoint> spacepoint_ptr(dh,i);
      
      larlite::spacepoint lite_spacepoint(spacepoint_ptr->ID(),
					  spacepoint_ptr->XYZ()[0],
					  spacepoint_ptr->XYZ()[1],
					  spacepoint_ptr->XYZ()[2],
					  spacepoint_ptr->ErrXYZ()[0],
					  spacepoint_ptr->ErrXYZ()[1],
					  spacepoint_ptr->ErrXYZ()[2],
					  spacepoint_ptr->Chisq() );
      
      // Store address map for downstream association
      fPtrIndex_sps[spacepoint_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_spacepoint);
    }
  }
  
  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::Track> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_track*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      art::Ptr<::recob::Track> track_ptr(dh,i);
      
      larlite::track track_lite;
      track_lite.set_track_id ( track_ptr->ID()   );
      // Direction & points
      for(size_t i=0; i<track_ptr->NumberTrajectoryPoints(); i++) {
	track_lite.add_vertex     (track_ptr->LocationAtPoint(i));
	track_lite.add_direction  (track_ptr->DirectionAtPoint(i));
      }
      // Covariance
      for(size_t i=0; i<track_ptr->NumberCovariance(); i++)
	track_lite.add_covariance (track_ptr->CovarianceAtPoint(i));
      // Momentum
      for(size_t i=0; i<track_ptr->NumberFitMomentum(); i++)
	track_lite.add_momentum   (track_ptr->MomentumAtPoint(i));
      
      // Store address map for downstream association
      fPtrIndex_track[track_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(track_lite);
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::Shower> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_shower*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      const art::Ptr<::recob::Shower> shower_ptr(dh,i);
      
      larlite::shower lite_shower;
      
      lite_shower.set_id(shower_ptr->ID());
      //light_shower.set_total_charge(shower_ptr->TotalCharge());
      lite_shower.set_direction(shower_ptr->Direction());
      lite_shower.set_direction_err(shower_ptr->DirectionErr());
      //light_shower.set_max_width(shower_ptr->MaxWidthX(),shower_ptr->MaxWidthY());
      //light_shower.set_distance_max_width(shower_ptr->DistanceMaxWidth());
      
      fPtrIndex_shower[shower_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_shower);
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::Vertex> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_vertex*)lite_dh;
    Double_t xyz[3]={0};
    
    for(size_t i=0; i < dh->size(); ++i) {
      
      const art::Ptr<::recob::Vertex> vtx_ptr(dh,i);
      
      vtx_ptr->XYZ(xyz);
      
      // get vertex info
      larlite::vertex lite_vtx(xyz,
			       vtx_ptr->ID());
      
      fPtrIndex_vertex[vtx_ptr] = std::make_pair(lite_data->size(),name_index);
      
      // Store data
      lite_data->push_back(lite_vtx);
      
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::anab::Calorimetry> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_calorimetry*)lite_dh;
    for(size_t i=0; i < dh->size(); ++i) {
      
      const art::Ptr<::anab::Calorimetry> calo_ptr(dh,i);
      
      larlite::calorimetry lite_calo;
      
      lite_calo.set_dedx(calo_ptr->dEdx());
      lite_calo.set_dqdx(calo_ptr->dQdx());
      lite_calo.set_residual_range(calo_ptr->ResidualRange());
      lite_calo.set_deadwire_range(calo_ptr->DeadWireResRC());
      lite_calo.set_kinetic_energy(calo_ptr->KineticEnergy());
      lite_calo.set_range(calo_ptr->Range());
      lite_calo.set_track_pitch(calo_ptr->TrkPitchVec());
      
      fPtrIndex_calo[calo_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_calo);
      
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::anab::ParticleID> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_partid*)lite_dh;
    
    for(size_t i=0; i<dh->size(); ++i) {
      
      art::Ptr<::anab::ParticleID> partid_ptr(dh,i);
      
      larlite::partid lite_partid( partid_ptr->Pdg(),
				   partid_ptr->Ndf(),
				   partid_ptr->MinChi2(),
				   partid_ptr->DeltaChi2(),
				   partid_ptr->Chi2Proton(),
				   partid_ptr->Chi2Kaon(),
				   partid_ptr->Chi2Pion(),
				   partid_ptr->Chi2Muon(),
				   partid_ptr->MissingE(),
				   partid_ptr->MissingEavg(),
				   partid_ptr->PIDA() );
      
      fPtrIndex_partid[partid_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_partid);
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::PFParticle> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_pfpart*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      art::Ptr<::recob::PFParticle> pfpart_ptr(dh,i);
      
      larlite::pfpart lite_pfpart(pfpart_ptr->PdgCode(),
				  pfpart_ptr->Self(),
				  pfpart_ptr->Parent(),
				  pfpart_ptr->Daughters());
      
      fPtrIndex_pfpart[pfpart_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_pfpart);
    }
  }

  template <class T>
  void ScanData(art::Handle<std::vector<T> > const &dh,
		::larlite::event_base* lite_dh)
  { throw cet::exception(__PRETTY_FUNCTION__)<<"Not implemented!"; }

  //
  // Getter for associated data product pointer 
  //
  template <> const std::map<art::Ptr< ::simb::MCTruth>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { return fPtrIndex_mctruth; }

  template <> const std::map<art::Ptr< ::simb::GTruth>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { return fPtrIndex_gtruth; }

  template <> const std::map<art::Ptr< ::simb::MCFlux>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { return fPtrIndex_mcflux; }

  template <> const std::map<art::Ptr< ::simb::MCParticle>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { return fPtrIndex_mcpart; }

  template <> const std::map<art::Ptr< ::sim::SimChannel>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { return fPtrIndex_simch; }

  template <> const std::map<art::Ptr< ::sim::MCShower>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { return fPtrIndex_mcshower; }

  template <> const std::map<art::Ptr< ::recob::OpHit>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { return fPtrIndex_ophit; }

  template <> const std::map<art::Ptr< ::recob::OpFlash>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { return fPtrIndex_opflash; }

  template <> const std::map<art::Ptr< ::recob::Hit>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { return fPtrIndex_hit; }

  template <> const std::map<art::Ptr< ::recob::Wire>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { return fPtrIndex_wire; }

  template <> const std::map<art::Ptr< ::recob::Cluster>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { return fPtrIndex_cluster; }

  template <> const std::map<art::Ptr< ::recob::Track>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { return fPtrIndex_track; }

  template <> const std::map<art::Ptr< ::recob::Shower>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { return fPtrIndex_shower; }

  template <> const std::map<art::Ptr< ::recob::Vertex>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { return fPtrIndex_vertex; }

  template <> const std::map<art::Ptr< ::recob::SpacePoint>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { return fPtrIndex_sps; }

  template <> const std::map<art::Ptr< ::recob::EndPoint2D>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { return fPtrIndex_end2d; }

  template <> const std::map<art::Ptr< ::anab::CosmicTag>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { return fPtrIndex_cosmictag; }

  template <> const std::map<art::Ptr< ::anab::Calorimetry>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { return fPtrIndex_calo; }

  template <> const std::map<art::Ptr< ::anab::ParticleID>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { return fPtrIndex_partid; }

  template <> const std::map<art::Ptr< ::recob::PFParticle>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { return fPtrIndex_pfpart; }

  template <class T>
  const std::map<art::Ptr<T>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { throw cet::exception(__PRETTY_FUNCTION__)<<"Not implemented for a specified data product type..."; }

  //
  // Type identifier functions
  //
  // simb
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::simb::GTruth> () const
  { return ::larlite::data::kGTruth; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::simb::MCTruth> () const
  { return ::larlite::data::kMCTruth; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::simb::MCParticle> () const
  { return ::larlite::data::kMCParticle; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::simb::MCFlux> () const
  { return ::larlite::data::kMCFlux; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::simb::MCTrajectory> () const
  { return ::larlite::data::kMCTrajectory; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::simb::MCNeutrino> () const
  { return ::larlite::data::kMCNeutrino; }
  // sim
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::sim::SimChannel> () const
  { return ::larlite::data::kSimChannel; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::sim::MCShower> () const
  { return ::larlite::data::kMCShower; }
  // raw
  // recob
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::Wire> () const
  { return ::larlite::data::kWire; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::Hit> () const
  { return ::larlite::data::kHit; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::Cluster> () const
  { return ::larlite::data::kCluster; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::SpacePoint> () const
  { return ::larlite::data::kSpacePoint; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::OpHit> () const
  { return ::larlite::data::kOpHit; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::OpFlash> () const
  { return ::larlite::data::kOpFlash; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::Seed> () const
  { return ::larlite::data::kSeed; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::Track> () const
  { return ::larlite::data::kTrack; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::Shower> () const
  { return ::larlite::data::kShower; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::Vertex> () const
  { return ::larlite::data::kVertex; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::EndPoint2D> () const
  { return ::larlite::data::kEndPoint2D; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::PFParticle> () const
  { return ::larlite::data::kPFParticle; }
  // anab
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::anab::CosmicTag> () const
  { return ::larlite::data::kCosmicTag; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::anab::Calorimetry> () const
  { return ::larlite::data::kCalorimetry; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::anab::ParticleID> () const
  { return ::larlite::data::kParticleID; }

  //
  // LocateLiteProduct implementation
  //
  template <class T>
  bool ScannerAlgo::LocateLiteProduct(art::Ptr<T> const ptr,
				      std::pair<size_t,size_t> &loc) const
  { 
    auto ptr_map = GetPtrMap<T>();
    auto id_iter = ptr_map.find(ptr);
    if(id_iter == ptr_map.end()) return false; 
    
    loc.first = (*id_iter).second.first;
    loc.second = (*id_iter).second.second;
    return true;
  }

  //
  // ScanAssociation implementation 
  //
  template <class T,class U>
  void ScannerAlgo::ScanAssociation(art::Event const& e,
				    art::Handle<std::vector<T> > &dh,
				    ::larlite::event_base* lite_dh) const
  { 
    art::FindManyP<U> ptr_coll_v(dh, e, lite_dh->name());
    auto ass_type = LiteDataType<U>();

    try{
      if(!ptr_coll_v.size()) return;
      const std::vector<art::Ptr<U> > ptr_coll = ptr_coll_v.at(0);
    }catch( art::Exception const& e){
      return;
    }
    // Instantiate association container. length = # of producers for associated data type
    std::vector< ::larlite::AssSet_t > ass_set_v(fModuleLabel_v[(size_t)(ass_type)].size());

    // Return if there's no data products stored for associated data type
    if(!(ass_set_v.size())) return;

    std::pair<size_t,size_t> lite_location;

    // Loop over origin data product vector, find associated objects and store association info
    for(size_t i=0; i<dh->size(); ++i) {

      const std::vector<art::Ptr<U> > ptr_coll = ptr_coll_v.at(i);

      // Association vector: one per associated data product producers
      std::vector<larlite::AssUnit_t> ass_unit_v(ass_set_v.size(),::larlite::AssUnit_t());

      for(auto& art_ptr : ptr_coll) {

	if(!LocateLiteProduct(art_ptr,lite_location)) continue;

	ass_unit_v[lite_location.second].push_back(lite_location.first);
      }
      for(size_t i=0; i<ass_set_v.size(); ++i)

	ass_set_v[i].push_back(ass_unit_v[i]);
      
    } // end looping over origin data products
	
    // Store associations in larlite data products
    for(size_t i=0; i<ass_set_v.size(); ++i) {

      // Loop over in one association set, store if there's any
      for(auto const& ass_unit : ass_set_v[i]) {

	if(ass_unit.size()) {

	  auto ass_name = fModuleLabel_v[(size_t)ass_type][i];

	  larlite::product_id ass_id(ass_type,ass_name);
	  
	  lite_dh->set_association(ass_id,ass_set_v[i]);

	  break;
	}
      }// end looping over association set
    }// end looping over a vector of association set
  }

  template <> void ScannerAlgo::ScanAssociation <::recob::Cluster,::recob::Cluster>(art::Event const& e,
										art::Handle<std::vector<::recob::Cluster> > &dh,
										::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::recob::EndPoint2D,::recob::EndPoint2D>(art::Event const& e,
										      art::Handle<std::vector<::recob::EndPoint2D> > &dh,
										      ::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::recob::Vertex,::recob::Vertex>(art::Event const& e,
										      art::Handle<std::vector<::recob::Vertex> > &dh,
										      ::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::anab::CosmicTag,::anab::CosmicTag>(art::Event const& e,
										      art::Handle<std::vector<::anab::CosmicTag> > &dh,
										      ::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::recob::SpacePoint,::recob::SpacePoint>(art::Event const& e,
										      art::Handle<std::vector<::recob::SpacePoint> > &dh,
										      ::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::recob::Track,::recob::Track>(art::Event const& e,
									    art::Handle<std::vector<::recob::Track> > &dh,
									    ::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::recob::Shower,::recob::Shower>(art::Event const& e,
									      art::Handle<std::vector<::recob::Shower> > &dh,
									      ::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::anab::Calorimetry,::anab::Calorimetry>(art::Event const& e,
										      art::Handle<std::vector<::anab::Calorimetry> > &dh,
										      ::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::anab::ParticleID,::anab::ParticleID>(art::Event const& e,
										    art::Handle<std::vector<::anab::ParticleID> > &dh,
										    ::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::simb::MCTruth,::simb::MCTruth>(art::Event const& e,
									      art::Handle<std::vector<::simb::MCTruth> > &dh,
									      ::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::simb::MCParticle,::simb::MCParticle>(art::Event const& e,
										    art::Handle<std::vector<::simb::MCParticle> > &dh,
										    ::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::recob::PFParticle,::recob::PFParticle>(art::Event const& e,
										      art::Handle<std::vector<::recob::PFParticle> > &dh,
										      ::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  //
  // LiteDataType
  //
  template <class T>
  const ::larlite::data::DataType_t ScannerAlgo::LiteDataType() const
  { throw cet::exception(__PRETTY_FUNCTION__)<<"Unsupported data type conversion attempted!";
    return ::larlite::data::kUndefined;
  }

  //
  // ProductID
  //
  template <class T>
  const ::larlite::product_id ScannerAlgo::ProductID(size_t name_index) const
  { auto data_type = LiteDataType<T>();
    if(fModuleLabel_v[(size_t)data_type].size() <= name_index)
      throw cet::exception(__PRETTY_FUNCTION__)
	<< "Length of registered products for data type " << ::larlite::data::kDATA_TREE_NAME[data_type].c_str()
	<< " is " << fModuleLabel_v[(size_t)data_type].size()
	<< " while you requested " << name_index;
    return ::larlite::product_id(data_type,fModuleLabel_v[(size_t)(data_type)][name_index]);
  }
}
#endif
