set( LARLITE_LIB -L$ENV{LARLITE_LIBDIR} -lLArLite_Base -lLArLite_DataFormat -lLArLite_LArUtil)

include_directories ( $ENV{LARLITE_COREDIR} )

link_libraries( ${LARLITE_LIB}) # All execs need this.

add_subdirectory(job)

art_make( EXCLUDE LLMetaMaker_service.cc
	  BASENAME_ONLY
	  LIBRARY_NAME     uboonecode_uboone_LiteMaker
	  LIB_LIBRARIES
		   ubraw_RawData
	  	   larcorealg_Geometry
		   ubcore_Geometry
		   ubobj_MuCS
		   ubobj_Trigger
		   larcore_Geometry_Geometry_service
	   	   larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
		   larsim_MCCheater_ParticleInventoryService_service
		   lardataobj_RecoBase
		   nusimdata::SimulationBase
		   lardataobj_MCBase
		   larcoreobj_SummaryData
		   lardataobj_RawData
		   lardata_RecoObjects
		   nusimdata::SimulationBase
		   art::Framework_Core
		   art::Framework_Principal

		   art::Framework_Services_Registry
		   art_root_io::tfile_support ROOT::Core
		   art_root_io::TFileService_service
		   art::Persistency_Common canvas::canvas
		   art::Persistency_Provenance canvas::canvas
		   art::Utilities canvas::canvas
		   messagefacility::MF_MessageLogger

		   fhiclcpp::fhiclcpp
		   cetlib::cetlib cetlib_except
		   ROOT::Geom
		   ROOT::XMLIO
		   ROOT::Gdml


          MODULE_LIBRARIES
		   uboonecode_uboone_LiteMaker
		   ubraw_RawData
	      	   ubevt::TFileMetadataMicroBooNE_service
	      	   ubevt::FileCatalogMetadataMicroBooNE_service
		   LLMetaMaker_service
		   ubobj_MuCS
		   ubobj_Trigger
                   lardataalg_DetectorInfo
		   larcorealg_Geometry
		   ubcore_Geometry
		   larcore_Geometry_Geometry_service
		   larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
		   nusimdata::SimulationBase
		   lardata_Utilities
		   lardataobj_MCBase
		   larcoreobj_SummaryData
		   lardataobj_RawData
		   lardataobj_RecoBase
		   lardata_RecoObjects
                   larevt_Filters
                   larevt_CalibrationDBI_IOVData
                   larevt_CalibrationDBI_Providers
		   nusimdata::SimulationBase
		   art::Framework_Core
		   art::Framework_Principal

		   art::Framework_Services_Registry
		   art_root_io::tfile_support ROOT::Core
		   art_root_io::TFileService_service
		   art::Persistency_Common canvas::canvas
		   art::Persistency_Provenance canvas::canvas
		   art::Utilities canvas::canvas
		   messagefacility::MF_MessageLogger

		   fhiclcpp::fhiclcpp
		   cetlib::cetlib cetlib_except
		   ROOT::Geom
		   ROOT::XMLIO
		   ROOT::Gdml

		)

simple_plugin(LLMetaMaker "service"
	      ubevt::TFileMetadataMicroBooNE_service
	      ubevt::FileCatalogMetadataMicroBooNE_service
              art::Framework_Principal
	      art::Utilities canvas::canvas
	      messagefacility::MF_MessageLogger

	      fhiclcpp::fhiclcpp
	      cetlib::cetlib cetlib_except

              BASENAME_ONLY
)

# install_headers()
install_fhicl()
install_source()
