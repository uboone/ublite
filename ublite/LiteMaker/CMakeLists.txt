cet_make_library(
  SOURCE
  ScannerAlgo.cxx
  LIBRARIES
  PUBLIC
  ubcore::Geometry
  ubobj::Trigger
  ubobj::MuCS
  larsim::EventWeight_Base
  lardata::LArPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::MCBase
  nusimdata::SimulationBase
  art::Framework_Principal
  canvas::canvas
  messagefacility::MF_MessageLogger
  larlite::DataFormat
)

cet_build_plugin(
  LiteScanner art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ublite::LiteMaker
  ubcore::Geometry
  ubobj::Trigger
  ubobj::MuCS
  larsim::EventWeight_Base
  larevt::ChannelStatusService
  lardata::Utilities
  lardata::LArPropertiesService
  lardata::DetectorPropertiesService
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::MCBase
  lardataobj::AnalysisBase
  nusimdata::SimulationBase
)

cet_build_plugin(
  UtilScanner art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubcore::Geometry
  lardata::Utilities
  lardata::LArPropertiesService
  lardata::DetectorPropertiesService
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  larlite::DataFormat
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  LLMetaMaker art::service
  LIBRARIES
  PUBLIC
  ubevt::Utilities
)

add_subdirectory(job)

# install_headers()
install_fhicl()
install_source()
