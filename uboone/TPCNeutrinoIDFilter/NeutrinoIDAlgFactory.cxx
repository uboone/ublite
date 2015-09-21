/**
 *  @file   NeutrinoIDAlgFactory.cxx
 * 
 *  @brief  A small module to instantiate various algorithms inheriting from NeutrinoIDAlgBase
 * 
 */

// Includes for algorithms
#include "TPCNeutrinoIDFilter/NeutrinoIDAlgFactory.h"
#include "TPCNeutrinoIDFilter/TrackPairPlusVertexAlg.h"


//------------------------------------------------------------------------------------------------------------------------------------------

namespace neutrinoid
{
    
std::unique_ptr< NeutrinoIDAlgBase > NeutrinoIDAlgFactory::MakeNeutrinoIDAlg(fhicl::ParameterSet const& p)
{
    std::string algName = p.get<std::string>("NeutrinoIDAlgName");
    
    std::unique_ptr< NeutrinoIDAlgBase > ptr;
    
    if(algName.compare("TrackPairPlusVertexAlg")==0)
    {
        std::unique_ptr< NeutrinoIDAlgBase > new_ptr(new TrackPairPlusVertexAlg(p));
        ptr.swap(new_ptr);
    }
    else{
        std::cout << "Algname is ... " << algName << std::endl;
        throw std::runtime_error("ERROR in NeutrinoIDAlgFactory: No registered Neutrino ID with that name.");
    }
    
    return std::move(ptr);
}

} // namespace neutrinoid
