/* automatic test */

#include "f_calculator.hpp"
#include <string>

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#error "MPI is not supported yet!"
#endif

namespace {
  template<class Tpsys>
  void CheckNearlist(Tpsys& sys) {
    const int num = getNumberOfParticleLocal();
    for(int i = 0; i < num; i++) {
      
    }
  }

  void CheckBondedlist() {
  

  }

  void CheckSorted() {
  

  }

  /*
    print out kinetic and configurational temperature.
  */
  void CheckTemperature() {
  
  }

  /*
    check Groot & Warren result
  */
  void CheckCompressibility() {
  
  }

  /*
    check Noguchi result
  */
  void CheckDiffusionConstant() {
  
  }

}

int main(int argc, char* argv[]) {
  PS::Initialize(argc, argv);
  const std::string mode(argv[1]);

  switch(mode) {
  case "nearlist":
    
    break;
  case "bondedlist":
    
    break;
  case "io":
    
    break;
  case "driftkick":

    break;
    
  case "observer":
    
    break;
  case "testall":
    std::cerr << "check all content\n";

    break;
  default:
    std::cerr << "I don't know such mode\n";
    std::cerr << "mode name: " << mode << std::endl;
  }
  
  PS::Finalize();
}
