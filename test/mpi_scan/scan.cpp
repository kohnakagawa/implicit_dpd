#include <iostream>
#include "particle_simulator.hpp"

int main(int argc, char* argv[]) {
  PS::Initialize(argc, argv);

  const PS::S32 rank = PS::Comm::getRank() + 3;
  PS::S32 val = 0;
  PS::Comm::exScan(&rank, &val, 1);
  std::cout << rank << " " << val << std::endl;

  PS::Finalize();
}
