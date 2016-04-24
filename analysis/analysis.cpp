#include "ptcl_buffer.hpp"
#include "rg.hpp"
#include "ade.hpp"
#include "locstress.hpp"

constexpr char Parameter::atom_type[21];

TrjAnalysis<FPDPD, Particle>* trjana_factory(const std::string& mode_name,
					     const std::string& cur_dir,
					     const char* trj_fname,
					     const char* run_fname) {
  if (mode_name == "Rg") {
    return new TrjAnalysisRg(cur_dir, trj_fname, run_fname);
  } else if (mode_name == "ADE") {
    return new TrjAnalysisAde(cur_dir, trj_fname, run_fname);
  } else if (mode_name == "Lstress") {
    return new TrjAnalysisLocStress(cur_dir, trj_fname, run_fname);
  } else {
    std::cerr << "Unknown mode " << mode_name << std::endl;
    std::exit(1);
  }
}

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: ./analysis dir_name analysis_mode ... \n";
    std::cerr << "argv[1] is target directory name " << std::endl;
    std::exit(1);
  }
  
  const std::string cur_dir = argv[1];
  const std::string mode_name = argv[2];
  const std::string trj_fname = cur_dir + "/traject.xyz";
  const std::string param_fname = cur_dir + "/run_param.txt";
  
  std::unique_ptr<TrjAnalysis<FPDPD, Particle> > ptr_analysis(trjana_factory(mode_name, cur_dir, trj_fname.c_str(), param_fname.c_str()));
  ptr_analysis->Initialize();
  ptr_analysis->DoAnalysis();
}
