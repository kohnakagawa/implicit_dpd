#include "ptcl_buffer.hpp"
#include "rg.hpp"
#include "ade.hpp"
#include "locstress.hpp"
#include "tube_rad.hpp"
#include "density_map.hpp"
#include "Reo.hpp"
#include "thickness.hpp"

// static members in Parameter
constexpr char Parameter::atom_type[21];
PS::F64vec Parameter::box_leng, Parameter::ibox_leng;
PS::U32 Parameter::time, Parameter::all_time, Parameter::step_mic, Parameter::step_mac;

namespace {
  const std::string suppording_mode = "Rg/Ade/LocStress/Rtube/Density/Reo/Thickness";
  
  TrjAnalysis<FPDPD, Particle>* trjana_factory(const std::string& mode_name,
					       const std::string& cur_dir,
					       const char* trj_fname,
					       const char* run_fname) {
    if (mode_name == "Rg") {
      return new TrjAnalysisRg(cur_dir, trj_fname, run_fname);
    } else if (mode_name == "Ade") {
      return new TrjAnalysisAde(cur_dir, trj_fname, run_fname);
    } else if (mode_name == "LocStress") {
      return new TrjAnalysisLocStress(cur_dir, trj_fname, run_fname);
    } else if (mode_name == "Rtube") {
      return new TrjAnalysisTubeRad(cur_dir, trj_fname, run_fname);
    } else if (mode_name == "Density") {
      return new TrjAnalysisDensity(cur_dir, trj_fname, run_fname);
    } else if (mode_name == "Reo") {
      return new TrjAnalysisReo(cur_dir, trj_fname, run_fname);
    } else if (mode_name == "Thickness") {
      return new TrjAnalysisThickness(cur_dir, trj_fname, run_fname);
    } else {
      std::cerr << "Unknown mode " << mode_name << std::endl;
      std::cerr << "Only suport " << suppording_mode << std::endl;
      std::exit(1);
    }
  }
}


int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: ./post_ana dir_name analysis_mode \n";
    std::cerr << "argv[1] is target directory name (should include traject.xyz and run_param.txt) \n";
    std::cerr << "argv[2] is execution mode " << suppording_mode << "." << std::endl;
    std::exit(1);
  }
  
  const std::string cur_dir   = argv[1], mode_name = argv[2];
  const std::string trj_fname = "/traject.xyz", param_fname = "/run_param.txt";
  
  std::unique_ptr<TrjAnalysis<FPDPD, Particle> > ptr_analysis(trjana_factory(mode_name, cur_dir, trj_fname.c_str(), param_fname.c_str()));
  ptr_analysis->Initialize();
  ptr_analysis->DoAnalysis();
}
