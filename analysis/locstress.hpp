#pragma once
#include "trjanalysis.hpp"
#include "adjust_axis.hpp"
// #include "mdstresslib/include/mds_stressgrid.h"
#include "parameter.hpp"
#include "io_util.hpp"

class TrjAnalysisLocStress : public TrjAnalysis<FPDPD, Particle> {
  void DebugDump(std::vector<Particle>& ptcls, const size_t time) {
    const std::string name_f = cur_dir_ + "/transformed.xyz";
    FILE* fp = io_util::xfopen(name_f.c_str(), "a");
    fprintf(fp, "%u\n", ptcls.size());
    fprintf(fp, "time %u\n", time);
    for (size_t i = 0; i < ptcls.size(); i++) {
      const auto prop = ptcls[i].prop;
      fprintf(fp, "%c %.15g %.15g %.15g\n",
	      Parameter::atom_type[prop], ptcls[i].r.x, ptcls[i].r.y, ptcls[i].r.z);
    }
    fclose(fp);
  }
  
public:
  TrjAnalysisLocStress(const std::string cur_dir,
		       const char* trj_fname,
		       const char* param_fname) : TrjAnalysis(cur_dir, trj_fname, param_fname) {}
  ~TrjAnalysisLocStress() override {}

  void DoAnalysis() override {
    const std::string name_fout = cur_dir_ + "/local_stress.txt";
    
    size_t time = 0;
    AxisAdjuster<Particle> aadjuster;
    while (true) {
      if (ptr_connector->ReadOneFrame(fin_trj, time)) break;
      std::cout << "time = " << time << std::endl;
      const auto num_patch = gen_connected_image(ptr_connector.get(), false);
      std::cout << "# of patch is " << num_patch << std::endl;
      
      // change base axis
      aadjuster.CreateMomInertia(ptr_connector->ptcls());
      aadjuster.DoTransform(ptr_connector->ptcls());
      
      DebugDump(ptr_connector->ptcls(), time);
    }
  }
};
