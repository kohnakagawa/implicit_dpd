#pragma once

#include <cstring>

namespace io_util {
  inline FILE* xfopen (const char* __restrict filename,
		       const char* __restrict  mode)
  {
    FILE* f = fopen(filename, mode);
    if (f == nullptr) {
      std::cerr << filename << ": Cannot open file\n.";
      std::exit(1);
    }
    return f;
  }
  
  /*
    xyz file format.
    <number of atoms>
    comment line
    <element> <X> <Y> <Z>
    ...
  */

  template<class FP>
  inline void WriteXYZForm(const FP* ptcl, const PS::U32 num, const PS::U32 time, FILE* fp) {
    fprintf(fp, "%u\n", num);
    fprintf(fp, "time %u\n", time);
    for(PS::U32 i = 0; i < num; i++)
      ptcl[i].writeAscii(fp);
  }
  
  template<class FP>
  inline void WriteXYZForm(const PS::ParticleSystem<FP>& sys, const PS::U32 num, const PS::U32 time, FILE* fp) {
    WriteXYZForm(&sys[0], num, time, fp);
  }

  template<class FP>
  inline void WriteXYZFormMPI(const PS::ParticleSystem<FP>& sys,
			      const PS::U32 time, const PS::S32 buf_size,
			      FILE* fp) {
    const PS::S32 num_proc = PS::Comm::getNumberOfProc();
    static PS::ReallocatableArray<FP> ptcl_buf(buf_size);
    static PS::ReallocatableArray<PS::S32> n_ptcl(num_proc);
    static PS::ReallocatableArray<PS::S32> n_ptcl_disp(num_proc + 1);
    PS::S32 n_loc = sys.getNumberOfParticleLocal();
    PS::Comm::allGather(&n_loc, 1, n_ptcl.getPointer());
    n_ptcl_disp[0] = 0;
    for(PS::S32 i = 0; i < num_proc; i++)
      n_ptcl_disp[i + 1] = n_ptcl_disp[i] + n_ptcl[i];

    assert(buf_size >= n_ptcl_disp[num_proc]);
    
    PS::Comm::allGatherV(sys.getParticlePointer(), n_loc,
			 ptcl_buf.getPointer(), n_ptcl.getPointer(), n_ptcl_disp.getPointer());
    if(PS::Comm::getRank() == 0)
      WriteXYZForm(&(ptcl_buf[0]), n_ptcl_disp[num_proc], time, fp);
  }
  
  template<class FP>
  inline void ReadXYZForm(FP* ptcl,
			  PS::U32& num,
			  PS::U32& cur_time,
			  FILE* fp) {
    char comment[4] = {'\0'}, buf = '0';
    char cmp_comm[] = {"time"};
    fscanf(fp, "%u\n", &num);
    fscanf(fp, "%s %u\n", comment, &cur_time);
    if(strcmp(comment, cmp_comm) == 0) {
      std::cout << "Restart configuration is successfully loaded.\n";
    } else {
      std::cerr << "The file format of init_config.xyz may be old.\n";
      std::cerr << "xyz file format\n";
      std::cerr << "<number of atoms>\n";
      std::cerr << "time %lf\n";
      std::cerr << "<element> <X> <Y> <Z> other mol inform\n";
      PS::Abort();
    }
    
    for(PS::U32 i = 0; i < num; i++)
      ptcl[i].readAscii(fp);
    fscanf(fp, "%c\n", &buf);
    if(!feof(fp)) {
      std::cerr << "# of lines is not equal to the information specified in file header.\n";
      std::exit(1);
    }
  }

  template<class FP>
  inline void ReadXYZForm(PS::ParticleSystem<FP>& sys, PS::U32& num, PS::U32& cur_time, FILE* fp) {
    ReadXYZForm(&sys[0], num, cur_time, fp);
  }
};
