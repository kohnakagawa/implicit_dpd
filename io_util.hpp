#pragma once

namespace io_util {
  inline FILE* xfopen (const char* __restrict filename,
		       const char* __restrict  mode)
  {
    FILE* f = fopen(filename, mode);
    if (f == NULL) {
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
  inline void WriteXYZForm(const FP* ptcl, const PS::U32 num, FILE* fp) {
    fprintf(fp, "%u\n", num);
    fprintf(fp, "MDPD configuration.\n");
    for(PS::U32 i = 0; i < num; i++)
      ptcl[i].writeAscii(fp);
  }
  
  template<class FP>
  inline void WriteXYZForm(const PS::ParticleSystem<FP>& sys, const PS::U32 num, FILE* fp) {
    WriteXYZForm(&sys[0], num, fp);
  }
  
  template<class FP>
  inline void ReadXYZForm(FP* ptcl, PS::U32& num, FILE* fp) {
    char comment[256] = {'\0'}, buf = '0';
    fscanf(fp, "%u\n", &num);
    fgets(comment, 256, fp);
    for(PS::U32 i = 0; i < num; i++)
      ptcl[i].readAscii(fp);
    fscanf(fp, "%c\n", &buf);
    if(!feof(fp)) {
      std::cerr << "# of lines is not equal to the information specified in file header.\n";
      std::exit(1);
    }
  }

  template<class FP>
  inline void ReadXYZForm(PS::ParticleSystem<FP>& sys, PS::U32& num, FILE* fp) {
    ReadXYZForm(&sys[0], num, fp);
  }
};
