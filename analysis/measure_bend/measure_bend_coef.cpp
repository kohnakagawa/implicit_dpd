#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <complex.h>
#include <fftw3.h>
#include <parallel/algorithm>

class BendCalculator final {
  struct set_XY {
    double q[2];
    double normq, coefbend, normhq;
  };

  double* hei = nullptr;
  fftw_complex* hei_q = nullptr;
  set_XY* set_qhei_q = nullptr, *set_qhei_qmean = nullptr;
  std::vector<int>	 hei_elem;
  int Nx, Ny, Nxy;
  double grid_leng[2], scL[2];
  fftw_plan	plan;

  std::ofstream fspect, fbendcf, foutmean;
  std::ifstream fin;
  const std::string cur_dir;

  int GenHash(const double* a) {
    int indx = static_cast<int>(a[0] / grid_leng[0]);
    int indy = static_cast<int>(a[1] / grid_leng[1]);
    if (indx == Nx) indx--;
    if (indy == Ny) indy--;
    const int hash = Ny * indx + indy;
    assert(hash >= 0 && hash < Nxy);
    return hash;
  }

public:
  explicit BendCalculator(const std::string cur_dir_) : cur_dir(cur_dir_) {}

  ~BendCalculator() {
    fftw_destroy_plan(plan);
    delete [] hei;
    fftw_free(hei_q);
    delete [] set_qhei_q;
    delete [] set_qhei_qmean;
  }

  void FileOpen() {
    const std::string fname = cur_dir + "/memb_height.txt";
    fin.open(fname.c_str());

    std::string line[6];
    for (int i = 0; i < 6; i++) {
      fin >> line[i];
    }
    Nx = std::stoi(line[1]);
    Ny = std::stoi(line[2]);
    Nxy = Nx * Ny;
    grid_leng[0] = std::stof(line[4]);
    grid_leng[1] = std::stof(line[5]);

    scL[0] = grid_leng[0] * Nx;
    scL[1] = grid_leng[1] * Ny;

    const std::string fout_name[] = {
      cur_dir + "/spectrum.txt",
      cur_dir + "/meanspect.txt",
      cur_dir + "/bendcoef.txt",
    };

    fspect.open(fout_name[0].c_str());
    foutmean.open(fout_name[1].c_str());
    fbendcf.open(fout_name[2].c_str());

    hei_elem.resize(Nxy, 0.0);
  }

  void LoadHeight() {
    double grid_pos[2], loc_hei;
    for (int i = 0; i < Nxy; i++) {
      fin >> grid_pos[0] >> grid_pos[1] >> loc_hei;
      hei[GenHash(grid_pos)] = loc_hei;
    }
  }

  bool IsEof() { return fin.eof(); }

  void SortbyNormq(const bool is_count) {
    for (int ix = 0; ix < Nx; ix++) {
      for (int iy = 0; iy < (Ny / 2 + 1); iy++) {
        const int         hash     = iy + (Ny / 2 + 1) * ix;
        const double	p_real	   = hei_q[hash][0]; //real part
        const double	p_imgn	   = hei_q[hash][1]; //imaginary part
        const double	norm_hei_q = p_real*p_real + p_imgn*p_imgn;
        const double	qx	   = (ix < Nx/2 + 1) ? 2.0 * M_PI * ix / (Nx * grid_leng[0]) : 2.0 * M_PI * (ix - Nx) / (Nx * grid_leng[0]);
        const double	qy	   = 2.0 * M_PI * iy / (Nx * grid_leng[1]);
        const double	norm_q	   = std::sqrt(qx * qx + qy * qy);
        const double      norm_q_2   = norm_q * norm_q;
        const double      norm_q_4   = norm_q_2 * norm_q_2;
        set_qhei_q[hash].q[0]        = qx;
        set_qhei_q[hash].q[1]        = qy;

        set_qhei_q[hash].normq	   = norm_q;
        set_qhei_q[hash].coefbend  = 1.0 / (norm_q_4 * norm_hei_q);
        set_qhei_q[hash].normhq	   = norm_hei_q;

        if (is_count) {
          set_qhei_qmean[hash].coefbend += 1.0 / (norm_q_4 * norm_hei_q);
          set_qhei_qmean[hash].normhq   += norm_hei_q;
        }
      }
    }
    __gnu_parallel::sort(set_qhei_q,
                         set_qhei_q + (Ny / 2 + 1 ) * Nx,
                         [](const set_XY& riL, const set_XY& riR) {
                           return riL.normq < riR.normq;
                         });
  }

  void Setfftw2d() { plan = fftw_plan_dft_r2c_2d(Nx, Ny, hei, hei_q, FFTW_ESTIMATE); }

  void Exefftw2d() { fftw_execute(plan); }

  void Normalize() {
    const double sqLxy = std::sqrt(scL[0] * scL[1]), invNxy = 1.0 / Nxy, coef = sqLxy * invNxy;
    for(int i = 0; i < Nx * ( Ny / 2 + 1 ); i++) {
      hei_q[i][0] *= coef;
      hei_q[i][1] *= coef;
    }
  }

  void MemAlloc() {
    hei            = new double [Nxy];
    hei_q          = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Ny / 2 + 1) * Nx);
    set_qhei_q		 = new set_XY [(Ny / 2 + 1) * Nx];
    set_qhei_qmean = new set_XY [(Ny / 2 + 1) * Nx];
  }

  void Initialize() {
    for(int ix = 0; ix < Nx; ix++) {
      for(int iy = 0; iy < Ny / 2 + 1; iy++) {
        const int       hash       = iy + (Ny / 2 + 1) * ix;
        const double	qx	   = (ix < Nx / 2 + 1) ? 2. * M_PI * ix / (Nx * grid_leng[0]) : 2.0 * M_PI * (ix - Nx) / (Nx * grid_leng[0]);
        const double	qy	   = 2.0 * M_PI * iy / (Nx * grid_leng[1]);
        const double	norm_q	   = std::sqrt(qx * qx + qy * qy);
        set_qhei_q[hash].q[0]        = qx;
        set_qhei_q[hash].q[1]        = qy;
        set_qhei_q[hash].normq = norm_q;

        set_qhei_qmean[hash].q[0]     = qx;
        set_qhei_qmean[hash].q[1]     = qy;
        set_qhei_qmean[hash].normq    = norm_q;
        set_qhei_qmean[hash].coefbend = 0.;
        set_qhei_qmean[hash].normhq   = 0.;
      }
    }
    __gnu_parallel::sort(set_qhei_q,
                         set_qhei_q + (Ny/2 + 1) * Nx,
                         [](const set_XY& riL, const set_XY& riR) {
                           return riL.normq < riR.normq;
                         });
  }

  void WriteQnorm() {
    for (int i = 0; i < Nx * (Ny / 2 + 1); i++) {
      if (set_qhei_q[i].normq != 0.0) {
        fspect << set_qhei_q[i].normq << " ";
        fbendcf << set_qhei_q[i].normq << " ";
      }
    }
    fspect << std::endl;
    fbendcf << std::endl;
  }

  void AllClear() {
    for(int i = 0; i < Nxy; i++) {
      hei[i] = 0.0;
      hei_elem[i] = 0;
    }
  }

  void WriteResult() {
    for (int i = 0; i < Nx * (Ny / 2 + 1); i++) {
      if (set_qhei_q[i].normq != 0.0) {
        fspect << set_qhei_q[i].normhq << " ";
        fbendcf << set_qhei_q[i].coefbend << " ";
      }
    }
    fspect << std::endl;
    fbendcf << std::endl;
  }

  void TakeMeanTime(const int beg_line, const int all_line) {
    int count = all_line - beg_line;
     for(int i = 0; i < Nx * (Ny / 2 + 1); i++) {
       set_qhei_qmean[i].coefbend /= count;
       set_qhei_qmean[i].normhq /= count;
     }
  }

  void WriteResultMean() {
    for (int i = 0; i < Nx * (Ny / 2 + 1); i++) {
      if (set_qhei_q[i].normq != 0) {
        foutmean << set_qhei_qmean[i].normq    << " "
                 << set_qhei_qmean[i].coefbend << " "
                 << set_qhei_qmean[i].normhq   << std::endl;
      }
    }
    std::cout << "spectrum element number " << Nx * (Ny / 2 + 1) << std::endl;
  }
};

int main(int argc, char *argv[]) {
  if (argc != 3) {
    std::cout << "usage: ";
    std::cout << argv[0] << " data_dir beg_time_point" << std::endl;
    std::exit(1);
  }

  const std::string cur_dir = argv[1];
  const int skip_line = std::atoi(argv[2]);
  BendCalculator bendcalc(cur_dir);

  bendcalc.FileOpen();
  bendcalc.MemAlloc();
  bendcalc.Initialize();
  bendcalc.WriteQnorm();
  bendcalc.AllClear();
  bendcalc.Setfftw2d();

  int cnt = 0;
  while (1) {
    bendcalc.AllClear();
    bendcalc.LoadHeight();
    if (bendcalc.IsEof()) break;
    bendcalc.Exefftw2d();
    bendcalc.Normalize();
    bendcalc.SortbyNormq(cnt >= skip_line);
    if (cnt >= skip_line) bendcalc.WriteResult();
    cnt++;
  }

  bendcalc.TakeMeanTime(skip_line, cnt);
  bendcalc.WriteResultMean();

  std::cout << "./output_dir/spectrum.txt  -> spectrum time evolution" << std::endl;
  std::cout << "./output_dir/meanspect.txt -> mean value of spectrum " << std::endl;
  std::cout << "./output_dir/bendcoef.txt  -> bend coef              " << std::endl;
}
