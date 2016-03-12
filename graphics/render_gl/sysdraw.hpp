#pragma once
#include <vector>
#include <bitset>
#include <fstream>
#include <array>
#include <GL/freeglut.h>
#include <memory>
#include <boost/algorithm/string.hpp>

class Jpegout;

class DrawSys {
  enum {
    SEED_N = 4,
  };

  template<int val, int n>
  struct power{
    static const int ret = val * power<val, n - 1>::ret;
  };
  template<int val>
  struct power<val, 0>{
    static const int ret = val;
  };

  struct particle{
    float r[3];
    int prop;
    float rad;
    int idx;
    
    void readFromXYZForm(std::ifstream& fin) {
#if 0
      std::string line;
      std::getline(fin, line);
      fin >> line;
      std::vector<std::string> v;
      boost::algorithm::split(v, line, boost::algorithm::is_space());
#else
      std::vector<std::string> v(17);
      for (int i = 0; i < 17; i++)
	fin >> v[i];
#endif
      
      r[0] = std::stof(v[1]);
      r[1] = std::stof(v[2]);
      r[2] = std::stof(v[3]);
      idx  = std::stoi(v[4]);
      prop = std::stoi(v[5]);
      if (idx == 10)
	rad = 0.01;
      else
	rad = 0.001;
    }
  };
  
  std::vector<std::array<GLfloat, 3> > p_color; 
  std::vector<std::array<GLfloat, 4> > lightpos;
  
  std::bitset<15>  draw_crit;
  std::vector<int> draw_crit_mask;
  int draw_crit_max, draw_crit_base;

  std::array<std::array<int, 2>, 12>  cubeedge;
  std::array<GLfloat, 3> nv;
  GLfloat cut_plane = 0.5, vertex[8][3];
  
  std::string cur_dir;
  
  float box_size[3] = {0.0, 0.0, 0.0}, inv_box_size[3] = {0.0, 0.0, 0.0}, invL = 0.0;
  float *p_fovy = nullptr, *p_perscenter = nullptr, *p_center2eye = nullptr, *p_base_z = nullptr;
  
  void *font = GLUT_BITMAP_TIMES_ROMAN_24;
  
  void RenderString2D(const char*, float, float);
protected:
  int cur_time = 0, time_step = 0, all_time = 0;
  bool swt_but = true, cut_but = false, cut_adv = false, crit_out = false, contined = true;
  int pN = 0;
  std::vector<particle> Particle;
  std::ifstream fin;
  
  std::unique_ptr<Jpegout> jpgout;
  
  void Drawxyz();
  void DrawCubic();
  void DrawAxis(float   , float, const float[][3]);
  void DrawSubAxis(float, float, const float[3]);
  
  void RenderCurTime();

  void RenderString3D(const char *,const float[3]);

  void RenderSphere(const particle& prtcl);
  void RenderPoint(const particle& prtcl);
  bool IsDrawnObject(const particle& prtcl);
  
  void ChangeNormalVector(int i);
  void ChangeCrossSection();
  
  void Dump2Jpg();
  void ChangeLookDirection(const int i);
  
public:
  DrawSys(const std::string& cur_dir_, const bool crit_out_, const int b_time);
  virtual ~DrawSys();
  void SetParams();
  
  void SetCallBackFunc() const;
  void SetColor(const GLfloat*);
  void SetLightPos(const GLfloat*);
  void GetMouseInfo(float*, float*, float*, float*);

  void InitCube();
  void InitWindow(int argc, char* argv[], const char* win_name) const;
  void InitColor() const;
  bool InitGlew() const;
  
  void FileOpen();
  
  particle ParseDataLine(const std::string& line) const ;
  void LoadParticleDat();
  
  virtual void Timer(int) = 0;
  virtual void Display() = 0;
  virtual void KeyBoard(unsigned char, int ,int) = 0;
  virtual void PrintDrawInfo() = 0;
  
  void Resize(int,int) const;
  void ChgDrawObject();
  
  void Execute(){
    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);
    glutMainLoop();
  }
  
  int Pow_n(int, int) const;
};

class AnimeDraw : public DrawSys{
public:
  AnimeDraw(const std::string& cur_dir_, const bool crit_out_, const int b_time) : DrawSys(cur_dir_, crit_out_, b_time){}
  ~AnimeDraw(){}
  void Timer(int);
  void Display();
  void KeyBoard(unsigned char, int, int);
  void PrintDrawInfo();
};

namespace callbacks{
  extern std::unique_ptr<DrawSys> drawsys;
  extern void wrap_display();
  extern void wrap_timer(int value);
  extern void wrap_resize(int w, int h);
  extern void wrap_keyboard(unsigned char key, int x, int y);
}
