#include <iostream>
#include <cstdlib>
#include <sstream>
#include <GL/glew.h>
#include <limits>
#include "sysdraw.hpp"
#include "jpegout.hpp"
#include "mousehandle.hpp"

void callbacks::wrap_display() {
  callbacks::drawsys->Display();
}
void callbacks::wrap_timer(int value) {
  callbacks::drawsys->Timer(value);
}
void callbacks::wrap_resize(int w, int h) {
  callbacks::drawsys->Resize(w, h);
}
void callbacks::wrap_keyboard(unsigned char key, int x, int y) {
  callbacks::drawsys->KeyBoard(key, x, y);
}

DrawSys::DrawSys(const std::string& cur_dir_, 
		 const bool crit_out_, 
		 const int b_time) : jpgout(new Jpegout()) {
  draw_crit_max  = power<2, SEED_N>::ret;
  draw_crit_base = draw_crit_max - 1;
  
  nv[0] = 1.0; nv[1] = 0.0, nv[2] = 0.0;
  
  cur_dir  = cur_dir_;
  crit_out = crit_out_;
}

DrawSys::~DrawSys() {}

void DrawSys::SetParams() {
  static_assert(SEED_N <= 15, "the number of max particle seed is 15");
  p_color.reserve(SEED_N);
  
  draw_crit.set();
  
  for(int i = 0; i < SEED_N; i++) {
    const int base_b = Pow_n(2, i);
    draw_crit_mask.push_back(base_b);
  }
}

void DrawSys::SetCallBackFunc() const {
  glutDisplayFunc(callbacks::wrap_display); 
  glutTimerFunc(50, callbacks::wrap_timer, 0);
  glutReshapeFunc(callbacks::wrap_resize);
  glutMouseFunc(callbacks::wrap_mclick);
  glutMotionFunc(callbacks::wrap_mmotion);
  glutMouseWheelFunc(callbacks::wrap_mwheel);
  glutKeyboardFunc(callbacks::wrap_keyboard);
}

void DrawSys::SetColor(const GLfloat* col) {
  const std::array<GLfloat, 3> buf = {col[0], col[1], col[2]};
  p_color.push_back(buf);
}

void DrawSys::SetLightPos(const GLfloat* pos) {
  const std::array<GLfloat ,4> buf = {pos[0], pos[1], pos[2], pos[3]};
  lightpos.push_back(buf);
}

void DrawSys::GetMouseInfo(float* fov_, float* perscent, float* cent2eye, float* base_z) {
  p_fovy       = fov_;
  p_perscenter = perscent;
  p_center2eye = cent2eye;
  p_base_z     = base_z;
};

void DrawSys::InitWindow(int argc, char* argv[], const char* win_name) const {
  glutInitWindowPosition(800, 400);
  glutInitWindowSize(1300, 1100);
  glutInit(&argc, argv); 
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
  glutCreateWindow(win_name);
}

void DrawSys::InitCube() {
  //down square
  vertex[0][0] = -0.5; vertex[0][1] = -0.5; vertex[0][2] = -0.5;
  vertex[1][0] = 0.5;  vertex[1][1] = -0.5; vertex[1][2] = -0.5;
  vertex[2][0] = 0.5;  vertex[2][1] = 0.5;  vertex[2][2] = -0.5;
  vertex[3][0] = -0.5; vertex[3][1] = 0.5;  vertex[3][2] = -0.5;

  //up square
  vertex[4][0] = -0.5; vertex[4][1] = -0.5; vertex[4][2] = 0.5;
  vertex[5][0] = 0.5;  vertex[5][1] = -0.5; vertex[5][2] = 0.5;
  vertex[6][0] = 0.5;  vertex[6][1] = 0.5;  vertex[6][2] = 0.5;
  vertex[7][0] = -0.5; vertex[7][1] = 0.5;  vertex[7][2] = 0.5;
  
  //edge connection
  cubeedge[0][0] = 0; cubeedge[0][1] = 1;
  cubeedge[1][0] = 1; cubeedge[1][1] = 2;
  cubeedge[2][0] = 2; cubeedge[2][1] = 3;
  cubeedge[3][0] = 3; cubeedge[3][1] = 0;
  cubeedge[4][0] = 4; cubeedge[4][1] = 5;
  cubeedge[5][0] = 5; cubeedge[5][1] = 6;
  cubeedge[6][0] = 6; cubeedge[6][1] = 7;
  cubeedge[7][0] = 7; cubeedge[7][1] = 4;
  cubeedge[8][0] = 0; cubeedge[8][1] = 4;
  cubeedge[9][0] = 1; cubeedge[9][1] = 5;
  cubeedge[10][0] = 2; cubeedge[10][1] = 6;
  cubeedge[11][0] = 3; cubeedge[11][1] = 7;
}

void DrawSys::InitColor() const {
  glClearColor(1.0, 1.0, 1.0, 1.0); //(Red,Green,Blue,A) Aは透明度 ウィンドウを塗りつぶす色を指定
  glEnable(GL_DEPTH_TEST); //デプスバッファ使用

  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);

  glEnable(GL_LIGHTING); //光源使用
  glEnable(GL_LIGHT0); //光源0を設定
}

#define CHECK_FILE_OPEN(fin, fname)					\
  do {									\
    if (!fin) {								\
      std::cerr << "File I/O error!" << std::endl;			\
      std::cerr << fname << " No such file" << std::endl;		\
      std::cerr << __FILE__ << " " << __LINE__ << std::endl;		\
      std::exit(1);							\
    }									\
  } while(false)							\

void DrawSys::FileOpen() {
  std::string str = cur_dir + "/traject.xyz";
  fin.open(str.c_str());
  CHECK_FILE_OPEN(fin, str);
  
  str = cur_dir + "/run_param.txt";
  std::ifstream fin_run(str.c_str());
  if (fin_run) {
    while (true) {
      std::string line;
      std::getline(fin_run, line);
      if (line.empty()) break;
      std::vector<std::string> v;
      boost::algorithm::split(v, line, boost::algorithm::is_space());
      if (v[0] == "core_amp_id") {
	core_amp_id.resize(v.size() - 1, 0);
	for (size_t i = 0; i < core_amp_id.size(); i++) {
	  core_amp_id[i] = std::stoi(v[i + 1]);
	}
      }
      if (v[0] == "influ_rad") {
	influ_rad = std::stof(v[1]);
      }
    }
  }
}

#undef CHECK_FILE_OPEN

void DrawSys::LoadParticleDat() {
  static bool is_first_call = true;
  
  // xyz format
  fin >> pN;
  std::string comment;
  fin >> comment >> cur_time;

  if (fin.eof()) {
    contined = false;
    return;
  }

  if (comment != "time") {
    std::cerr << "XYZ format may be old.\n";
    std::exit(1);
  }

  Particle.resize(pN);
  for (int i = 0; i < pN; i++) {
    Particle[i].readFromXYZForm(fin, core_amp_id, influ_rad);
    
    if (box_size[0] < Particle[i].r[0]) box_size[0] = Particle[i].r[0];
    if (box_size[1] < Particle[i].r[1]) box_size[1] = Particle[i].r[1];
    if (box_size[2] < Particle[i].r[2]) box_size[2] = Particle[i].r[2];
  }

  if (is_first_call) {
    const auto min_len = *std::min_element(box_size, box_size + 3);
    invL = 1.0 / min_len;
    pN_init = pN;
    is_first_call = false;
  }
  
  box_size[0] *= invL; box_size[1] *= invL; box_size[2] *= invL;
  inv_box_size[0] = 1.0 / box_size[0]; inv_box_size[1] = 1.0 / box_size[1]; inv_box_size[2] = 1.0 / box_size[2];
  
  for (int i = 0; i < pN; i++) {
    Particle[i].r[0] = Particle[i].r[0] * invL - box_size[0] * 0.5;
    Particle[i].r[1] = Particle[i].r[1] * invL - box_size[1] * 0.5;
    Particle[i].r[2] = Particle[i].r[2] * invL - box_size[2] * 0.5;

    Particle[i].rad *= invL;
  }
}

void DrawSys::Drawxyz() {
  const float axis_pos[3] = {0.0, 0.0, 0.0};
  glPushMatrix();
  glTranslatef(axis_pos[0], axis_pos[1], axis_pos[2]);
  glBegin(GL_LINES);

  glColor3f(1.0, 0.0, 0.0);
  glVertex3f(1.5, 0.0, 0.0);
  glVertex3f(0.0, 0.0, 0.0);

  glColor3f(0.0, 1.0, 0.0);
  glVertex3f(0.0, 1.5, 0.0);
  glVertex3f(0.0, 0.0, 0.0);

  glColor3f(0.0 , 0.0, 1.0);
  glVertex3f(0.0, 0.0, 1.5);
  glVertex3f(0.0, 0.0, 0.0);
  glEnd();
  glPopMatrix();
}

void DrawSys::DrawCubic() {
  const GLfloat color[3] = {0.0, 0.0, 0.0};
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);
  glBegin(GL_LINES);
  for (int i = 0; i < 12; ++i) {
    glVertex3fv(vertex[cubeedge[i][0]]);
    glVertex3fv(vertex[cubeedge[i][1]]);
  }
  glEnd();
}

void DrawSys::DrawAxis(float d, float s,const float col[][3]) {
  const float origin[3] = {-1.0, 0.0, 0.0};

  glPushMatrix();
  glTranslatef(origin[0], origin[1], origin[2]);
  DrawSubAxis(d,s,col[0]);
  glPopMatrix();

  glPushMatrix();
  glTranslatef(origin[0], origin[1], origin[2]);
  glRotatef(90.0, 1.0, 0.0, 0.0);
  DrawSubAxis(d,s,col[1]);
  glPopMatrix();

  glPushMatrix();
  glTranslatef(origin[0], origin[1], origin[2]);
  glRotatef(-90.0, 0.0, 0.0, 1.0);
  DrawSubAxis(d,s,col[2]);
  glPopMatrix();

  //Render Axis name
  /*const float n_pos[3][3] = {{1.23,0.0,0.0},{0.0,1.23,0.0},{0.0,0.0,1.23}};
  const float black[3] = {0.0, 0.0, 0.0};
  glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,black);
  RenderString3D("X",n_pos[0]);
  RenderString3D("Y",n_pos[1]);
  RenderString3D("Z",n_pos[2]);*/
}

void DrawSys::DrawSubAxis(float d, float s, const float col[3]) {
  glBegin(GL_QUAD_STRIP);
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, col);
  for (float i=0; i<=6.0; i+=1.0) {
    const float t = i * 2.0 * M_PI / 6.0;
    glNormal3f(std::cos(t)    , 0.0, std::sin(t)    );
    glVertex3f((d * std::cos(t) ), 0.0, (d * std::sin(t) ) );
    glVertex3f((d * std::cos(t) ), s  , (d * std::sin(t) ) );
  }
  glEnd();
  glTranslatef(0.0, s, 0.0);
  glRotatef(-90.0, 1.0, 0.0, 0.0);
  glutSolidCone(2.0 * d, 4.0 * d, 5, 5);
}

void DrawSys::Resize(int w,int h) const {
  glViewport(0, 0, w, h);
  
  glMatrixMode(GL_PROJECTION); //透視変換行列設定
  glLoadIdentity(); //変換行列初期化
  
  gluPerspective(*p_fovy, (float)w / (float)h, 0.1, 200.0); //透視投影　カメラの画角、アスペクト比を指定

  glMatrixMode(GL_MODELVIEW); //モデルビュー変換行列設定
  glLoadIdentity();

  //光源と視点が一緒に動く場合
  const GLfloat light0pos[4] = {lightpos[0][0],
				lightpos[0][1],
				lightpos[0][2],
				lightpos[0][3]};

  glLightfv(GL_LIGHT0, GL_POSITION, light0pos); 
  //

  gluLookAt(p_perscenter[0] + p_center2eye[0],
	    p_perscenter[1] + p_center2eye[1],
	    p_perscenter[2] + p_center2eye[2],
	    p_perscenter[0],
	    p_perscenter[1],
	    p_perscenter[2],
	    p_base_z[0],
	    p_base_z[1],
	    p_base_z[2]);
}

void DrawSys::ChgDrawObject() {
  draw_crit_base++;
  if (draw_crit_base == draw_crit_max) draw_crit_base = 1;
  
  for (int i=0; i<SEED_N; i++)
    draw_crit[i] = (draw_crit_base & draw_crit_mask[i]) == draw_crit_mask[i];

}

int DrawSys::Pow_n(int x, int n) const {
  int a = 1;
  while(n > 0) {
    if( (n & 1) == 0) {
      x *= x;
      n >>= 1;
    } else {
      a *= x;
      n--;
    }
  }
  return a;
}

void DrawSys::RenderCurTime() {
  std::stringstream ss;
  ss << "time = " << cur_time ;
  glColor3d(0.0, 0.0, 0.0);
  RenderString2D(ss.str().c_str(), 0.8, 0.8);
  glPopAttrib();
}

void DrawSys::RenderString2D(const char *str,float x,float y) {
  glWindowPos2f(x,y);

  while(*str) {
    glutBitmapCharacter(font, *str);
    ++str;
  }
}

void DrawSys::RenderString3D(const char *str,const float r[3]) {
  glRasterPos3f(r[0], r[1], r[2]);
  glutBitmapString(font, reinterpret_cast<const unsigned char*>(str));
}

void DrawSys::RenderSphere(const particle& prtcl) {
  int prop = prtcl.prop;
  if (prtcl.idx >= pN_init) prop = 2;
#if 1  
  if (((prtcl.idx / 4) % 2) && (prop == 0)) prop = 2;
#endif
  GLfloat color[3] = {p_color[prop][0], p_color[prop][1], p_color[prop][2]};
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);
  glPushMatrix();
  glTranslated(prtcl.r[0], prtcl.r[1], prtcl.r[2]);
  // glutSolidSphere(prtcl.rad, 10, 10);
  glutSolidSphere(0.02, 10, 10);
  glPopMatrix();
}

void DrawSys::RenderPoint(const particle& prtcl) {
  
}

bool DrawSys::IsDrawnObject(const particle& prtcl) {
  if (!cut_but) {
    return draw_crit[prtcl.prop];
  } else {
    const float in_prod = nv[0] * prtcl.r[0] * inv_box_size[0] + nv[1] * prtcl.r[1] * inv_box_size[1] + nv[2] * prtcl.r[2] * inv_box_size[2];
    const bool ret = (in_prod < cut_plane) && draw_crit[prtcl.prop];
    return ret;
  }
}

void DrawSys::ChangeNormalVector(const int i) {
  nv[0] = nv[1] = nv[2] = 0.0;
  switch (i) {
  case 0:
    nv[0] = 1.0;
    break;
  case 1:
    nv[1] = 1.0;
    break;
  case 2:
    nv[2] = 1.0;
    break;
  case 3:
    nv[0] = nv[1] = nv[2] = 1.0 / std::sqrt(3.f);
    break;
  default:
    break;
  }
}

void DrawSys::ChangeCrossSection() {
#if 0
  static float sign = 1.0f;
  if (cut_plane >= 0.2f) {
    sign = -1.0f;
  } else if(cut_plane <= -0.2f) {
    sign = 1.0f;
  }
#else
  const float sign = 1.0f;
  if (cut_plane >= 0.5 * box_size[0]) {
    cut_plane = 0.0;
  }
#endif
  if (cut_adv) cut_plane += sign * 0.03f;
}

void DrawSys::Dump2Jpg() {
  if (crit_out) {
    jpgout->GetWindowInfo();
    jpgout->SnapijgImage();
    jpgout->SaveImgJpeg(cur_time, all_time, cur_dir, time_step);
  }
}

void DrawSys::ChangeLookDirection(const int i) {
  p_center2eye[0] = p_center2eye[1] = p_center2eye[2] = 0.5;
  p_center2eye[i] = 6.0;
  p_base_z[0] = p_base_z[1] = p_base_z[2] = 0.0;
  if (i == 2)
    p_base_z[1] = 1.0;
  else
    p_base_z[2] = 1.0;
}

bool DrawSys::InitGlew() const {
  const GLenum err = glewInit();
  if (err != GLEW_OK) {
    std::cerr << glewGetErrorString(err) << std::endl;
    return false;
  }
  return true;
}

void AnimeDraw::Timer(int value) {
  if (!contined) {
    if (!crit_out)
      swt_but = false;
    else
      std::exit(0);
  }
  glutPostRedisplay();
  glutTimerFunc(50, callbacks::wrap_timer, 0);
}

void AnimeDraw::Display() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); //カラーバッファ,デプスバッファ指定
  /*光源と視点が一緒に動かず固定されている場合
    const GLfloat light0pos[4] = {lightpos[0][0],
    lightpos[0][1],
    lightpos[0][2],
    lightpos[0][3]};
    glLightfv(GL_LIGHT0, GL_POSITION, light0pos);*/
    
  if (swt_but) LoadParticleDat();
  
  for (int i = 0; i < pN; i++)
    if (IsDrawnObject(Particle[i])) RenderSphere(Particle[i]);
  
  ChangeCrossSection();
  //DrawCubic();
  //Drawxyz();

  const float col[3][3] = {{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}};
  // DrawAxis(0.02, 0.3, col);
  RenderCurTime();
  Dump2Jpg();
  glutSwapBuffers();
}

void AnimeDraw::KeyBoard(unsigned char key, int x, int y) {
  const int wid = glutGet(GLUT_WINDOW_WIDTH), hei = glutGet(GLUT_WINDOW_HEIGHT);
  switch (key) {
    case 'b':
      swt_but = true;
      break;
    case 'B':
      cut_adv = true;
      break;
    case 's':
      swt_but = false;
      break;
    case 'S':
      cut_adv = false;
      break;
    case 'q':
      swt_but = true;
      cut_but = false;
      break;
    case 'Q':
      swt_but = false;
      cut_but = false;
      break;
    case 'x':
      swt_but  = false;
      cut_adv  = cut_but = true;
      ChangeNormalVector(0);
      Resize(wid, hei);
      break;
    case 'y': 
      swt_but  = false;
      cut_adv  = cut_but = true;
      ChangeNormalVector(1);
      Resize(wid, hei);
      break;
    case 'z':
      swt_but = false;
      cut_adv = cut_but = true;
      ChangeNormalVector(2);
      Resize(wid, hei);
      break;
    case 'c':
      swt_but = false;
      cut_adv = cut_but = true;
      ChangeNormalVector(3);
      Resize(wid, hei);
      break;
    case 'X':
      ChangeLookDirection(0);
      Resize(wid, hei);
      break;
    case 'Y':
      ChangeLookDirection(1);
      Resize(wid, hei);
      break;
    case 'Z':
      ChangeLookDirection(2);
      Resize(wid, hei);
      break;
    case 'h':
      swt_but = false;
      ChgDrawObject();
      Resize(wid, hei);
      break;
    default:
      break;		  
  }
}

void AnimeDraw::PrintDrawInfo() {
  std::cout << "b restart drawing\n"       ;
  std::cout << "s stop drawing\n"          ;
  std::cout << "c start cutting(skew)\n"   ;
  std::cout << "x start cutting(x axis)\n" ;
  std::cout << "y start cutting(y axis)\n" ;
  std::cout << "z start cutting(z axis)\n" ;
  std::cout << "B restart cutting\n"       ;
  std::cout << "S stop cutting\n"          ;
  std::cout << "q quit cutting mode\n"     ;
  std::cout << "X look along x axis\n"     ;
  std::cout << "Y look along y axis\n"     ;
  std::cout << "Z look along z axis\n"     ;
  std::cout << "h change visible object.\n";
}
