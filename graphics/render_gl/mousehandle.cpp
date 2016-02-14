#include "sysdraw.hpp"
#include "mousehandle.hpp"
#include <iostream>

void callbacks::wrap_mclick(int but, int state, int x, int y) {
  callbacks::mousehandle->MouseClick(but, state, x, y);
}
void callbacks::wrap_mmotion(int x, int y) {
  callbacks::mousehandle->MouseMotion(x, y);
}
void callbacks::wrap_mwheel(int wheel_n, int direct, int x, int y) {
  callbacks::mousehandle->MouseWheel(wheel_n, direct, x, y);
}

MouseHandle::MouseHandle(float center_x, float center_y, float center_z, float eDist, float fov) {
  eyeDistance  = eDist;
    
  persCenter_[0] = center_x;
  persCenter_[1] = center_y;
  persCenter_[2] = center_z;
    
  fovy_ = fov;

  phi   = M_PI * 0.35f;
  theta = M_PI * 0.3f;
    
  center2eye_[0] = eyeDistance * std::cos(phi) * std::sin(theta);
  center2eye_[1] = eyeDistance * std::sin(phi) * std::sin(theta);
  center2eye_[2] = eyeDistance * std::cos(theta);

  ebase_z_[0] = cos(phi) * std::sin(theta - 0.5f * M_PI);
  ebase_z_[1] = sin(phi) * std::sin(theta - 0.5f * M_PI);
  ebase_z_[2] = cos(theta - 0.5f * M_PI);
}

void MouseHandle::RotPersVect(float* vec, float* rot_axis, float theta) {
  quaternion quate_vec(0.0, vec[0], vec[1], vec[2]);
  quaternion quate_axis(theta, rot_axis);
  quaternion quate_temp = quate_vec * quate_axis;
  quate_axis.conj();
  quate_vec = quate_axis * quate_temp;
  quate_vec.ret_array3(vec);
}

void MouseHandle::MouseClick(int button, int state, int x, int y) {
  switch (button) {
  case GLUT_LEFT_BUTTON:
    but = 1;
    if (state == GLUT_DOWN) {
	x_bef = x;
	y_bef = y;
      }
    break;
  default:
    break;
  }
  if (state == GLUT_UP) but = 0;
}

void MouseHandle::MouseMotion(int x, int y) {
  const int  width = glutGet(GLUT_WINDOW_WIDTH);
  const int  height = glutGet(GLUT_WINDOW_HEIGHT);
  if (but == 1) {
    const float dx = (float) (x - x_bef);
    const float dy = (float) (y - y_bef);
    const float theta_u = dx * 0.005;
    const float theta_w = -dy * 0.005;
      
    float vec_u[] = {ebase_z_[0], ebase_z_[1], ebase_z_[2]};
    float vec_w[] = {0.0, 0.0, 0.0};
    CrossProduct(center2eye_, vec_u, vec_w);

    RotPersVect(center2eye_, vec_u, theta_u);
    RotPersVect(center2eye_, vec_w, theta_w);
    RotPersVect(ebase_z_,    vec_w, theta_w);

    callbacks::wrap_resize(width, height);
  }
}

void MouseHandle::MouseWheel(int wheel_n, int direct, int x, int y) {
  const int  width  = glutGet(GLUT_WINDOW_WIDTH);
  const int  height = glutGet(GLUT_WINDOW_HEIGHT);
  std::cout << wheel_n << " " << direct << std::endl;
  fovy_ -= direct * 1.5;
  callbacks::wrap_resize(width, height);
}
