#pragma once
#include <memory>
#include <cmath>

class MouseHandle { 
  struct quaternion {
    float t, x, y, z;
    inline quaternion() {
      t = x = y = z = 0.0;
    }
    inline explicit quaternion(const float a[]) : t(0.0), x(a[0]), y(a[1]), z(a[2]) {}
    inline quaternion(const float t_, const float x_, const float y_, const float z_) : t(t_), x(x_), y(y_), z(z_) {}
    inline quaternion(const float radian, float axis[]) {
      const float inv_dr = 1.0f / std::sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
      axis[0] *= inv_dr;
      axis[1] *= inv_dr;
      axis[2] *= inv_dr;

      const float sss = std::sin(0.5f * radian);
      
      t = std::cos(0.5f * radian);
      x = sss * axis[0];
      y = sss * axis[1];
      z = sss * axis[2];
    }
    
    inline const quaternion operator*(const quaternion& obj) const {
      const quaternion ret(t * obj.t - x * obj.x - y * obj.y - z * obj.z,
			   t * obj.x + x * obj.t + y * obj.z - z * obj.y,
			   t * obj.y + y * obj.t + z * obj.x - x * obj.z,
			   t * obj.z + z * obj.t + x * obj.y - y * obj.x);
      return ret;
    }
    
    inline void conj() {
      x *= -1.0;
      y *= -1.0;
      z *= -1.0;
    }
    
    inline void ret_array3(float vec[3]) {
      vec[0] = x; vec[1] = y; vec[2] = z;
    }
  };
  
  float phi, theta, fovy_;
  float eyeDistance;
  float persCenter_[3];
  float center2eye_[3];  
  float ebase_z_[3];
  
  int x_bef = 0, y_bef = 0;
  int but = 0;
  
  void CrossProduct(const float* a, const float* b, float* c) const {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
  }
public:
  MouseHandle(float, float, float, float, float);

  float* fovy()       { return &fovy_; }
  float* persCenter() { return persCenter_; }
  float* center2eye() { return center2eye_; }
  float* ebase_z()    { return ebase_z_; }

  void RotPersVect(float*, float*, float);
  void MouseClick(int, int, int, int);
  void MouseMotion(int, int);
  void MouseWheel(int, int, int, int);
  
};

namespace callbacks{
  extern std::unique_ptr<MouseHandle> mousehandle;
  extern void wrap_mclick(int but, int state, int x, int y);
  extern void wrap_mmotion(int x, int y);
  extern void wrap_mwheel(int wheel_n, int direct, int x, int y);
}
