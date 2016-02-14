/*
How to use

Please insert the following sentences before calling glutSwapBuffers()

GetWindowInfo();
SnapijgImage();
SaveImgJpeg(current_time,all_time_steps,current_dir);

The picture whose size is larger than 2048 pixel is declined.
*/

#include "jpegout.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <GL/freeglut.h>

void Jpegout::InitJpegOjbects(jpeg_compress_struct& cinfo, jpeg_error_mgr& jerr) const {
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
}

void Jpegout::SetImageQuality(jpeg_compress_struct& cinfo) const {
  cinfo.image_width		= b_wid;
  cinfo.image_height		= b_hei;
  cinfo.input_components	= RGB_ELEM;
  cinfo.in_color_space		= JCS_RGB;
  jpeg_set_defaults(&cinfo);
  jpeg_set_quality(&cinfo, IMG_QUALITY, TRUE);
}

FILE* Jpegout::RetOutputFilePtr(const int jpeg_time, const std::string& cur_dir, jpeg_compress_struct& cinfo) const {
  std::stringstream ss;
  ss << cur_dir << "/time" << std::setw(numb_digit) << std::setfill('0') << jpeg_time << ".jpeg";
  return fopen(ss.str().c_str(), "wb");
}

void Jpegout::WriteImage2File(jpeg_compress_struct& cinfo) {
  jpeg_start_compress(&cinfo, TRUE);
  for (int i=0; i< b_hei; i++)
    img[i] = ptr_buffer + (b_hei - i) * RGB_ELEM * b_wid;
  jpeg_write_scanlines(&cinfo, &*img.begin(), b_hei);
  jpeg_finish_compress(&cinfo);
}

void Jpegout::CleanUp(jpeg_compress_struct& cinfo, FILE* outfile) const {
  jpeg_destroy_compress(&cinfo);
  fclose(outfile);
}

void Jpegout::GetWindowInfo() {
  const int wid = glutGet(GLUT_WINDOW_WIDTH);
  const int hei = glutGet(GLUT_WINDOW_HEIGHT);  
  if (wid == b_wid && hei == b_hei) return;

  b_wid = wid;
  b_hei = hei;
  ijg_buffer.clear();
  ijg_buffer.resize(RGB_ELEM * wid * hei);
  img.resize(hei);
}

void Jpegout::SnapijgImage() {
  //http://www.khronos.org/opengles/sdk/1.1/docs/man/glPixelStorei.xml
  glReadBuffer(GL_FRONT); 
  ptr_buffer = reinterpret_cast<JSAMPLE*>(&*ijg_buffer.begin());
  glReadPixels(0, 0, b_wid, b_hei, GL_RGB, GL_UNSIGNED_BYTE, ptr_buffer);
}

void Jpegout::SaveImgJpeg(const int current_time, const int all_time_steps, const std::string& current_dir, const int time_step) {
  jpeg_compress_struct cinfo;
  jpeg_error_mgr jerr;
  
  InitJpegOjbects(cinfo, jerr);
  FILE* outfile = RetOutputFilePtr(jpeg_time, current_dir, cinfo);
  jpeg_stdio_dest(&cinfo, outfile);
  SetImageQuality(cinfo);
  WriteImage2File(cinfo);
  CleanUp(cinfo, outfile);
  
  jpeg_time++;
  if (jpeg_time > MAX_TIME) std::cerr << "jpeg_time exceeds predefined maximum time\n";
}
