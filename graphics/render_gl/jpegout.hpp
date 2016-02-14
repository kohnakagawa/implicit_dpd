#pragma once

#include <vector>
#include <cmath>
#include <fstream>
#include <jpeglib.h>

class Jpegout {
  enum {
    MAX_TIME = 1000000,
    RGB_ELEM = 3,
    IMG_QUALITY = 90, //0 ~ 100
  };
  
  int b_wid = -1, b_hei = -1;
  std::vector<float> ijg_buffer;
  int jpeg_time = 0;
  static constexpr int numb_digit = static_cast<int>(log10(MAX_TIME) ) + 1;

  std::vector<JSAMPROW> img;
  JSAMPLE* ptr_buffer = nullptr;
  
  void InitJpegOjbects(jpeg_compress_struct& cinfo, jpeg_error_mgr& jerr) const ;
  void SetImageQuality(jpeg_compress_struct& cinfo) const ;
  FILE* RetOutputFilePtr(const int jpeg_time, const std::string& cur_dir, jpeg_compress_struct& cinfo) const ;
  void WriteImage2File(jpeg_compress_struct& cinfo);
  void CleanUp(jpeg_compress_struct& cinfo, FILE* outfile) const ;
public:
  Jpegout() {}
  void GetWindowInfo();
  void SnapijgImage();
  void SaveImgJpeg(const int current_time, const int all_time_steps, const std::string& cur_dir, const int time_step);
};



