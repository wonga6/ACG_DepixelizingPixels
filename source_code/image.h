#ifndef _IMAGE_H_
#define _IMAGE_H_

#include <cassert>
#include <string>
#include <iostream>
#include <limits>

class ImageColor {
public:
  ImageColor(int r_=255, int g_=255, int b_=255) : r(r_),g(g_),b(b_) {}
  bool isWhite() const { return r==255 && g==255 && b==255; }
  bool isBlack() const { return r==0 && g==0 && b==0; }
  int r,g,b;
};


// ====================================================================
// IMAGE CLASS
// can be saved and loaded from the standard file format .ppm 

class Image {
public:
  // ========================
  // CONSTRUCTOR & DESTRUCTOR
  Image() : width(0), height(0), data(NULL) {}
  void Allocate(int w, int h) {
    width = w;
    height = h;
    delete [] data;
    if (width == 0 && height == 0) {
      data = NULL;
    } else {
      assert (width > 0 && height > 0);
      data = new ImageColor[width*height]; 
    }
  }
  ~Image() {
    delete [] data; 
  }

  Image(const Image &image) { 
    data=NULL;
    copy_helper(image); }
  const Image& operator=(const Image &image) { 
    if (this != &image)
      copy_helper(image);
    return *this; }

  void copy_helper(const Image &image) {
    Allocate (image.Width(), image.Height());
    for (int i = 0; i < image.Width(); i++) {
      for (int j = 0; j < image.Height(); j++) {
        this->SetPixel(i,j,image.GetPixel(i,j));
      }
    }
  }

  // =========
  // ACCESSORS
  int Width() const { return width; }
  int Height() const { return height; }
  const ImageColor& GetPixel(int x, int y) const {
    assert(x >= 0);
    assert(x < width);
    assert(y >= 0);
    assert(y < height);
    return data[y*width + x]; }
  ImageColor& GetPixel(int x, int y) {
    assert(x >= 0);
    assert(x < width);
    assert(y >= 0);
    assert(y < height);
    return data[y*width + x]; }

  // =========
  // MODIFIERS
  void SetAllPixels(const ImageColor &value) {
    for (int i = 0; i < width*height; i++) {
      data[i] = value; } }
  void SetPixel(int x, int y, const ImageColor &value) {
    assert(x >= 0 && x < width);
    assert(y >= 0 && y < height);
    data[y*width + x] = value; }

  // ===========
  // LOAD & SAVE
  bool Load(const std::string &filename);
  bool Save(const std::string &filename) const; 
  
private:
  // ==============
  // REPRESENTATION
  int width;
  int height;
  ImageColor *data;
};

// ====================================================================
#endif
