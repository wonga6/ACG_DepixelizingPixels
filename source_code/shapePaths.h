#ifndef __shapePaths_h__
#define __shapePaths_h__

class ShapePaths {
 public:
 ShapePaths(const Color &col) : c(col) {}
  const Color& getColor() const { return c; } 
  const std::vector<std::pair<float,float> >& getSubCurve(int i) const {
    return pts[i];
  }
  const int numSubCurves() const {
    return pts.size();
  }
  void addSubCurve(std::vector<std::pair<float,float> >& subCurve) {
    pts.push_back(subCurve);
  }
  void combineCurves(){
  }
 private:
  std::vector<std::vector<std::pair<float,float> > > pts;
  Color c;
};

#endif
