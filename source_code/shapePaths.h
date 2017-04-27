#ifndef __shapePaths_h__
#define __shapePaths_h__

#include <cmath>

float myDistance(const std::pair<float,float>& a, const std::pair<float,float>& b){
  float distX = b.first - a.first;
  float distY = b.second - a.second;
  float dist = sqrt(distX*distX + distY*distY);
  return dist;
}

class ShapePaths {
 public:
 ShapePaths(const Color &col) : c(col) {}
  const Color& getColor() const { return c; } 
  const std::vector<std::pair<float,float> >& getSubCurve(int i) const {
    if (i >= pts.size()){
      return pts[0];
    }
    return pts[i];
  }
  const int numSubCurves() const {
    return pts.size();
  }
  void addSubCurve(std::vector<std::pair<float,float> >& subCurve) {
    pts.push_back(subCurve);
  }
  void combineCurves(){
    while (pts.size() > 1){
      int removeCurve = pts.size() - 1;
      int parentCurve = 0;
      int ptInRemove = 0;
      int ptInParent = 0;
      float currDist = myDistance(pts[removeCurve][ptInRemove],
				pts[parentCurve][ptInParent]);
      for (int i = 0; i < removeCurve; i++){
	for (int j = 0; j < pts[i].size(); j++){
	  for (int k = 0; k < pts[removeCurve].size(); k++){
	    float dist = myDistance(pts[removeCurve][k],pts[i][j]);
	    if (dist < currDist){
	      parentCurve = i;
	      ptInParent = j;
	      ptInRemove = k;
	      currDist = dist;
	    }
	  }
	}
      }
      std::vector<std::pair<float,float> > newVec;
      ptInParent;
      for (int i = 0; i < ptInParent+1; i++){
	newVec.push_back(pts[parentCurve][i]);
      }
      for (int i = ptInRemove; i < pts[removeCurve].size(); i++){
	newVec.push_back(pts[removeCurve][i]);
      }
      for (int i = 0; i < ptInRemove+1; i++){
	newVec.push_back(pts[removeCurve][i]);
      }
      for (int i = ptInParent; i < pts[parentCurve].size(); i++){
	newVec.push_back(pts[parentCurve][i]);
      }
      pts[parentCurve] = newVec;
      /*
      pts[removeCurve].insert(removeCurve+ptInRemove,
			      pts[removeCurve][ptInRemove]);
      pts[parentCurve].insert(parentCurve+ptInParent,
			      pts[parentCurve][ptInRemove]);
      pts[parentCurve].insert(parentCurve+ptInParent+1,
			      pts[removeCurve].begin(),
			      removeCurve+ptInRemove+1);
      pts[parentCurve].insert(parentCurve+ptInParent+1,
			      removeCurve+ptInRemove+1,
			      pts[removeCurve].end());
      */
      pts.pop_back();
    }
  }
 private:
  std::vector<std::vector<std::pair<float,float> > > pts;
  Color c;
};

#endif
