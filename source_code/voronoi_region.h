#ifndef __voronoi_region_h__
#define __voronoi_region_h__

#include <list>
#include <map>
#include <list>
#include <set>
#include <utility>
#include <cmath>

#include "image.h"

bool floatEquals(float a, float b){
  float EPSILON = 0.0001;
  return abs(a - b) < EPSILON;
}

// class to store the edge in a Voronoi diagram
class VoronoiEdge {
	public:
		VoronoiEdge(const std::pair<float, float> &start,
			const std::pair<float, float> &end) : startPoint(start),
			endPoint(end) { }

		std::pair<float, float> startPoint;
		std::pair<float, float> endPoint;

		bool operator<(const VoronoiEdge &ve) const {
		  if (startPoint.first < ve.startPoint.first)
		    return true;
		  else if (floatEquals(startPoint.first, ve.startPoint.first)){
		    if (startPoint.second < ve.startPoint.second)
		      return true;
		    else if (floatEquals(startPoint.second,ve.startPoint.second)){
		      if (endPoint.first < ve.endPoint.first)
			return true;
		      else if (floatEquals(endPoint.first, ve.endPoint.first) &&
			       endPoint.second < ve.endPoint.second)
			return true;
		    }
		  }
		  return false;
		}
};

bool operator==(const VoronoiEdge& a, const VoronoiEdge& b) {

	if(floatEquals(a.startPoint.first, b.startPoint.first) && 
	   floatEquals(a.startPoint.second,b.startPoint.second) &&
	   floatEquals(a.endPoint.first, b.endPoint.first) && 
	   floatEquals(a.endPoint.second, b.endPoint.second)) {
	  return true;
	}
	
	return false;
}



// class to store the points that define the shape of the region
// and the ImageColor

class VoronoiRegion {
	public:
		// CONSTRUCTORS
		VoronoiRegion(const ImageColor &col) : c(col) {}

		// MODIFIER FUNCTIONS
		void addPoint(float x, float y) { pts.push_back(std::make_pair(x, y)); }
		void setImageColor(const ImageColor &clr) { c = clr; }

		void storeX (int xCor) { x = xCor; }
		void storeY (int yCor) { y = yCor; }

		// ACCESSOR FUNCTIONS
		const int getX() const { return x; }
		const int getY() const { return y; }
		const ImageColor& getImageColor() const { return c; }
		const std::list<std::pair<float, float> >& getPts() const { return pts; }
		const std::set<VoronoiEdge>& getEdges() { return edges; }

		// OTHER FUNCTIONs

		// returns true if point is in pts
		bool hasPoint(float x, float y) {
			float EPSILON = 0.0001;

			for(std::list<std::pair<float, float> >::iterator itr = pts.begin();
				itr != pts.end(); itr++) {
				if((x > itr->first - EPSILON && x < itr->first + EPSILON) &&
					y > itr->second - EPSILON && y < itr->second + EPSILON) {
					return true;
				}
			}

			return false;
		}

		// removes point (x,y) from pts
		void removePoint(float x, float y) {
			float EPSILON = 0.0001;

			for(std::list<std::pair<float, float> >::iterator itr = pts.begin();
				itr != pts.end(); itr++) {
				if((x > itr->first - EPSILON && x < itr->first + EPSILON) &&
					y > itr->second - EPSILON && y < itr->second + EPSILON) {
					pts.erase(itr);
					break;
				}
			}
		}

		// add the edges between the points of the list - points in clockwise order
		void addEdges() {
			std::list<std::pair<float, float> >::iterator itr = pts.begin();

			std::cout << "All pts:" << std::endl;
			for(; itr != pts.end(); itr++) {
				std::cout << "(" << itr->first << "," << itr->second << ")" << std::endl;
			}
			
			std::cout << "----------------------" << std::endl;
			itr = pts.begin();
			// add the edges between the points
			while(itr != pts.end()) {
				std::pair<float, float> start = *itr;
				itr++;

				if(itr == pts.end()) break;

				std::pair<float, float> end = *itr;

				VoronoiEdge tmp(start, end);
				std::pair<std::set<VoronoiEdge>::iterator, bool> result = edges.insert(tmp);
				std::cout << "("<< start.first << " " << start.second << ")  "  << "(" << end.first << " " << end.second <<") " << result.second << std::endl;
			}

			// connect the last point in the list and the first point
			VoronoiEdge connect(pts.back(), pts.front());
			edges.insert(connect);

			std::cout << std::endl;
		}

		// erase an edge from the set
		void removeEdge(const VoronoiEdge &ve) {
		  std::set<VoronoiEdge>::iterator itr = edges.find(ve);
		  
		  // if not in the set test for the opposite direction
		  if(itr == edges.end()) {
		    VoronoiEdge opp(ve.endPoint, ve.startPoint);
		    itr = edges.find(opp);
		    
		    if(itr == edges.end()) {
		      return;
		    }
		  }
		  edges.erase(*itr);
		}

		// determine if set has an edge (or its opposite)
		bool hasEdge(const VoronoiEdge &ve) {
			std::set<VoronoiEdge>::iterator itr = edges.find(ve);

			// if not in the set test for the opposite direction
			if(itr == edges.end()) {
				VoronoiEdge opp(ve.endPoint, ve.startPoint);
				itr = edges.find(opp);

				if(itr == edges.end()) {
					return false;
				}
			}

			return true;
		}

		// PRINT FUNCTIONS
		void print() const {
			std::cout << "Place in sim graph: " << x << " " << y << std::endl;
			std::cout << "ImageColor: " << c.r << " " << c.g << " " << c.b << std::endl;
			for(std::list<std::pair<float, float> >::const_iterator itr = pts.begin(); 
				itr != pts.end(); itr++) {
				std::cout << itr->first << " " << itr->second << std::endl;
			}
		}

		void printEdges() const {
			std::cout << "Place in sim graph: " << x << " " << y << std::endl;
			std::cout << "ImageColor: " << c.r << " " << c.g << " " << c.b << std::endl;

			for(std::set<VoronoiEdge>::const_iterator itr = edges.begin(); itr != edges.end();
				itr++) {
				std::cout << "( " << itr->startPoint.first << "," << itr->startPoint.second << " )";
				std::cout << " -----> ( " << itr->endPoint.first << "," << itr->endPoint.second << " )";
				std::cout << std::endl;
			}
		}

	private:
		std::list<std::pair<float, float> > pts;
		ImageColor c;

		std::set<VoronoiEdge> edges;

		// store coordinates on graph
		int x;
		int y;
};

#endif
