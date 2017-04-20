#ifndef __voronoi_region_h__
#define __voronoi_region_h__

#include <list>
#include <map>
#include <list>
#include <utility>

#include "image.h"

// class to store the points that define the shape of the region
// and the color

class VoronoiRegion {
	public:
		// CONSTRUCTORS
		VoronoiRegion(const Color &col) : c(col) {}

		// MODIFIER FUNCTIONS
		void addPoint(float x, float y) { pts.push_back(std::make_pair(x, y)); }
		void setColor(const Color &clr) { c = clr; }

		void storeX (int xCor) { x = xCor; }
		void storeY (int yCor) { y = yCor; }

		// ACCESSOR FUNCTIONS
		const Color& getColor() const { return c; }
		const std::list<std::pair<float, float> >& getPts() const { return pts; }

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

		// PRINT FUNCTIONS
		void print() const {
			std::cout << "Place in sim graph: " << x << " " << y << std::endl;
			std::cout << "Color: " << c.r << " " << c.g << " " << c.b << std::endl;
			for(std::list<std::pair<float, float> >::const_iterator itr = pts.begin(); 
				itr != pts.end(); itr++) {
				std::cout << itr->first << " " << itr->second << std::endl;
			}
		}

	private:
		std::list<std::pair<float, float> > pts;
		Color c;

		// store coordinates on graph
		int x;
		int y;
};

#endif
