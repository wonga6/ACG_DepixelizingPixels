#ifndef __node_h__
#define __node_h__

#include "edge.h"

class Node {
	public:
		// CONSTRUCTORS
		Node() : xCoor(0), yCoor(0) {}
		Node(int x, int y, Edge e) : xCoor(x), yCoor(y), edg(e) {}

		// MODIFER FUNCTIONS
		void setXCoor(int x) { xCoor = x; }
		void setYCoor(int y) { yCoor = y; }
		void setEdge(const Edge& e) { edg = e; }

		// ACCESSOR FUNCTIONS
		int getXCoor() const { return xCoor; }
		int getYCoor() const { return yCoor; }
		Edge getEdge() const { return edg; }

	private:
		int xCoor;
		int yCoor;
		Edge edg;
};

#endif
