#ifndef __node_h__
#define __node_h__

#include "edge.h"

class Node {
	public:
		Node() : xCoor(0), yCoor(0) {}
		Node(int x, int y, Edge e) xCoor(x), yCoor(y), edg(e) {}
	private:
		int xCoor;
		int yCoor;
		Edge edg;
};

#endif
