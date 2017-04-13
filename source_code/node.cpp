#include "node.h"

// set the edge weight of a node's edge with the ending (endX, endY)
void Node::setEdgeWeight(int endX, int endY, float weight) {
	for(unsigned int i = 0; i < edges.size(); i++) {
		if(edges[i].getEnd().first == endX && edges[i].getEnd().second == endY) {
			edges[i].setWeight(weight);
		}
	}
}

// gets the iterator to an edge with the ending (x, y) in the edges vector
std::vector<Edge>::iterator Node::getEdge(int x, int y) {
	for(std::vector<Edge>::iterator itr = edges.begin(); itr != edges.end(); itr++) {
		if(itr->getEnd().first == x && itr->getEnd().second == y) {
			return itr;
		}
	}

	return edges.end();
}

// return if the node has an edge with the ending (endX, endY)
bool Node::hasEdge(int endX, int endY) const {
	for(unsigned int i = 0; i < edges.size(); i++) {
		if(edges[i].getEnd().first == endX && edges[i].getEnd().second == endY) {
			return true;
		}
	}

	return false;
}
