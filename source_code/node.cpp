#include "node.h"

const Node& Node::operator=(const Node& n) {
	xCoor = n.getXCoor();
	yCoor = n.getYCoor();
	edges = n.getEdges();

	return *this;
}

// set the edge weight of a node's edge with the ending (endX, endY)
void Node::setEdgeWeight(int endX, int endY, float weight) {
	for(unsigned int i = 0; i < edges.size(); i++) {
		if(edges[i].getEnd()->getXCoor() == endX && edges[i].getEnd()->getYCoor() == endY) {
			edges[i].addWeight(weight);
			break;
		}
	}
}

// gets the iterator to an edge with the ending (x, y) in the edges vector
std::vector<Edge>::iterator Node::getEdge(int x, int y) {
	for(std::vector<Edge>::iterator itr = edges.begin(); itr != edges.end(); itr++) {
		if(itr->getEnd()->getXCoor() == x && itr->getEnd()->getYCoor() == y) {
			return itr;
		}
	}

	return edges.end();
}

// return if the node has an edge with the ending (endX, endY)
bool Node::hasEdge(int endX, int endY) const {
	for(unsigned int i = 0; i < edges.size(); i++) {
		if(edges[i].getEnd()->getXCoor() == endX && edges[i].getEnd()->getYCoor() == endY) {
			return true;
		}
	}

	return false;
}

void Node::removeEdge2(int x, int y) {
	for(std::vector<Edge>::iterator itr = edges.begin(); itr != edges.end(); itr++) {
		if(itr->getEnd()->getXCoor() == x && itr->getEnd()->getYCoor() == y) {
			edges.erase(itr);
			break;
		}
	}
}

// two nodes are the same if they have the same
// x, y coordinate in the graph since Nodes shouldn't overlap
bool operator==(const Node& a, const Node& b) {
	return (a.getXCoor() == b.getXCoor() &&
		    a.getYCoor() == b.getYCoor());
}