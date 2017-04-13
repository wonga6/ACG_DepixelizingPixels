#include "edge.h"

Edge::Edge(int sX, int sY, int eX, int eY, float w) {
	startPoint = std::make_pair(sX, sY);
	endPoint = std::make_pair(eX, eY);
}

const Edge& Edge::operator=(const Edge& e) {
	startPoint = e.getStart();
	endPoint = e.getEnd();
	weight = e.getWeight();
	return *this;
}

bool operator==(const Edge& a, const Edge& b) {
	return (a.getStart() == b.getStart() &&
		a.getEnd() == b.getEnd() &&
		a.getWeight() == b.getWeight());
}
