#include "node.h"

Edge::Edge(Node *start, Node* end, float w) {
	startPoint = start;
	endPoint = end;
	weight = w;
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
