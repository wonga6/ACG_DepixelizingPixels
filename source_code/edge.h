#ifndef __edge_h__
#define __edge_h__

class Edge {
	public:
		Edge() : startPoint(0), endPoint(0), weight(0) {}
		Edge(int s, int e, float w) : startPoint(s), endPoint(e), weight(w) {}

		int getStart() { return startPoint; }
		int getEnd() { return endPoint; }
		float getWeight() { return weight; }

	private:
		int startPoint;
		int endPoint;
		float weight;
};

#endif
