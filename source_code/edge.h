#ifndef __edge_h__
#define __edge_h__

class Edge {
	public:
		// CONSTRUCTORS
		Edge() : startPoint(0), endPoint(0), weight(0) {}
		Edge(int s, int e, float w) : startPoint(s), endPoint(e), weight(w) {}

		// OPERATORS
		const Edge& operator=(const Edge& e);

		// MODIFIER FUNCTIONS
		void setStart(int s) { startPoint = s; }
		void setEnd(int e) { endPoint = e; }
		void setWeight(float w) { weight = w; }

		// ACCESSOR FUNCTIONS
		int getStart() const { return startPoint; }
		int getEnd() const { return endPoint; }
		float getWeight() const { return weight; }

	private:
		int startPoint;
		int endPoint;
		float weight;
};

bool operator==(const Edge& a, const Edge& b);

#endif
