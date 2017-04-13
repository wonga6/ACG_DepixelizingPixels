#ifndef __edge_h__
#define __edge_h__

#include <utility>

class Edge {
	public:
		// CONSTRUCTORS
		Edge() : weight(0) {}
		Edge(int sX, int sY, int eX, int eY, float w);

		// OPERATORS
		const Edge& operator=(const Edge& e);

		// MODIFIER FUNCTIONS
		void setStart(int x, int y) { startPoint = std::make_pair(x, y); }
		void setEnd(int x, int y) { endPoint = std::make_pair(x, y); }
		void setWeight(float w) { weight = w; }

		// ACCESSOR FUNCTIONS
		std::pair<int, int> getStart() const { return startPoint; }
		std::pair<int, int> getEnd() const { return endPoint; }
		float getWeight() const { return weight; }

	private:
		std::pair<int, int> startPoint;
		std::pair<int, int> endPoint;
		float weight;
};

bool operator==(const Edge& a, const Edge& b);

#endif
