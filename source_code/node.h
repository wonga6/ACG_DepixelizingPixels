#ifndef __node_h__
#define __node_h__

#include <vector>

// ===========================================================================================
// NODE CLASS

/* Each Node has 4 edges to Nodes - N, NE, E, SE
			N   NE
			|  /
		   Node - E
		       \
		       SE
	Each Node stores its x, y coordinates in the
	vector<vector<Node> > graph

	Each Node stores a vector of all its Edges
*/

// forward declaration of Edge class
class Edge;

class Node {
	public:
		// CONSTRUCTORS
		Node() : xCoor(0), yCoor(0) {}
		Node(int x, int y) : xCoor(x), yCoor(y) {}

		// OPERATORS
		const Node& operator=(const Node& e);

		// MODIFER FUNCTIONS
		void setXCoor(int x) { xCoor = x; }
		void setYCoor(int y) { yCoor = y; }
		void addEdge(const Edge& e) { edges.push_back(e); }
		void removeEdge(std::vector<Edge>::iterator place) { edges.erase(place); }
		void setEdgeWeight(int endX, int endY, float weight);

		// ACCESSOR FUNCTIONS
		int getXCoor() const { return xCoor; }
		int getYCoor() const { return yCoor; }
		int numEdges() const { return edges.size(); }
		std::vector<Edge>::iterator getEdge(int x, int y);
		const std::vector<Edge>& getEdges() const { return edges; }

		bool hasEdge(int endX, int endY) const;

	private:
		int xCoor;
		int yCoor;
		std::vector<Edge> edges;
};

// ===========================================================================================
// EDGE CLASS

/*
	Each Edge stores:
		-it's weight
		-it's start and end Nodes
*/

class Edge {
	public:
		// CONSTRUCTORS
		Edge() : weight(0) {}
		Edge(const Node &start, const Node &end, float w);

		// OPERATORS
		const Edge& operator=(const Edge& e);

		// MODIFIER FUNCTIONS
		void setStart(const Node &n) { startPoint = n; }
		void setEnd(const Node &n) { endPoint = n; }
		void addWeight(float w) { weight += w; }

		// ACCESSOR FUNCTIONS
		Node getStart() const { return startPoint; }
		Node getEnd() const { return endPoint; }
		float getWeight() const { return weight; }

	private:
		Node startPoint;
		Node endPoint;
		float weight;
};

// ===========================================================================================

// COMPARISION OPERATORS
bool operator==(const Edge& a, const Edge& b);
bool operator==(const Node& a, const Node& b);

#endif
