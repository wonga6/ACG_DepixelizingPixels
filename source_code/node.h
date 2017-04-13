#ifndef __node_h__
#define __node_h__

#include "edge.h"
#include <vector>

class Node {
	public:
		// CONSTRUCTORS
		Node() : xCoor(0), yCoor(0) {}
		Node(int x, int y) : xCoor(x), yCoor(y) {}

		// MODIFER FUNCTIONS
		void setXCoor(int x) { xCoor = x; }
		void setYCoor(int y) { yCoor = y; }
		void addEdge(const Edge& e) { edges.push_back(e); }
		void removeEdge(std::vector<Edge>::iterator place) { edges.erase(place); }
		void setEdgeWeight(int endX, int endY, float weight) {
			for(unsigned int i = 0; i < edges.size(); i++) {
				if(edges[i].getEnd().first == endX && edges[i].getEnd().second == endY) {
					edges[i].setWeight(weight);
				}
			}
		}

		// ACCESSOR FUNCTIONS
		int getXCoor() const { return xCoor; }
		int getYCoor() const { return yCoor; }
		int numEdges() const { return edges.size(); }
		std::vector<Edge>::iterator getEdge(int x, int y) {
			for(std::vector<Edge>::iterator itr = edges.begin(); itr != edges.end(); itr++) {
				if(itr->getEnd().first == x && itr->getEnd().second == y) {
					return itr;
				}
			}

			return edges.end();
		}
		const std::vector<Edge>& getEdges() { return edges; }

		bool hasEdge(int endX, int endY) const {
			for(unsigned int i = 0; i < edges.size(); i++) {
				if(edges[i].getEnd().first == endX && edges[i].getEnd().second == endY) {
					return true;
				}
			}

			return false;
		}

	private:
		int xCoor;
		int yCoor;
		std::vector<Edge> edges;
};

#endif
