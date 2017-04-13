#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include "node.h"
#include "image.h"

void calcYUV(const Color &c, float &y, float &u, float &v);
bool compareYUV(const Color &current, const Color &connected);
void findNode(const std::vector<Node> &graph, Node &n, int x, int Y);
void islandsHeuristic(std::vector<Node> &graph);

int main(int argc, char* argv[]) {

	// command line args: input, output file names
	if(argc != 3) {
		std::cerr << "Wrong number of command line args:"
			  << " need input and output file name" << std::endl;
		exit(0);
	}
	
	std::string inputFile = argv[1];
	std::string outputFile = argv[2];

	Image image;
	image.Load(inputFile);

	std::vector<Node> graph;

	// fill the graph with nodes and the edges between them
	for(int i = 0; i < image.Width(); i++) {
		for(int j = 0; j < image.Height(); j++) {
			// in 4 directions connect edges
			Node n(i, j);

			// north
			if(j-1 >= 0) {
				Edge e(i, j, i, j-1, 0);
				n.addEdge(e);
			}

			// north east
			if(i+1 < image.Width() && j-1 >= 0) {
				Edge e(i, j, i+1, j-1, 0);
				n.addEdge(e);
			}

			// east
			if(i+1 < image.Width()) {
				Edge e(i, j, i+1, j, 0);
				n.addEdge(e);
			}

			// south east
			if(i+1 < image.Width() && j+1 < image.Height()) {
				Edge e(i, j, i+1, j+1, 0);
				n.addEdge(e);
			}

			graph.push_back(n);
		}
	}

	// search through graph for edges connecting pixels with dissimilar colors
	for(unsigned int i = 0; i < graph.size(); i++) {
		Color current = image.GetPixel(graph[i].getXCoor(), graph[i].getYCoor());

		// iterate through all edges of the node
		std::vector<Edge> edges = graph[i].getEdges();
		for(std::vector<Edge>::iterator itr = edges.begin(); itr != edges.end(); ) {

			std::pair<int, int> otherPixel = itr->getEnd();
			Color otherColor = image.GetPixel(otherPixel.first, otherPixel.second);

			// compare color of pixel and the pixel connected by the edge
			if(compareYUV(current, otherColor)) {
				itr++;
			}
			else {
				graph[i].removeEdge(itr);
				itr = edges.erase(itr);
			}
		}
	}

	// remove diagonals between 2x2 square of pixels that is fully connected
	for(unsigned int i = 0; i < graph.size(); i++) {
		Node current = graph[i];

		// get the nodes current is connected to
		Node above, right;
		findNode(graph, above, current.getXCoor(), current.getYCoor() - 1);
		findNode(graph, right, current.getXCoor() + 1, current.getYCoor());

		// check if fully connected by edges
		bool connected = current.hasEdge(current.getXCoor(), current.getYCoor() - 1);
		connected = current.hasEdge(current.getXCoor() + 1, current.getYCoor() - 1);
		connected = current.hasEdge(current.getXCoor() + 1, current.getYCoor());
		connected = above.hasEdge(above.getXCoor() + 1, above.getYCoor());
		connected = above.hasEdge(above.getXCoor() + 1, above.getYCoor() + 1);
		connected = right.hasEdge(right.getXCoor(), right.getYCoor() - 1);

		// if fully connected remove diagonals
		if(connected) {
			current.removeEdge(current.getEdge(current.getXCoor() + 1, current.getYCoor() - 1));
			above.removeEdge(above.getEdge(above.getXCoor() + 1, above.getYCoor() + 1));
		}

	}

	islandsHeuristic(graph);


	return 0;
}

// find a node with (startX, startY) and (endX, endY)
void findNode(const std::vector<Node> &graph, Node &n, int x, int y) {
	for(unsigned int i = 0; i < graph.size(); i++) {
		if(graph[i].getXCoor() == x && graph[i].getYCoor() == y) {
			n = graph[i];
			return;
		}
	}
}

// calculate the YUV of the RGB color
void calcYUV(const Color &c, float &y, float &u, float &v) {
	float w_r = 0.299;
	float w_g = 0.587;
	float w_b = 0.114;
	float u_max = 0.436;
	float v_max = 0.615;

	y = w_r*c.r + w_g*c.g + w_b*c.b;
	u = u_max * ((c.b - y) / (1 - w_b));
	v = v_max * ((c.r - y) / (1 - w_r));
 }

// returns true if the colors are similar, false if dissimilar
bool compareYUV(Color current, Color connected) {

	float currY, currU, currV;
	float connY, connU, connV;

	calcYUV(current, currY, currU, currV);
	calcYUV(connected, connY, connU, connV);

	float diffY = connY - currY;
	float diffU = connU - currU;
	float diffV = connV - currV;

	if(diffY > (48.0f / 255) || diffU > (7.0f / 255) || 
		diffV > (6.0f) / 255) {
		return false;
	}


	return true;
}

// weigh edges connecting islands
void islandsHeuristic(std::vector<Node> &graph) {
	for(unsigned int i = 0; i < graph.size(); i++) {
		Node current = graph[i];

		// get the nodes current is connected to
		Node above, right;
		findNode(graph, above, current.getXCoor(), current.getYCoor() - 1);
		findNode(graph, right, current.getXCoor() + 1, current.getYCoor());

		// check if box has border edges
		bool connectedLeftEdge= current.hasEdge(current.getXCoor(), current.getYCoor() - 1);
		bool connectedBottomEdge = current.hasEdge(current.getXCoor() + 1, current.getYCoor());
		bool connectedTopEdge = above.hasEdge(above.getXCoor() + 1, above.getYCoor());
		bool connectedRightEdge = right.hasEdge(right.getXCoor(), right.getYCoor() - 1);

		// if box has border edges, then has more than just diagonal edges
		if(connectedRightEdge || connectedTopEdge || connectedLeftEdge || connectedBottomEdge) {
			continue;
		}

		// if diagonals are the only edges of the box:

		// weigh north east diagonal if exists
		bool connectedDiagonalNE = current.hasEdge(current.getXCoor() + 1, current.getYCoor() - 1);
		if(connectedDiagonalNE) {
			current.setEdgeWeight(current.getXCoor() + 1, current.getYCoor() - 1, 5);
		}


		// weigh south east diagonal if exists
		bool connectedDiagonalSE = current.hasEdge(current.getXCoor() + 1, current.getYCoor() + 1);
		if(connectedDiagonalSE) {
			current.setEdgeWeight(current.getXCoor() + 1, current.getYCoor() + 1, 5);
		}

	}
}