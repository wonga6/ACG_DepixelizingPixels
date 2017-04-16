#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include "node.h"
#include "image.h"

typedef std::vector<std::vector<Node> > Graph;

// ===========================================================================================

void createGraph(Graph& similarity_graph, const Image &image);

// functions for disconnecting nodes of dissimilar colors
void disconnectDissimilar(Graph& similarity_graph, const Image &image);
void calcYUV(const Color &c, float &y, float &u, float &v);
bool compareYUV(const Color &current, const Color &connected);

// functions for removing diagonals in 2x2 squares of pixels
void removeDiagonals(Graph &similarity_graph);

// function for weighing diagonal edges between islands
void islandsHeuristic(Graph &graph, int x, int y);

void sparsePixelHeuristic(Graph &graph, const Image &image, int x, int y);

// ===========================================================================================
int main(int argc, char* argv[]) {

	// command line args: input, output file names
	if(argc != 3) {
		std::cerr << "Wrong number of command line args:"
			  << " need input and output file name" << std::endl;
		exit(0);
	}
	
	std::string inputFile = argv[1];
	std::string outputFile = argv[2];

	// load the image
	Image image;
	image.Load(inputFile);

	// create the graph
	Graph similarity_graph;
	createGraph(similarity_graph, image);

	// start of algorithm
	disconnectDissimilar(similarity_graph, image);
	removeDiagonals(similarity_graph);

	for(unsigned int x = 0; x < similarity_graph.size(); x++) {
		for(unsigned int y = 0; y < similarity_graph[x].size(); y++) {
			islandsHeuristic(similarity_graph, x, y);
		}
	}


	return 0;
}

// create the similarity graph from the image
void createGraph(Graph& similarity_graph, const Image &image) {

	// add all the nodes to the graph

	// create all the rows
	for(int x = 0; x < image.Width(); x++) {
		// create all the columns
		std::vector<Node> column;
		for(int y = 0; y < image.Height(); y++) {
			column.push_back(Node(x, y));
		}
		similarity_graph.push_back(column);
	}


	// add all the edges between the nodes (connect all nodes at first)
	for(unsigned int x = 0; x < similarity_graph.size(); x++) {
		for(unsigned int y = 0; y < similarity_graph[x].size(); y++) {
			// add the North edge
			if(y-1 >= 0) {
				Edge e(similarity_graph[x][y], similarity_graph[x][y-1], 0);
				similarity_graph[x][y].addEdge(e);
			}

			// add the North-east edge
			if(x+1 < similarity_graph.size() && y-1 >= 0) {
				Edge e(similarity_graph[x][y], similarity_graph[x+1][y-1], 0);
				similarity_graph[x][y].addEdge(e);
			}

			// add the East edge
			if(x+1 < similarity_graph.size()) {
				Edge e(similarity_graph[x][y], similarity_graph[x+1][y], 0);
				similarity_graph[x][y].addEdge(e);
			}

			// add the South-east edge
			if(x+1 < similarity_graph.size() && y+1 < similarity_graph[x].size()) {
				Edge e(similarity_graph[x][y], similarity_graph[x+1][y+1], 0);
				similarity_graph[x][y].addEdge(e);
			}
		}
	}
}

// ===========================================================================================

// disconnect edges of nodes with dissimilar colors
void disconnectDissimilar(Graph& similarity_graph, const Image &image) {
	// go through all the nodes
	for(unsigned int x = 0; x < similarity_graph.size(); x++) {
		for(unsigned int y = 0; y < similarity_graph[x].size(); y++) {
			Color current = image.GetPixel(x, y);

			// check North node
			if(y-1 >= 0) {
				Color north = image.GetPixel(x, y-1);
				
				// if colors aren't similar
				if(!compareYUV(current, north)) {
					// get the edge
					std::vector<Edge>::iterator itr = similarity_graph[x][y].getEdge(x, y-1);
					if(itr != similarity_graph[x][y].getEdges().end()) {
						// disconnect the nodes
						similarity_graph[x][y].removeEdge(itr);
					}
				}
			}
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
bool compareYUV(const Color &current, const Color &connected) {

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

// ===========================================================================================

// remove diagonals between 2x2 square of pixels that is fully connected

/*  A fully connected 2x2 block will have the following Node structure
		 N----NE
		 | \ / |
		 | / \ |
	    Node---E
*/
void removeDiagonals(Graph &graph) {

	// remove diagonals from 2x2 squares of pixels that are fully connected
	for(unsigned int x = 0; x < graph.size(); x++) {
		for(unsigned int y = 0; y < graph.size(); y++) {
			Node current = graph[x][y];
			bool connected = true;

			// if fully connected node will have a North and a East edge
			connected = current.hasEdge(x, y-1);
			connected = current.hasEdge(x+1, y);

			// if fully connected node above will have a East edge
			if(y-1 >= 0) {
				connected = graph[x][y-1].hasEdge(x+1, y-1);
			}

			// if fully connected node to the right will have a North edge
			if(x+1 < graph.size()) {
				connected = graph[x+1][y].hasEdge(x+1, y-1);
			}

			// if fully connected disconnect the diagonals if they exist
			if(connected) {
				// remove diagonal to North-east
				std::vector<Edge>::iterator itr = current.getEdge(x+1, y-1);
				if(itr != current.getEdges().end()) {
					graph[x][y].removeEdge(itr);
				}

				// remove above node's diagonal to the right node (N to E)
				if(y-1 >= 0) {
					itr = graph[x][y-1].getEdge(x+1, y);
					if(itr != graph[x][y-1].getEdges().end()) {
						graph[x][y-1].removeEdge(itr);
					}
				}
			}
		}
	}
}

// ===========================================================================================

// weigh edges connecting islands
// check the endpoints of the node's diagonals
// if endpoint has only one edge it's an island
void islandsHeuristic(Graph &graph, int x, int y) {
	
	// check for diagonal going North-east
	if(x+1 < graph.size() && y-1 >= 0) {
		Node northEast = graph[x+1][y-1];
		if(northEast.getEdges().size() == 1) {
			graph[x][y].setEdgeWeight(x+1, y-1, 5);
		}
	}

	// check the diagonal going South-east
	if(x+1 < graph.size() && y+1 < graph[x].size()) {
		Node southEast = graph[x+1][y+1];
		if(southEast.getEdges().size() == 1) {
			graph[x][y].setEdgeWeight(x+1, y+1, 5);
		}
	}
}

// weight of the edge is the difference in the number of nodes with the same color in the 
// 8x8 window
void sparsePixelHeuristic(Graph &graph, const Image &image, int x, int y) {

	// check that there are diagonals forming a cross
	std::vector<Edge>::const_iterator itr = graph[x][y].getEdges().end();
	if(x+1 < graph.size() && y-1 >= 0) {
		itr = graph[x][y].getEdge(x+1, y-1);
		itr = graph[x][y-1].getEdge(x+1, y);
	}

	if(itr == graph[x][y].getEdges().end()) return;

	Color current = image.GetPixel(x, y);
	Color other = image.GetPixel(x, y-1);

	int currentColorCount = 0;
	int otherColorCount = 0;

	// go through the 8x8 square around the diagonals
	for(int i = x-4; i < x+4; i++) {
		for(int j = y-4; j < y+4; j++) {
			// make sure not to walk off the grid
			if(i < 0 || j < 0 || i > graph.size() || j > graph[i].size()) continue;

			// compare the colors
			Color pix = image.GetPixel(i, j);

			if(pix.r == current.r && pix.g == current.g && pix.b == current.b) {
				currentColorCount++;
			}
			else if(pix.r == other.r && pix.g == other.g && pix.b == other.b) {
				otherColorCount++;
			}
		}
	}

	// weigh the edge with the smaller color count
	if(currentColorCount < otherColorCount) {
		graph[x][y].setEdgeWeight(x+1, y-1, otherColorCount - currentColorCount);
	}
	else {
		graph[x][y-1].setEdgeWeight(x+1, y, currentColorCount - otherColorCount);
	}
}