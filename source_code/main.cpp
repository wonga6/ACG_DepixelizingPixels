#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <list>

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
void removeDiagonals2x2(Graph &similarity_graph);

// functions for weighing diagonals
void islandsHeuristic(Graph &graph, int x, int y);

int curveHeuristicHelp(Graph &graph,
		       int startNodex, int startNodey,
		       int currNodex, int currNodey,
		       int prevNodex, int prevNodey);
void curveHeuristic(Graph &graph, int nodex, int nodey);

void sparsePixelHeuristic(Graph &graph, const Image &image, int x, int y);

void removeDiagonals(Graph& graph);

// functions for creating Voronoi diagram
void makeVoronoi(const Graph &graph, std::list<Node> &nodes, std::list<Edge> &edges);

// ===========================================================================================
// FUNCTIONS FOR DEBUGGING
void printGraph(const Graph &graph); 

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

	std::cout << "Finished making graph" << std::endl;
	//printGraph(similarity_graph);

	// start of algorithm
	disconnectDissimilar(similarity_graph, image);
	//printGraph(similarity_graph);

	removeDiagonals2x2(similarity_graph);
	//printGraph(similarity_graph);

	// go through the graph and run the 3 hueristics on the cross diagonals
	for(unsigned int x = 0; x < similarity_graph.size(); x++) {
		for(unsigned int y = 0; y < similarity_graph[x].size(); y++) {

			// check that there are diagonals forming a cross
			bool cross = false;
			int n = 0;
			if(x+1 < similarity_graph.size() && y+1 < similarity_graph[x].size()) {
				cross = similarity_graph[x][y].hasEdge(x+1, y+1);
				n += cross;

				cross = similarity_graph[x][y+1].hasEdge(x+1, y);
				n += cross;
			}

			if(n != 2) continue;

			// if there's a cross perform the heuristics on it
			curveHeuristic(similarity_graph, x, y);
			sparsePixelHeuristic(similarity_graph, image, x, y);
			islandsHeuristic(similarity_graph, x, y);
		}
	}

	std::cout << "Finished heuristics" << std::endl;
	//printGraph(similarity_graph);

	// remove the diagonals based on the aggregated weight
	removeDiagonals(similarity_graph);
	printGraph(similarity_graph);


	// create Voronoi diagram
	std::list<Node> VNodes;
	std::list<Edge> VEdges;

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

	std::cout << "Finished putting nodes in graph" << std::endl;

	// add all the edges between the nodes (connect all nodes at first)
	for(unsigned int x = 0; x < similarity_graph.size(); x++) {
		for(unsigned int y = 0; y < similarity_graph[x].size(); y++) {
			// add the North edge
			if(y > 0) {
				Edge e(&similarity_graph[x][y], &similarity_graph[x][y-1], 0);
				similarity_graph[x][y].addEdge(e);
			}

			// add the North-east edge
			if(x+1 < similarity_graph.size() && y > 0) {
				Edge e(&similarity_graph[x][y], &similarity_graph[x+1][y-1], 0);
				similarity_graph[x][y].addEdge(e);
			}

			// add the East edge
			if(x+1 < similarity_graph.size()) {
				Edge e(&similarity_graph[x][y], &similarity_graph[x+1][y], 0);
				similarity_graph[x][y].addEdge(e);
			}

			// add the South-east edge
			if(x+1 < similarity_graph.size() && y+1 < similarity_graph[x].size()) {
				Edge e(&similarity_graph[x][y], &similarity_graph[x+1][y+1], 0);
				similarity_graph[x][y].addEdge(e);
			}
		}
	}

	std::cout << "Finished setting diagonals" << std::endl;
}

// ===========================================================================================

// disconnect edges of nodes with dissimilar colors
void disconnectDissimilar(Graph& similarity_graph, const Image &image) {
	// go through all the nodes
	for(unsigned int x = 0; x < similarity_graph.size(); x++) {
		for(unsigned int y = 0; y < similarity_graph[x].size(); y++) {
			Color current = image.GetPixel(x, y);

			// check North node
			if(y > 0) {
				Color north = image.GetPixel(x, y-1);
				
				// if colors aren't similar
				if(!compareYUV(current, north)) {
					//std::cout << "similar" << std::endl;
					// get the edge
					std::vector<Edge>::iterator itr = similarity_graph[x][y].getEdge(x, y-1);
					if(itr != similarity_graph[x][y].getEdges().end()) {
						// disconnect the nodes
						similarity_graph[x][y].removeEdge(itr);
					}
				}
			}

			// check North-east node
			if(x+1 < similarity_graph.size() && y > 0) {
				Color northEast = image.GetPixel(x+1, y-1);
				
				// if colors aren't similar
				if(!compareYUV(current, northEast)) {
					//std::cout << "similar" << std::endl;
					// get the edge
					std::vector<Edge>::iterator itr = similarity_graph[x][y].getEdge(x+1, y-1);
					if(itr != similarity_graph[x][y].getEdges().end()) {
						// disconnect the nodes
						similarity_graph[x][y].removeEdge(itr);
					}
				}
			}

			// check East node
			if(x+1 < similarity_graph.size()) {
				Color east = image.GetPixel(x+1, y);
				
				// if colors aren't similar
				if(!compareYUV(current, east)) {
					//std::cout << "similar" << std::endl;
					// get the edge
					std::vector<Edge>::iterator itr = similarity_graph[x][y].getEdge(x+1, y);
					if(itr != similarity_graph[x][y].getEdges().end()) {
						// disconnect the nodes
						similarity_graph[x][y].removeEdge(itr);
					}
				}
			}

			// check South-east node
			if(x+1 < similarity_graph.size() && y+1 < similarity_graph[x].size()) {
				Color southEast = image.GetPixel(x+1, y+1);
				
				// if colors aren't similar
				if(!compareYUV(current, southEast)) {
					//std::cout << "similar" << std::endl;
					// get the edge
					std::vector<Edge>::iterator itr = similarity_graph[x][y].getEdge(x+1, y+1);
					if(itr != similarity_graph[x][y].getEdges().end()) {
						// disconnect the nodes
						similarity_graph[x][y].removeEdge(itr);
					}
				}
			}
		}
	}

	std::cout << "Finished disconnecting dissimilar nodes w/ pixels of dissimilar colors" << std::endl;
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
void removeDiagonals2x2(Graph &graph) {

	// remove diagonals from 2x2 squares of pixels that are fully connected
	for(unsigned int x = 0; x < graph.size(); x++) {
		for(unsigned int y = 0; y < graph.size(); y++) {
			Node current = graph[x][y];
			int n = 0;
			bool connected = false;

			// if fully connected node above will have a East edge
			// and current node will have a North edge
			if(y > 0) {
				connected = current.hasEdge(x, y-1);
				n += connected;

				connected = graph[x][y-1].hasEdge(x+1, y-1);
				n += connected;
			}

			// if fully connected node to the right will have a North edge
			// and current node will have an East edge
			if(x+1 < graph.size()) {
				connected = current.hasEdge(x+1, y);
				n += connected;

				connected = graph[x+1][y].hasEdge(x+1, y-1);
				n += connected;
			}

			// if fully connected disconnect the diagonals if they exist
			if(n == 4) {
				if(y > 0) {
					graph[x][y].removeEdge2(x+1, y-1);
					graph[x][y-1].removeEdge2(x+1, y);
				}
			}
		}
	}

	std::cout << "Finished removing diagonsals from 2x2 squares" << std::endl;
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

	Color current = image.GetPixel(x, y);
	Color other = image.GetPixel(x, y+1);

	int currentColorCount = 0;
	int otherColorCount = 0;

	// go through the 8x8 square around the diagonals
	for(int i = x-4; i < x+4; i++) {
		for(int j = y-4; j < y+4; j++) {
			// make sure not to walk off the grid
			if(i < 0 || j < 0 || i >= graph.size() || j >= graph[i].size()) continue;

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
		graph[x][y].setEdgeWeight(x+1, y+1, otherColorCount - currentColorCount);
	}
	else {
		graph[x][y+1].setEdgeWeight(x+1, y, currentColorCount - otherColorCount);
	}
}

int curveHeuristicHelp(Graph &graph,
		       int startNodex, int startNodey,
		       int currNodex, int currNodey,
		       int prevNodex, int prevNodey) {
  if (startNodex == currNodex && startNodey == currNodey)
    return 0;
  int connectX = -1;
  int connectY = -1;
  int connections = 0;
  if (currNodex != 0 && currNodey != 0)
    if ((prevNodex != currNodex-1 || prevNodey != currNodey-1) &&
	graph[currNodex-1][currNodey-1].hasEdge(currNodex, currNodey)){
      connectX = currNodex-1;
      connectY = currNodey-1;
      connections++;
    }
  if (currNodey != 0)
    if ((prevNodex != currNodex || prevNodey != currNodey-1) &&
	graph[currNodex][currNodey].hasEdge(currNodex, currNodey-1)){
      connectX = currNodex;
      connectY = currNodey-1;
      connections++;
    }
  if (currNodex != 0)
    if ((prevNodex != currNodex-1 || prevNodey != currNodey) &&
	graph[currNodex-1][currNodey].hasEdge(currNodex, currNodey)){
      connectX = currNodex-1;
      connectY = currNodey;
      connections++;
    }
  if (currNodex != 0 && currNodey+1 != graph[currNodex].size())
    if ((prevNodex != currNodex-1 || prevNodey != currNodey+1) &&
	graph[currNodex-1][currNodey+1].hasEdge(currNodex, currNodey)){
      connectX = currNodex-1;
      connectY = currNodey+1;
      connections++;
    }
  if (currNodey+1 != graph[currNodex].size())
    if ((prevNodex != currNodex || prevNodey != currNodey+1) &&
	graph[currNodex][currNodey+1].hasEdge(currNodex, currNodey)){
      connectX = currNodex;
      connectY = currNodey+1;
      connections++;
    }
  if (currNodex+1 != graph.size() && currNodey != 0)
    if ((prevNodex != currNodex+1 || prevNodey != currNodey-1) &&
	graph[currNodex][currNodey].hasEdge(currNodex+1, currNodey-1)){
      connectX = currNodex+1;
      connectY = currNodey-1;
      connections++;
    }
  if (currNodex+1 != graph.size())
    if ((prevNodex != currNodex+1 || prevNodey != currNodey) &&
	graph[currNodex][currNodey].hasEdge(currNodex+1, currNodey)){
      connectX = currNodex+1;
      connectY = currNodey;
      connections++;
    }
  if (currNodex+1 != graph.size() && currNodey+1 != graph[currNodex].size())
    if ((prevNodex != currNodex+1 || prevNodey != currNodey+1) &&
	graph[currNodex][currNodey].hasEdge(currNodex+1, currNodey+1)){
      connectX = currNodex+1;
      connectY = currNodey+1;
      connections++;
    }
  if (connections==1){
    return 1 + curveHeuristicHelp(graph, startNodex,startNodey,
				  connectX,connectY,
				  currNodex,currNodey);
  }
  return 0;
}

void curveHeuristic(Graph &graph, int nodex, int nodey){
  int diag1 = 0;
  int diag2 = 0;
  diag1 += curveHeuristicHelp(graph, nodex+1, nodey+1, nodex, nodey,
			      nodex+1,nodey+1);
  diag1 += curveHeuristicHelp(graph, nodex, nodey, nodex+1, nodey+1,
			      nodex, nodey);
  diag2 += curveHeuristicHelp(graph, nodex+1, nodey, nodex, nodey+1,
			      nodex+1, nodey);
  diag2 += curveHeuristicHelp(graph, nodex, nodey+1, nodex+1, nodey,
			      nodex, nodey+1);
  graph[nodex][nodey].setEdgeWeight(nodex+1,nodey+1,diag1);
  graph[nodex+1][nodey].setEdgeWeight(nodex,nodey+1,diag2);
} 

// remove the diagonal with the smaller weight
void removeDiagonals(Graph& graph) {
	// go through the graph and run the 3 hueristics on the cross diagonals
	for(unsigned int x = 0; x < graph.size(); x++) {
		for(unsigned int y = 0; y < graph[x].size(); y++) {

			// check that there are diagonals forming a cross
			bool cross = false;
			int n = 0;
			if(x+1 < graph.size() && y+1 < graph[x].size()) {
				cross = graph[x][y].hasEdge(x+1, y+1);
				n += cross;

				cross = graph[x][y+1].hasEdge(x+1, y);
				n += cross;
			}

			if(n != 2) continue;

			// get weights of crosses
			std::vector<Edge>::iterator seEdgeItr = graph[x][y].getEdge(x+1, y+1);
			std::vector<Edge>::iterator otherEdgeItr = graph[x][y+1].getEdge(x+1, y);

			if(seEdgeItr->getWeight() > otherEdgeItr->getWeight()) {
				graph[x][y+1].removeEdge(otherEdgeItr);
			}
			else {
				graph[x][y].removeEdge(seEdgeItr);
			}
		}
	}	
}

// ===========================================================================================

void makeVoronoi(const Graph &graph, std::list<Node> &nodes, std::list<Edge> &edges) {
	// go through the similarity graph
	for(unsigned int x = 0; x < graph.size(); x++) {
		for(unsigned int y = 0; y < graph[x].size(); y++) {

		}
	}
}




void printGraph(const Graph &graph) {
	std::cout << graph.size() << std::endl;

	for(unsigned int x = 0; x < graph.size(); x++) {
		for(unsigned int y = 0; y < graph.size(); y++) {
			std::cout << x << " " << y << std::endl;

			std::vector<Edge> edges = graph[x][y].getEdges();
			for(unsigned int i = 0; i < edges.size(); i++) {
				std::cout << "    " << edges[i].getEnd()->getXCoor() << " " << edges[i].getEnd()->getYCoor() << " "
						  << edges[i].getWeight() << std::endl;
			}
		}
	}
}