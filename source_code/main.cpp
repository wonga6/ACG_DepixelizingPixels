#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <list>
#include <list>
#include <utility>

#include "node.h"
#include "image.h"
#include "voronoi_region.h"

typedef std::vector<std::vector<Node> > Graph;
typedef std::vector<std::vector<VoronoiRegion> > V_Graph;
typedef std::vector<std::pair<float, float> > Spline;

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
void makeVoronoi(const Graph &graph, const Image &image, V_Graph &voronoi);
bool edgeExists(int startX, int startY, int endX, int endY, const Graph& graph);

// function for simplifying the Voronoi diagram
void simplifyVoronoi(V_Graph &voronoi);
bool isBoundary(float x, float y, int boundX, int boundY);
void addEdgesVoronoi(V_Graph &voronoi);

// functions to create shapes
void createShapes(std::vector<std::vector<VoronoiRegion> > &shapes, const Graph &similarity_graph,
	const V_Graph &voronoi);
void inShape(const VoronoiRegion &vr, const VoronoiRegion &connected,
	std::vector<std::vector<VoronoiRegion> > &shapes);
void removeSimilar(std::vector<std::vector<VoronoiRegion> > &shapes);

// functions for B-splines
void makeBSplines(std::vector<Spline> &bsplines, const V_Graph &voronoi);
std::pair<float, float> cubicBSpline(const std::vector<std::pair<float, float> > &points, float t);

// ===========================================================================================
// FUNCTIONS FOR DEBUGGING

void printGraph(const Graph &graph); 
void printVoronoi(const V_Graph &voronoi);
void printVoronoiEdges(const V_Graph &voronoi);
void printShapes(const std::vector<std::vector<VoronoiRegion> > &shapes);
void printShapeEdges(const std::vector<std::vector<VoronoiRegion> > &shapes);

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
	//printGraph(similarity_graph);


	// create Voronoi diagram
	V_Graph voronoi;
	makeVoronoi(similarity_graph, image, voronoi);
	//printVoronoi(voronoi);

	std::cout << std::endl << "Simplified Voronoi" << std::endl;

	simplifyVoronoi(voronoi);
	//printVoronoi(voronoi);

	std::cout << "Voronoi Edges" << std::endl;
	addEdgesVoronoi(voronoi);
	//printVoronoiEdges(voronoi);

	std::cout << "Create Shapes" << std::endl;
	std::vector<std::vector<VoronoiRegion> > shapes;
	createShapes(shapes, similarity_graph, voronoi);
	//printShapes(shapes);

	std::cout << "Remove Similar" << std::endl;
	removeSimilar(shapes);
	printShapeEdges(shapes);

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

// determine if there is an edge between the two nodes
bool edgeExists(int startX, int startY, int endX, int endY, const Graph& graph) {
	return graph[startX][startY].hasEdge(endX, endY);
}

// make the Voronoi diagram - represented as a 2d vector of vector of pairs(that are the points)
void makeVoronoi(const Graph &graph, const Image &image, V_Graph &voronoi) {
	// go through the similarity graph
	for(unsigned int x = 0; x < graph.size(); x++) {
		std::vector<VoronoiRegion> tmp;

		// assumes that the point (x,y) is in the corner of a square area
		for(unsigned int y = 0; y < graph[x].size(); y++) {
			VoronoiRegion shape(image.GetPixel(x, y));

			float xCenter = x + 0.5;
			float yCenter = y + 0.5;

			// Top left
			if(x > 0 && y > 0) {
				// if there's an edge going to the top left
				if(graph[x-1][y-1].hasEdge(x, y)) {
					shape.addPoint(xCenter-0.75, yCenter-0.25);
					shape.addPoint(xCenter-0.25, yCenter-0.75);
				}
				else if(graph[x-1][y].hasEdge(x, y-1)) {
					shape.addPoint(xCenter-0.25, yCenter-0.25);
				}
				else {
					shape.addPoint(xCenter-0.5, yCenter-0.5);
				}
			}
			else {
				shape.addPoint(xCenter-0.5, yCenter-0.5);
			}

			// Top
			shape.addPoint(xCenter, yCenter-0.5);

			// Top Right
			if(x+1 < graph.size() && y > 0) {
				if(graph[x][y].hasEdge(x+1, y-1)) {
					shape.addPoint(xCenter+0.25, yCenter-0.75);
					shape.addPoint(xCenter+0.75, yCenter-0.25);
				}
				else if(graph[x][y-1].hasEdge(x+1, y)) {
					shape.addPoint(xCenter+0.25, yCenter-0.25);
				}
				else {
					shape.addPoint(xCenter+0.5, yCenter-0.5);
				}
			}
			else if(x < graph.size()) {
				shape.addPoint(xCenter+0.5, yCenter-0.5);
			}

			// Right
			shape.addPoint(xCenter+0.5, yCenter);

			// Bottom Right
			if(x+1 < graph.size() && y+1 < graph[x].size()) {
				if(graph[x][y].hasEdge(x+1, y+1)) {
					shape.addPoint(xCenter+0.75, yCenter+0.25);
					shape.addPoint(xCenter+0.25, yCenter+0.75);
				}
				else if(graph[x][y+1].hasEdge(x+1, y)) {
					shape.addPoint(xCenter+0.25, yCenter+0.25);
				}
				else {
					shape.addPoint(xCenter+0.5, yCenter+0.5);
				}
			}
			else {
				shape.addPoint(xCenter+0.5, yCenter+0.5);
			}

			// Bottom
			shape.addPoint(xCenter, yCenter+0.5);

			// Bottom Left
			if(x > 0 && y+1 < graph[x].size()) {
				if(graph[x-1][y+1].hasEdge(x, y)) {
					shape.addPoint(xCenter-0.25, yCenter+0.75);
					shape.addPoint(xCenter-0.75, yCenter+0.25);
				}
				else if(graph[x-1][y].hasEdge(x, y+1)) {
					shape.addPoint(xCenter-0.25, yCenter+0.25);
				}
				else {
					shape.addPoint(xCenter-0.5, yCenter+0.5);
				}
			}
			else if(y < graph[x].size()) {
				shape.addPoint(xCenter-0.5, yCenter+0.5);
			}

			// Left
			shape.addPoint(xCenter-0.5, yCenter);

			// add the color of the shape
			shape.setColor(image.GetPixel(x,y));

			// store coordinates for debugging
			shape.storeX(x);
			shape.storeY(y);

			tmp.push_back(shape);
		}

		voronoi.push_back(tmp);
	}
}

// simplify the Voronoi diagram
void simplifyVoronoi(V_Graph &voronoi) {
	// go through shapes in Voronoi diagrma
	for(unsigned int x = 0; x < voronoi.size(); x++) {
		for(unsigned int y = 0; y < voronoi[x].size(); y++) {

			// go through all the pts in a shape
			std::list<std::pair<float, float> > pts = voronoi[x][y].getPts();
			for(std::list<std::pair<float, float> >::iterator itr = pts.begin(); 
				itr != pts.end(); itr++) {

				std::vector<std::pair<int, int> > indices;
				indices.push_back(std::make_pair(x, y));

				// check the 8 surrounding squares for the pt
				
				// TOP LEFT
				if(x > 0 && y > 0) {
					if(voronoi[x-1][y-1].hasPoint(itr->first, itr->second)) {
						indices.push_back(std::make_pair(x-1, y-1));
					}
				}

				// TOP
				if(y > 0) {
					if(voronoi[x][y-1].hasPoint(itr->first, itr->second)) {
						indices.push_back(std::make_pair(x, y-1));
					}
				}

				// TOP RIGHT
				if(x+1 < voronoi.size() && y > 0) {
					if(voronoi[x+1][y-1].hasPoint(itr->first, itr->second)) {
						indices.push_back(std::make_pair(x+1, y-1));
					}
				}

				// RIGHT
				if(x+1 < voronoi.size()) {
					if(voronoi[x+1][y].hasPoint(itr->first, itr->second)) {
						indices.push_back(std::make_pair(x+1, y));
					}
				}

				// BOTTOM RIGHT
				if(x+1 < voronoi.size() && y+1 < voronoi[x].size()) {
					if(voronoi[x+1][y+1].hasPoint(itr->first, itr->second)) {
						indices.push_back(std::make_pair(x+1, y+1));
					}
				}

				// BOTTOM
				if(y+1 < voronoi[x].size()) {
					if(voronoi[x][y+1].hasPoint(itr->first, itr->second)) {
						indices.push_back(std::make_pair(x, y+1));
					}
				}

				// BOTTOM LEFT
				if(x > 0 && y+1 < voronoi[x].size()) {
					if(voronoi[x-1][y+1].hasPoint(itr->first, itr->second)) {
						indices.push_back(std::make_pair(x-1, y+1));
					}
				}

				// LEFT
				if(x > 0) {
					if(voronoi[x-1][y].hasPoint(itr->first, itr->second)) {
						indices.push_back(std::make_pair(x-1, y));
					}
				}

				// skip points with valence 2 on the boundary
				if(isBoundary(itr->first, itr->second, voronoi.size(), voronoi[x].size())) {
					continue;
				}

				// if valance (same) is 2, erase pt from wherever it occurs in diagram
				if(indices.size() == 2) {
					for(unsigned int i = 0; i < indices.size(); i++) {
						voronoi[indices[i].first][indices[i].second].removePoint(itr->first, itr->second);
					}
				}
			}
		}
	}
}

// determine if point (x,y) is on the bounds of the grid
bool isBoundary(float x, float y, int boundX, int boundY) {
	float EPSILON = 0.0001;

	if((x > -EPSILON && x < EPSILON ) || (y > -EPSILON && y < EPSILON) ||
		(x > boundX - EPSILON && x < boundX + EPSILON) ||
		(y > boundY - EPSILON && y < boundY + EPSILON)) {
		return true;
	}

	return false;
}

// add edges between nodes of a shape
void addEdgesVoronoi(V_Graph &voronoi) {
	for(unsigned int x = 0; x < voronoi.size(); x++) {
		for(unsigned int y = 0; y < voronoi[x].size(); y++) {
			voronoi[x][y].addEdges();
		}
	}
}

// ===========================================================================================
// FUNCTIONS TO CREATE SHAPES

// main create shape function
void createShapes(std::vector<std::vector<VoronoiRegion> > &shapes, const Graph &similarity_graph,
	const V_Graph &voronoi) {

	// go through voronoi diagram & similarity graph
	for(unsigned int x = 0; x < voronoi.size(); x++) {
		for(unsigned int y = 0; y < voronoi[x].size(); y++) {

			bool hasConnection = false;

			// check the Bottom Left
			if(x > 0 && y+1 < voronoi[x].size()) {
				if(similarity_graph[x-1][y+1].hasEdge(x, y)) {
					VoronoiRegion vr = voronoi[x][y];
					VoronoiRegion connected = voronoi[x-1][y+1];
					inShape(vr, connected, shapes);
					hasConnection = true;
				}
			}

			// check the Left
			if(x > 0) {
				if(similarity_graph[x-1][y].hasEdge(x, y)) {
					VoronoiRegion vr = voronoi[x][y];
					VoronoiRegion connected = voronoi[x-1][y];
					inShape(vr, connected, shapes);
					hasConnection = true;
				}
			}

			// check the Top Left
			if(y > 0 && x > 0) {
				if(similarity_graph[x-1][y-1].hasEdge(x, y)) {
					VoronoiRegion vr = voronoi[x][y];
					VoronoiRegion connected = voronoi[x-1][y-1];
					inShape(vr, connected, shapes);
					hasConnection = true;	
				}
			}

			// check the Top
			if(y > 0) {
				if(similarity_graph[x][y].hasEdge(x, y-1)) {
					VoronoiRegion vr = voronoi[x][y];
					VoronoiRegion connected = voronoi[x][y-1];
					inShape(vr, connected, shapes);
					hasConnection = true;
				}
			}

			if(!hasConnection) {
				std::vector<VoronoiRegion> tmp(1, voronoi[x][y]);
				shapes.push_back(tmp);
			}

		}
	}

}

// add Voronoi shape to shapes vector
void inShape(const VoronoiRegion &vr, const VoronoiRegion &connected,
	std::vector<std::vector<VoronoiRegion> > &shapes) {

	for(unsigned int i = 0; i < shapes.size(); i++) {
		for(unsigned int j = 0; j < shapes[i].size(); j++) {

			// if the coonnected shape is in the vector add shape to vector
			if(shapes[i][j].getX() == connected.getX() && shapes[i][j].getY() == connected.getY()) {
				shapes[i].push_back(vr);
				return;
			}

		}
	}

	// else create a new spot in the vector
	std::vector<VoronoiRegion> tmp(1, vr);
	shapes.push_back(tmp);
}

// remove edges that are shared between shapes as they're not needed
void removeSimilar(std::vector<std::vector<VoronoiRegion> > &shapes) {

	// loop over all the shapes
	for(unsigned int i = 0; i < shapes.size(); i++) {
		for(unsigned int j = 0; j < shapes[i].size(); j++) {

			// compare the edges in a shape to the ones in the shapes of the vector
			for(unsigned int k = 0; k < shapes[i].size(); k++) {
				if(k == j) continue;

				// remove eedges that are shared between the sets
				std::set<VoronoiEdge> edgesJ = shapes[i][j].getEdges();
				std::set<VoronoiEdge>::iterator marker = edgesJ.begin();
				for(std::set<VoronoiEdge>::iterator itr = edgesJ.begin(); itr != edgesJ.end() 
					&& marker != edgesJ.end(); ) {

					if(shapes[i][k].hasEdge(*itr)) {
						marker++;
						shapes[i][k].removeEdge(*itr);
						itr = marker;
					}
					else {
						itr++;
						marker++;
					}
				}
			}

		}
	}

}

// ===========================================================================================
void makeBSplines(std::vector<Spline> &bsplines, const V_Graph &voronoi) {
	
	// go through all the shapes in the voronoi diagram
	for(unsigned int x = 0; x < voronoi.size(); x++) {
		for(unsigned int y = 0; y < voronoi[x].size(); y++) {

		}
	}
}

// cubic b splines equation
std::pair<float, float> cubicBSpline(const std::vector<std::pair<float, float> > &points, float t) {
	assert(points.size() == 4);

	// p_i-3
	float x_3 = (1 - t) * (1 - t) * (1 - t);
	x_3 = (x_3 / 6) * points[0].first;

	float y_3 = (1 - t) * (1 - t) * (1 - t);
	y_3 = (y_3 / 6) * points[0].second;

	// p_i-2
	float x_2 = (3 * (t * t * t)) - (6 * (t * t)) + 4;
	x_2 = (x_2 / 6) * points[1].first;

	float y_2 = (3 * (t * t * t)) - (6 * (t * t)) + 4;
	y_2 = (y_2 / 6) * points[1].second;

	// p_i-1
	float x_1 = (-3 * (t * t * t)) + (3 * (t * t )) + (3 * t) + 1;
	x_1 = (x_1 / 6) * points[2].first;

	float y_1 = (-3 * (t * t * t)) + (3 * (t * t )) + (3 * t) + 1;
	y_1 = (y_1 / 6) * points[2].second;

	// p_i
	float x_0 = ((t * t * t) / 6) * points[3].first;

	float y_0 = ((t * t * t) / 6) * points[3].second;

	// sum point
	float x = x_3 + x_2 + x_1 + x_0;
	float y = y_3 + y_2 + y_1 + y_0;

	std::pair<float, float> bSplinePoint = std::make_pair(x, y);
	return bSplinePoint;
}

// ===========================================================================================
// DEBUGGING FUNCTIONS

// print the nodes and their connected nodes in the similarity graph
void printGraph(const Graph &graph) {
	std::cout << "Print Similarity Graph" << std::endl;
	std::cout << "Graph size: " << graph.size() << std::endl;

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

// print the points of the regions of the voronoi diagram
void printVoronoi(const V_Graph &voronoi) {
	std::cout << "Print Voronoi Diagram" << std::endl;

	for(unsigned int x = 0; x < voronoi.size(); x++) {
		for(unsigned int y = 0; y < voronoi[x].size(); y++) {
			voronoi[x][y].print();
			std::cout << std::endl;
		}
	}
}

// print the edges of each region in the Voronoi diagram
void printVoronoiEdges(const V_Graph &voronoi) {
	std::cout << "Print the edges of each shape of the Voronoi diagram" << std::endl;

	for(unsigned int x = 0; x < voronoi.size(); x++) {
		for(unsigned int y = 0; y < voronoi[x].size(); y++) {
			voronoi[x][y].printEdges();
			std::cout << std::endl;
		}
	}
}

// print the shape groups
void printShapes(const std::vector<std::vector<VoronoiRegion> > &shapes) {

	std::cout << "==========================" << std::endl;

	for(unsigned int i = 0; i < shapes.size(); i++) {
		for(unsigned int j = 0; j < shapes[i].size(); j++) {
			Color c = shapes[i][j].getColor();
			std::cout << "Color: " << c.r << " " << c.g << " " << c.b << std::endl;
			std::cout << shapes[i][j].getX() << " " << shapes[i][j].getY() << std::endl;
		}

		std::cout << std::endl;
	}

	std::cout << "==========================" << std::endl;
}

// print the edges of the shape groups
void printShapeEdges(const std::vector<std::vector<VoronoiRegion> > &shapes) {

	std::cout << "==========================" << std::endl;

	for(unsigned int i = 0; i < shapes.size(); i++) {
		for(unsigned int j = 0; j < shapes[i].size(); j++) {
			shapes[i][j].printEdges();
		}

		std::cout << std::endl;
	}

	std::cout << "==========================" << std::endl;
}