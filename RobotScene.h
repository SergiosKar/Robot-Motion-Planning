#ifndef ROBOT_SCENE_H
#define ROBOT_SCENE_H

#include <VVRScene/scene.h>
#include <VVRScene/canvas.h>
#include <VVRScene/utils.h>
#include <GeoLib.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <ctime>

#include  <list>
using namespace std;




// A directed graph using adjacency list representation
class Graph
{
	int V;    // No. of vertices in graph
	list<int> *adj; // Pointer to an array containing adjacency lists

	// A recursive function used by printAllPaths()
	void printAllPathsUtil(int, int, bool[], int[], int &);

public:
	Graph(int V);   // Constructor
	void addEdge(int u, int v);
	void printAllPaths(int s, int d);
};




struct Tri
{
	C2DPoint *v1;
	C2DPoint *v2;
	C2DPoint *v3;
	float area;

	Tri(C2DPoint *v1, C2DPoint *v2, C2DPoint *v3) : v1(v1), v2(v2), v3(v3) {
		area = to_C2D().GetArea();
	}

	C2DTriangle to_C2D() const { return C2DTriangle(*v1, *v2, *v3); }

	vvr::Triangle2D to_vvr(vvr::Colour col = vvr::Colour::black, bool filled = false) const {
		vvr::Triangle2D t(v1->x, v1->y, v2->x, v2->y, v3->x, v3->y, col);
		t.setSolidRender(filled);
		return t;
	}

	bool operator < (const Tri& other) const {
		if (area != other.area) return (area < other.area);
		else if (v1 != other.v1) return v1 < other.v1;
		else if (v2 != other.v2) return v2 < other.v2;
		else if (v3 != other.v3) return v3 < other.v3;
	}
};

class RobotScene : public vvr::Scene
{
public:
	RobotScene();
	
    const char* getName() const override {
        return "UNIVERSITY OF PATRAS - VVR GROUP - COMPUTATIONAL GEOMETRY LAB";
    }

protected:
    void draw() override;
    void reset() override;
    void keyEvent(unsigned char key, bool up, int modif) override;
    
	
	//void MinkowskiSum(C2DPolygon robot, C2DPolygon &obstacle);



private:
    vvr::Canvas2D m_canvas;

	C2DPolygon bound;
	std::vector<C2DPoint> boundpoints;

    std::vector<C2DPoint> robotpts;
	C2DCircle robot1;
	C2DLine robot2;
	C2DPolygon robot3;
	
	std::vector<C2DPolygon> obstpolygons;
	std::vector<C2DPolygon> cspace_obsts;


    int m_style_flag;
    float m_lw_canvas;
    float m_sz_pt;
	float m_lw_tris;
	

	int robottype;
	C2DPointSet m_pts;
	std::vector<C2DLine> VoronoiLineset;
	std::vector<C2DPoint> VoronoiNodes;
	C2DPointSet tr_pts;
	std::vector<Tri> m_tris;
	std::vector<C2DLine> m_path;

	

private:
	bool readFile();

	void VoronoiViaDelauny();
	void processPoint(C2DPoint *p);
	
	
	C2DPolygon Polygon_Intersection(C2DPolygon p1, C2DPolygon p2);
	C2DPolygon HalfPlane_Intesection(std::vector<C2DPolygon> H);
		

	std::vector<C2DLine> RobotScene::path3DOF(C2DPoint start, C2DPoint end, std::vector<C2DPolygon> &cs_obsts);
};




void MinkowskiSumCircle(C2DCircle rb, std::vector<C2DPolygon> obstacles, std::vector<C2DPolygon> &cs_obsts);
void MinkowskiSumLine(C2DLine rb, std::vector<C2DPolygon> obstacles, std::vector<C2DPolygon> &cs_obsts, vvr::Canvas2D &canvas);
void MinkowskiSumRect(C2DPolygon rb, std::vector<C2DPolygon> obstacles, std::vector<C2DPolygon> &cs_obsts, vvr::Canvas2D &canvas);



C2DCircle GetCircumCircle(const C2DTriangle &t);

bool IsDelaunay(const C2DTriangle &t, const C2DPointSet &pset);

bool FindAdjacentTriangle(std::vector<Tri> &tris, C2DPoint *p1, C2DPoint *p2, unsigned *tri_adj_index, C2DPoint **opposite_vertex);

void FindViolations(std::vector<Tri> &tris, const C2DPointSet &ptset, std::vector<unsigned> &violations);

void ShowViolations(std::vector<Tri> &tris, const std::vector<unsigned> &violations, vvr::Canvas2D &canvas, const vvr::Colour &col);

void TriangulationRecurrentFunction(std::vector<Tri> &tris, vvr::Canvas2D &canvas, const C2DPointSet &pts, C2DPoint *p1new, C2DPoint *p2new, C2DPoint *p3new);


void FindAllAdjacent(std::vector<Tri> &tris, Tri triangle, int index, vvr::Canvas2D &m_canvas, std::vector<C2DLine>  &vorlines);

int NodeExists(C2DPoint node, std::vector<C2DPoint> &nodevector);

bool LineExists(C2DLine line, std::vector<C2DLine> &vorlines);

std::vector<C2DPoint> findPath(C2DPoint start, C2DPoint end, std::vector<C2DPoint> &nodevector, std::vector<C2DLine> &VoronoiLineset, vvr::Canvas2D &canvas);

std::vector<C2DLine> VisibilityGraph(C2DPoint rbpt, std::vector<C2DPolygon> cs_obsts, C2DPolygon bound);

std::vector<C2DLine> pathViaVisibility(C2DPoint start, C2DPoint end, std::vector<C2DPolygon> cs_obsts, vvr::Canvas2D &m_canvas, C2DPolygon bound);

C2DPolygon Polygon_Union(C2DPolygon p1, C2DPolygon p2);

#endif // ROBOT_SCENE
