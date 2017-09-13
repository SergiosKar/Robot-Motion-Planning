#include "Robotscene.h"
#include <VVRScene/utils.h>
#include <VVRScene/canvas.h>
#include <VVRScene/mesh.h>
#include <VVRScene/settings.h>
#include <algorithm>
#include <iostream>
#include <MathGeoLib.h>
#include <fstream>
#include <string>
#include <ctime>
#include <stdio.h>

#define FLAG_SHOW_CANVAS 1
#define PI 3.14159265359

using namespace std;
using namespace vvr;
using namespace math;



RobotScene::RobotScene()
{
    m_bg_col = Colour(0x44, 0x44, 0x44);
    m_create_menus = false;
	m_lw_canvas = 2; m_perspective_proj = true;
	m_hide_log = false;
	m_hide_sliders = false;
	m_fullscreen = true;
 
    m_sz_pt = 8;
	m_lw_tris = 1;

	C2DPoint *A = new C2DPoint(-450, -330);
	C2DPoint *B = new C2DPoint(+450, -330);
	C2DPoint *C = new C2DPoint(+450, +330);
	C2DPoint *D = new C2DPoint(-450, +330);

	

	boundpoints.push_back(*A);
	boundpoints.push_back(*B);
	boundpoints.push_back(*C);
	boundpoints.push_back(*D);

	tr_pts.Add(A);
	tr_pts.Add(B);
	tr_pts.Add(C);
	tr_pts.Add(D);

	bound = C2DPolygon(&boundpoints[0],boundpoints.size(),true);

	///descritize bound///////
	C2DPoint Aa = C2DPoint(0, -328);
	C2DPoint Ba = C2DPoint(+448, 0);
	C2DPoint Ca = C2DPoint(-448, 0);
	C2DPoint Da = C2DPoint(0, +328);

	C2DPoint E = C2DPoint(-450 / 2, -328);
	C2DPoint F = C2DPoint(+450 / 2, -328);
	C2DPoint G = C2DPoint(+450 / 2, +328);
	C2DPoint H = C2DPoint(-450 / 2, +328);

	C2DPoint J = C2DPoint(+448, -328 / 2);
	C2DPoint K = C2DPoint(+448, +328 / 2);
	C2DPoint L = C2DPoint(-448, +328 / 2);
	C2DPoint M = C2DPoint(-448, -328 / 2);

	C2DPoint N = C2DPoint(-450 / 2, -328);
	C2DPoint O = C2DPoint(+450 / 2, -328);
	C2DPoint P = C2DPoint(+450 / 2, +328);
	C2DPoint Q = C2DPoint(-450 / 2, +328);




	m_pts.AddCopy(Aa);
	m_pts.AddCopy(Ba);
	m_pts.AddCopy(Ca);
	m_pts.AddCopy(Da);
	m_pts.AddCopy(E);
	m_pts.AddCopy(F);
	m_pts.AddCopy(G);
	m_pts.AddCopy(H);
	m_pts.AddCopy(J);
	m_pts.AddCopy(K);
	m_pts.AddCopy(L);
	m_pts.AddCopy(M);
	m_pts.AddCopy(N);
	m_pts.AddCopy(O);
	m_pts.AddCopy(P);
	m_pts.AddCopy(Q);
///////////////////////////////////
	m_tris.push_back(Tri(A, D, C));
	m_tris.push_back(Tri(A, B, C));


	reset();
}

void RobotScene::reset()
{
    Scene::reset();
    m_style_flag = 0;
    m_style_flag |= FLAG_SHOW_CANVAS;
 	
	readFile();
}



void RobotScene::keyEvent(unsigned char key, bool up, int modif)
{
	Scene::keyEvent(key, up, modif);
	key = tolower(key);

	switch (key){

	case ('c') : {
	}
			
	case('v') : {
		
		if (robottype == 0);
			//pathViaVisibility(robot1.GetCentre(),)
		else if (robottype == 1);
			//VisibilityGraph(robot2.GetMidPoint(), VisibilityLines, cspace_obsts);
		else if (robottype == 2);
			//pathViaVisibility(robot3.GetCentroid(), C2DPoint(-50, -300), cspace_obsts, m_canvas);

		robot2.RotateToRight(90, robot2.GetMidPoint());
		m_path=path3DOF(robot2.GetMidPoint(), C2DPoint(-300, -300),  cspace_obsts);

					
	}
						

			

	 
	}
			
	
}

void RobotScene::draw()
{
    enterPixelMode();



    //! Draw violations and anything else added to canvas
    if (m_style_flag & FLAG_SHOW_CANVAS) {
        Shape::DEF_LINE_WIDTH = m_lw_canvas;
        m_canvas.draw();
    }


	////draw robots//////////

	
	m_canvas.add(robot1, Colour::red, true);
	m_canvas.add(robot2, Colour::red);
	vvr::draw(robot3, Colour::red, true);


	///draw bound
	for (int i = 0; i < 4; i++){
		if (i==3)
			m_canvas.add(C2DLine(boundpoints[i], boundpoints[0]));
		else
			m_canvas.add(C2DLine(boundpoints[i], boundpoints[i + 1]));

	}
	
	//! Draw points
	/*
	Shape::DEF_POINT_SIZE = 8;
	vvr::draw(m_pts, Colour::red);
	*/

	/*
	//! Draw triangles
	Shape::DEF_LINE_WIDTH = m_lw_tris;

	vector<Tri>::const_iterator tri_iter;
	for (tri_iter = m_tris.begin(); tri_iter != m_tris.end(); ++tri_iter) {
		tri_iter->to_vvr(Colour::magenta).draw();
	}
	*/
	
	
    //! Draw obstacles
    Shape::DEF_POINT_SIZE = m_sz_pt;
	for (int i = 0; i < obstpolygons.size(); i++){
		vvr::draw(obstpolygons.at(i), Colour::green, false);
	}

	
	//! Draw C-SPace
	Shape::DEF_POINT_SIZE = m_sz_pt;
	for (int i = 0; i < cspace_obsts.size(); i++){
		vvr::draw(cspace_obsts.at(i), Colour::green, false);
	}
	
	/*
	//! Draw voronoi
	Shape::DEF_LINE_WIDTH = m_lw_tris;
	for (int i = 0; i < VoronoiLineset.size(); i++){
		m_canvas.add(VoronoiLineset.at(i), Colour::orange);
	}
	
	Shape::DEF_LINE_WIDTH = m_lw_tris;
	for (int i = 0; i < VoronoiNodes.size(); i++){
		m_canvas.add(VoronoiNodes.at(i), Colour::blue);
	}*/
	
	///draw path
	for (int i = 0; i < m_path.size(); i++){
		m_canvas.add(m_path[i], Colour::blue);
		

	}
	std::cout << "path" << m_path.size();
	/*
	//draw Visibility Graph
	Shape::DEF_LINE_WIDTH = m_lw_tris;
	for (int i = 0; i < VisibilityLines.size(); i++){
		m_canvas.add(VisibilityLines.at(i), Colour::blue);
	}*/

	m_canvas.add(C2DPoint(-300, -300), Colour::white);

    returnFromPixelMode();
}

bool RobotScene::readFile()
{
	
	robotpts.clear();
	obstpolygons.clear();

	fstream input("../resources/coordinates.txt");
	if (!input.is_open()){
		std::cout << "File not opened" << endl;
	}
	else{
		std::cout << "Loading shape " << endl;
	}


	//////////robot type///////////////////
	char robotch;
	input >> robotch;
	robottype = robotch - '0';

	if (robottype == 0){
		robot1 = C2DCircle(C2DPoint(-200 - 10 / sqrt(2), 100 - 10 / sqrt(2)), 10);
	}
	else if (robottype == 1){
		robot2 = C2DLine(C2DPoint(300, 0), C2DPoint(300, 280));
	}
	else if (robottype == 2){
		C2DPointSet rb3pt;
		rb3pt.AddCopy(C2DPoint(320, 320));
		rb3pt.AddCopy(C2DPoint(300, 320));
		rb3pt.AddCopy(C2DPoint(300, 300));
		rb3pt.AddCopy(C2DPoint(320, 300));

		robot3 = C2DPolygon(rb3pt);
	}
	
	//C2DPoint((-200 - 10) / sqrt(2), (100 - 10) / sqrt(2)), 10);
	

	char ch;  //  comma that separates two numbers
	C2DPoint tmp;
	char dot;

	
	char numofobst;
	input >> numofobst;
	int num = numofobst-'0';/////number of obstacles
	
	
	///////////obstacles creation///////////////////

	C2DPolygon poly;
	C2DPointSet obst;
	
	char ptsnum;

	for (int i = 0; i < num; i++){
		
		input >> ptsnum;
		int numpts = ptsnum - '0';

		for (int j = 0; j < numpts; j++){

			input >> tmp.x >> ch >> tmp.y;

			obst.AddCopy(tmp);
			//m_pts.AddCopy(tmp);
			std::cout << ch << endl;
		}
		
			
		poly = C2DPolygon(obst,true);
	/*	//obstpolygons.push_back(poly);
		
		///////////////descretisize obstacles more////
		for (int i = 0; i < poly.GetLineCount(); i++){
			if (poly.GetLine(i)->GetLength()>50){
				//obst.AddCopy(poly.GetLine(i)->GetMidPoint());
			}
		}
		/////////////////////////////////////////////////////
		//poly = C2DPolygon(obst, true);*/
		obstpolygons.push_back(poly);
		



		obst.DeleteAll();
		poly.Clear();
		
			

		

	}
	///////////////////////check functions/////////////////////////////////

	/*
	C2DPolygon rand_poly;
	rand_poly.CreateRandom(bound, 20, 20);
	

	rand_poly.GetPointsCopy(m_pts);
	*/
	
	


	/*C2DPolygon convexhull;
	convexhull.CreateConvexHull(m_pts);
	for (int i = 0; i < m_pts.size();i++)*/

	
	//VoronoiViaDelauny();
	/*
	vector<unsigned> violations;
	FindViolations(m_tris, tr_pts, violations);
	ShowViolations(m_tris, violations, m_canvas,Colour::magenta);
	*/


	

	

	//VoronoiViaDelauny();
	


	
	
	///////////////////////////////////////////////////////////////////////////////
	//input.close();


	
	return true;
}



void MinkowskiSumCircle(C2DCircle rb, std::vector<C2DPolygon> obstacles, std::vector<C2DPolygon> &cspace_obsts){

	

	//check for non convex
	
	///to polygwno paei clockwise

	std::cout << "C Space" << endl;

	///descritize circular robot////

	double radius = rb.GetRadius();
	
	std::vector<C2DPoint> circlular_robot_points;

	
	for (int i = -radius; i < radius; i++){
		double y;
		y = sqrt(radius * radius - i * i);
		C2DPoint point(i, y);
		circlular_robot_points.push_back(point);
	}
	
	for (int i = radius; i > -radius; i--){
		double y;
		y = sqrt(radius * radius - i * i);
		C2DPoint point(i, -y);
		circlular_robot_points.push_back(point);
	}

	std::vector<C2DPoint> current_sum_points;
	////minkowski for each obstacle
	
	for (int i = 0; i < obstacles.size(); i++){
		C2DPointSet obstacle_points;
		obstacles[i].GetPointsCopy(obstacle_points);

		
		for (int r = 0; r < circlular_robot_points.size(); r++){
			
			for (int j = 0; j < obstacle_points.size(); j++){
				C2DPoint obstacle_point(*obstacle_points.GetAt(j));
				C2DPoint sum_point(circlular_robot_points[r].x + obstacle_point.x, circlular_robot_points[r].y + obstacle_point.y);
				current_sum_points.push_back(sum_point);
			}
		}
		
		//compute the convex hull of all the current_sum_points
		C2DPoint *newpts = new C2DPoint[current_sum_points.size()];
		for (int k = 0; k < current_sum_points.size(); k++)
			newpts[k] = current_sum_points[k];
		C2DPolygon cloud_polygon;
		cloud_polygon.Create(newpts, current_sum_points.size());
		
		C2DPolygon current_minkowski_sum;
		current_minkowski_sum.CreateConvexHull(cloud_polygon);
		delete[] newpts;
		
		
		
		std::cout << obstacles[i].GetLineCount();
		current_sum_points.clear();

		cspace_obsts.push_back(current_minkowski_sum);

	}
	
	for (int w = cspace_obsts.size() - 1; w >= 0; w--){

		for (int q = cspace_obsts.size() - 1; q >= 0; q--){
			if (w != q)
			{
				if (cspace_obsts[w].Overlaps(cspace_obsts[q])){
					cspace_obsts[q] = Polygon_Union(cspace_obsts[q], cspace_obsts[w]);
					cspace_obsts.erase(cspace_obsts.begin() + w);
				}
			}
		}
	}
	



	
}

void MinkowskiSumLine(C2DLine rb, std::vector<C2DPolygon> obstacles,std::vector<C2DPolygon> &cspace_obsts,vvr::Canvas2D &canvas){
	

	std::cout << "C Space" << endl;
	std::vector<C2DPolygon> temp_cspace;
	std::vector<C2DPoint> current_sum_points;

	
	C2DPoint rma(rb.GetMidPoint());
	////move to center//////////////

	///question: move centroid to (0,0) or one edge to (0,0)
	rb.Move(C2DVector(rma, C2DPoint(0, 0)));
	
	C2DPoint rp1(rb.GetPointFrom());
	C2DPoint rp2(rb.GetPointTo());
	C2DPoint rm(rb.GetMidPoint());

	////-R/////
	rp1.x = -rp1.x;
	rp1.y = -rp1.y;
	rp2.x = -rp2.x;
	rp2.y = -rp2.y;

	/////////Minkowski sum for each of the obstacles
	for (int i = 0; i < obstacles.size(); i++){
		
		C2DPointSet obstacle_points;
		obstacles[i].GetPointsCopy(obstacle_points);

		//for each point of the obstacle find the p+(-r) sum
		for (int j = 0; j < obstacle_points.size(); j++){
			C2DPoint obstacle_point(*obstacle_points.GetAt(j));
			C2DPoint sum_point1(rp1.x + obstacle_point.x, rp1.y + obstacle_point.y);
			C2DPoint sum_point2(rp2.x + obstacle_point.x, rp2.y + obstacle_point.y);
			current_sum_points.push_back(sum_point1);
			current_sum_points.push_back(sum_point2);
		}

		//compute the convex hull of all the current_sum_points
		C2DPoint *newpts = new C2DPoint[current_sum_points.size()];
		for (int k = 0; k < current_sum_points.size(); k++)
			newpts[k] = current_sum_points[k];
		C2DPolygon cloud_polygon;
		cloud_polygon.Create(newpts, current_sum_points.size());
		
		C2DPolygon current_minkowski_sum;

		current_minkowski_sum.CreateConvexHull(cloud_polygon);
		delete[] newpts;
		
		

		temp_cspace.push_back(current_minkowski_sum);
		
		current_sum_points.clear();
	}

	
	
	bool addflag=false;
	for (int w = temp_cspace.size() - 1; w >= 0; w--){
		
		for (int q = temp_cspace.size() - 1; q >= 0; q--){
			if (w!=q)
			{
				if (temp_cspace[w].Overlaps(temp_cspace[q])){
					cspace_obsts.push_back(Polygon_Union(temp_cspace[q], temp_cspace[w]));
					
					temp_cspace.erase(temp_cspace.begin() + w);
					addflag = true;

				}
			}
		}
		if (addflag == false)
			cspace_obsts.push_back(temp_cspace[w]);
	}
	
	std::cout << "obst num" << cspace_obsts.size();

	

	
}

void MinkowskiSumRect(C2DPolygon rb, std::vector<C2DPolygon> obstacles, std::vector<C2DPolygon> &cspace_obsts, vvr::Canvas2D &canvas){


	std::vector<C2DPoint> current_sum_points;

	//extract the robot points and find its center
	C2DPointSet rbptset;
	
	//move to center//////////
	rb.Move(C2DVector(rb.GetCentroid(), C2DPoint(0, 0)));
	rb.GetPointsCopy(rbptset);
	

	//compute -R(0,0)
	rbptset[0].x = -rbptset[0].x;
	rbptset[0].y = -rbptset[0].y;
	rbptset[1].x = -rbptset[1].x;
	rbptset[1].y = -rbptset[1].y;
	rbptset[2].x = -rbptset[2].x;
	rbptset[2].y = -rbptset[2].y;
	rbptset[3].x = -rbptset[3].x;
	rbptset[3].y = -rbptset[3].y;

	/////////Minkowski sum for each of the obstacles
	for (int i = 0; i < obstacles.size(); i++){
		C2DPointSet obstacle_points;
		obstacles[i].GetPointsCopy(obstacle_points);

		//for each point of the polygon robot
		for (int r = 0; r < rbptset.size(); r++){
			//for each point of the obstacle find the p+(-r) sum
			for (int j = 0; j < obstacle_points.size(); j++){
				C2DPoint obstacle_point(*obstacle_points.GetAt(j));
				C2DPoint sum_point(rbptset[r].x + obstacle_point.x, rbptset[r].y + obstacle_point.y);
				current_sum_points.push_back(sum_point);
			}
		}

		//compute the convex hull of all the current_sum_points
		C2DPoint *newpts = new C2DPoint[current_sum_points.size()];
		for (int k = 0; k < current_sum_points.size(); k++)
			newpts[k] = current_sum_points[k];
		C2DPolygon cloud_polygon;
		cloud_polygon.Create(newpts, current_sum_points.size());
		//define the current minkowski_sum and add it to the vector of sums
		C2DPolygon current_minkowski_sum;
		current_minkowski_sum.CreateConvexHull(cloud_polygon);
		delete[] newpts;
	

		cspace_obsts.push_back(current_minkowski_sum);
		
		current_sum_points.clear();
		std::cout << "size" <<obstacles[i].GetLineCount();
	}

	for (int w = cspace_obsts.size() - 1; w >= 0; w--){

		for (int q = cspace_obsts.size() - 1; q >= 0; q--){
			if (w != q)
			{
				if (cspace_obsts[w].Overlaps(cspace_obsts[q])){
					cspace_obsts[q] = Polygon_Union(cspace_obsts[q], cspace_obsts[w]);
					cspace_obsts.erase(cspace_obsts.begin() + w);
				}
			}
		}
	}
}

void RobotScene::VoronoiViaDelauny(){


	/////run voronoi on cSpace///////
	for (int i = 0; i < cspace_obsts.size(); i++){

		for (int j = 0; j < cspace_obsts[i].GetLineRectCount(); j++){
			m_pts.AddCopy(*cspace_obsts[i].GetPoint(j));
		}

	}
	
	//////triangulation///////

	for (int i = 0; i < m_pts.size(); i++){
		processPoint(&m_pts[i]);

	}

	
	////extract voronoi///////

	std::vector<Tri> vor_tris = m_tris;

	for (int i = 0; i < vor_tris.size(); i++){
		//m_canvas.add(m_tris[i].to_C2D().GetCircumCentre(), Colour::blue);
		FindAllAdjacent(vor_tris, vor_tris[i], i, m_canvas,VoronoiLineset);
	}

	
	//////remove extra lines////////
	
	C2DLine vline;
	C2DLine boundline;
	bool erasionflag;////due to exception
	for (int i = VoronoiLineset.size()-1; i>=0; i--){
		vline = VoronoiLineset[i];
		for (int j = 0; j < cspace_obsts.size(); j++){
			if (cspace_obsts[j].Crosses(vline) || cspace_obsts[j].Contains(vline))
				VoronoiLineset.erase(VoronoiLineset.begin() + i);
				
		
		}
		if (bound.Crosses(vline))
				VoronoiLineset.erase(VoronoiLineset.begin() + i);
		//if (!bound.Contains(vline))
			//VoronoiLineset.erase(VoronoiLineset.begin() + i);

		
	}
	

	
	/*
	///remove duplicates lines
	C2DLine dline1;
	C2DLine dline2;
	for (int i = VoronoiLineset.size() - 1; i >= 0; i--){

		dline1 = VoronoiLineset[i];
		for (int j = VoronoiLineset.size() - 1; j >= 0; j--){
			dline2 = VoronoiLineset[j];

			if (dline1.GetPointFrom().x == dline2.GetPointTo().x && 
				dline1.GetPointFrom().y == dline2.GetPointTo().y &&
				dline1.GetPointTo().x == dline2.GetPointFrom().x &&
				dline1.GetPointTo().y == dline2.GetPointFrom().y 
				)
				VoronoiLineset.erase(VoronoiLineset.begin() + i);
				//m_canvas.add(VoronoiLineset[i], Colour::white);
			
		}
		
	
	}*/

	vector<unsigned> violations;
	FindViolations(m_tris, tr_pts, violations);
	ShowViolations(m_tris, violations, m_canvas, Colour::magenta);
	
	
	////create node vector///////
	/*
	C2DPoint v_node;
	C2DPoint v_node2;
	for (int i = 0; i < VoronoiLineset.size(); i++){
		
		v_node = VoronoiLineset[i].GetPointFrom();
		v_node2 = VoronoiLineset[i].GetPointTo();
		
		if (NodeExists(v_node, VoronoiNodes) == -1)
			VoronoiNodes.push_back(v_node);

		if (NodeExists(v_node2, VoronoiNodes) == -1)
			VoronoiNodes.push_back(v_node2);

		
	}*/
	
	std::cout << "size" << VoronoiLineset.size();

	


	/*
	std::cout << "size" << VoronoiNodes.size();
	
	///////create graph////////////

	Graph g(VoronoiNodes.size());
	
	for (int i = 0; i < VoronoiLineset.size(); i++){

		v_node = VoronoiLineset[i].GetPointFrom();
		v_node2 = VoronoiLineset[i].GetPointTo();
		g.addEdge(NodeExists(v_node, VoronoiNodes), NodeExists(v_node2, VoronoiNodes));
		//g.addEdge(NodeExists(v_node2, VoronoiNodes), NodeExists(v_node, VoronoiNodes));


	}std::cout << "size" << VoronoiNodes.size();
	int s = 2, d = 10;
	cout << "Following are all different paths from " << s
		<< " to " << d << endl;
	g.printAllPaths(s, d);
	
	*/

	

	

}


void RobotScene::processPoint(C2DPoint *p){


	
	//! Find enclosing triangle.
	unsigned i_enclosing;
	unsigned count_enclosing = 0;
	for (int i = 0; i < m_tris.size(); i++) {
		if (m_tris[i].to_C2D().Contains(*p)) {
			count_enclosing++;
			i_enclosing = i;
		}
	}

	//! [Check 2]
	//! If no enclosing triangle was found.
	//! Or if more than one were found.

	if (count_enclosing != 1) {
		delete p;
		return;
	}

	

	
	vector<Tri> tris_new;
	Tri &tri_enclosing = m_tris[i_enclosing];
	tris_new.push_back(Tri(p, tri_enclosing.v1, tri_enclosing.v2));
	tris_new.push_back(Tri(p, tri_enclosing.v2, tri_enclosing.v3));
	tris_new.push_back(Tri(p, tri_enclosing.v3, tri_enclosing.v1));

	//! [Check 3]
	//! Check if any of the 3 triangles are colinear. (Use GeoLib's function)

	if (tris_new[0].to_C2D().Collinear() ||
		tris_new[1].to_C2D().Collinear() ||
		tris_new[2].to_C2D().Collinear())
	{
		delete p;
		return;
	}

	//! HERE: We have a valid point, and we can proceed
	tr_pts.Add(p);

	m_canvas.clear();
	m_tris.erase(m_tris.begin() + i_enclosing);



	for (int i = 0; i < 3; i++)
	{
		bool did_flip = false;

		if (!IsDelaunay(tris_new[i].to_C2D(), tr_pts))
		{
			C2DPoint *p1 = tris_new[i].v1;
			C2DPoint *p2 = tris_new[i].v2;
			C2DPoint *p3 = tris_new[i].v3;

			C2DPoint *v_opposite = NULL;
			unsigned tri_adjacent_index = NULL;

			bool adj_exists = FindAdjacentTriangle(m_tris, p2, p3, &tri_adjacent_index, &v_opposite);

			if (adj_exists)
			{
				//...
				m_tris.erase(m_tris.begin() + tri_adjacent_index);

				m_tris.push_back(Tri(p, p2, v_opposite));
				//m_canvas.add(new Triangle2D(Tri(p, p2, v_opposite).to_vvr(vvr::Colour::green, true)));


				//m_canvas.add(GetCircumCircle(Tri(p, p2, v_opposite).to_C2D()), vvr::Colour::blue, false);

				if (!IsDelaunay(Tri(p, p2, v_opposite).to_C2D(), tr_pts))
					TriangulationRecurrentFunction(m_tris, m_canvas, tr_pts, p, p2, v_opposite);


				m_tris.push_back(Tri(p, p3, v_opposite));
				//m_canvas.add(new Triangle2D(Tri(p, p3, v_opposite).to_vvr(vvr::Colour::green, true)));


				//m_canvas.add(GetCircumCircle(Tri(p, p3, v_opposite).to_C2D()), vvr::Colour::blue, false);//kuklos delaugnhy


				if (!IsDelaunay(Tri(p, p3, v_opposite).to_C2D(), tr_pts))
					TriangulationRecurrentFunction(m_tris, m_canvas, tr_pts, p, p3, v_opposite);

				did_flip = true;
				//...
				//...
				//printf("adjacent exist:%d  \n", i);
			}
		}

		if (!did_flip)
		{
			//...
			C2DPoint *p2 = tris_new[i].v2;
			C2DPoint *p3 = tris_new[i].v3;

			m_tris.push_back(Tri(p, p2, p3));
			//	printf("not fliped: %d \n ", i);
			//...
			//...
		}
	}

	//! Visualize the violations.

	
	
}

void FindAllAdjacent(std::vector<Tri> &tris, Tri triangle, int index, vvr::Canvas2D &m_canvas, std::vector<C2DLine> &vorlines){

	unsigned *tri_adj_index = NULL;

	C2DPoint p1 = triangle.to_C2D().GetPoint1();
	C2DPoint p2 = triangle.to_C2D().GetPoint2();
	C2DPoint p3 = triangle.to_C2D().GetPoint3();

	//vtris.erase(vtris.begin() + index);

	for (int i = 0; i < tris.size(); i++){

		C2DLine l;
		C2DPoint  v1 = tris[i].to_C2D().GetPoint1();
		C2DPoint  v2 = tris[i].to_C2D().GetPoint2();
		C2DPoint  v3 = tris[i].to_C2D().GetPoint3();


		if (v1 == p1 && v2 == p2 && i != index) {
			l = C2DLine(tris[index].to_C2D().GetCircumCentre(), tris[i].to_C2D().GetCircumCentre());
			if (LineExists(l,vorlines)==false)
				vorlines.push_back(l);
			std::cout << i << endl;
		}
		if (v1 == p2 && v2 == p1 && i != index) {
			l = C2DLine(tris[index].to_C2D().GetCircumCentre(), tris[i].to_C2D().GetCircumCentre());
			if (LineExists(l, vorlines) == false)
			vorlines.push_back(l); std::cout << i << endl;
		}
		if (v2 == p1 && v3 == p2 && i != index) {
			l = C2DLine(tris[index].to_C2D().GetCircumCentre(), tris[i].to_C2D().GetCircumCentre());
			if (LineExists(l, vorlines) == false)
			vorlines.push_back(l); std::cout << i << endl;
		}
		if (v2 == p2 && v3 == p1 && i != index) {
			l = C2DLine(tris[index].to_C2D().GetCircumCentre(), tris[i].to_C2D().GetCircumCentre());
			if (LineExists(l, vorlines) == false)
			vorlines.push_back(l); std::cout << i << endl;
		}
		if (v3 == p1 && v1 == p2 && i != index) {
			l = C2DLine(tris[index].to_C2D().GetCircumCentre(), tris[i].to_C2D().GetCircumCentre());
			if (LineExists(l, vorlines) == false)
			vorlines.push_back(l); std::cout << i << endl;
		}
		if (v3 == p2 && v1 == p1 && i != index) {
			l = C2DLine(tris[index].to_C2D().GetCircumCentre(), tris[i].to_C2D().GetCircumCentre());
			if (LineExists(l, vorlines) == false)
			vorlines.push_back(l); std::cout << i << endl;
		}

		if (v1 == p2 && v2 == p3 && i != index) {
			l = C2DLine(tris[index].to_C2D().GetCircumCentre(), tris[i].to_C2D().GetCircumCentre());
			if (LineExists(l, vorlines) == false)
			vorlines.push_back(l); std::cout << i << endl;
		}
		if (v1 == p3 && v2 == p2 && i != index) {
			l = C2DLine(tris[index].to_C2D().GetCircumCentre(), tris[i].to_C2D().GetCircumCentre());
			if (LineExists(l, vorlines) == false)
			vorlines.push_back(l); std::cout << i << endl;
		}
		if (v2 == p2 && v3 == p3 && i != index) {
			l = C2DLine(tris[index].to_C2D().GetCircumCentre(), tris[i].to_C2D().GetCircumCentre());
			if (LineExists(l, vorlines) == false)
			vorlines.push_back(l); std::cout << i << endl;
		}
		if (v2 == p3 && v3 == p2 && i != index) {
			l = C2DLine(tris[index].to_C2D().GetCircumCentre(), tris[i].to_C2D().GetCircumCentre());
			if (LineExists(l, vorlines) == false)
			vorlines.push_back(l); std::cout << i << endl;
		}
		if (v3 == p2 && v1 == p3 && i != index) {
			l = C2DLine(tris[index].to_C2D().GetCircumCentre(), tris[i].to_C2D().GetCircumCentre());
			if (LineExists(l, vorlines) == false)
			vorlines.push_back(l); std::cout << i << endl;
		}
		if (v3 == p3 && v1 == p2 && i != index) {
			l = C2DLine(tris[index].to_C2D().GetCircumCentre(), tris[i].to_C2D().GetCircumCentre());
			if (LineExists(l, vorlines) == false)
			vorlines.push_back(l); std::cout << i << endl;
		}

		if (v1 == p1 && v2 == p3 && i != index) {
			l = C2DLine(tris[index].to_C2D().GetCircumCentre(), tris[i].to_C2D().GetCircumCentre());
			if (LineExists(l, vorlines) == false)
			vorlines.push_back(l); std::cout << i << endl;
		}
		if (v1 == p3 && v2 == p1 && i != index) {
			l = C2DLine(tris[index].to_C2D().GetCircumCentre(), tris[i].to_C2D().GetCircumCentre());
			if (LineExists(l, vorlines) == false)
			vorlines.push_back(l); std::cout << i << endl;
		}
		if (v2 == p1 && v3 == p3 && i != index) {
			l = C2DLine(tris[index].to_C2D().GetCircumCentre(), tris[i].to_C2D().GetCircumCentre());
			if (LineExists(l, vorlines) == false)
			vorlines.push_back(l); std::cout << i << endl;
		}
		if (v2 == p3 && v3 == p1 && i != index) {
			l = C2DLine(tris[index].to_C2D().GetCircumCentre(), tris[i].to_C2D().GetCircumCentre());
			if (LineExists(l, vorlines) == false)
			vorlines.push_back(l); std::cout << i << endl;
		}
		if (v3 == p1 && v1 == p3 && i != index) {
			l = C2DLine(tris[index].to_C2D().GetCircumCentre(), tris[i].to_C2D().GetCircumCentre());
			if (LineExists(l, vorlines) == false)
			vorlines.push_back(l); std::cout << i << endl;
		}
		if (v3 == p3 && v1 == p1 && i != index) {
			l = C2DLine(tris[index].to_C2D().GetCircumCentre(), tris[i].to_C2D().GetCircumCentre());
			if (LineExists(l, vorlines) == false)
			vorlines.push_back(l); std::cout << i << endl;
		}
	}


}

C2DCircle GetCircumCircle(const C2DTriangle &t)
{
	

	C2DCircle circle;
	circle.SetCircumscribed(
		t.GetPoint1(),
		t.GetPoint2(),
		t.GetPoint3()
		);
	return circle;
}

bool IsDelaunay(const C2DTriangle &t, const C2DPointSet &pset)
{
	

	for (int i = 0; i < pset.size(); i++) {
		C2DCircle c = GetCircumCircle(t);
		c.SetRadius(c.GetRadius()*0.99);
		if (c.Contains(*pset.GetAt(i))) {
			return false;
		}
	}

	return true;
}

void FindViolations(vector<Tri> &tris, const C2DPointSet &ptset, vector<unsigned> &violations)
{
	

	violations.clear();

	for (int i = 0; i < tris.size(); i++)
	{
		Tri &tri = tris[i];
		C2DTriangle t(*tri.v1, *tri.v2, *tri.v3);
		if (!IsDelaunay(t, ptset)) {
			violations.push_back(i);
		}
	}
}

void ShowViolations(vector<Tri> &tris, const vector<unsigned> &violations, Canvas2D &canvas, const Colour &col)
{
	for (int i = 0; i < violations.size(); i++) {
		Tri &tri = tris[violations[i]];
		C2DTriangle t(*tri.v1, *tri.v2, *tri.v3);
		canvas.add(GetCircumCircle(t), col, false);
	}
}

bool FindAdjacentTriangle(vector<Tri> &tris, C2DPoint *p1, C2DPoint *p2, unsigned *tri_adj_index, C2DPoint **opp_ver)
{
	

	for (int i = 0; i < tris.size(); i++)
	{
		C2DPoint  *v1 = tris[i].v1;
		C2DPoint  *v2 = tris[i].v2;
		C2DPoint  *v3 = tris[i].v3;

		if (v1 == p1 && v2 == p2) { *opp_ver = v3; *tri_adj_index = i; return true; }
		if (v1 == p2 && v2 == p1) { *opp_ver = v3; *tri_adj_index = i; return true; }
		if (v2 == p1 && v3 == p2) { *opp_ver = v1; *tri_adj_index = i; return true; }
		if (v2 == p2 && v3 == p1) { *opp_ver = v1; *tri_adj_index = i; return true; }
		if (v3 == p1 && v1 == p2) { *opp_ver = v2; *tri_adj_index = i; return true; }
		if (v3 == p2 && v1 == p1) { *opp_ver = v2; *tri_adj_index = i; return true; }
	}

	return false;
}

void TriangulationRecurrentFunction(vector<Tri> &tris, Canvas2D &canvas, const C2DPointSet &pts, C2DPoint *p1new, C2DPoint *p2new, C2DPoint *p3new){

	printf("%d", tris.size());


	bool did_flip = false;

	Tri tr = Tri(p1new, p2new, p3new);

	if (!IsDelaunay(tr.to_C2D(), pts))
	{
		C2DPoint *p1 = p1new;
		C2DPoint *p2 = p2new;
		C2DPoint *p3 = p3new;

		C2DPoint *v_opposite = NULL;
		unsigned tri_adjacent_index = NULL;

		bool adj_exists = FindAdjacentTriangle(tris, p2, p3, &tri_adjacent_index, &v_opposite);

		if (adj_exists)
		{
			//...
			tris.erase(tris.begin() + tris.size() - 1);///vazw to teleutaio na to elejw
			///h na tou pernaw to i kai na i++ se kathe anadromh
			/// h sunarthsh pou n t vriskei san find adjacent
			tris.erase(tris.begin() + tri_adjacent_index);

			tris.push_back(Tri(p1, p2, v_opposite));
			//canvas.add(new Triangle2D(Tri(p1, p2, v_opposite).to_vvr(vvr::Colour::green, true)));
			if (!IsDelaunay(Tri(p1, p2, v_opposite).to_C2D(), pts))
				TriangulationRecurrentFunction(tris, canvas, pts, p1, p2, v_opposite);


			tris.push_back(Tri(p1, p3, v_opposite));
			//canvas.add(new Triangle2D(Tri(p1, p3, v_opposite).to_vvr(vvr::Colour::green, true)));
			if (!IsDelaunay(Tri(p1, p3, v_opposite).to_C2D(), pts))
				TriangulationRecurrentFunction(tris, canvas, pts, p1, p3, v_opposite);

			did_flip = true;
			//...
			//...
			//	printf("adjacent exist:%d  \n", i);
		}
	}




}

int NodeExists(C2DPoint node, std::vector<C2DPoint> &nodevector){

	for (int i = 0; i < nodevector.size(); i++){

		if (node.x == nodevector.at(i).x && node.y == nodevector.at(i).y)
			return i;

		
	}
	return -1;

}

bool LineExists(C2DLine line, std::vector<C2DLine> &vorlines){

	for (int i = 0; i < vorlines.size(); i++){
		if (line.GetPointFrom() == vorlines[i].GetPointFrom()
			&& line.GetPointTo() == vorlines[i].GetPointTo())
			return true;

		if (line.GetPointFrom() == vorlines[i].GetPointTo()
			&& line.GetPointTo() == vorlines[i].GetPointFrom())
			return true;
	}
	return  false;

}

std::vector<C2DPoint>  findPath(C2DPoint start, C2DPoint end, std::vector<C2DPoint> &vornodes, std::vector<C2DLine> &vorlines,vvr::Canvas2D &canvas){
	

	///fix loops and deadends
	std::vector<C2DPoint> path;
	
	//////find closest node from start//////

	


	C2DCircle cir;
	bool stflag = false;

	for (int i = 0; i < 10000; i++){
		cir = C2DCircle(start, i);
		for (int n = 0; n < vorlines.size(); n++){
			if (cir.Contains(vorlines[n].GetPointFrom())){
				start = C2DPoint(vorlines[n].GetPointFrom());
				stflag = true;
				
				break;

			}
			if (cir.Contains(vorlines[n].GetPointTo())){
				start = C2DPoint(vorlines[n].GetPointTo());
				stflag = true;
				break;

			}

		}
		if (stflag == true)
			break;
	}

	//////find closest node from end//////

	
	
	 stflag = false;
	
	for (int i = 0; i < 10000; i++){
		cir = C2DCircle(end, i);
		for (int n = 0; n < vorlines.size(); n++){
			if (cir.Contains(vorlines[n].GetPointFrom())){
				end = C2DPoint(vorlines[n].GetPointFrom());
				stflag = true;
				break;

			}
			if (cir.Contains(vorlines[n].GetPointTo())){
				end = C2DPoint(vorlines[n].GetPointTo());
				stflag = true;
				break;

			}

		}
		if (stflag == true)
			break;
	}
	
	//////////find path/////////////

	bool pathflag=false;
	path.push_back(start);
	C2DPoint prevPoint;;
	C2DPoint curPoint;
	curPoint = C2DPoint(start);

	for(int i=0;i<20;i++){
	//while(curPoint!=end){


		for (int i = 0; i < vorlines.size(); i++){
			if (curPoint == vorlines[i].GetPointFrom() && vorlines[i].GetPointTo() != prevPoint ){

				
				pathflag = true;
				prevPoint = C2DPoint(vorlines[i].GetPointFrom());
				curPoint = C2DPoint(vorlines[i].GetPointTo());
				path.push_back(curPoint);
				break;

			}
			if (curPoint == vorlines[i].GetPointTo() && vorlines[i].GetPointFrom() != prevPoint){

				
				pathflag = true;
				prevPoint = C2DPoint(vorlines[i].GetPointTo());
				curPoint = C2DPoint(vorlines[i].GetPointFrom());
				path.push_back(curPoint);
				break;

			}

			//curPoint = C2DPoint(prevPoint);
		
		}
		if (pathflag == false)
			break;



	}

	
	return path;
}

std::vector<C2DLine> pathViaVisibility(C2DPoint start, C2DPoint end, std::vector<C2DPolygon> cs_obsts, vvr::Canvas2D &m_canvas, C2DPolygon bound){


	std::vector<C2DLine> currentpath;

	std::vector<C2DLine> vislines;
	//while (start!=end)
	for (int q = 0; q < 10;q++)
	{

		//if (start == end)
			//break;
		
		bool crossflag = false;

		///check if goes directly to end
		for (int i = 0; i < cs_obsts.size(); i++){
			if (cs_obsts[i].Crosses(C2DLine(start, end))){ 
				crossflag = true; break;
			}

		}
		if (crossflag == false){
			
			currentpath.push_back(C2DLine(start, end));
			start = C2DPoint(end);
			vislines.clear();
		}
		else{
			vislines = VisibilityGraph(start, cs_obsts,bound);

			std::vector<double> fns;
			///f(n)= g(n)+h(n) ,g(n)= length,h(n)= distance from end
			for (int i = 0; i < vislines.size(); i++){
				fns.push_back(vislines[i].GetPointTo().Distance(end) );

			}
			double min = 10000;
			int minindex;
			///go to min
			for (int i = 0; i < fns.size(); i++){
				if (fns[i] < min){
					min = fns[i];
					minindex = i;
				}
			}
			start = C2DPoint(vislines[minindex].GetPointTo());
			if (currentpath.size() != 0){
				if (start == currentpath.back().GetPointTo())
					break;
			}

			currentpath.push_back(vislines[minindex]);
			vislines.clear();
			
			
			
		}
	}
	return currentpath;
	

}

std::vector<C2DLine> VisibilityGraph(C2DPoint rbpt, std::vector<C2DPolygon> cs_obsts, C2DPolygon bound){

	
	std::vector<C2DLine> visLines;
	C2DLine newline;
	bool crossflag = false;
	C2DPointSet intersects;
	for (int i = 0; i < cs_obsts.size(); i++){

		for (int j = 0; j < cs_obsts[i].GetPointsCount(); j++){
			
			////line from robot////////////
			newline = C2DLine(rbpt, *cs_obsts[i].GetPoint(j));
			
			//check if crosses obstacle
			for (int k = 0; k < cs_obsts.size(); k++){
				cs_obsts[k].Crosses(newline, &intersects);
				std::cout <<"inter "<< intersects.size()<<endl;
				


				if (intersects.size()>0 ){
					intersects.DeleteAll();
					crossflag = true;
					//break;
				}
				intersects.DeleteAll();
						
							
			}
			///////check if is inside obstacle
			if (cs_obsts[i].Contains(newline.GetMidPoint())){

				crossflag = true;
			
			}
			////checks if overlaps line of polygon
			for (int m = 0; m < cs_obsts[i].GetLineCount(); m++){
				if (newline.GetPointFrom() == cs_obsts[i].GetLine(m)->GetPointFrom() &&
					newline.GetPointTo() == cs_obsts[i].GetLine(m)->GetPointTo()){

					crossflag = false;
					break;
				}
				if (newline.GetPointFrom() == cs_obsts[i].GetLine(m)->GetPointTo() &&
					newline.GetPointTo() == cs_obsts[i].GetLine(m)->GetPointFrom()){

					crossflag = false;
					break;
				}			
			
			}
			////checks if crosses bound
			if (bound.Crosses(newline))
				crossflag = true;


			if (crossflag == false)
				visLines.push_back(newline);
			crossflag = false;
			
		}
	}



	return visLines;
	//////between obstacles
	/*
	for (int i = 0; i < cs_obsts.size(); i++){
		for (int j = 0; j < cs_obsts.size(); j++){
				
			for (int k = 0; k < cs_obsts[i].GetPointsCount(); k++){
				for (int l = 0; cs_obsts[j].GetPointsCount(); l++){
					
					newline = C2DLine(*cs_obsts[i].GetPoint(k), *cs_obsts[j].GetPoint(l));
					//check if crosses obstacle
					for (int w = 0; w < cs_obsts.size(); w++){
						cs_obsts[w].Crosses(newline, &intersects);
						std::cout << "inter " << intersects.size() << endl;
						if (intersects.size()>0){
							intersects.DeleteAll();
							crossflag = true;
							break;
						}
						intersects.DeleteAll();

					}
					if (crossflag == false)
						visLines.push_back(newline);
					crossflag = false;

				
				
				}
			}
		
		}
	}

	*/
}

std::vector<C2DLine> RobotScene::path3DOF(C2DPoint start, C2DPoint end, std::vector<C2DPolygon> &cs_obsts){

	std::vector<C2DLine> cur_path;
	std::vector<C2DLine> path;


	
	MinkowskiSumLine(robot2, obstpolygons, cs_obsts, m_canvas);
	return path;

	


}

C2DPolygon Polygon_Union(C2DPolygon p1, C2DPolygon p2){

	
	C2DPolygon p3;

	C2DPointSet inters;
	C2DPointSet p1pts;
	C2DPointSet p2pts;
	p1.GetPointsCopy(p1pts);
	p2.GetPointsCopy(p2pts);

	C2DPointSet p3pts;

	C2DHoledPolyBaseSet holpol;

	p1.GetUnion(p2, holpol);

	p3=C2DPolygon(*holpol[0].GetRim());
	



	std::cout << "p1 " << p1pts.size();
	std::cout << "p2" << p2pts.size();
	std::cout << "polygon unionn" << inters.size();
	//p3 = C2DPolygon(inters, true);


	return p3;
	



}




int main(int argc, char* argv[])
{
	try
	{
		return vvr::mainLoop(argc, argv, new RobotScene);
	}
	catch (std::string exc)
	{
		cerr << exc << endl;
		return 1;
	}
	catch (...)
	{
		std::cerr << "Unknown exception" << std::endl;
		return 1;
	}
}
