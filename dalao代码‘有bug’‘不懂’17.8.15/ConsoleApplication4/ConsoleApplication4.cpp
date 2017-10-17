#include <igl/readOFF.h>
#include <igl/writeOBJ.h>
#include <igl/jet.h>
#include <igl/viewer/Viewer.h>
#include "tutorial_shared_path.h"
#include "ray_mesh_intersect.h"
#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <stdlib.h>
#include <algorithm>
#include <queue>
#define INF 0x3f3f3f3f
#define pi 3.1415926535898

const double eps = 1e-8;
typedef Eigen::Vector3d Point;

double MINY = INF;
double arf;

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd N_faces;
Eigen::VectorXd Z;
Eigen::MatrixXd C;
Eigen::Vector3i dir;
Point insPoint;
std::set<int> treset;
std::set<int> pset;
std::set<int> P;

//define the point struct with coordinates to sort
struct pPoint {

	double x, y, z;
	bool operator < (const pPoint & b) const {
		return y > b.y;
	}

	pPoint(double b, double c, double d) {
		x = b;
		y = c;
		z = d;
	}

	pPoint(Point a) {
		x = a(0);
		y = a(1);
		z = a(2);
	}

	double getDis(pPoint bb, pPoint & ss) {
		double x1 = x, x2 = bb.x, y1 = y, y2 = bb.y, z1 = z, z2 = bb.z;
		double tana = std::tan(arf);
		double s = std::sqrt((x1 - x2)*(x1 - x2) + (z1 - z2)*(z1 - z2));
		/*double a = (s + (y1 - y2)*std::tan(arf)) / 2;
		double b = (s + (y2 - y1)*std::tan(arf)) / 2;
		double y0 = ((y1 + y2)*std::tan(arf) - s) / (2 * std::tan(arf));
		double xx = x1 + a / (a + b)*(x2 - x1);
		double zz = z1 + a / (a + b)*(z2 - z1);
		ss.x = xx, ss.y = y0, ss.z = zz;*/
		double y0 = (y1 + y2) / 2 - s / (2 * std::tan(arf));
		double b = (y2 - y0)*tana;
		double a = s - b;
		double x0, z0;
		x0 = (a*x2 + b*x1) / (a + b);
		z0 = (a*z2 + b*z1) / (a + b);
		ss.x = x0;
		ss.z = z0;
		return y1 - y0;
	}
};

void addEdge(Point & a, Point & b, igl::viewer::Viewer & viewer) {
	//std::cout << "add edge here" << std::endl;
	viewer.data.add_edges(a, b, Eigen::RowVector3d(1, 0, 0));
}

bool isIntersect(Point m, Eigen::Vector3d vl, Point p, Point & O) {
	double v1 = vl(0), v2 = vl(1), v3 = vl(2);
	double a = p(0), b = p(1), c = p(2), u = m(0), v = m(1), w = m(2);
	if (v1 == 0) return false;
	else {
		double t = (v2*(v - b) - v3*(w - c) - v1*(u - a)) / v1;
		O(0) = u + t;
		O(1) = v;
		O(2) = w;
		if (t > 0)
			return true;
		return false;
	}
}

double triangleArea(Eigen::Vector3d A, Eigen::Vector3d B, Eigen::Vector3d C) {
	Eigen::Vector3d AB, AC, AP;
	AB = B - A;
	AC = C - A;
	AP(0) = AB(1)*AC(2) - AC(1)*AB(2);
	AP(1) = AC(0)*AB(2) - AB(0)*AC(2);
	AP(2) = AB(0)*AC(1) - AC(0)*AB(1);
	double m = sqrt(AP(0)*AP(0) + AP(1)*AP(1) + AP(2)*AP(2));
	return m / 2.0;
}

bool isInTriangle(Eigen::Vector3d tA, Eigen::Vector3d tB, Eigen::Vector3d tC, Eigen::Vector3d tP) {
	double a1 = triangleArea(tA, tB, tC);
	double a2 = triangleArea(tA, tB, tP);
	double a3 = triangleArea(tA, tP, tC);
	double a4 = triangleArea(tP, tB, tC);
	if (abs(a2 + a3 + a4 - a1) < eps) return true;
	return false;
}

void checkPoint() {
	std::vector<std::vector<int>> hash;
	for (int i = 0; i < V.rows(); i++) {
		int index = i;
		if (treset.find(index) != treset.end()) continue;
		std::vector<int> temp;
		for (int j = 0; j < F.rows(); j++) {
			int a = F.row(j)(0);
			int b = F.row(j)(1);
			int c = F.row(j)(2);
			if (a == index) temp.push_back(b), temp.push_back(c);
			else if (b == index) temp.push_back(a), temp.push_back(c);
			else if (c == index) temp.push_back(a), temp.push_back(b);
		}
		hash.push_back(temp);
	}
	for (auto it = hash.begin(); it != hash.end(); it++) {
		int point = it - hash.end();
		double minz = 0x3f3f3f3f;
		for (auto ait = (*it).begin(); ait != (*it).end(); ait++) {
			minz = std::min(V.row(*ait)(1), minz);
		}
		if (minz = V.row(point)(1)) pset.insert(point);
	}
	std::cout << pset.size() << std::endl;
}

void findPoint() {

	//right

	using namespace std;
	bool * vis = new bool[V.rows()];
	memset(vis, 1, sizeof(vis));

	for (int i = 0; i < N_faces.rows(); i++) {
		//if the traingle face down
		if (N_faces(i, 1) < 0) {
			//find the lowest point;
			double miny = min(V(F(i, 0), 1), min(V(F(i, 1), 1), V(F(i, 2), 1)));
			//if (fabs(miny - MINY) < eps) continue;

			/*//cout << V(F(i, 0), 1) << " " << V(F(i, 1), 1) << " " << V(F(i, 2), 1) << endl;
			for (int j = 0; j < 3; j++) {
			if (fabs(V(F(i, j), 1) - miny) > eps) {
			auto it = P.find(F(i, j));
			if (it != P.end()) {
			P.erase(it);
			}
			}
			else
			{
			P.insert(F(i, j));
			}
			}*/

			for (int j = 0; j < 3; j++) {
				//if the current point is the lowest,erase it
				if (fabs(V(F(i, j), 1) - miny) > eps) {
					if (vis[F(i, j)]) {
						vis[F(i, j)] = 0;
					}
				}
			}
			//the point that has not been erased is the lowest point global or local
			for (int i = 0; i < V.rows(); i++) {
				if (vis[i]) {
					P.insert(i);
				}
			}
		}
	}

	delete[] vis;

	cout << P.size() << endl;

}

void getMiny() {

	using namespace std;

	for (int i = 0; i < V.rows(); i++) {
		MINY = min(V(i, 1), MINY);
	}

}

void outputNormals() {
	//using namespace std;
	double dx, dy, dz;
	dx = dir(0), dy = dir(1), dz = dir(2);
	Z = N_faces.col(0);
	int rows = N_faces.rows();
	//cout << rows << endl;
	//int cnt = 0;
	for (int i = 0; i < rows; i++) {
		double ax, ay, az;
		ax = N_faces.row(i)(0);
		ay = N_faces.row(i)(1);
		az = N_faces.row(i)(2);
		double cs = ay*ay / ax*ax + ay*ay + az*az;
		if (ay<0 && cs > sin(arf)*sin(arf)) {
			Z(i) = 0.7;
			//std::cout << az << std::endl;
			treset.insert(i);
		}
		else
		{
			Z(i) = 0.002;
		}
	}
	//cout << cnt << endl; 
	//cout << Z << endl;
}

void init() {

	using namespace std;

	double lx = INF, ly = INF, lz = INF, hx = -INF, hy = -INF, hz = -INF;
	for (int i = 0; i < V.rows(); i++) {
		lx = min(V(i, 0), lx);
		ly = min(V(i, 1), ly);
		lz = min(V(i, 2), lz);
		hx = max(V(i, 0), hx);
		hy = max(V(i, 1), hy);
		hz = max(V(i, 2), hz);
	}

	insPoint = Eigen::Vector3d::Random();
	bool flag = false;
	Point O;
	int cnt;
	while (!flag) {
		//cout << O << endl;
		insPoint.setRandom();
		insPoint(0) = abs(insPoint(0))*hx;
		insPoint(1) = abs(insPoint(1))*hy;
		insPoint(2) = abs(insPoint(2))*hz;
		cnt = 0;
		for (int i = 0; i < F.rows(); i++) {
			if (isIntersect(insPoint, N_faces.row(i), V.row(F(i, 0)), O) && isInTriangle(V.row(F(i, 0)), V.row(F(i, 1)), V.row(F(i, 2)), O)) {
				cnt++;
			}
		}
		if (cnt % 2 == 1) {
			cout << "The coordinate of the point inside the MESH is " << O << endl;
			return;
		}
	}
}

double closetPoint(Point & source) {

	srand((unsigned)time(NULL));
	//set the rand times
	int times = 1000;

	double xz, xz_y;
	Point direc;
	igl::Hit hit;
	igl::Hit temp;
	bool ans = false;
	while (times--) {
		//generate an angle on the plane of x and z
		xz = rand() / double(RAND_MAX) * 2 * pi;
		//generate an angle between the ray and the y axis
		xz_y = rand() / double(RAND_MAX) * arf;
		direc(0) = cos(xz);
		//face down
		direc(1) = -sin(xz_y);
		direc(2) = sin(xz);
		if (igl::ray_mesh_intersect(source, direc, V, F, temp)) {
			if (!ans || temp.t < hit.t) {
				hit = temp;
				ans = 1;
			}
		}

		//std::cout << "stop here " << times << std::endl;
	}

	if (ans) {
		return hit.t;
	}
	else {
		return source(1) - MINY;
	}

}

void solve(igl::viewer::Viewer & viewer) {

	using namespace std;

	set<pPoint> q;


	for (auto it = treset.begin(); it != treset.end(); it++) {
		P.insert(F(*it, 0));
		P.insert(F(*it, 1));
		P.insert(F(*it, 2));
	}

	for (auto it = P.begin(); it != P.end(); it++) {
		//viewer.data.add_points(V.row(*it), Eigen::RowVector3d(1, 0, 0));
		q.insert(pPoint(V.row(*it)));
	}

	while (!q.empty()) {
		pPoint p = pPoint((*q.begin()).x, (*q.begin()).y, (*q.begin()).z);
		//c is the intersection of the cones of cp and cj;
		pPoint c(0, 0, 0), tmp(0, 0, 0);
		q.erase(q.begin());
		double dis = INF;
		auto minin = q.begin();
		for (auto j = q.begin(); j != q.end(); j++) {
			double temp = p.getDis(*j, tmp);
			//printf("(%.lf, %.lf, %.lf)----------(%.lf,%.lf,%.lf),(%.lf,%.lf,%.lf)\n", p.x, p.y, p.z, (*j).x,(*j).y,(*j).z , tmp.x, tmp.y, tmp.z);
			if (temp < dis) {
				minin = j;
				dis = temp;
				c = tmp;
				cout << "c" << c.x << " " << c.y << " " << c.z << endl;
			}
		}
		if (minin == q.end()) break;
		pPoint ap = { (*minin).x,(*minin).y,(*minin).z };
		cout << "c" << c.x << " " << c.y << " " << c.z << endl;
		q.erase(minin);
		double s = 0;
		//s is the distance between the point which is the intersection of p and mesh,means m in the article, and the source p;
		s = closetPoint(Point(p.x, p.y, p.z));
		//if the intersection c is closer than the intersection of the source and the mesh,insert the merged point c into the set;
		//besides,add an adge between p and c,ap and c;IMA
		if (dis<s) {
			q.insert(c);
			cout << c.y << endl;
			viewer.data.add_points(Point(c.x, c.y, c.z), Eigen::RowVector3d(1, 0, 0));
			//cout << c.x << " " << c.y << " " << c.z << endl;
			//addEdge(Point(p.x, p.y, p.z), Point(c.x, c.y, c.z), viewer);
			//addEdge(Point(ap.x, ap.y, ap.z), Point(c.x, c.y, c.z), viewer);	
		}

		std::cout << "stop here for solve(),q.size()==" << q.size() << std::endl;
	}
	std::cout << "solve finished here" << std::endl;
}

int main(int argc, char *argv[])
{
	// Load a mesh in OFF format
	igl::readOFF(TUTORIAL_SHARED_PATH "/bunny.off", V, F);
	//igl::writeOBJ("cube.obj", V, F);
	igl::per_face_normals(V, F, N_faces);
	dir << 0, 1, 0;

	arf = pi / 6;

	//init();
	outputNormals();

	// Plot the mesh
	igl::viewer::Viewer viewer;

	viewer.data.set_mesh(V, F);

	igl::jet(Z, true, C);

	viewer.data.set_colors(C);

	getMiny();
	//checkPoint();
	findPoint();

	solve(viewer);

	/*for (auto it = P.begin(); it != P.end(); it++) {
	viewer.data.add_points(V.row(*it), Eigen::RowVector3d(1, 0, 0));
	}*/

	/*Eigen::MatrixXd tmp;
	tmp << 0, 0, 1;
	viewer.data.add_points(tmp, Eigen::RowVector3d(1, 0, 0));

	std::cout << "-------------------------------------------------------" << std::endl;*/
	std::cout << "solve finished" << std::endl;
	viewer.launch();
}