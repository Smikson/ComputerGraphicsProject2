#include <windows.h>
#include <Windows.h>
#include <gl/Gl.h>
#include <gl/glu.h>
#include <gl/glut.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <istream>
#include <cstdlib>
#include <string>

using namespace std;


// Function to take the max of two numbers
//<<<<<<<<<<<<<<<<<<<<<<< max >>>>>>>>>>>>>>>>>>>>
double max(double a, double b) {
	return a > b? a : b;
}

// Create Point/Vector struct (simply to hold 3d x/y/z)
//<<<<<<<<<<<<<<<<<<<<<<< Point >>>>>>>>>>>>>>>>>>>>
struct Point {
	double x;
	double y;
	double z;
};

struct PointList {
	Point pts[100];
	int size;
};

// Set up a point as a copy of another point
//<<<<<<<<<<<<<<<<<<<<<<< copy >>>>>>>>>>>>>>>>>>>>
void copy(Point &dest, Point ori) {
	dest.x = ori.x;
	dest.y = ori.y;
	dest.z = ori.z;
}

// Does the dot product of two point/vector variables
//<<<<<<<<<<<<<<<<<<<<<<< dot >>>>>>>>>>>>>>>>>>>>
double dot(Point p1, Point p2) {
	return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}

// Does the cross product of two point/vector variables
//<<<<<<<<<<<<<<<<<<<<<<< dot >>>>>>>>>>>>>>>>>>>>
Point cross(Point p1, Point p2) {
	// Calculate the cross product
	Point ret;
	ret.x = p1.y * p2.z - p1.z * p2.y;
	ret.y = p1.z * p2.x - p1.x * p2.z;
	ret.z = p1.x * p2.y - p1.y * p2.x;
	
	// Return the result
	return ret;
}

// Calculates the magnitude of a point/vector
//<<<<<<<<<<<<<<<<<<<<<<< magnitude >>>>>>>>>>>>>>>>>>>>
double magnitude(Point p) {
	return sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
}

// Normalizes a point/vector
//<<<<<<<<<<<<<<<<<<<<<<< normalize >>>>>>>>>>>>>>>>>>>>
void normalize(Point &p) {
	// Calculate the magnitude of the point/vector
	double mag = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
	
	// If the magnitude is 0, do nothing
	if (mag == 0) {
		return;
	}
	
	// Otherwise, divide each component by the magnitude
	p.x /= 1.0*mag;
	p.y /= 1.0*mag;
	p.z /= 1.0*mag;
}

// Overload the addition and subtraction for points/vectors to do operations component-wise
Point operator-(Point p1, Point p2) {
	Point res;
	res.x = p1.x - p2.x;
	res.y = p1.y - p2.y;
	res.z = p1.z - p2.z;
	return res;
}

Point operator+(Point p1, Point p2) {
	Point res;
	res.x = p1.x + p2.x;
	res.y = p1.y + p2.y;
	res.z = p1.z + p2.z;
	return res;
}

// Camera Class
//<<<<<<<<<<<<<<<<<<<<<<< Camera >>>>>>>>>>>>>>>>>>>>
class Camera {
	private:
		Point eye, look, up, u, v, n;
		double viewAngle, aspect, nearDist, farDist; // view volume shape
		void setModelviewMatrix(); // tell OpenGL where the camera is
	
	public:
		Camera(); // constructor
		void set(Point Eye, Point Look, Point Up); // like gluLookAt()
		void roll(double angle); // roll it
		void pitch(double angle); // increase pitch
		void yaw(double angle); // yaw it
		void slide(double delU, double delV, double delN); // slide it
		void setShape(double vAng, double asp, double nearD, double farD);
		void getShape(double &vAng, double &asp, double &nearD, double &farD);
		Point getEye();
};

// Base constructor, set all values to 0 (not a useful camera)
Camera::Camera() {
	eye.x = 0; eye.y = 0; eye.z = 0;
	look.x = 0; look.y = 0; look.z = 0;
	up.x = 0; up.y = 0; up.z = 0;
	u.x = 0; u.y = 0; u.z = 0;
	v.x = 0; v.y = 0; v.z = 0;
	n.x = 0; n.y = 0; n.z = 0;
	viewAngle = 0;
	aspect = 0;
	nearDist = 0;
	farDist = 0;
}

// Implementation of functions for camera class
// setModelviewMatrx (tell OpenGL where the camera is)
void Camera::setModelviewMatrix() {
	// Load modelview matrix with existing camera values
	double m[16];
	
	// construct a 3d vector for the eye of the camera
	Point eVec;
	eVec.x = eye.x; eVec.y = eye.y; eVec.z = eye.z;
	m[0] = u.x; m[4] = u.y; m[8] = u.z; m[12] = -dot(eVec, u);
	m[1] = v.x; m[5] = v.y; m[9] = v.z; m[13] = -dot(eVec, v);
	m[2] = n.x; m[6] = n.y; m[10] = n.z;m[14] = -dot(eVec, n);
	m[3] = 0;	m[7] = 0;	m[11] = 0;	m[15] = 1.0;
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixd(m); // Load OpenGL's modelview matrix
}

// set (functions like a gluLookAt(), sets up the camera)
void Camera::set(Point Eye, Point Look, Point Up) {
	// Create a modelview matrix and send it to OpenGL
	copy(eye, Eye); copy(look, Look); copy(up, Up);
	
	// Make n, normalize it
	n.x = eye.x - look.x; n.y = eye.y - look.y; n.z = eye.z - look.z;
	normalize(n);
	
	// make u = up X n, normalize it
	Point upCrossN;
	upCrossN = cross(up, n);
	copy(u, upCrossN);
	normalize(u);
	
	// make v = n X u
	Point nCrossU;
	nCrossU = cross(n, u);
	copy(v, nCrossU);
	
	// Tell OpenGL
	setModelviewMatrix();
}

// setShape (sets the private variables and utilizes gluPerspective to set up the size/shape for where the camera is looking)
void Camera::setShape(double vAng, double asp, double nearD, double farD) {
	// Set the variables
	viewAngle = vAng; aspect = asp; nearDist = nearD; farDist = farD;
	
	// pass to gluPerspective
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(viewAngle, aspect, nearDist, farDist);
}

// getShape - sets the input variables equal to the class variables
void Camera::getShape(double &vAng, double &asp, double &nearD, double &farD) {
	vAng = viewAngle; asp = aspect; nearD = nearDist; farD = farDist;
}

// slide (moves camera without rotating it by altering u, v, or n for forward/back, up/down, or left/right respectively)
void Camera::slide(double delU, double delV, double delN) {
	// Adjust eye and look by the specified amounts
	eye.x += delU * u.x + delV * v.x + delN * n.x;
	eye.y += delU * u.y + delV * v.y + delN * n.y;
	eye.z += delU * u.z + delV * v.z + delN * n.z;
	look.x += delU * u.x + delV * v.x + delN * n.x;
	look.y += delU * u.y + delV * v.y + delN * n.y;
	look.z += delU * u.z + delV * v.z + delN * n.z;
	setModelviewMatrix();
}

// roll (rotates camera along its "n axis")
void Camera::roll(double angle) {
	// roll the camera through angle degrees
	double cs = cos(3.14159265/180.0 * angle);
	double sn = sin(3.14159265/180.0 * angle);
	
	// Copy u, store temporarily in t, then perform calculations
	Point t;
	copy(t, u);
	u.x = cs*t.x - sn*v.x; u.y = cs*t.y - sn*v.y; u.z = cs*t.z - sn*v.z;
	v.x = sn*t.x + cs*v.x; v.y = sn*t.y + cs*v.y; v.z = sn*t.z + cs*v.z;
	setModelviewMatrix();
}

// pitch (rotates camera along its "u axis")
void Camera::pitch(double angle) {
	// pitch the camera through angle degrees
	double cs = cos(3.14159265/180.0 * angle);
	double sn = sin(3.14159265/180.0 * angle);
	
	// Copy u, store temporarily in t, then perform calculations
	Point t;
	copy(t, n);
	n.x = cs*t.x - sn*v.x; n.y = cs*t.y - sn*v.y; n.z = cs*t.z - sn*v.z;
	v.x = sn*t.x + cs*v.x; v.y = sn*t.y + cs*v.y; v.z = sn*t.z + cs*v.z;
	setModelviewMatrix();
}

// yaw (rotates camera along its "v axis")
void Camera::yaw(double angle) {
	// roll the camera through angle degrees
	double cs = cos(3.14159265/180.0 * angle);
	double sn = sin(3.14159265/180.0 * angle);
	
	// Copy u, store temporarily in t, then perform calculations
	Point t;
	copy(t, u);
	u.x = cs*t.x - sn*n.x; u.y = cs*t.y - sn*n.y; u.z = cs*t.z - sn*n.z;
	n.x = sn*t.x + cs*n.x; n.y = sn*t.y + cs*n.y; n.z = sn*t.z + cs*n.z;
	setModelviewMatrix();
}

// Returns a copy of the eye variable for reference of where the camera is
Point Camera::getEye() {
	Point ret;
	copy(ret, eye);
	return ret;
}


// Create a global camera object to be used by the program
Camera cam;
PointList dino;

// Sets the initial dino coordinates
//<<<<<<<<<<<<<<<<<<<<<<<< setDino >>>>>>>>>>>>>>>>>
void setDino() {
	// Dino consists of 9 8-vertex rectangular prisms according to the long list of points
	dino.size = 72;
	double xs[] = {1,2,2,1,1,2,2,1,1.2,1.3,1.3,1.2,1.2,1.3,1.3,1.2,1.7,1.8,1.8,1.7,1.7,1.8,1.8,1.7,1.2,1.3,1.3,1.2,1.2,1.3,1.3,1.2,1.7,1.8,1.8,1.7,1.7,1.8,1.8,1.7,1.45,1.55,1.55,1.45,1.45,1.55,1.55,1.45,1.3,1.7,1.7,1.3,1.3,1.7,1.7,1.3,1.35,1.45,1.45,1.35,1.35,1.45,1.45,1.35,1.55,1.65,1.65,1.55,1.55,1.65,1.65,1.55};
	double ys[] = {1,1,2,2,1,1,2,2,.25,.25,1,1,.25,.25,1,1,.25,.25,1,1,.25,.25,1,1,.25,.25,1,1,.25,.25,1,1,.25,.25,1,1,.25,.25,1,1,2,2,2.5,2.5,2,2,2.5,2.5,2.5,2.5,2.85,2.85,2.5,2.5,2.85,2.85,2.65,2.65,2.75,2.75,2.65,2.65,2.75,2.75,2.65,2.65,2.75,2.75,2.65,2.65,2.75,2.75};
	double zs[] = {.5,.5,.5,.5,2.5,2.5,2.5,2.5,.9,.9,.9,.9,1,1,1,1,.9,.9,.9,.9,1,1,1,1,2,2,2,2,2.1,2.1,2.1,2.1,2,2,2,2,2.1,2.1,2.1,2.1,2.2,2.2,2.2,2.2,2.35,2.35,2.35,2.35,2.15,2.15,2.15,2.15,2.75,2.75,2.75,2.75,2.75,2.75,2.75,2.75,2.8,2.8,2.8,2.8,2.75,2.75,2.75,2.75,2.8,2.8,2.8,2.8};
	for (int i = 0; i < 72; i++) {
		Point p;
		p.x = xs[i];
		p.y = ys[i];
		p.z = zs[i];
		
		dino.pts[i] = p;
	}
}

// Move dino left/right/front/back
//<<<<<<<<<<<<<<<<<<<<<<<< translateDino >>>>>>>>>>>>>>>>>
void translateDino(double left, double right, double front, double back) {
	// Adjust all points in the dino according to what is given
	for (int i = 0; i < dino.size; i++) {
		dino.pts[i].x += right - left;
		dino.pts[i].z += front - back;
	}
}

// Arrows function to interactively move camera
//<<<<<<<<<<<<<<<<<<<<<<<< arrows >>>>>>>>>>>>>>>>>
void arrows(int key, int x, int y) {
	switch(key) {
		// Arrow Controls for camera
		case GLUT_KEY_UP: cam.slide(0, 0, -0.2); break;
		case GLUT_KEY_DOWN: cam.slide(0, 0, 0.2); break;
		case GLUT_KEY_LEFT: cam.slide(-0.2, 0, 0); break;
		case GLUT_KEY_RIGHT: cam.slide(0.2, 0, 0); break;
	}
	
	glutPostRedisplay();
}

// Keyboard function to interactively move camera
//<<<<<<<<<<<<<<<<<<<<<<<< moveCam >>>>>>>>>>>>>>>>>
void moveCam(unsigned char key, int x, int y) {
	switch(key) {
		// Controls for camera
		case 'q': cam.roll(-1.0); break;
		case 'w': cam.roll(1.0); break;
		case 'e': cam.pitch(1.0); break;
		case 'd': cam.pitch(-1.0); break;
		case 'a': cam.yaw(1.0); break;
		case 's': cam.yaw(-1.0); break;
		// Controls for dino, move left, right, up, down
		case 'j': translateDino(-0.2,0,0,0); break;
		case 'l': translateDino(0,-0.2,0,0); break;
		case 'i': translateDino(0,0,0.2,0); break;
		case 'k': translateDino(0,0,0,0.2); break;
	}
	
	glutPostRedisplay();
}

void setLightColor(PointList plane, double baseRed, double baseGreen, double baseBlue) {
	// Constants from the book
	double par = 0.24725;
	double pag = 0.1995;
	double pab = 0.0745;
	
	double pdr = 0.75164;
	double pdg = 0.60648;
	double pdb = 0.22648;
	
	double psr = 0.628281;
	double psg = 0.555802;
	double psb = 0.366065;
	
	double f = 51.2;
	
	// Keep a running sum of the intensity for each vertex (to be averaged)
	double Ir = 0;
	double Ig = 0;
	double Ib = 0;
	for (int i = 0; i < plane.size; i++) {
		// Calculate the intensity of the current vertex
		double Is = 1.0;
		Point eye = cam.getEye();
		Point vertex;
		vertex.x = plane.pts[i].x; vertex.y = plane.pts[i].y; vertex.z = plane.pts[i].z;
		Point prevVertex;
		prevVertex.x = plane.pts[(i-1)%plane.size].x; prevVertex.y = plane.pts[(i-1)%plane.size].y; prevVertex.z = plane.pts[(i-1)%plane.size].z;
		Point nextVertex;
		nextVertex.x = plane.pts[(i+1)%plane.size].x; nextVertex.y = plane.pts[(i+1)%plane.size].y; nextVertex.z = plane.pts[(i+1)%plane.size].z;
		Point light;
		light.x = 2; light.y = 10; light.z = 8;
		
		Point v = eye - vertex;
		Point s = light - vertex;
		Point a = prevVertex - vertex;
		Point b = nextVertex - vertex;
		Point m1 = cross(a, b);
		Point m2 = cross(b, a);
		
		// Calculate lambert
		double lambert = 0;
		if (dot(m1, v) > 0) {
			lambert = max(0, dot(m1, s) / (magnitude(s) * magnitude(m1)));
		}
		if (dot(m2, v) > 0) {
			lambert = max(0, dot(m2, s) / (magnitude(s) * magnitude(m2)));
		}
		
		// Calculate phong
		Point h = s + v;
		double phong = 0;
		if (dot(m1, h) > 0) {
			phong = max(0, dot(m1, h) / (magnitude(h) * magnitude(m1)));
		}
		if (dot(m2, h) > 0) {
			phong = max(0, dot(m2, h) / (magnitude(h) * magnitude(m2)));
		}
		
		// Put them all together, add to sum
		Ir += Is*par + Is*pdr*lambert + Is*psr*pow(phong, f);
		Ig += Is*pag + Is*pdg*lambert + Is*psg*pow(phong, f);
		Ib += Is*pab + Is*pdb*lambert + Is*psb*pow(phong, f);
	}
	
	// Set the color to the average r/g/b calculated (divide by number of vertexes in plane)
	glColor3d(Ir/plane.size * baseRed, Ig/plane.size * baseGreen, Ib/plane.size * baseBlue);
}

// Draws a rectangular prism with lines (assumes stored as first 'n' make the "base", next 'n' make the "top")
//<<<<<<<<<<<<<<<<<<<<<<<< drawPrismLines >>>>>>>>>>>>>>>>>
void drawPrismLines(PointList prism, int n) {
	// Draw the "base" plane
	glBegin(GL_LINES);
	for (int i = 0; i < n; i++) {
		glVertex3f(prism.pts[i].x, prism.pts[i].y, prism.pts[i].z);
		glVertex3f(prism.pts[(i+1)%n].x, prism.pts[(i+1)%n].y, prism.pts[(i+1)%n].z);
	}
	
	// Draw the "top" plane
	for (int i = n; i < 2*n; i++) {
		glVertex3f(prism.pts[i].x, prism.pts[i].y, prism.pts[i].z);
		glVertex3f(prism.pts[(i+1)%n + n].x, prism.pts[(i+1)%n + n].y, prism.pts[(i+1)%n + n].z);
	}
	
	// Draw each line between the top and the bottom: the side "walls"
	for (int i = 0; i < n; i++) {
		glVertex3f(prism.pts[i].x, prism.pts[i].y, prism.pts[i].z);
		glVertex3f(prism.pts[i+n].x, prism.pts[i+n].y, prism.pts[i+n].z);
	}
	glEnd();
}

// Draws a filled rectangular prism (assumes stored as first 'n' make the "base", next 'n' make the "top")
//<<<<<<<<<<<<<<<<<<<<<<<< drawPrismSolid >>>>>>>>>>>>>>>>>
void drawPrismSolid(PointList prism, int n, double red, double green, double blue) {
	// Draw the "base" plane
	// Create a PointList to set the color
	PointList base;
	base.size = n;
	for (int i = 0; i < n; i++) {
		Point p;
		p.x = prism.pts[i].x;
		p.y = prism.pts[i].y;
		p.z = prism.pts[i].z;
		base.pts[i] = p;
	}
	setLightColor(base, red, green, blue);
	
	// Draw the plane
	glBegin(GL_POLYGON);
	for (int i = 0; i < n; i++) {
		glVertex3f(prism.pts[i].x, prism.pts[i].y, prism.pts[i].z);
	}
	glEnd();
	
	// Draw the "top" plane
	// Create a PointList to set the color
	PointList top;
	top.size = n;
	for (int i = n; i < 2*n; i++) {
		Point p;
		p.x = prism.pts[i].x;
		p.y = prism.pts[i].y;
		p.z = prism.pts[i].z;
		top.pts[i-n] = p;
	}
	setLightColor(top, red, green, blue);
	
	// Draw the plane
	glBegin(GL_POLYGON);
	for (int i = n; i < 2*n; i++) {
		glVertex3f(prism.pts[i].x, prism.pts[i].y, prism.pts[i].z);
	}
	glEnd();
	
	// Draw each plane between all the points: the side "walls"
	for (int i = 0; i < n; i++) {
		// Create a PointList for the wall so that it can be passed to the light function
		PointList wall;
		wall.size = 4;
		Point p1; p1.x = prism.pts[i].x; p1.y = prism.pts[i].y; p1.z = prism.pts[i].z; wall.pts[0] = p1;
		Point p2; p2.x = prism.pts[(i + 1) % n].x; p2.y = prism.pts[(i + 1) % n].y; p2.z = prism.pts[((i + 1) % n) + n].z; wall.pts[1] = p2;
		Point p3; p3.x = prism.pts[((i + 1) % n) + n].x; p3.y = prism.pts[((i + 1) % n) + n].y; p3.z = prism.pts[((i + 1) % n) + n].z; wall.pts[2] = p3;
		Point p4; p4.x = prism.pts[i+n].x; p4.y = prism.pts[i+n].y; p4.z = prism.pts[i+n].z; wall.pts[3] = p4;
		setLightColor(wall, red, green, blue);
		
		// Draw the plane look at the same points from each half of the stored prism
		glBegin(GL_POLYGON);
		glVertex3f(prism.pts[i].x, prism.pts[i].y, prism.pts[i].z);
		glVertex3f(prism.pts[(i + 1) % n].x, prism.pts[(i + 1) % n].y, prism.pts[(i + 1) % n].z);
		glVertex3f(prism.pts[((i + 1) % n) + n].x, prism.pts[((i + 1) % n) + n].y, prism.pts[((i + 1) % n) + n].z);
		glVertex3f(prism.pts[i+n].x, prism.pts[i+n].y, prism.pts[i+n].z);
		glEnd();
	}
}

// Draws the dino based on what is stored in the global variable (lines to form wire-like surface commented out, can be used for debugging)
//<<<<<<<<<<<<<<<<<<<<<<<< drawDino >>>>>>>>>>>>>>>>>
void drawDino() {
	// Draw all 9 rectangular prisms with the correct color
	// The body
	PointList body;
	body.size = 8;
	for (int i = 0; i < 8; i++) {
		body.pts[i] = dino.pts[i];
	}
	
	drawPrismSolid(body, 4, 0, 1, 0);
	//glColor3f(0.0f, 0.0f, 0.0f);
	//drawPrismLines(body, 4);
	
	// The 4 legs
	PointList leg1;
	leg1.size = 8;
	for (int i = 0; i < 8; i++) {
		leg1.pts[i] = dino.pts[i+8];
	}
	
	drawPrismSolid(leg1, 4, 0, 0, 1);
	//glColor3f(0.0f, 0.0f, 0.0f);
	//drawPrismLines(leg1, 4);
	
	PointList leg2;
	leg2.size = 8;
	for (int i = 0; i < 8; i++) {
		leg2.pts[i] = dino.pts[i+16];
	}
	
	drawPrismSolid(leg2, 4, 0, 0, 1);
	//glColor3f(0.0f, 0.0f, 0.0f);
	//drawPrismLines(leg2, 4);
	
	PointList leg3;
	leg3.size = 8;
	for (int i = 0; i < 8; i++) {
		leg3.pts[i] = dino.pts[i+24];
	}
	
	drawPrismSolid(leg3, 4, 0, 0, 1);
	//glColor3f(0.0f, 0.0f, 0.0f);
	//drawPrismLines(leg3, 4);
	
	PointList leg4;
	leg4.size = 8;
	for (int i = 0; i < 8; i++) {
		leg4.pts[i] = dino.pts[i+32];
	}
	
	drawPrismSolid(leg4, 4, 0, 0, 1);
	//glColor3f(0.0f, 0.0f, 0.0f);
	//drawPrismLines(leg4, 4);
	
	// The neck
	PointList neck;
	neck.size = 8;
	for (int i = 0; i < 8; i++) {
		neck.pts[i] = dino.pts[i+40];
	}
	
	drawPrismSolid(neck, 4, 0, 0, 1);
	//glColor3f(0.0f, 0.0f, 0.0f);
	//drawPrismLines(neck, 4);
	
	// The head
	PointList head;
	head.size = 8;
	for (int i = 0; i < 8; i++) {
		head.pts[i] = dino.pts[i+48];
	}
	
	drawPrismSolid(head, 4, 0, 1, 0);
	//glColor3f(0.0f, 0.0f, 0.0f);
	//drawPrismLines(head, 4);
	
	// The eyes
	PointList eye1;
	eye1.size = 8;
	for (int i = 0; i < 8; i++) {
		eye1.pts[i] = dino.pts[i+56];
	}
	
	drawPrismSolid(eye1, 4, 1, 0, 0);
	//glColor3f(0.0f, 0.0f, 0.0f);
	//drawPrismLines(eye1, 4);
	
	PointList eye2;
	eye2.size = 8;
	for (int i = 0; i < 8; i++) {
		eye2.pts[i] = dino.pts[i+64];
	}
	
	drawPrismSolid(eye2, 4, 1, 0, 0);
	//glColor3f(0.0f, 0.0f, 0.0f);
	//drawPrismLines(eye2, 4);
}

// Draws a tree passed. Assumed to have 32 points. (lines to form wire-like surface commented out, can be used for debugging)
void drawTree(PointList tree) {
	// First layer (the bark)
	PointList layer1;
	layer1.size = 8;
	for (int i = 0; i < 8; i++) {
		layer1.pts[i] = tree.pts[i];
	}
	
	drawPrismSolid(layer1, 4, 0.71, 0.396, 0.114);
	//glColor3f(0.0f, 0.0f, 0.0f);
	//drawPrismLines(layer1, 4);
	
	// Second layer (.75 out)
	PointList layer2;
	layer2.size = 8;
	for (int i = 0; i < 8; i++) {
		layer2.pts[i] = tree.pts[i+8];
	}
	
	drawPrismSolid(layer2, 4, 0, 1, 0);
	//glColor3f(0.0f, 0.0f, 0.0f);
	//drawPrismLines(layer2, 4);
	
	// Third layer (1.5 out)
	PointList layer3;
	layer3.size = 8;
	for (int i = 0; i < 8; i++) {
		layer3.pts[i] = tree.pts[i+16];
	}
	
	drawPrismSolid(layer3, 4, 0, 1, 0);
	//glColor3f(0.0f, 0.0f, 0.0f);
	//drawPrismLines(layer3, 4);
	
	// Fourth layer (.75 out)
	PointList layer4;
	layer4.size = 8;
	for (int i = 0; i < 8; i++) {
		layer4.pts[i] = tree.pts[i+24];
	}
	
	drawPrismSolid(layer4, 4, 0, 1, 0);
	//glColor3f(0.0f, 0.0f, 0.0f);
	//drawPrismLines(layer4, 4);
}

// Draws the scene
//<<<<<<<<<<<<<<<<<<<<<<<< drawScene >>>>>>>>>>>>>>>>>
void drawScene(void) {
	// Clear the screen
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	
	// Draw the light as a big white dot so that if the user zooms out they can see it
	glColor3f(1.0f, 1.0f, 1.0f);
	glBegin(GL_POINTS);
	glVertex3f(2, 10, 8);
	glEnd();
	
	// Draw the dino at its current coordinates
	drawDino();
	
	// Create a tree
	PointList tree;
	tree.size = 32;
	
	// Do the first layer (bark)
	double xs[] = {1,1.75,1.75,1,1,1.75,1.75,1};
	double ys[] = {.25,.25,3.25,3.25,.25,.25,3.25,3.25};
	double zs[] = {4,4,4,4,4.75,4.75,4.75,4.75};
	
	PointList layer1;
	layer1.size = 8;
	for (int i = 0; i < 8; i++) {
		Point p;
		p.x = xs[i];
		p.y = ys[i];
		p.z = zs[i];
		
		layer1.pts[i] = p;
		tree.pts[i] = p;
	}
	
	// Do the second layer (.75 extra all directions in green)
	double xs2[] = {.25,2.5,2.5,.25,.25,2.5,2.5,.25};
	double ys2[] = {3.25,3.25,4,4,3.25,3.25,4,4};
	double zs2[] = {3.25,3.25,3.25,3.25,5.5,5.5,5.5,5.5};
	
	PointList layer2;
	layer2.size = 8;
	for (int i = 0; i < 8; i++) {
		Point p;
		p.x = xs2[i];
		p.y = ys2[i];
		p.z = zs2[i];
		
		layer2.pts[i] = p;
		tree.pts[i+8] = p;
	}
	
	// Do the third layer (1.5 extra all directions (except y which just keep going up by .75) in green)
	double xs3[] = {-.5,3.25,3.25,-.5,-.5,3.25,3.25,-.5};
	double ys3[] = {4,4,4.75,4.75,4,4,4.75,4.75};
	double zs3[] = {2.5,2.5,2.5,2.5,6.25,6.25,6.25,6.25};
	
	PointList layer3;
	layer3.size = 8;
	for (int i = 0; i < 8; i++) {
		Point p;
		p.x = xs3[i];
		p.y = ys3[i];
		p.z = zs3[i];
		
		layer3.pts[i] = p;
		tree.pts[i+16] = p;
	}
	
	// Do the fourth layer (.75 extra all directions in green)
	double xs4[] = {.25,2.5,2.5,.25,.25,2.5,2.5,.25};
	double ys4[] = {4.75,4.75,5.5,5.5,4.75,4.75,5.5,5.5};
	double zs4[] = {3.25,3.25,3.25,3.25,5.5,5.5,5.5,5.5};
	
	PointList layer4;
	layer4.size = 8;
	for (int i = 0; i < 8; i++) {
		Point p;
		p.x = xs4[i];
		p.y = ys4[i];
		p.z = zs4[i];
		
		layer4.pts[i] = p;
		tree.pts[i+24] = p;
	}
	
	// Draw the tree created
	drawTree(tree);
	
	// Copy and draw the tree to make another smaller tree farther away
	PointList tree2;
	tree2.size = 32;
	for (int i = 0; i < 32; i++) {
		Point p;
		p.x = tree.pts[i].x * .75 - 4;
		p.y = tree.pts[i].y * .75 + .25;
		p.z = tree.pts[i].z * .75 - 4;
		
		tree2.pts[i] = p;
	}
	drawTree(tree2);
	
	// Copy and draw the tree again to start behind the dino
	PointList tree3;
	tree3.size = 32;
	for (int i = 0; i < 32; i++) {
		Point p;
		p.x = tree.pts[i].x;
		p.y = tree.pts[i].y * .75 + .25;
		p.z = tree.pts[i].z * .75 - 4;
		
		tree3.pts[i] = p;
	}
	drawTree(tree3);
	
	// Send all output to display
	glFlush();
	
	// Display the screen just made
	glutSwapBuffers();
}

// Runs the code setting the GL functions to the appropriate from above
//<<<<<<<<<<<<<<<<<<<<<<<< main >>>>>>>>>>>>>>>>>>>>>>
int main(int argc, char **argv) {
	setDino();
	
	glutInit(&argc, argv);          // initialize the toolkit
	glutInitDisplayMode(GLUT_DEPTH | GLUT_RGB); // set the display mode
	glutInitWindowSize(640,480);     // set the window size
	glutInitWindowPosition(100, 150); // set the window position on the screen
	glutCreateWindow("Move the camera and dino around the scene!"); // open the screen window(with its exciting title)
	glutSpecialFunc(arrows);
	glutKeyboardFunc(moveCam);
	glutDisplayFunc(drawScene);     // register the redraw function
	glViewport(0,0,640,480);
	Point eye, look, up;
	eye.x = 12; eye.y = 5.5; eye.z = 12;
	look.x = 0; look.y = 0; look.z = 0;
	up.x = 0; up.y = 1; up.z = 0;
	cam.setShape(30.0, 64.0/48.0, 0.5, 100.0);
	cam.set(eye, look, up);
	glClearColor(0.1, 0.1, 0.1, 0.0);   // set the bg color to white
	glPointSize(16.0);		            // set the point size to 16 by 16 pixels
	glutMainLoop(); 		     // go into a perpetual loop
}


