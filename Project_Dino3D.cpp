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


// Create Point/Vector struct (simply to hold 3d x/y/z)
//<<<<<<<<<<<<<<<<<<<<<<< Point >>>>>>>>>>>>>>>>>>>>
struct Point {
	double x;
	double y;
	double z;
};

struct PointList {
	Point pts[50];
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


// Create a global camera object to be used by the program
Camera cam;

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
	}
	
	glutPostRedisplay();
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
void drawPrismSolid(PointList prism, int n) {
	// Draw the "base" plane
	glBegin(GL_POLYGON);
	for (int i = 0; i < n; i++) {
		glVertex3f(prism.pts[i].x, prism.pts[i].y, prism.pts[i].z);
	}
	glEnd();
	
	// Draw the "top" plane
	glBegin(GL_POLYGON);
	for (int i = n; i < 2*n; i++) {
		glVertex3f(prism.pts[i].x, prism.pts[i].y, prism.pts[i].z);
	}
	glEnd();
	
	// Draw each plane between all the points: the side "walls"
	for (int i = 0; i < n; i++) {
		glBegin(GL_POLYGON);
		
		// Draw the plane look at the same points from each half of the stored prism
		glVertex3f(prism.pts[i].x, prism.pts[i].y, prism.pts[i].z);
		glVertex3f(prism.pts[(i + 1) % n].x, prism.pts[(i + 1) % n].y, prism.pts[(i + 1) % n].z);
		glVertex3f(prism.pts[((i + 1) % n) + n].x, prism.pts[((i + 1) % n) + n].y, prism.pts[((i + 1) % n) + n].z);
		glVertex3f(prism.pts[i+n].x, prism.pts[i+n].y, prism.pts[i+n].z);
		
		glEnd();
	}
}

// Draws the dodecahedron
//<<<<<<<<<<<<<<<<<<<<<<<< draw >>>>>>>>>>>>>>>>>
void draw(void) {
	// Clear the screen
	glClear(GL_COLOR_BUFFER_BIT);
	
	// Draw Dodecahedron -- Change this to a dinosour!
	glColor3f(0.0f, 1.0f, 0.0f);
	//glutSolidDodecahedron();
	glColor3f(0.0f, 0.0f, 0.0f);
	//glutWireDodecahedron();
	
	// Create a test rectangular prism
	double xs[] = {1,2,2,1,1,2,2,1};
	double ys[] = {1,1,2,2,1,1,2,2};
	double zs[] = {.5,.5,.5,.5,2.5,2.5,2.5,2.5};
	
	PointList test;
	test.size = 8;
	for (int i = 0; i < 8; i++) {
		Point p;
		p.x = xs[i];
		p.y = ys[i];
		p.z = zs[i];
		
		test.pts[i] = p;
	}
	
	
	glColor3f(0.0f, 1.0f, 0.0f);
	drawPrismSolid(test, 4);
	glColor3f(0.0f, 0.0f, 0.0f);
	drawPrismLines(test, 4);
	
	// Send all output to display
	glFlush();
	
	// Display the screen just made
	glutSwapBuffers();
}

// Runs the code setting the GL functions to the appropriate from above
//<<<<<<<<<<<<<<<<<<<<<<<< main >>>>>>>>>>>>>>>>>>>>>>
int main(int argc, char **argv) {
	glutInit(&argc, argv);          // initialize the toolkit
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB); // set the display mode
	glutInitWindowSize(640,480);     // set the window size
	glutInitWindowPosition(100, 150); // set the window position on the screen
	glutCreateWindow("Move the camera and dino around the scene!"); // open the screen window(with its exciting title)
	glutSpecialFunc(arrows);
	glutKeyboardFunc(moveCam);
	glutDisplayFunc(draw);     // register the redraw function
	glViewport(0,0,640,480);
	Point eye, look, up;
	eye.x = 4; eye.y = 4; eye.z = 4;
	look.x = 0; look.y = 0; look.z = 0;
	up.x = 0; up.y = 1; up.z = 0;
	cam.setShape(30.0, 64.0/48.0, 0.5, 100.0);
	cam.set(eye, look, up);
	glClearColor(1.0, 1.0, 1.0, 0.0);   // set the bg color to white
	glutMainLoop(); 		     // go into a perpetual loop
}


