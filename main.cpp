#include <windows.h>   // use as needed for your system
#include <iostream>
#include <fstream>
#include <gl/Gl.h>
#include <gl/glu.h>
#include <gl/glut.h>
using namespace std;

//<<<<<<<<<<<<<<<<<<<<<<< myInit >>>>>>>>>>>>>>>>>>>>
 void myInit(void)
 {
    glClearColor(0, 0, 0, 0.0);      // set the bg color to a bright white
    glColor3f(255, 255, 255);           // set the drawing color to black 
 	glPointSize(8.0);		            //set the point size to 4 by 4 pixels
	glMatrixMode(GL_PROJECTION);// set up appropriate matrices- to be explained 
	glLoadIdentity();// to be explained
	gluOrtho2D(0.0, 640.0, 0.0, 480.0);// to be explained
}

//<<<<<<<<<<<<<<<<<<<<<<<< myDisplay >>>>>>>>>>>>>>>>>
// the redraw function
void myDisplay(void)
{
	fstream inStream;
	
	glClear(GL_COLOR_BUFFER_BIT); // clear the screen
	GLint numpolys, numLines, x ,y;
	
	for(int r = 0; r<5; r++)
	{
		for (int c = 0; c<5; c++)
		{
			inStream.open("dino.txt", ios ::in); // open the file
			if(inStream.fail())
				return;
			inStream >> numpolys; // read the number of polylines
			for(int j = 0; j < numpolys; j++) // read each polyline   
			{
				inStream >> numLines;
				glBegin(GL_LINE_STRIP); // draw the next polyline
				for (int i = 0; i < numLines; i++)
				{
					inStream >> x >> y; // read the next x, y pair
					x=(int)r*128+x/5;
					y=(int)c*96+y/5;
					glVertex2i(x, y);
				}
				glEnd();			
			}
			inStream.close();
		}
	}
	
	glFlush();
	
}
//<<<<<<<<<<<<<<<<<<<<<<<< main >>>>>>>>>>>>>>>>>>>>>>
int main(int argc, char **argv)
{
	glutInit(&argc, argv);          // initialize the toolkit
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB); // set the display mode
	glutInitWindowSize(640,480);     // set the window size
	glutInitWindowPosition(100, 150); // set the window position on the screen
	glutCreateWindow("my first attempt"); // open the screen window(with its exciting title)
	glutDisplayFunc(myDisplay);     // register the redraw function
	myInit();                   
	glutMainLoop(); 		     // go into a perpetual loop
	
	return 0;
}
