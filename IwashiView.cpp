/*
 * ShishamoView.cpp
 *
 *  Created on: 2013/02/15
 *      Author: mnsaru
 */

using namespace std;

#include <stdio.h>
#include <iostream>
using std::ios;
using std::cin;
using std::cout;
using std::endl;

#include <math.h>

#include <vector>
#include <sys/types.h>
#include <dirent.h>
#include <sys/dir.h>
#include <sys/stat.h>

#include <GL/glui.h>
//#include <GL/glut.h>
//#include <GL/freeglut.h>

#include "ReadCSV.h"
#include "OffsetHilbert.h"
#include "ExtractSurface.h"

//#define filepath "./data/"
#define filepath "/var/run/media/nagaso/sotoHD/iwana/2/"
//#define filepath "/var/run/media/nagaso/sotoHD/yamame/1/"
#define data_width 2
#define data_length 10000
#define sampling_rate 0.00000001

//iwana1
//#define x_steps 158
//#define y_steps 41
//iwana2
#define x_steps 157
#define y_steps 35
//iwana3
//#define x_steps 137
//#define y_steps 29
//iwana4
//#define x_steps 158
//#define y_steps 41
//iwana5
//#define x_steps 154
//#define y_steps 33
//iwana11
//#define x_steps 141
//#define y_steps 31

//yamame1
//#define x_steps 131
//#define y_steps 32
//yamame2
//#define x_steps 158
//#define y_steps 41
//yamame3
//#define x_steps 158
//#define y_steps 41
//yamame4
//#define x_steps 158
//#define y_steps 41

#define filenum x_steps * y_steps

#define sonicvelo 1500 //onsoku m/s
#define ini_ignore 4000 //for coloring


vector<direct *> entries;

double dataview[filenum][data_length][data_width];
int point_status_3D[filenum][10000];
int mountain[filenum][10000][2];
int viewchanged_flag = 0;
double colorlange;
int reftimebottom;
#define headerchangenum 1 //the x axe number one before header name changes
#define thresholdVal spinnerValf

//-------varieties for Delaunay view----
#include "Delaunay2D.h"
std::set<Tercel::Vector> vertices;
std::set<Tercel::Triangle> triangles;

//-------varieties for glut-----------
int WinID[2]; //ウィンドウID
int WindowNum = 0;
int WinFlag[2];
const char *WindowName[]={"main_window", "sub_window"};

//-------varieties for glui-----------
float rotary[16] = {
	1.0, 0.0, 0.0, 0.0,
	0.0, 1.0, 0.0, 0.0,
	0.0, 0.0, 1.0, 0.0,
	0.0, 0.0, 0.0, 1.0
};

float spinnerValf;
float spinnerSurfVal;
int x_begin;
int x_end;
int y_begin;
int y_end;
double z_begin;
double z_end;
float z_front;
float z_back;
float trans_ary[2] = {0.0, 0.0};
float trans_aryZ[] = {0.0};
int rot_deg = 0;
float rot_axs_x = 0;
float rot_axs_y = 0;
float rot_axs_z = 0;
int sub_x;
int sub_y;
float trans_ary_sub[] = {0.0};
float trans_aryZ_sub[] = {0.0};
float brightness;
int pointsize = 1;
int radiobuttonVal;
double peak_threshold;

double getmaxvol()
{
	double volmax = 0;

	for(int i = 0; i < filenum; i++)
		for( int j = ini_ignore; j < data_length; j++)
		{
			if( volmax < dataview[i][j][1] )
				volmax = dataview[i][j][1];
		}


	return volmax;
}
double gettimemax()
{
	double volmax = 0;

	for(int i = 0; i < filenum; i++)
		for( int j = 0; j < data_length; j++)
		{
			if( volmax < dataview[i][j][0] )
				volmax = dataview[i][j][0];
		}


	return volmax;
}
double gettimemin()
{
	double volmin = 0;

	for(int i = 0; i < filenum; i++)
		for( int j = 0; j < data_length; j++)
		{
			if( volmin > dataview[i][j][1] )
				volmin = dataview[i][j][1];
		}


	return volmin;
}

int gettimebottom()
{
	double valmax = 0;
	int timebottom = 0;

	for( int i = 0; i < data_length; i++)
		if( dataview[0][i][1] > valmax)
		{
			valmax = dataview[0][i][1];
			timebottom = i;
		}

	return timebottom;
}

double makeDistance( double time )
{
	return time * sonicvelo * 1000 / 2;
}

void readOnly(int xs, int ys, int count)
{
	bool headerflag = 0;
	string filename = entries[count]->d_name;
	//cout << "writing file name..." << filename << endl;
	ReadCSV test(filename, xs, ys, headerflag);


	for(int k = 0; k < data_length; k++){
		dataview[count][k][0] = test.data[k][0];
		dataview[count][k][1] = test.data[k][1];
	}
}

double matchingAmpFunc(double data0[data_length], double (*data1)[data_width])
{
	double max_0, max_1;
	double ratio;

	for(int i = 0; i < data_length; i++)
	{
		if(max_0 < fabs(data0[i]))
			max_0 = fabs(data0[i]);
		if(max_1 < fabs(data1[i][1]))
			max_1 = fabs(data1[i][1]);
	}

	ratio = max_1 / max_0;
	return ratio;
}

void readAndHilbert(int xs, int ys, int count)
{
	bool headerflag = 0;
	string filename = entries[count]->d_name;
	//cout << "writing file name..." << filename << endl;
	ReadCSV test(filename, xs, ys, headerflag);

	double temp[data_length];
	for(int k = 0; k < data_length; k++)
		temp[k] = test.data[k][1];

	OffsetHilbert test2(temp);


	//double ratio = matchingAmpFunc(temp, test2.data2);

	for(int k = 0; k < data_length; k++){
		dataview[count][k][0] = test.data[k][0];
		dataview[count][k][1] = test2.data2[k][1];
		//dataview[count][k][1] = test2.data2[k][1] / ratio;
		if(xs == 0 && ys == 0)
			cout << dataview[count][k][1] << endl;
	}

	if(xs == 0 && ys == 0){
			cout << "testing OffsetHilbert..."
					<< endl;
			reftimebottom = gettimebottom(); //get the distance to the mounting board from 0,0 data signal
	}
}

void calcImpedance()
{
	double imp; //impedance
	double Sig_ini = 1.0; //strength of initial signal
	double impWater = 1.0 * sonicvelo;

	for( int i = 0; i < filenum; i++)
		for( int j = 0; j < data_length; j++)
		{
			imp = -1 * (dataview[i][j][1] + Sig_ini) / (dataview[i][j][1] - Sig_ini) * impWater;
			dataview[i][j][1] = imp;

			if(i == 1000)
				cout << "imp " << imp << endl;
		}
}

void extractSurface(int x, int y, int count, int point_status[])
{
	double temp[data_length];
	for( int i = 0; i < data_length; i++)
		temp[i] = dataview[count][i][1];

	ExtractSurface ES(reftimebottom, temp, spinnerSurfVal);

	for( int i = 0; i < data_length; i++)
		point_status[i] = 0;

	for( int i = 0; i <data_length; i++)
		point_status[i] = ES.point_status[i];

	if(x == sub_x && y == sub_y){
		peak_threshold = ES.threshold;
	}

	for(int i = 0; i < 10000; i++)
	{
		for(int j = 0; j < 2; j++)
		{
			if(i < ES.mountain_num)
				mountain[count][i][j] = ES.mountain_place[i][j];
			else
				mountain[count][i][j] = 0;
		}
	}
}


void makeArray()
{
	int count = 0;

	for( int i = 0; i < y_steps; i++)
		for( int j = 0; j < x_steps; j++)
		{
			//readOnly(i, j, count);
			readAndHilbert(j, i, count);
			//extractSurface(i, j, count);
			count++;
		}
	z_begin = makeDistance(0);
	z_end =  makeDistance(gettimemax());
	//calcImpedance();

	cout << "z min/max..." << z_begin << "/" << z_end << endl;
	cout << "last count filenum = " << count << endl;

	colorlange = getmaxvol();
	cout << "max_vol..." << colorlange << endl;
}


//----------------------------functions for indicate main_window-----------------------------------------
void setColor(double incol, double lancol)
{
	float R, G, B;
	//double geta = 0.3;
	//monocolor
	//lancol = 0.02;
	lancol *= brightness;

	R = incol / lancol;
	G = incol / lancol;
	B = incol / lancol;

	if(R >= 1.0)
		R = 1.0;
	if(G >= 1.0)
		G = 1.0;
	if(B >= 1.0)
		B = 1.0;


	glColor3f(R, G, B);
}

void Draw2d()
{
	int multnum_y_axis = 100;
	int peak_count = 0;
	int mountain_count = 0;
	bool mountain_inout_flag = 0;
	double y_min, y_max;
	int filecount = sub_y * x_steps + sub_x;

	//extractSurface(sub_x, sub_y, filecount, peak_place);
	glPointSize(1.0);
	glBegin(GL_POINTS);
	if(viewchanged_flag == 0)// drawing for pointview
	{
		for(int i = 0; i < data_length; i++){
			GLdouble x_point = makeDistance(dataview[filecount][i][0]);
			GLdouble y_point = dataview[filecount][i][1] / colorlange * multnum_y_axis;
			double distance = makeDistance(dataview[filecount][i][0]);

			if(	distance > z_front && distance < z_back )
			{
				glColor3f(1.0, 1.0, 1.0);
				glVertex2d(x_point, y_point);

				if(y_point > y_max)
					y_max = distance;
				if(y_point < y_min)
					y_min = distance;
			}
		}
	}
	else //for surfaceview and one point view
	{
		for(int i = 0; i < data_length; i++){
			GLdouble x_point = makeDistance(dataview[filecount][i][0]);
			GLdouble y_point = dataview[filecount][i][1] / colorlange * multnum_y_axis;
			double distance = makeDistance(dataview[filecount][i][0]);

			if(	distance > z_front && distance < z_back )
			{
				if(point_status_3D[filecount][i] == 0)
					glColor3f(1.0, 1.0, 1.0);
				else if( point_status_3D[filecount][i] == 1)
					glColor3f(0.0, 1.0, 0.0);
				else
					glColor3f(1.0, 0.0, 0.0);

				glVertex2d(x_point, y_point);

				if(y_point > y_max)
					y_max = distance;
				if(y_point < y_min)
					y_min = distance;
			}
		}
	}

	glEnd();
	//draw threshold line of point view
	glBegin(GL_LINES);
	glColor3f(1.0, 0.0, 0.0);
	glVertex2d(z_front, thresholdVal * multnum_y_axis);
	glVertex2d(z_back, thresholdVal * multnum_y_axis);
	glEnd();
	//draw threshold line of peak detection
	glBegin(GL_LINES);
	glColor3f(0.0, 0.0, 1.0);
	glVertex2d(z_front, peak_threshold * multnum_y_axis);
	glVertex2d(z_back, peak_threshold * multnum_y_axis);
	glEnd();

}

void Draw3d()
{

	if(radiobuttonVal != 3)//drawing not for delaunay view
	{
		glPointSize(pointsize);
		glBegin(GL_POINTS);
		for( int i = 0; i < filenum; i++)
		{
			GLdouble y_point = int( i / x_steps);
			GLdouble x_point = i % x_steps;

			if( radiobuttonVal == 0) //point view
			{
				for(int j = 0; j < data_length; j++)
				{
					GLdouble z_point = makeDistance(dataview[i][j][0]);

					if( fabs(dataview[i][j][1] / colorlange) > thresholdVal
							&& x_point > x_begin
							&& x_point < x_end
							&& y_point > y_begin
							&& y_point < y_end
							&& z_point > z_front
							&& z_point < z_back
						)
					{
						setColor( dataview[i][j][1], colorlange);
						if(x_point == sub_x && y_point == sub_y)
							glColor3f(1.0, 0.0, 0.0);

						glVertex3d(x_point, y_point, z_point);
					}
				}
			}
			else if( radiobuttonVal == 1) //surface view
			{
				for(int j = 0; j < data_length; j++)
				{
					if(point_status_3D[i][j] == 2)
					{
						GLdouble z_point = makeDistance(dataview[i][j][0]);
						if( x_point > x_begin
							&& x_point < x_end
							&& y_point > y_begin
							&& y_point < y_end
							&& z_point > z_front
							&& z_point < z_back
							)
						{
							setColor( dataview[i][j][1], colorlange);
							if(x_point == sub_x && y_point == sub_y)
								glColor3f(1.0, 0.0, 0.0);

							glVertex3d(x_point, y_point, z_point);
						}
					}

				}
			}
			else if( radiobuttonVal == 2) //one point view
			{
				int top_point_count = 0;
				for(int j = 0; j < data_length; j++)
				{
					if(top_point_count == 0 && point_status_3D[i][j] == 2 )
					{
						GLdouble z_point = makeDistance(dataview[i][j][0]);
						if( x_point > x_begin
							&& x_point < x_end
							&& y_point > y_begin
							&& y_point < y_end
							&& z_point > z_front
							&& z_point < z_back
							)
						{
							setColor( dataview[i][j][1], colorlange);
							if(x_point == sub_x && y_point == sub_y)
								glColor3f(1.0, 0.0, 0.0);

							glVertex3d(x_point, y_point, z_point);
						}
						top_point_count++;
					}
				}
			}
		}
		glEnd();
	}
	else if( radiobuttonVal == 3) // for delaunay view
	{
		double zpoint[filenum] = {0.0};// memo z cord. of each surface point
		glPointSize(pointsize);

		for( int i = 0; i < filenum; i++)
		{
			GLdouble y_point = int( i / x_steps);
			GLdouble x_point = i % x_steps;
			int top_point_count = 0;

			if( x_point > x_begin
					&& x_point < x_end
					&& y_point > y_begin
					&& y_point < y_end
				)
			{
				for( int j  = 0; j < data_length; j++)
				{
					if( top_point_count == 0 && point_status_3D[i][j] == 2)
					{
						Tercel::Vector v;
						v.x = x_point;
						v.y = y_point;
						vertices.insert(v);

						zpoint[i] = makeDistance(dataview[i][j][0]);
						top_point_count++;
					}
				}
			}
		}
		Tercel::Delaunay2d::getDelaunayTriangles( vertices, &triangles); // Delaunay triangulation invoke

		typedef std::set<Tercel::Triangle> TriangleSet;
		for( TriangleSet::iterator it = triangles.begin(); it != triangles.end(); ++it) //draw triangulated triangles
		{
			Tercel::Triangle t = *it;
			int x1 = t.p1->x;int x2 = t.p2->x;int x3 = t.p3->x;
			int y1 = t.p1->y;int y2 = t.p2->y;int y3 = t.p3->y;
			int filecount1 = y1 * x_steps + x1;
			int filecount2 = y2 * x_steps + x2;
			int filecount3 = y3 * x_steps + x3;
			GLdouble z1 = zpoint[filecount1];
			GLdouble z2 = zpoint[filecount2];
			GLdouble z3 = zpoint[filecount3];

			if(z1 != 0 && z2 != 0 && z3 != 0)
			{

				/*glBegin(GL_TRIANGLES);

				glVertex2d(t.p1->x, t.p1->y);
				glVertex2d(t.p2->x, t.p2->y);
				glVertex2d(t.p3->x, t.p3->y);

				glEnd();*/

				glBegin(GL_LINE_STRIP);
				glVertex3d(x1, y1, z1);
				glVertex3d(x2, y2, z2);
				glVertex3d(x3, y3, z3);
				glVertex3d(x1, y1, z1);

				glEnd();
			}
		}
	}


	//cout << z_begin << "  " << z_end << endl;
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT);


	glPushMatrix();
	gluLookAt( 0.0, 0.0, 300.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	//glRotatef( 90.0, 0.0, 0.0, 1.0);
	glRotatef( rot_deg, rot_axs_x, rot_axs_y, rot_axs_z);
	glTranslated( -1 * x_steps / 2, -1 * y_steps / 2, 0.0);

	glTranslatef(trans_ary[0],trans_ary[1], trans_aryZ[0]);
	glMultMatrixf( rotary );

	//Draw2d();
	Draw3d();
	glPopMatrix();

	glutSwapBuffers();
	glutPostRedisplay();
}


void gluiCallback(int num)
{
	exit(0);
}
void funcReset(int num)
{
	for(int i = 0; i < 16; i++){
		rotary[i] = 0.0;
		if(i == 0 || i == 5 || i == 10 || i == 15)
			rotary[i] = 1.0;
	}
}
void funcResetTrans(int num)
{
	trans_ary[0] = 0.0; trans_ary[1] = 0.0;
	trans_aryZ[0] = 0.0;
}
void funcResetSub(int num)
{
	sub_x = 99;
	sub_y = 99;
}

void init()
{
	glClearColor(0.0, 0.0, 0.0, 1.0);
}

void idle()
{
	for(int loop = 0; loop < WindowNum; ++loop){
		if(WinFlag[loop]>0)
		{
			//printf("idle, loop=%d, WinFlag=%d\n", loop, WinFlag[loop]);
			glutSetWindow(WinID[loop]);
			glutPostRedisplay(); //redraw (calls display() fnc)
		}
	}
}


void keyboard(unsigned char key, int x, int y)
{
	switch (key) {
	case 'q':
	case 'Q':
	case '\033': // esc
		exit(0);
	default:
		break;
	}
}

void resize( GLsizei w, GLsizei h)
{
	GLfloat fAspect;
	//configure region for drawing picture in window
	glViewport(0, 0, w, h);

	fAspect = (GLfloat)w / (GLfloat)h;

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective( 30.0f, fAspect, 1.0, 600.0);//sawaranai

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void Init()
{
	glutInitWindowPosition(0,0);
	glutInitWindowSize(1000,1000);
}

void callBacks()
{
	//call back func
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(resize);
	init();
}


//-----------------------functions to set up sub window-------------------------
void display2()
{
	glClear(GL_COLOR_BUFFER_BIT);

	glPushMatrix();
	gluLookAt( 0.0, 0.0, 50.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	//glRotatef( 90.0, 0.0, 0.0, 1.0);
	//glRotatef( 180.0, 0.0, 1.0, 0.0);
	glTranslated(0.0, 0.0, 0.0);

	glTranslatef(trans_ary_sub[0], 0.0, trans_aryZ_sub[0]);

	Draw2d();

	glPopMatrix();

	glutSwapBuffers();
	glutPostRedisplay();
}

void resize2( GLsizei w, GLsizei h)
{
	GLfloat fAspect;
	//configure region for drawing picture in window
	glViewport(0, 0, w, h);

	fAspect = (GLfloat)w / (GLfloat)h;

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective( 30.0f, fAspect, 1.0, 600.0);//sawaranai

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void Init2()
{
	glutInitWindowPosition(1000,0);
	glutInitWindowSize(800,500);
}
void callBacks2()
{
	glutDisplayFunc(display2);
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(resize2);

	//init();
}




char* targetdir;

int dotOmit(struct dirent *darray)
{
	if( strcmp( darray->d_name,".") && strcmp( darray->d_name, ".."))
	{
		return 1;
	}

	return 0;
}

int timesort(const dirent **v1, const dirent **v2)
{
	struct dirent *d1, *d2;
	struct stat s1, s2;
	char fnm[256];

	d1 = *(struct dirent **)v1;
	d2 = *(struct dirent **)v2;

	//cout << targetdir << endl;
	//cout << fnm << endl;

	strcpy(fnm, targetdir);
	strcat(fnm, d1->d_name);
	stat(fnm, &s1);

	strcpy(fnm, targetdir);
	strcat(fnm, d2->d_name);
	stat(fnm, &s2);

	return (int)(s1.st_mtime - s2.st_mtime);
}

int fileFilter(const direct *dir)
{
	const char * s = dir->d_name;
	int csvlen = strlen(s) - 4;
	if ( csvlen >= 0 )
		if ( strncmp(s+csvlen,".csv",4) == 0)
			return 1;
	return 0;
}

void getDirEntries(const string& dirname, vector<direct *> &entries)
{
	direct** darray;

	targetdir = const_cast<char*>(dirname.c_str());
	//cout << dirname << "  " << targetdir << endl;
	int entryCount = scandir(const_cast<char*>(dirname.c_str()), &darray, fileFilter, timesort);

	for(int k = 0; k < entryCount; k++)
	{
		entries.push_back(darray[k]);
	}
}

void getFileList()
{
	string filename = filepath;
	//vector<direct *> entries;
	getDirEntries(filename, entries);

	cout << "making file list... " << endl;

	for(int k = 0; k < entries.size(); k++){
		cout << k << ".\t" << entries[k]->d_name << endl;
	}

}

void viewchange_callback(int val)
{
	//cout << "val is" << val << " " << radiobuttonVal << endl;
	if( radiobuttonVal == 0){
		cout << "point view" << endl;
		viewchanged_flag = 0;
	}
	else if( radiobuttonVal == 1)
	{
		cout << "surface view" << endl;
		viewchanged_flag = 1;
	}
	else if ( radiobuttonVal == 2)
	{
		cout << "one point view" << endl;
		viewchanged_flag = 2;
	}
	else
	{
		cout << "delaunay view" << endl;
		viewchanged_flag = 3;
	}


	for(int i = 0; i < filenum; i++)
	{
		int y_point = int( i / x_steps);
		int x_point = i % x_steps;
		int point_status[10000];
		extractSurface( x_point, y_point, i, point_status);
		for(int j = 0; j < 10000; j++)
			point_status_3D[i][j] = point_status[j];
	}
	//viewchanged_flag = 1;
	cout << "peak_place list has changed" << endl;


}

int main(int argc, char *argv[])
{
	getFileList();

	//3dview---------------------------------------
	makeArray();




	//initiation GLUT
	glutInit(&argc, argv);//initialize OpenGL environment
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);//double buffering

	//funcs for main_window
	Init();
	WinID[WindowNum] = glutCreateWindow(WindowName[WindowNum]);
	callBacks();
	WinFlag[WindowNum] = 1;
	WindowNum = WindowNum + 1;

	//funcs for sub_window
	Init2();
	WinID[WindowNum] = glutCreateWindow(WindowName[WindowNum]);
	callBacks2();
	WinFlag[WindowNum] = 1;
	WindowNum=WindowNum+1;

	glutIdleFunc(idle);

	//GLUI window and add controls

	//controls for 3D view
	GLUI *glui = GLUI_Master.create_glui("control", 0);

	GLUI_Rotation *view_rot = glui->add_rotation("Rotation", rotary);
		glui->add_button("rot_reset", 0, funcReset);
	GLUI_EditText *segment_edittext_rot_deg =
			glui->add_edittext( "rot_deg", GLUI_EDITTEXT_INT, &rot_deg);
	GLUI_EditText *segment_edittext_rot_axs_x =
				glui->add_edittext( "rot_axs_x", GLUI_EDITTEXT_FLOAT, &rot_axs_x);
	GLUI_EditText *segment_edittext_rot_axs_y =
				glui->add_edittext( "rot_axs_y", GLUI_EDITTEXT_FLOAT, &rot_axs_y);
	GLUI_EditText *segment_edittext_rot_axs_z =
				glui->add_edittext( "rot_axs_z", GLUI_EDITTEXT_FLOAT, &rot_axs_z);

	GLUI_Translation *translation_xy = glui->add_translation( "TranslationXY", GLUI_TRANSLATION_XY, trans_ary);
	GLUI_Translation *translation_z = glui->add_translation( "TranslationZ", GLUI_TRANSLATION_Z, trans_aryZ);
	glui->add_button("trans_reset", 0, funcResetTrans);


	GLUI_Spinner *segment_spinner =
			glui->add_spinner( "threshold", GLUI_SPINNER_FLOAT, &spinnerValf);
	segment_spinner->set_float_limits(0.0, 1.0, GLUI_LIMIT_CLAMP);
	segment_spinner->set_float_val(0.006);

	GLUI_EditText *segment_edittext_xmin =
			glui->add_edittext( "x_min", GLUI_EDITTEXT_INT, &x_begin);
	segment_edittext_xmin->set_int_limits( 0, x_steps, GLUI_LIMIT_CLAMP);
	segment_edittext_xmin->set_int_val(0);
	GLUI_EditText *segment_edittext_xmax =
			glui->add_edittext( "x_max", GLUI_EDITTEXT_INT, &x_end);
	segment_edittext_xmax->set_int_limits( 0, x_steps, GLUI_LIMIT_CLAMP);
	segment_edittext_xmax->set_int_val(x_steps);
	GLUI_EditText *segment_edittext_ymin =
			glui->add_edittext( "y_min", GLUI_EDITTEXT_INT, &y_begin);
	segment_edittext_ymin->set_int_limits( 0, y_steps, GLUI_LIMIT_CLAMP);
	segment_edittext_ymin->set_int_val(0);
	GLUI_EditText *segment_edittext_ymax =
			glui->add_edittext( "y_max", GLUI_EDITTEXT_INT, &y_end);
	segment_edittext_ymax->set_int_limits( 0, y_steps, GLUI_LIMIT_CLAMP);
	segment_edittext_ymax->set_int_val(y_steps);
	GLUI_EditText *segment_edittext_zmin =
			glui->add_edittext( "z_min", GLUI_EDITTEXT_FLOAT, &z_front);
	segment_edittext_zmin->set_float_limits( 0.0, 100.0, GLUI_LIMIT_CLAMP);
	segment_edittext_zmin->set_float_val(10.0);
	//z_begin = segment_edittext_zmin->get_float_val();
	GLUI_EditText *segment_edittext_zmax =
			glui->add_edittext( "z_max", GLUI_EDITTEXT_FLOAT, &z_back);
	segment_edittext_zmax->set_float_limits( 0.0, 100.0, GLUI_LIMIT_CLAMP);
	segment_edittext_zmax->set_float_val(35);
	GLUI_EditText *segment_edittext_brightness =
				glui->add_edittext( "brightness", GLUI_EDITTEXT_FLOAT, &brightness);
		segment_edittext_brightness->set_float_limits( 0.0, 1.0, GLUI_LIMIT_CLAMP);
		segment_edittext_brightness->set_float_val(0.04);
	GLUI_EditText *segment_edittext_point =
				glui->add_edittext( "point size", GLUI_EDITTEXT_INT, &pointsize);
		segment_edittext_point->set_int_limits( 0, 10000, GLUI_LIMIT_CLAMP);
		segment_edittext_point->set_int_val(1);
	GLUI_Spinner *segment_spinner_surface =
			glui->add_spinner( "surface_extraction", GLUI_SPINNER_FLOAT, &spinnerSurfVal, radiobuttonVal, viewchange_callback);
		segment_spinner_surface->set_float_limits(0.0, 10000.0, GLUI_LIMIT_CLAMP);
		segment_spinner_surface->set_float_val(5.0);

	GLUI_RadioGroup *radio_button =
			glui->add_radiogroup(&radiobuttonVal, radiobuttonVal, viewchange_callback);
		glui->add_radiobutton_to_group(radio_button, "point view");
		glui->add_radiobutton_to_group(radio_button, "surface view");
		glui->add_radiobutton_to_group(radio_button, "one point view");
		glui->add_radiobutton_to_group(radio_button, "delaunay view");
	//controls for 2D view
	glui->add_separator();
	GLUI_EditText *segment_edittext_sub_x =
			glui->add_edittext( "sub_x", GLUI_EDITTEXT_INT, &sub_x);
	segment_edittext_sub_x->set_int_limits(0.0, x_steps, GLUI_LIMIT_CLAMP);
	segment_edittext_sub_x->set_int_val(70);
	GLUI_EditText *segment_edittext_sub_y =
				glui->add_edittext( "sub_y", GLUI_EDITTEXT_INT, &sub_y);
		segment_edittext_sub_y->set_int_limits(0.0, y_steps, GLUI_LIMIT_CLAMP);
		segment_edittext_sub_y->set_int_val(20);
	GLUI_Translation *translation_x_sub = glui->add_translation( "TranslationX", GLUI_TRANSLATION_X, trans_ary_sub);
		translation_x_sub->set_speed( 0.5 );
	GLUI_Translation *translation_z_sub = glui->add_translation( "TranslationZ", GLUI_TRANSLATION_Z, trans_aryZ_sub);
		translation_z_sub->set_speed( 0.5 );
	glui->add_button("reset_sub", 0, funcResetSub);


	glui->add_button("Exit", 0, gluiCallback);


	glutMainLoop();



	/*
	//2dview--------------------------------
	glutInitWindowPosition(100,100);
	glutInitWindowSize(1000,250);
	glutInit(&argc, argv);//initialize OpenGL environment
	glutInitDisplayMode(GLUT_DOUBLE);//double buffering
	glutInitDisplayMode(GLUT_RGBA);
	glutCreateWindow("test");
	glutDisplayFunc(display2d);
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(resize2d);
	init();
	glutMainLoop();
	*/

	return 0;
}


