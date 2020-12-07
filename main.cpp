//A program for SPH method.
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "SPH.h"
using namespace std;


SPH::Sph* sph0 = new SPH::Sph();
#define NUM_THREADS 4

bool Compute_code = true;
bool Read_file = !Compute_code;

static double fovy = 50, aspect = 1, zFar = 1.0 + SPH::h, zNear = -1.0;//投影参数
static double eyex = 3 * SPH::l, eyey = 2 * SPH::m, eyez = 0.5 * SPH::h;
static double centerx = SPH::l / 2, centery = SPH::m / 2, centerz = SPH::h / 2;
static double upx = 0, upy = 0.0, upz = 1.0;//三维观察

void init(void) {
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);

	gluPerspective(fovy, aspect, zFar, zNear);
	gluLookAt(eyex, eyey, eyez, centerx, centery, centerz, upx, upy, upz);
}

static void display(void)
{
	init();

	if(Compute_code){ sph0->step(); }
	else if (Read_file) { sph0->OnlyReadFileStep(); }

	if (sph0->getstep() > sph0->getMaxstep()) {
		system("pause");
		exit(0);
	}

	glutSwapBuffers();
}

static void idle(void)
{
	glutPostRedisplay();
}

int main(int argc, char *argv[])
{
	omp_set_num_threads(NUM_THREADS);
	glutInit(&argc, argv);
	glutInitWindowSize(600, 600);
	glutInitWindowPosition(10, 10);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_STENCIL);

	glutCreateWindow("SPH method");

	init();

	glutDisplayFunc(display);
	glutIdleFunc(idle);

	glutMainLoop();

	return EXIT_SUCCESS;
}

