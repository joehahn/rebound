#ifdef OPENGL
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#ifdef _APPLE
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif 
#include "zpr.h"
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "tree.h"
#include "opengl.h"

GLuint DListSPHERE;
#ifndef COLLISIONS_NONE
int display_spheres = 1;
#else
int display_spheres = 0;
#endif
int display_init_done = 0;
int display_pause_sim = 0;
int display_pause = 0;
int display_tree = 0;

void displayKey(unsigned char key, int x, int y){
	switch(key){
		case 'q':
			exit(0);
			break;
		case ' ':
			display_pause_sim=!display_pause_sim;
			if (display_pause_sim){
				printf("Pause.\n");
				glutIdleFunc(NULL);
			}else{
				printf("Resume.\n");
				glutIdleFunc(iterate);
			}
			break;
		case 's':
			display_spheres = !display_spheres;
			break;
		case 'r':
			zprReset(0.7/boxsize_max);
			break;
		case 't':
			display_tree = !display_tree;
			break;
		case 'd':
			display_pause = !display_pause;
			break;
	}
	display();
}

#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)
void display_cell(struct cell* node){
	if (node == NULL) return;
	glTranslatef(node->x,node->y,node->z);
	glutWireCube(node->w);
	glTranslatef(-node->x,-node->y,-node->z);
	for (int i=0;i<8;i++) {
		display_cell(node->oct[i]);
	}
}
void display_entire_tree(){
	for(int i=0;i<root_nx;i++){
	for(int j=0;j<root_ny;j++){
	for(int k=0;k<root_nz;k++){
		int index = (k*root_ny+j)*root_nx+i;
		display_cell(root[index]);
	}
	}
	}
}
#endif

void display(){
	if (display_pause) return;
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	if (display_spheres){
		glDisable(GL_BLEND);                    
		glDepthMask(GL_TRUE);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		GLfloat lightpos[] = {0, boxsize_max, boxsize_max, 0.f};
		glLightfv(GL_LIGHT0, GL_POSITION, lightpos);
	}else{
		glEnable(GL_BLEND);                    
		glDepthMask(GL_FALSE);
		glDisable(GL_DEPTH_TEST);
		glDisable(GL_LIGHTING);
		glDisable(GL_LIGHT0);
	}
	glTranslatef(0,0,-boxsize_max);
	glPointSize(5.);
	glEnable(GL_POINT_SMOOTH);
	glVertexPointer(3, GL_DOUBLE, sizeof(struct particle), particles);
	for (int i=-nghostx;i<=nghostx;i++){
	for (int j=-nghosty;j<=nghosty;j++){
	for (int k=-nghostz;k<=nghostz;k++){
		struct ghostbox gb = get_ghostbox(i,j,k);
		glTranslatef(gb.shiftx,gb.shifty,gb.shiftz);
		if (display_spheres){
			// Drawing Spheres
			glColor4f(1.0,1.0,1.0,1.0);
#ifndef COLLISIONS_NONE
#ifdef _APPLE
			for (int i=0;i<N;i++){
				struct particle p = particles[i];
				glTranslatef(p.x,p.y,p.z);
				glScalef(p.r,p.r,p.r);
				glCallList(DListSPHERE);
				glScalef(1./p.r,1./p.r,1./p.r);
				glTranslatef(-p.x,-p.y,-p.z);
			}
#endif
#endif
		}else{
			// Drawing Points
			glEnableClientState(GL_VERTEX_ARRAY);
			glColor4f(1.0,0.0,0.0,0.9);
			glDrawArrays(GL_POINTS, 0, N_active_first);
			glColor4f(1.0,1.0,0.0,0.6);
			glDrawArrays(GL_POINTS, N_active_first, N_active_last-N_active_first);
			glColor4f(1.0,1.0,0.0,0.9);
			glDrawArrays(GL_POINTS, N_active_last, N-N_active_last);
			glDisableClientState(GL_VERTEX_ARRAY);
		}
		// Drawing Tree
		glColor4f(1.0,0.0,0.0,0.4);
#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)
		if (display_tree){
			glColor4f(1.0,0.0,0.0,0.4);
			display_entire_tree();
		}
#endif
		glTranslatef(-gb.shiftx,-gb.shifty,-gb.shiftz);
	}
	}
	}
	glColor4f(1.0,0.0,0.0,0.4);
	glScalef(boxsize_x,boxsize_y,boxsize_z);
	glutWireCube(1);
	glScalef(1./boxsize_x,1./boxsize_y,1./boxsize_z);
	glutSwapBuffers();
	glTranslatef(0,0,boxsize_max);
}

void init_display(int argc, char* argv[]){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );
	glutInitWindowSize(700,700);
	glutCreateWindow("nbody");
	zprInit(0.7/boxsize_max);
	glutDisplayFunc(display);
	glutIdleFunc(iterate);
	glutKeyboardFunc(displayKey);
	glEnable(GL_BLEND);                    
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);  
	
	// Sphere
#ifdef _APPLE
	DListSPHERE = glGenLists(1);
	GLUquadricObj *sphere;
	glNewList(DListSPHERE, GL_COMPILE);
	sphere = gluNewQuadric();
	gluSphere(sphere, 1.f, 20, 20);
	gluDeleteQuadric(sphere);
	glEndList();
#endif
  	
	// Light

	glCullFace(GL_BACK);
	glShadeModel ( GL_SMOOTH );
	glEnable( GL_NORMALIZE );
	glEnable(GL_COLOR_MATERIAL);
	static GLfloat light[] = {0.7f, 0.7f, 0.7f, 1.f};
	static GLfloat lightspec[] = {0.2f, 0.2f, 0.2f, 1.f};
	static GLfloat lmodel_ambient[] = { 0.15, 0.14, 0.13, 1.0 };

	glLightfv(GL_LIGHT0, GL_DIFFUSE, light );
	glLightfv(GL_LIGHT0, GL_SPECULAR, lightspec );
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);

	static GLfloat sphere_mat[] = {0.8f, 0.8f, 0.8f, 1.f};
	static GLfloat sphere_spec[] = {1.0f, 1.0f, 1.0f, 1.f};
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, sphere_mat);
	glMaterialfv(GL_FRONT, GL_SPECULAR, sphere_spec);
	glMaterialf(GL_FRONT, GL_SHININESS, 80);

	display_init_done =1; 

	glutMainLoop();
}

#endif