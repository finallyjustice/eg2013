#include "sph_header.h"
#include "sph_data.h"
#include "sph_timer.h"
#include "sph_system.h"
#include "sph_bitmap.h"
#include "MarchingCube.h"
#include <GL\glew.h>
#include <GL\glut.h>

using namespace std;

#pragma comment(lib, "glew32.lib") 

unsigned long *screenData;
int frameNum;
int create_video=0;

SPHSystem *sph;
MarchingCube *mc;

Timer *sph_timer;
char *window_title;

GLuint v;
GLuint f;
GLuint p;

GLfloat LightAmbient[]= { 1.0f, 1.0f, 1.0f, 1.0f };
GLfloat LightDiffuse[]= { 0.1f, 0.1f, 1.0f, 1.0f };
GLfloat LightPosition[]= { -20.0f, -20.0f, -20.0f, 1.0f };

int create_ray=0;
void make_mesh(int n)
{
	/*location <0,3,-35> \n \
	//look_at <0,0,0> \n \*/
	char filename[] = "Povray_Shots/000000.pov";
	char header[] = "#include \"colors.inc\" \n \
					camera { \n \
					location <0,4,-25> \n \
					look_at <0,-7,1> \n \
					} \n \
					light_source {<13, 13, -20> color White} \n \
					mesh {\n";
					
	int i = 18;

	while (n > 0) {
		filename[i--] += n % 10;
		n /= 10;
	}

	FILE *fp;

	fp=fopen(filename, "w");

	fprintf(fp, header);

	for(uint count=0; count<mc->num_triangle; count++)
	{
		float t0x=mc->t1[count].x;
		float t0y=mc->t1[count].y;
		float t0z=mc->t1[count].z;

		float t1x=mc->t2[count].x;
		float t1y=mc->t2[count].y;
		float t1z=mc->t2[count].z;

		float t2x=mc->t3[count].x;
		float t2y=mc->t3[count].y;
		float t2z=mc->t3[count].z;

		float nt0x=mc->n1[count].x;
		float nt0y=mc->n1[count].y;
		float nt0z=mc->n1[count].z;

		float nt1x=mc->n2[count].x;
		float nt1y=mc->n2[count].y;
		float nt1z=mc->n2[count].z;

		float nt2x=mc->n3[count].x;
		float nt2y=mc->n3[count].y;
		float nt2z=mc->n3[count].z;

		fprintf(fp, "smooth_triangle { <%f, %f, %f>, <%f, %f, %f>, <%f, %f, %f>, <%f, %f, %f>, <%f, %f, %f>, <%f, %f, %f> }  \n", t0x, t0y, t0z, nt0x, nt0y, nt0z, t1x, t1y, t1z, nt1x, nt1y, nt1z, t2x, t2y, t2z, nt2x, nt2y, nt2z);
	}

	fprintf(fp, "pigment { color red 0.3 green 0.3 blue 0.8 transmit 0.2} finish { ambient 0.3 diffuse 0.5 } } plane { <0,-1,0>, 10 pigment {checker color White, color Grey}}\n");

	fclose(fp);
}

uint global=0;
void screenShot(int n)
{
	Bitmap b;
	char filename[] = "Screen_Shots/000000.bmp";
	int i = 18;

	while (n > 0) {
		filename[i--] += n % 10;
		n /= 10;
	}

	glReadPixels(0, 0, window_width, window_height, GL_RGBA, GL_UNSIGNED_BYTE, screenData);
	b.makeBMP(screenData, window_width, window_height);
	b.Save(filename);

	global++;

	if(global == 1000)
	{
		exit(0);
	}
}

void set_shaders()
{
	char *vs=NULL;
	char *fs=NULL;

	vs=(char *)malloc(sizeof(char)*10000);
	fs=(char *)malloc(sizeof(char)*10000);
	memset(vs, 0, sizeof(char)*10000);
	memset(fs, 0, sizeof(char)*10000);

	FILE *fp;
	char c;
	int count;

	fp=fopen("shader/shader.vs", "r");
	count=0;
	while((c=fgetc(fp)) != EOF)
	{
		vs[count]=c;
		count++;
	}
	fclose(fp);

	fp=fopen("shader/shader.fs", "r");
	count=0;
	while((c=fgetc(fp)) != EOF)
	{
		fs[count]=c;
		count++;
	}
	fclose(fp);

	v=glCreateShader(GL_VERTEX_SHADER);
	f=glCreateShader(GL_FRAGMENT_SHADER);

	const char *vv;
	const char *ff;
	vv=vs;
	ff=fs;

	glShaderSource(v, 1, &vv, NULL);
	glShaderSource(f, 1, &ff, NULL);

	int success;

	glCompileShader(v);
	glGetShaderiv(v, GL_COMPILE_STATUS, &success);
	if(!success)
	{
		char info_log[5000];
		glGetShaderInfoLog(v, 5000, NULL, info_log);
		printf("Error in vertex shader compilation!\n");
		printf("Info Log: %s\n", info_log);
	}

	glCompileShader(f);
	glGetShaderiv(f, GL_COMPILE_STATUS, &success);
	if(!success)
	{
		char info_log[5000];
		glGetShaderInfoLog(f, 5000, NULL, info_log);
		printf("Error in fragment shader compilation!\n");
		printf("Info Log: %s\n", info_log);
	}

	p=glCreateProgram();
	glAttachShader(p, v);
	glAttachShader(p, f);
	glLinkProgram(p);
	glUseProgram(p);

	free(vs);
	free(fs);
}

void draw_box(float ox, float oy, float oz, float width, float height, float length)
{
    glLineWidth(1.0f);
	glColor3f(1.0f, 0.0f, 0.0f);

    glBegin(GL_LINES);   
		
        glVertex3f(ox, oy, oz);
        glVertex3f(ox+width, oy, oz);

        glVertex3f(ox, oy, oz);
        glVertex3f(ox, oy+height, oz);

        glVertex3f(ox, oy, oz);
        glVertex3f(ox, oy, oz+length);

        glVertex3f(ox+width, oy, oz);
        glVertex3f(ox+width, oy+height, oz);

        glVertex3f(ox+width, oy+height, oz);
        glVertex3f(ox, oy+height, oz);

        glVertex3f(ox, oy+height, oz+length);
        glVertex3f(ox, oy, oz+length);

        glVertex3f(ox, oy+height, oz+length);
        glVertex3f(ox, oy+height, oz);

        glVertex3f(ox+width, oy, oz);
        glVertex3f(ox+width, oy, oz+length);

        glVertex3f(ox, oy, oz+length);
        glVertex3f(ox+width, oy, oz+length);

        glVertex3f(ox+width, oy+height, oz);
        glVertex3f(ox+width, oy+height, oz+length);

        glVertex3f(ox+width, oy+height, oz+length);
        glVertex3f(ox+width, oy, oz+length);

        glVertex3f(ox, oy+height, oz+length);
        glVertex3f(ox+width, oy+height, oz+length);

    glEnd();
}

void init_sph_system()
{
	real_world_origin.x=-10.0f;
	real_world_origin.y=-10.0f;
	real_world_origin.z=-10.0f;

	real_world_side.x=20.0f;
	real_world_side.y=10.0f;
	real_world_side.z=10.0f;

	sph=new SPHSystem();
	sph->init_system();

	screenData = new unsigned long[(int)window_width * (int)window_height];
	frameNum = 0;

	sph_timer=new Timer();
	window_title=(char *)malloc(sizeof(char)*50);
}

void init()
{
	glewInit();

	glViewport(0, 0, window_width, window_height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(45.0, (float)window_width/window_height, 10.0f, 100.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0f, 0.0f, -3.0f);
}

void init_ratio()
{
	sim_ratio.x=real_world_side.x/sph->world_size.x;
	sim_ratio.y=real_world_side.y/sph->world_size.y;
	sim_ratio.z=real_world_side.z/sph->world_size.z;

	mc=new MarchingCube(sph->scalar_size.x, sph->scalar_size.y, sph->scalar_size.z, sph->scalar, sph->scalar_pos, real_world_origin, sim_ratio, sph->scalar_step, 0.019f);
}

void render_particles()
{
	glPointSize(1.0f);
	glColor3f(0.2f, 0.2f, 1.0f);

	for(uint i=0; i<sph->num_particle; i++)
	{
		if(sph->mem[i].level == 0)
		{
			glColor3f(1.0f, 0.2f, 0.2f);
			glBegin(GL_POINTS);
			glVertex4f(sph->mem[i].pos.x*sim_ratio.x+real_world_origin.x, 
						sph->mem[i].pos.y*sim_ratio.y+real_world_origin.y,
						sph->mem[i].pos.z*sim_ratio.z+real_world_origin.z,
						0.0f);
			glEnd();
		}
		else
		{
			glColor3f(1.0f, 1.0f, 0.2f);
			glBegin(GL_POINTS);
			glVertex4f(sph->mem[i].pos.x*sim_ratio.x+real_world_origin.x, 
						sph->mem[i].pos.y*sim_ratio.y+real_world_origin.y,
						sph->mem[i].pos.z*sim_ratio.z+real_world_origin.z,
						1.0f);
			glEnd();
		}
	}
}

void display_mesh()
{
	glBegin(GL_TRIANGLES);
	for(uint i=0; i<mc->num_triangle; i++)
	{
		glNormal3f(mc->n1[i].x, mc->n1[i].y, mc->n1[i].z);
		glVertex3f(mc->t1[i].x, mc->t1[i].y, mc->t1[i].z);

		glNormal3f(mc->n2[i].x, mc->n2[i].y, mc->n2[i].z);
		glVertex3f(mc->t2[i].x, mc->t2[i].y, mc->t2[i].z);

		glNormal3f(mc->n3[i].x, mc->n3[i].y, mc->n3[i].z);
		glVertex3f(mc->t3[i].x, mc->t3[i].y, mc->t3[i].z);
	}
	glEnd();
}

void display_func()
{
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_DEPTH_TEST);
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glPushMatrix();

	if(buttonState == 1)
	{
		xRot+=(xRotLength-xRot)*0.1f;
		yRot+=(yRotLength-yRot)*0.1f;
	}

	glTranslatef(xTrans, yTrans, zTrans);
    glRotatef(xRot, 1.0f, 0.0f, 0.0f);
    glRotatef(yRot, 0.0f, 1.0f, 0.0f);

	//sph->animation();

	/////////////////////////////////////

	if(sph->sys_running != 0)
	{
		//clock_t start_sim=clock();

		sph->build_table();
		if(sph->sys_render == 1)
		{
			sph->interp_scalar();
		}

		if(sph->enable_split == 1 && sph->tick == 0 || sph->enable_merge == 1 && sph->tick == 1)
		{
			//clock_t start_ana=clock();
			sph->comp_scalar();
			sph->find_surface();
			sph->comp_dist();
			//clock_t end_ana=clock();
			//double span_ana=static_cast<double>(end_ana-start_ana)/CLOCKS_PER_SEC*1000;
			//printf("ANA: %lf\n", span_ana);

			//clock_t start_ene=clock();
			//sph->comp_energy();
			//clock_t end_ene=clock();
			//double span_ene=static_cast<double>(end_ene-start_ene)/CLOCKS_PER_SEC*1000;
			//printf("ENE: %lf\n", span_ene);
		}

		sph->comp_dens_pres();
		sph->comp_force_adv();

		//clock_t start_dyn=clock();
		if(sph->enable_split == 1 && sph->tick == 0)
		{
			sph->comp_split();
		}

		if(sph->enable_merge == 1 && sph->tick == 1)
		{
			sph->comp_merge();
		}
		//clock_t end_dyn=clock();
		//double span_dyn=static_cast<double>(end_dyn-start_dyn)/CLOCKS_PER_SEC*1000;
		//printf("DYN: %lf\n", span_dyn);

		sph->advection();

		sph->tick=1-sph->tick;	

		//clock_t end_sim=clock();
		//double span_sim=static_cast<double>(end_sim-start_sim)/CLOCKS_PER_SEC*1000;
		//printf("SIM: %lf\n", span_sim);
	}
	//////////////////////////////////////

	if(sph->sys_render == 0)
	{
		glUseProgram(p);
		render_particles();
	}

	glUseProgram(0);
	//draw_box(real_world_origin.x, real_world_origin.y, real_world_origin.z, real_world_side.x, real_world_side.y, real_world_side.z);

	if(sph->sys_render == 1)
	{
		glEnable(GL_LIGHTING);
		glLightfv(GL_LIGHT1, GL_AMBIENT, LightAmbient);				
		glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse);			
		glLightfv(GL_LIGHT1, GL_POSITION,LightPosition);			
		glEnable(GL_LIGHT1);

		glColor3f(0.2f, 0.2f, 0.8f);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		mc->run();
		glUseProgram(p);
		render_particles();
		//display_mesh();
		glDisable(GL_LIGHTING);
	}

	if(create_ray == 1 && sph->sys_running == 1)
	{
		if (frameNum % 2 == 0)
		{
			make_mesh(frameNum / 2);
		}
	}

	glPopMatrix();

	if(create_video == 1 && sph->sys_running == 1)
	{
		if (frameNum % 2 == 0)
		{
			screenShot(frameNum / 2);
		}
	}

	if(create_ray == 1 || create_video == 1)
	{
		frameNum++;
	}

    glutSwapBuffers();
	
	sph_timer->update();
	memset(window_title, 0, 50);
	sprintf(window_title, "SPH System 3D. FPS: %f", sph_timer->get_fps());
	glutSetWindowTitle(window_title);
}

void idle_func()
{
	glutPostRedisplay();
}

void reshape_func(GLint width, GLint height)
{
	window_width=width;
	window_height=height;

	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(45.0, (float)width/height, 0.001, 100.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0f, 0.0f, -3.0f);
}

void keyboard_func(unsigned char key, int x, int y)
{
	if(key == ' ')
	{
		sph->sys_running=1-sph->sys_running;
	}

	if(key == 'w')
	{
		zTrans += 0.3f;
	}

	if(key == 's')
	{
		zTrans -= 0.3f;
	}

	if(key == 'a')
	{
		xTrans -= 0.3f;
	}

	if(key == 'd')
	{
		xTrans += 0.3f;
	}

	if(key == 'q')
	{
		yTrans -= 0.3f;
	}

	if(key == 'e')
	{
		yTrans += 0.3f;
	}

	if(key == '1')
	{
		sph->enable_split=1-sph->enable_split;
		printf("Split Status: %u\n", sph->enable_split);
	}

	if(key == '2')
	{
		sph->enable_merge=1-sph->enable_merge;
		printf("Merge Status: %u\n", sph->enable_merge);
	}

	if(key == '3')
	{
		create_video=1-create_video;
		printf("Video Status: %u\n", create_video);
	}

	if(key == '4')
	{
		sph->sys_render=1-sph->sys_render;
		printf("Render Status: %u\n", sph->sys_render);
	}

	if(key == '5')
	{
		create_ray=1-create_ray;
		printf("Ray Status: %d\n", create_ray);
	}

	glutPostRedisplay();
}

void mouse_func(int button, int state, int x, int y)
{
	if (state == GLUT_DOWN)
	{
        buttonState = 1;
	}
    else if (state == GLUT_UP)
	{
        buttonState = 0;
	}

    ox = x; oy = y;

    glutPostRedisplay();
}

void motion_func(int x, int y)
{
    float dx, dy;
    dx = (float)(x - ox);
    dy = (float)(y - oy);

	if (buttonState == 1) 
	{
		xRotLength += dy / 5.0f;
		yRotLength += dx / 5.0f;
	}

	ox = x; oy = y;

	glutPostRedisplay();
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(window_width, window_height);
    glutCreateWindow("SPH Fluid 3D");

	init_sph_system();
	init();
	init_ratio();
	set_shaders();
	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_NV);
	glEnable(GL_POINT_SPRITE_ARB);
	glTexEnvi(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
	glDepthMask(GL_TRUE);
	glEnable(GL_DEPTH_TEST);

	glShadeModel(GL_SMOOTH);						
	//glEnable(GL_LIGHTING);

    glutDisplayFunc(display_func);
	glutReshapeFunc(reshape_func);
	glutKeyboardFunc(keyboard_func);
	glutMouseFunc(mouse_func);
	glutMotionFunc(motion_func);
	glutIdleFunc(idle_func);

    glutMainLoop();

    return 0;
}
