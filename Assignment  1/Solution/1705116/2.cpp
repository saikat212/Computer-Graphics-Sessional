#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>

#include <windows.h>
#include <GL/glut.h>

#define pi (2*acos(0.0))
using namespace std ;

double cameraHeight;
double cameraAngle;
double changed_cameraHight = 4;
double changed_cameraAngle = 0.04;
int drawgrid;
int drawaxes;
double angle;

int MoveForward = 0;
int MoveBackward = 0;
int RotateLeft =0;
int RotateRight =0;


#define SLICES 50
#define STACKS 50

#define WHEEL_RADIUS 20
#define WHEEL_WIDTH 10
#define AXLE_WIDTH (0.3*WHEEL_WIDTH)

#define MOVE_ANGLE 3.0
#define MOVE_DISTANCE (WHEEL_RADIUS*MOVE_ANGLE*pi/180.0)
#define ANGLE_CHANGE_FOR_ROTATION 1.0

#define c1 0.7
#define c2 0.7
#define c3 0.7

#define number_of_point 100

double FIXED_CURRENT_LEN;


struct point
{
	double x,y,z;
};

class Point3D {

public:
    double x,y,z;

};

class Vector {
public:
    double x,y,z;

};

Point3D pos;
Vector u,r,l;

Point3D points[number_of_point];

Point3D wCenter;
Vector wFront;
double wheelAngleX;
double wheelRotationAngle;


Vector get_cross_Product(Vector a, Vector b)
{

   Vector result;
    result.x = a.y*b.z - b.y*a.z;
    result.y = a.z*b.x - b.z*a.x;
    result.z = a.x*b.y - b.x*a.y;
    return result;
}


Vector get_rotate_vector(Vector v , Vector refer_vect , double rotational_angle)
{

    Vector result_vect , perpendicular_vect;
    double rotational_angle_in_radian = rotational_angle*pi/180.0;

    perpendicular_vect = get_cross_Product(refer_vect,v);


    result_vect.x = v.x*cos(rotational_angle_in_radian) + perpendicular_vect.x * sin(rotational_angle_in_radian);
    result_vect.y = v.y*cos(rotational_angle_in_radian) + perpendicular_vect.y * sin(rotational_angle_in_radian);
    result_vect.z = v.z*cos(rotational_angle_in_radian) + perpendicular_vect.z * sin(rotational_angle_in_radian);

    return result_vect;

}


void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);

			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
		}glEnd();
	}
}


void drawGrid()
{
	int i;
	if(drawgrid==1)
	{
		glColor3f(0.6, 0.6, 0.6);	//grey
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -90, 0);
				glVertex3f(i*10,  90, 0);

				//lines parallel to X-axis
				glVertex3f(-90, i*10, 0);
				glVertex3f( 90, i*10, 0);
			}
		}glEnd();
	}
}

void Axle_YZ()
{
      glColor3f(c1,c2,c3);
    glBegin(GL_QUADS);
    {
        glVertex3f(0,AXLE_WIDTH,-WHEEL_RADIUS);
        glVertex3f(0,AXLE_WIDTH,WHEEL_RADIUS);
        glVertex3f(0,-AXLE_WIDTH,-WHEEL_RADIUS);
        glVertex3f(0,-AXLE_WIDTH,WHEEL_RADIUS);

    } glEnd();
}

void Axle_XY()
{
      glColor3f(c1,c2,c3);
    glBegin(GL_QUADS);
    {
        glVertex3f(-WHEEL_RADIUS,  AXLE_WIDTH, 0);
        glVertex3f( WHEEL_RADIUS, -AXLE_WIDTH, 0);
        glVertex3f( WHEEL_RADIUS,  AXLE_WIDTH, 0);
        glVertex3f(-WHEEL_RADIUS, -AXLE_WIDTH, 0);
    } glEnd();
}

void draw_XAxle()
{
    Axle_YZ();
    Axle_XY();
}

double x_coordinate(int i)
{
    double x_value;
    double total_angel = 2*pi;
    double ratio_angel = ((double)i/(double)SLICES);
    double theta = ratio_angel*total_angel;
    x_value = WHEEL_RADIUS*cos(theta);
    return x_value;
}

double y_coordinate(int i)
{
    double y_value;
    double total_angel = 2*pi;
    double ratio_angel = ((double)i/(double)SLICES);
    double theta = ratio_angel*total_angel;
    y_value = WHEEL_RADIUS*sin(theta);
    return y_value;
}

void point_generation()
{
   int i;
    for(i=0;i<=SLICES;i++){
        points[i].x = x_coordinate(i);
        points[i].y = y_coordinate(i);
    }
}

void segment_generation()
{
    int lower_half = -STACKS/2;
    int upper_half = STACKS/2;
    double width_per_stacks = WHEEL_WIDTH/STACKS;


    for(int k = lower_half ; k < upper_half; k++){
        for(int i=0;i<SLICES;i++){
            double color_value = ((double)i/(double)SLICES);

            glColor3f(color_value,color_value,color_value);
            glBegin(GL_QUADS);
            {

                double z2_val = (k+1)*WHEEL_WIDTH/STACKS ;
                glVertex3f(points[i+1].x,points[i+1].y,z2_val);
                glVertex3f(points[i].x,points[i].y,z2_val);

                double z1_val = k*width_per_stacks;
                glVertex3f(points[i].x,points[i].y,z1_val);
                glVertex3f(points[i+1].x,points[i+1].y,z1_val);
            }
            glEnd();
        }
    }
}
void draw_Wheel()
{
    glColor3f(c1,c2,c3);
    point_generation();
    segment_generation();

}

double ChangedValue(double v)
{
    return v*MOVE_DISTANCE;
}
void CenterChangeForForward()
{
   wCenter.x += ChangedValue(wFront.x);
   wCenter.y += ChangedValue(wFront.y);
   wCenter.z += ChangedValue(wFront.z);

}
void GoForward()
{
    wheelRotationAngle += MOVE_ANGLE;
}
void ForwardOperation()
{
     CenterChangeForForward();
     GoForward();

}

void CenterChangeForBackward()
{
    wCenter.x -= ChangedValue(wFront.x);
    wCenter.y -= ChangedValue(wFront.y);
    wCenter.z -= ChangedValue(wFront.z);
}
void GoBackward()
{
   wheelRotationAngle -= MOVE_ANGLE;
}
void BackwardOperation()
{
     CenterChangeForBackward();
     GoBackward();

}

void LeftRotate()
{
    wheelAngleX += ANGLE_CHANGE_FOR_ROTATION;

}

Vector get_Z_unit_vector()
{
    Vector reference_vect;

    reference_vect.x = 0;
    reference_vect.y = 0;
    reference_vect.z = 1;

    return reference_vect;
}

void ChangedDirectionForLeftRotate()
{
        Vector reference_vect = get_Z_unit_vector();
        wFront = get_rotate_vector(wFront, reference_vect,ANGLE_CHANGE_FOR_ROTATION);
}

void RotateLeftOperation()
{
    LeftRotate();
    ChangedDirectionForLeftRotate();
}


void RightRotate()
{
    wheelAngleX -= ANGLE_CHANGE_FOR_ROTATION;
}

void ChangedDirectionForRightRotate()
{
    Vector reference_vect = get_Z_unit_vector();
    wFront = get_rotate_vector(wFront, reference_vect, -ANGLE_CHANGE_FOR_ROTATION);
}

void RotateRightOperation()
{
    RightRotate();
    ChangedDirectionForRightRotate();
}

void keyboardListener(unsigned char key, int x,int y)
{
    switch(key)
    {
    case 'w':
        ForwardOperation();
        break;

    case 's':
        BackwardOperation();
        break;

    case 'a':
        RotateLeftOperation();
        break;

    case 'd':
        RotateRightOperation();
        break;

    default:
        break;
    }
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:
			cameraHeight -= changed_cameraHight;
			break;
		case GLUT_KEY_UP:
			cameraHeight += changed_cameraHight;
			break;
		case GLUT_KEY_RIGHT:
			cameraAngle += changed_cameraAngle;
			break;
		case GLUT_KEY_LEFT:
			cameraAngle -= changed_cameraAngle;
			break;
		case GLUT_KEY_PAGE_UP:
			break;
		case GLUT_KEY_PAGE_DOWN:
			break;
		case GLUT_KEY_INSERT:
			break;
		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}

void SetUpCamera()
{
    gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
}


void Center_setup()
{
    glTranslatef(wCenter.x, wCenter.y, wCenter.z);

}
void Z_Axis_Reference_Rotation()
{
    glRotatef(wheelAngleX, 0, 0, 1);
}

void Y_Axis_Reference_Rotation()
{
    glRotatef(wheelRotationAngle, 0, 1, 0);
}


void property_print()
{
    cout<<"CameraHeight: "<<cameraHeight<<endl;
    cout<<"CameraAngle: " << cameraAngle<<endl;

}

void wCenter_print()
{

    cout<<"wCenter.x: "<<wCenter.x<<endl;
    cout<<"wCenter.y: "<<wCenter.y<<endl;
    cout<<"wCenter.z: "<<wCenter.z<<endl;

}
void wFront_print()
{

    cout<<"wCenter.x: "<<wFront.x<<endl;
    cout<<"wCenter.y: "<<wFront.y<<endl;
    cout<<"wCenter.z: "<<wFront.z<<endl;

}
void property_print1()
{
    cout<<"wheelAngleX:"<<wheelAngleX<<endl;
    cout<<"wheelRotationAngle"<<wheelRotationAngle<<endl;
}

void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);
	//initialize the matrix
	glLoadIdentity();


	SetUpCamera();

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/

	drawAxes();
	drawGrid();

    Center_setup();
    Z_Axis_Reference_Rotation();
    Y_Axis_Reference_Rotation();

	draw_XAxle();
	glRotatef(90,1,0,0);
	draw_Wheel();


	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){

	glutPostRedisplay();
}

void SetInitialValues()
{
    drawgrid=1;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;

    wFront.x = 1;
    wFront.y = 0;
    wFront.z = 0;

    wCenter.x = 0;
	wCenter.y = 0;
	wCenter.z = WHEEL_RADIUS;


    wheelRotationAngle = 0;
	wheelAngleX = 0;

}

void init(){

	/************************
	/ set-up projection here
	************************/

    SetInitialValues();

    //clear the screen
	glClearColor(0,0,0,0);

	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(80,	1,	1,	1000.0);

}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
