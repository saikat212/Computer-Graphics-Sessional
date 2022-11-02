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
int drawgrid;
int drawaxes;
double angle;

#define MAXIMUM_LENGTH 20.0
#define LENGTH_CHANGE 0.7

#define CHANGED_ANGLE 0.8
#define POS_CHANGE 3

#define SLICES 50
#define STACKS 50


double FIXED_CURRENT_LEN;


struct point
{
    double x,y,z;
};

class Point3D
{

public:
    double x,y,z;

};

class Vector
{
public:
    double x,y,z;

};

Point3D pos;
Vector u,r,l;
Point3D points[100][100];
Point3D cylinder_points[100][100];

void x_ref_rotation(double theta)
{
    glRotatef(theta,1,0,0);
}
void y_ref_rotation(double theta)
{
    glRotatef(theta,0,1,0);
}
void z_ref_rotation(double theta)
{
    glRotatef(theta,0,0,1);
}


Vector get_cross_Product(Vector a, Vector b)
{

    Vector result;
    result.x = a.y*b.z - b.y*a.z;
    result.y = a.z*b.x - b.z*a.x;
    result.z = a.x*b.y - b.x*a.y;
    return result;
}


Vector get_rotate_vector(Vector v, Vector refer_vect, double rotational_angle)
{

    Vector result_vect, perpendicular_vect;
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
        glBegin(GL_LINES);
        {
            glVertex3f( 100,0,0);
            glVertex3f(-100,0,0);

            glVertex3f(0,-100,0);
            glVertex3f(0, 100,0);

            glVertex3f(0,0, 100);
            glVertex3f(0,0,-100);
        }
        glEnd();
    }
}


void drawGrid()
{
    int i;
    if(drawgrid==1)
    {
        glColor3f(0.6, 0.6, 0.6);	//grey
        glBegin(GL_LINES);
        {
            for(i=-8; i<=8; i++)
            {

                if(i==0)
                    continue;	//SKIP the MAIN axes

                //lines parallel to Y-axis
                glVertex3f(i*10, -90, 0);
                glVertex3f(i*10,  90, 0);

                //lines parallel to X-axis
                glVertex3f(-90, i*10, 0);
                glVertex3f( 90, i*10, 0);
            }
        }
        glEnd();
    }
}

void drawSquare(double a)
{
    //glColor3f(1.0,0.0,0.0);
    glBegin(GL_QUADS);
    {
        glVertex3f( a,a,0);
        glVertex3f( a,-a,0);
        glVertex3f(-a,-a,0);
        glVertex3f(-a, a,0);
    }
    glEnd();
}

double get_h_value(double radius,int i)
{
    double h_value;
    double factor = (pi*(0.5));
    double ratio_val = ((double)i/(double)STACKS);
    double theta = ratio_val*factor;
    h_value = radius*sin(theta);
    return h_value;
}

double get_r_value(double radius,int i)
{
    double r_value;
    double factor = (pi*(0.5));
    double ratio_val = ((double)i/(double)STACKS);
    double theta = ratio_val*factor;
    r_value = radius*cos(theta);
    return r_value;
}


double get_x_value(double radius,int i)
{
    double r_value;
    double factor = (pi*(0.5));
    double ratio_val = ((double)i/(double)SLICES);
    double theta = ratio_val*factor;
    r_value = radius*cos(theta);
    return r_value;
}

double get_y_value(double radius,int i)
{
    double r_value;
    double factor = (pi*(0.5));
    double ratio_val = ((double)i/(double)SLICES);
    double theta = ratio_val*factor;
    r_value = radius*sin(theta);
    return r_value;
}
void point_generate(double radius)
{
    int i,j;
    double h,r;
    double factor = (pi*(0.5));
    //generate points
    for(i=0; i<=STACKS; i++)
    {
        h = get_h_value(radius,i);
        r = get_r_value(radius,i);

        for(j=0; j<=SLICES; j++)
        {
            points[i][j].x = get_x_value(r,j);
            points[i][j].y = get_y_value(r,j);
            points[i][j].z = h;
        }
    }
}

void upper_part(int i, int j)
{
    glBegin(GL_QUADS);
    {

        glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);

        glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);

        glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);

        glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);

    }
    glEnd();
}

void draw_quads_using_generated_points()
{
    for(int i=0; i<STACKS; i++)
    {
        for(int j=0; j<SLICES; j++)
        {
            upper_part(i,j);
        }
    }
}

void draw_part_of_Sphere(double radius)
{

    point_generate(radius);
    draw_quads_using_generated_points();

}
double get_cylinder_h_value(double height,int i)
{
    double ratio_val = (double)i/(double)STACKS;
    return height*ratio_val;
}


void cylinder_points_generation(double radius, double height)
{
    int i,j;
    double h;
    for(i=0; i<=STACKS; i++)
    {
        h = get_cylinder_h_value(height,i);

        for(j=0; j<=SLICES; j++)
        {
            cylinder_points[i][j].x= get_x_value(radius,j);
            cylinder_points[i][j].y= get_y_value(radius,j);
            cylinder_points[i][j].z=h;
        }
    }
}

void draw_cylinder_quads_using_generated_points()
{

    for(int i=0; i<STACKS; i++)
    {
        for(int j=0; j<SLICES; j++)
        {
            glBegin(GL_QUADS);
            {
                glVertex3f(cylinder_points[i][j].x,cylinder_points[i][j].y,cylinder_points[i][j].z);
                glVertex3f(cylinder_points[i][j+1].x,cylinder_points[i][j+1].y,cylinder_points[i][j+1].z);
                glVertex3f(cylinder_points[i+1][j+1].x,cylinder_points[i+1][j+1].y,cylinder_points[i+1][j+1].z);
                glVertex3f(cylinder_points[i+1][j].x,cylinder_points[i+1][j].y,cylinder_points[i+1][j].z);
                //lower half
                glVertex3f(cylinder_points[i][j].x,cylinder_points[i][j].y,-cylinder_points[i][j].z);
                glVertex3f(cylinder_points[i][j+1].x,cylinder_points[i][j+1].y,-cylinder_points[i][j+1].z);
                glVertex3f(cylinder_points[i+1][j+1].x,cylinder_points[i+1][j+1].y,-cylinder_points[i+1][j+1].z);
                glVertex3f(cylinder_points[i+1][j].x,cylinder_points[i+1][j].y,-cylinder_points[i+1][j].z);
            }
            glEnd();
        }
    }
}
void draw_part_of_Cylinder(double radius, double height)
{
    cylinder_points_generation(radius,height);
    draw_cylinder_quads_using_generated_points();

}



void upper_hemi_sphere()
{
     glColor3f(1, 0, 0);

    glPushMatrix();
    glTranslatef(FIXED_CURRENT_LEN,FIXED_CURRENT_LEN,FIXED_CURRENT_LEN);
    draw_part_of_Sphere(MAXIMUM_LENGTH - FIXED_CURRENT_LEN);
    glPopMatrix();


    glPushMatrix();
    z_ref_rotation(90);
    glTranslatef(FIXED_CURRENT_LEN,FIXED_CURRENT_LEN,FIXED_CURRENT_LEN);
    draw_part_of_Sphere(MAXIMUM_LENGTH - FIXED_CURRENT_LEN);
    glPopMatrix();


    glPushMatrix();
    z_ref_rotation(180);
    glTranslatef(FIXED_CURRENT_LEN,FIXED_CURRENT_LEN,FIXED_CURRENT_LEN);
    draw_part_of_Sphere(MAXIMUM_LENGTH - FIXED_CURRENT_LEN);
    glPopMatrix();


    glPushMatrix();
    z_ref_rotation(270);
    glTranslatef(FIXED_CURRENT_LEN,FIXED_CURRENT_LEN,FIXED_CURRENT_LEN);
    draw_part_of_Sphere(MAXIMUM_LENGTH - FIXED_CURRENT_LEN);
    glPopMatrix();


}
void lower_hemi_sphere()
{

    glPushMatrix();
    glTranslatef(FIXED_CURRENT_LEN, FIXED_CURRENT_LEN, -FIXED_CURRENT_LEN);
    glRotatef(180,1,1,0);
    draw_part_of_Sphere(MAXIMUM_LENGTH - FIXED_CURRENT_LEN);
    glPopMatrix();


    glPushMatrix();
    glRotatef(90,0,0,1);
    glTranslatef(FIXED_CURRENT_LEN, FIXED_CURRENT_LEN, -FIXED_CURRENT_LEN);
    glRotatef(180,1,1,0);
    draw_part_of_Sphere(MAXIMUM_LENGTH - FIXED_CURRENT_LEN);
    glPopMatrix();


    glPushMatrix();
    glRotatef(180,0,0,1);
    glTranslatef(FIXED_CURRENT_LEN, FIXED_CURRENT_LEN, -FIXED_CURRENT_LEN);
    glRotatef(180,1,1,0);
    draw_part_of_Sphere(MAXIMUM_LENGTH - FIXED_CURRENT_LEN);
    glPopMatrix();


    glPushMatrix();
    glRotatef(270,0,0,1);
    glTranslatef(FIXED_CURRENT_LEN, FIXED_CURRENT_LEN, -FIXED_CURRENT_LEN);
    glRotatef(180,1,1,0);
    draw_part_of_Sphere(MAXIMUM_LENGTH - FIXED_CURRENT_LEN);
    glPopMatrix();

}

void draw_final_Sphere()
{

    upper_hemi_sphere();
    lower_hemi_sphere();

}

void cylinder_sides()
{
    glPushMatrix();
    glTranslatef(FIXED_CURRENT_LEN,FIXED_CURRENT_LEN,0);
    draw_part_of_Cylinder(MAXIMUM_LENGTH - FIXED_CURRENT_LEN, FIXED_CURRENT_LEN);
    glPopMatrix();


    glPushMatrix();
    z_ref_rotation(90);
    glTranslatef(FIXED_CURRENT_LEN,FIXED_CURRENT_LEN,0);
    draw_part_of_Cylinder(MAXIMUM_LENGTH - FIXED_CURRENT_LEN, FIXED_CURRENT_LEN);
    glPopMatrix();

    glPushMatrix();
    z_ref_rotation(180);
    glTranslatef(FIXED_CURRENT_LEN,FIXED_CURRENT_LEN,0);
    draw_part_of_Cylinder(MAXIMUM_LENGTH - FIXED_CURRENT_LEN, FIXED_CURRENT_LEN);
    glPopMatrix();

    glPushMatrix();
    z_ref_rotation(270);
    glTranslatef(FIXED_CURRENT_LEN,FIXED_CURRENT_LEN,0);
    draw_part_of_Cylinder(MAXIMUM_LENGTH - FIXED_CURRENT_LEN, FIXED_CURRENT_LEN);
    glPopMatrix();

}

void cylinder_upper_part()
{
    glPushMatrix();
    x_ref_rotation(45);
    z_ref_rotation(45);
    glTranslatef(FIXED_CURRENT_LEN,FIXED_CURRENT_LEN,0);
    glRotatef(90,1,1,0);
    draw_part_of_Cylinder(MAXIMUM_LENGTH - FIXED_CURRENT_LEN, FIXED_CURRENT_LEN);
    glPopMatrix();

    glPushMatrix();
    z_ref_rotation(90);
    x_ref_rotation(45);
    z_ref_rotation(45);

    glTranslatef(FIXED_CURRENT_LEN,FIXED_CURRENT_LEN,0);
    glRotatef(90,1,1,0);
    draw_part_of_Cylinder(MAXIMUM_LENGTH - FIXED_CURRENT_LEN, FIXED_CURRENT_LEN);
    glPopMatrix();

    glPushMatrix();
    z_ref_rotation(180);
    x_ref_rotation(45);
    z_ref_rotation(45);

    glTranslatef(FIXED_CURRENT_LEN,FIXED_CURRENT_LEN,0);
    glRotatef(90,1,1,0);
    draw_part_of_Cylinder(MAXIMUM_LENGTH - FIXED_CURRENT_LEN, FIXED_CURRENT_LEN);
    glPopMatrix();


    glPushMatrix();
    z_ref_rotation(270);
    x_ref_rotation(45);
    z_ref_rotation(45);

    glTranslatef(FIXED_CURRENT_LEN,FIXED_CURRENT_LEN,0);
    glRotatef(90,1,1,0);
    draw_part_of_Cylinder(MAXIMUM_LENGTH - FIXED_CURRENT_LEN, FIXED_CURRENT_LEN);
    glPopMatrix();

}

void cylinder_lower_part()
{
    glPushMatrix();
    glRotatef(-45,1,0,0);
    glRotatef(45,0,0,1);
    glTranslatef(FIXED_CURRENT_LEN,FIXED_CURRENT_LEN,0);
    glRotatef(90,1,1,0);
    draw_part_of_Cylinder(MAXIMUM_LENGTH - FIXED_CURRENT_LEN, FIXED_CURRENT_LEN);
    glPopMatrix();


    glPushMatrix();
    glRotatef(90,0,0,1);
    glRotatef(-45,1,0,0);
    glRotatef(45,0,0,1);
    glTranslatef(FIXED_CURRENT_LEN,FIXED_CURRENT_LEN,0);
    glRotatef(90,1,1,0);
    draw_part_of_Cylinder(MAXIMUM_LENGTH - FIXED_CURRENT_LEN, FIXED_CURRENT_LEN);
    glPopMatrix();

    glPushMatrix();
    glRotatef(180,0,0,1);
    glRotatef(-45,1,0,0);
    glRotatef(45,0,0,1);
    glTranslatef(FIXED_CURRENT_LEN,FIXED_CURRENT_LEN,0);
    glRotatef(90,1,1,0);
    draw_part_of_Cylinder(MAXIMUM_LENGTH - FIXED_CURRENT_LEN, FIXED_CURRENT_LEN);
    glPopMatrix();

    glPushMatrix();
    glRotatef(270,0,0,1);
    glRotatef(-45,1,0,0);
    glRotatef(45,0,0,1);
    glTranslatef(FIXED_CURRENT_LEN,FIXED_CURRENT_LEN,0);
    glRotatef(90,1,1,0);
    draw_part_of_Cylinder(MAXIMUM_LENGTH - FIXED_CURRENT_LEN, FIXED_CURRENT_LEN);
    glPopMatrix();
}

void draw_final_Cylinder()
{
    glColor3f(0, 1, 0);   //green color

    cylinder_sides();
    cylinder_upper_part();
    cylinder_lower_part();


}


void top_up()
{
    glPushMatrix();
    glTranslatef(0,0,MAXIMUM_LENGTH);
    drawSquare(FIXED_CURRENT_LEN);
    glPopMatrix();
}

void top_down()
{
    glPushMatrix();
    glTranslatef(0,0,-MAXIMUM_LENGTH);
    drawSquare(FIXED_CURRENT_LEN);
    glPopMatrix();
}

void x_positive_side()
{
    glPushMatrix();
    glTranslatef(MAXIMUM_LENGTH,0,0);
    //glRotatef(90,0,1,0);
    y_ref_rotation(90);
    drawSquare(FIXED_CURRENT_LEN);
    glPopMatrix();
}

void x_negative_side()
{
    glPushMatrix();
    glTranslatef(-MAXIMUM_LENGTH,0,0);
    //glRotatef(90,0,1,0);
    y_ref_rotation(90);
    drawSquare(FIXED_CURRENT_LEN);
    glPopMatrix();
}

void y_positive_side()
{
    glPushMatrix();
    glTranslatef(0,MAXIMUM_LENGTH,0);
    //glRotatef(90,1,0,0);
    x_ref_rotation(90);
    drawSquare(FIXED_CURRENT_LEN);
    glPopMatrix();

}

void y_negative_side()
{
    glPushMatrix();
    glTranslatef(0,-MAXIMUM_LENGTH,0);
    x_ref_rotation(90);
    drawSquare(FIXED_CURRENT_LEN);
    glPopMatrix();
}
void draw_Cube()
{
    glColor3f(1,1,1);

    top_up();
    top_down();
    x_positive_side();
    x_negative_side();
    y_positive_side();
    y_negative_side();

}



void draw_combination_shape()
{

    draw_Cube();
    draw_final_Sphere();
    draw_final_Cylinder();
}

void rotate_look_left()
{
    l = get_rotate_vector(l,u,CHANGED_ANGLE);
    r = get_rotate_vector(r,u,CHANGED_ANGLE);
}

void rotate_look_right()
{
    l = get_rotate_vector(l,u, -CHANGED_ANGLE);
    r = get_rotate_vector(r,u, -CHANGED_ANGLE);
}

void look_up()
{
    u = get_rotate_vector(u,r,CHANGED_ANGLE);
    l = get_rotate_vector(l,r,CHANGED_ANGLE);
}

void look_down()
{
    u = get_rotate_vector(u,r, -CHANGED_ANGLE);
    l = get_rotate_vector(l,r, -CHANGED_ANGLE);
}

void tilt_clockwise()
{
    u = get_rotate_vector(u,l,CHANGED_ANGLE);
    r = get_rotate_vector(r,l,CHANGED_ANGLE);
}

void tilt_counterclockwise()
{
    u = get_rotate_vector(u,l, -CHANGED_ANGLE);
    r = get_rotate_vector(r,l, -CHANGED_ANGLE);
}
void keyboardListener(unsigned char key, int x,int y)
{
    switch(key)
    {

    case '1':
        rotate_look_left();
        break;

    case '2':
        rotate_look_right();
        break;
    case '3':
        look_up();

        break;

    case '4':
        look_down();

        break;

    case '5':
        tilt_clockwise();

        break;

    case '6':
        tilt_counterclockwise();
        break;

    default:
        break;
    }
}


void MoveForward()
{
    pos.x += l.x*POS_CHANGE;
    pos.y +=  l.y*POS_CHANGE;
    pos.z += l.z*POS_CHANGE;

}
void MoveBackward()
{
    pos.x += - l.x*POS_CHANGE;
    pos.y += - l.y*POS_CHANGE;
    pos.z += - l.z*POS_CHANGE;
}
void MoveRight()
{
    pos.x +=  r.x*POS_CHANGE;
    pos.y +=  r.y*POS_CHANGE;
    pos.z +=  r.z*POS_CHANGE;
}

void MoveLeft()
{
    pos.x +=  -r.x*POS_CHANGE;
    pos.y +=  -r.y*POS_CHANGE;
    pos.z +=  -r.z*POS_CHANGE;

}
void PgUp_MoveUP()
{
    pos.x +=  u.x*POS_CHANGE;
    pos.y +=  u.y*POS_CHANGE;
    pos.z +=  u.z*POS_CHANGE;
}

void PgDn_MoveDown()
{
    pos.x +=  -u.x*POS_CHANGE;
    pos.y +=  -u.y*POS_CHANGE;
    pos.z +=  -u.z*POS_CHANGE;
}

void CubeToSphere()
{
    FIXED_CURRENT_LEN -= LENGTH_CHANGE;
    if(FIXED_CURRENT_LEN < 0)
    {
        FIXED_CURRENT_LEN = 0;
    }
}
void SphereToCube()
{
    FIXED_CURRENT_LEN += LENGTH_CHANGE;

    if(FIXED_CURRENT_LEN > MAXIMUM_LENGTH)
    {
        FIXED_CURRENT_LEN =  MAXIMUM_LENGTH;
    }
}



double x_coordinate(int i,double r)
{
    double x_value;
    double total_angel = 2*pi;
    double ratio_angel = ((double)i/(double)SLICES);
    double theta = ratio_angel*total_angel;
    x_value =r*cos(theta);
    return x_value;
}

double y_coordinate(int i,double r)
{
    double y_value;
    double total_angel = 2*pi;
    double ratio_angel = ((double)i/(double)SLICES);
    double theta = ratio_angel*total_angel;
    y_value = r*sin(theta);
    return y_value;
}

void position_print()
{
    cout<<"posx: "<<pos.x<<endl;
    cout<<"posy: "<<pos.y<<endl;
    cout<<"posz: "<<pos.z<<endl;
}

void u_print()
{
    cout<<"ux: "<<u.x<<endl;
    cout<<"uy: "<<u.y<<endl;
    cout<<"uz: "<<u.z<<endl;
}
void r_print()
{
    cout<<"rx: "<<r.x<<endl;
    cout<<"ry: "<<r.y<<endl;
    cout<<"rz: "<<r.z<<endl;
}
void l_print()
{
    cout<<"lx: "<<l.x<<endl;
    cout<<"ly: "<<l.y<<endl;
    cout<<"lz: "<<l.z<<endl;
}

void property_print()
{
    cout<<"CameraHeight: "<<cameraHeight<<endl;
    cout<<"CameraAngle: " << cameraAngle<<endl;

}

void specialKeyListener(int key, int x,int y)
{
    switch(key)
    {
    case GLUT_KEY_DOWN:
        MoveBackward();
        break;


    case GLUT_KEY_UP:
        MoveForward();
        break;

    case GLUT_KEY_RIGHT:
        MoveRight();
        break;

    case GLUT_KEY_LEFT:
        MoveLeft();

        break;

    case GLUT_KEY_PAGE_UP:
        PgUp_MoveUP();

        break;
    case GLUT_KEY_PAGE_DOWN:
        PgDn_MoveDown();
        break;

    case GLUT_KEY_INSERT:
        break;

    case GLUT_KEY_HOME:

        CubeToSphere();
        break;

    case GLUT_KEY_END:

        SphereToCube();
        break;

    default:
        break;
    }
}


void mouseListener(int button, int state, int x, int y) 	//x, y is the x-y of the screen (2D)
{
    switch(button)
    {
    case GLUT_LEFT_BUTTON:
        if(state == GLUT_DOWN) 		// 2 times?? in ONE click? -- solution is checking DOWN or UP
        {
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



void display()
{

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

    gluLookAt(pos.x,pos.y,pos.z, pos.x+l.x,pos.y+l.y,pos.z+l.z,   u.x,u.y,u.z);


    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/
    //add objects

    drawAxes();
    drawGrid();
    draw_combination_shape();
    //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}


void animate()
{
    angle+=0.05;
    //codes for any changes in Models, Camera
    glutPostRedisplay();
}

void pos_set()
{
    pos.x = 100;
    pos.y = 100;
    pos.z = 50;
}
void u_set()
{
    u.x = 0;
    u.y = 0;
    u.z = 1;
}

void l_set()
{

    l.x = -1.0/sqrt(2);
    l.y = -1.0/sqrt(2);
    l.z = 0;


}
void r_set()
{
    r.x = -1.0/sqrt(2);
    r.y = 1.0/sqrt(2);
    r.z = 0;
}
void init()
{
    //codes for initialization
    drawgrid=0;
    drawaxes=1;
    cameraHeight=150.0;
    cameraAngle=1.0;
    angle=0;

    /************************
    / set-up projection here
    ************************/

    pos_set();
    r_set();
    u_set();
    l_set();
    FIXED_CURRENT_LEN = MAXIMUM_LENGTH*0.50;

    //clear the screen
    glClearColor(0,0,0,0);

    //load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    //initialize the matrix
    glLoadIdentity();

    //give PERSPECTIVE parameters
    gluPerspective(80,	1,	1,	1000.0);
    //field of view in the Y (vertically)
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
}

int main(int argc, char **argv)
{
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
