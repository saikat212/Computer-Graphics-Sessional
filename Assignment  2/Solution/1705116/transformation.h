#include <iomanip>
#include <cmath>
#include <string>
#include <cstdlib>
#include <ctime>
#include <limits>

using namespace std;
#define INF numeric_limits<double>::infinity()

#include "point.h"
#define PI 2.0 * acos(0.0)

class Transformation
{
public:
    double matrix[4][4];

    void IdentityMatrixGeneration()
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                if (i == j)
                {
                    matrix[i][j] = 1.0;
                }
                else
                {
                    matrix[i][j] = 0.0;
                }
            }
        }
    }

    Transformation()
    {
        IdentityMatrixGeneration();
    }

    Point RodriguesFormula(Point x, Point a, double theta)
    {
        // x * cos(theta * PI / 180.0)
        double factor1 = cos(theta * PI / 180.0);
        Point p1 = vector_scalar_multiplication(x, factor1);

        // a * (a * x) * (1 - cos(theta * PI / 180.0))
        double factor2 = 1 - factor1;
        double a_dot_x = vector_dot_multiplication(a, x);
        double factor3 = a_dot_x * factor2;
        Point p2 = vector_scalar_multiplication(a, factor3);

        // (a ^ x) * sin(theta * PI / 180.0)
        double factor4 = sin(theta * PI / 180.0);
        Point a_cross_x = vector_cross_multiplication(a, x);
        Point p3 = vector_scalar_multiplication(a_cross_x, factor4);

        Point temp = AdditionOfPoint(p1, p2);
        Point result_point = AdditionOfPoint(temp, p3);

        return result_point;
    }

    void MatrixMultiplication(Transformation transformation)
    {

        int rows = 4;
        int cols = 4;
        double temp[rows][cols];

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                temp[i][j] = 0.0;

                for (int k = 0; k < rows; k++)
                {
                    double product = transformation.matrix[k][j] * matrix[i][k];
                    temp[i][j] += product;
                }
            }
        }
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                matrix[i][j] = temp[i][j];
            }
        }
    }

    void TranslationMatrixGeneration(double tx, double ty, double tz)
    {
        int i = 0;
        int j = 3;
        matrix[i][j] = tx;
        matrix[i + 1][j] = ty;
        matrix[i + 2][j] = tz;
    }

    void ScalingMatrixGeneration(double sx, double sy, double sz)
    {
        int i = 0;
        matrix[i][i] = sx;
        matrix[i + 1][i + 1] = sy;
        matrix[i + 2][i + 2] = sz;
    }

    void RotationMatrixGeneration(double angle, double ax, double ay, double az)
    {
        Point a;
        a.x = ax;
        a.y = ay;
        a.z = az;
        a.normalize();


        Point x_axis_point(1.0, 0.0, 0.0);
        Point y_axis_point(0.0, 1.0, 0.0);
        Point z_axis_point(0.0, 0.0, 1.0);

        Point c1 = RodriguesFormula(x_axis_point, a, angle);
        Point c2 = RodriguesFormula(y_axis_point, a, angle);
        Point c3 = RodriguesFormula(z_axis_point, a, angle);


        int i = 0;
        int j = 0;
        int m = 1;
        int n = 2;
        matrix[i][j] = c1.x;
        matrix[i + 1][j] = c1.y;
        matrix[i + 2][j] = c1.z;

        matrix[i][m] = c2.x;
        matrix[i + 1][m] = c2.y;
        matrix[i + 2][m] = c2.z;

        matrix[i][n] = c3.x;
        matrix[i + 1][n] = c3.y;
        matrix[i + 2][n] = c3.z;
    }

    void ViewMatrixGeneration(Point eye, Point look, Point up)
    {


        Point l = SubtractionOfPoint(look, eye);
        l.normalize();
        Point r = vector_cross_multiplication(l, up);
        r.normalize();
        Point u = vector_cross_multiplication(r, l);


        Transformation translationTransformation;
        translationTransformation.TranslationMatrixGeneration(-eye.x, -eye.y, -eye.z);


        matrix[0][0] = r.x;
        matrix[0][1] = r.y;
        matrix[0][2] = r.z;

        matrix[1][0] = u.x;
        matrix[1][1] = u.y;
        matrix[1][2] = u.z;

        matrix[2][0] = -l.x;
        matrix[2][1] = -l.y;
        matrix[2][2] = -l.z;


        MatrixMultiplication(translationTransformation);
    }

    void ProjectionMatrixGeneration(double fovY, double aspectRatio, double near, double far)
    {

        double fovX = aspectRatio * fovY;
        double factor1 = tan(fovX / 2.0 * PI / 180.0);
        double r = near * factor1;
        double factor2 = tan(fovY / 2.0 * PI / 180.0);
        double t = near * factor2;


        matrix[0][0] = near / r;
        matrix[1][1] = near / t;
        matrix[2][2] = -(far + near) / (far - near);
        matrix[3][3] = 0.0;
        matrix[2][3] = -(2.0 * far * near) / (far - near);
        matrix[3][2] = -1.0;
    }

    ~Transformation()
    {
        IdentityMatrixGeneration();
    }

    void showMatrix()
    {
        cout << "Transformation class matrix : " << endl;
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                cout << matrix[i][j] << "--";
            }
            cout << endl;
        }
    }
};

struct Color
{
    int redValue;
    int greenValue;
    int blueValue;
};

struct Triangle
{
    Point corners[3];
    Color rgb;
};

Transformation two_4_by_4_matrix_multiplication(Transformation m1, Transformation m2)
{

    Transformation temp;
    temp.MatrixMultiplication(m1);
    temp.MatrixMultiplication(m2);
    return temp;
}

Point one_4_by_4_and_one_4_by_1_matrix_multiplication(Transformation m, Point p)
{

    double temp[4];

    for (int i = 0; i < 4; i++)
    {
        temp[i] = 0.0;

        for (int j = 0; j < 4; j++)
        {
            temp[i] += m.matrix[i][j] * ((j == 0) ? p.x : ((j == 1) ? p.y : ((j == 2) ? p.z : p.w)));
        }
    }
    return Point(temp[0], temp[1], temp[2], temp[3]);
}

void print_point_in_file(Point p1, Point p2, Point p3)
{
}
