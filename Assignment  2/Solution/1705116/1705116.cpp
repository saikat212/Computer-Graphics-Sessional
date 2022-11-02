#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <cstdlib>
#include <stack>
#include <ctime>
#include "bitmap_image.hpp"
#include <limits>
#include "transformation.h"
#include "clipping.h"


using namespace std;
#define INF numeric_limits<double>::infinity()



int main(int argc, char **argv)
{
    ifstream input;
    ofstream output;
    string testCaseDir = "4";

    /* ####  WARNING:: set current directory path of test-cases */

    string currentPath = "C:/Users/SAIKAT/Desktop/SAIKAT/BUET_4_1/CSE 410 - Computer Graphics Sessional/my/test-cases/";

    input.open(currentPath + testCaseDir + "/scene.txt");
    if (!input.is_open())
    {
        exit(EXIT_FAILURE);
        cout<<"input<<scene.txt file open failed."<<endl;
    }

    double eyeX, eyeY, eyeZ;
    double lookX, lookY, lookZ;
    double upX, upY, upZ;
    double fovY, aspectRatio, near, far;

    input >> eyeX >> eyeY >> eyeZ;
    input >> lookX >> lookY >> lookZ;
    input >> upX >> upY >> upZ;
    input >> fovY >> aspectRatio >> near >> far;

    stack<Transformation> MatrixStack;
    Transformation bottomMatrix;
    MatrixStack.push(bottomMatrix);

    output.open(currentPath + testCaseDir + "/stage1.txt");

    if (!output.is_open())
    {
        exit(EXIT_FAILURE);
        cout<<"output >>stage1.txt file open failed."<<endl;
    }

    int triangle_counter =0;
    int push_counter = 0;

    string display_command;

    while (true)
    {
        input >> display_command;

        if (display_command == "triangle")
        {
            Point p1;
            Point p2;
            Point p3;

            input >> p1.x >> p1.y >> p1.z;
            input >> p2.x >> p2.y >> p2.z;
            input >> p3.x >> p3.y >> p3.z;

            Transformation top_matrix = MatrixStack.top();
            p1 = one_4_by_4_and_one_4_by_1_matrix_multiplication(top_matrix, p1);
            p1.scale();

            top_matrix = MatrixStack.top();
            p2 = one_4_by_4_and_one_4_by_1_matrix_multiplication(top_matrix, p2);
            p2.scale();

            top_matrix = MatrixStack.top();
            p3 = one_4_by_4_and_one_4_by_1_matrix_multiplication(top_matrix, p3);
            p3.scale();

            output << fixed << setprecision(7) << p1.x << ' ' << p1.y << ' ' << p1.z << endl;
            output << fixed << setprecision(7) << p2.x << ' ' << p2.y << ' ' << p2.z << endl;
            output << fixed << setprecision(7) << p3.x << ' ' << p3.y << ' ' << p3.z << endl;
            output << endl;

            triangle_counter++;
        }

        else if (display_command == "push")
        {

            Transformation temp;
            Transformation top_matrix = MatrixStack.top();
            temp = two_4_by_4_matrix_multiplication(temp, top_matrix);
            MatrixStack.push(temp);
            push_counter++;
        }

        else if (display_command == "pop")
        {
            if (push_counter == 0)
            {
                cout << display_command << ": pop on empty stack" << endl;
                cout<<"pop failed."<<endl;
                break;
            }

            MatrixStack.pop();
            push_counter--;
        }

        else if (display_command == "scale")
        {
            double sx, sy, sz;
            input >> sx >> sy >> sz;

            Transformation scalingTransformation;
            scalingTransformation.ScalingMatrixGeneration(sx, sy, sz);

            Transformation top_matrix = MatrixStack.top();
            Transformation temp = two_4_by_4_matrix_multiplication(top_matrix, scalingTransformation);

            MatrixStack.pop();
            MatrixStack.push(temp);
        }

        else if (display_command == "translate")
        {
            double tx, ty, tz;
            input >> tx >> ty >> tz;

            Transformation translationTransformation;
            translationTransformation.TranslationMatrixGeneration(tx, ty, tz);
            Transformation top_matrix = MatrixStack.top();
            Transformation temp = two_4_by_4_matrix_multiplication(top_matrix, translationTransformation);

            MatrixStack.pop();
            MatrixStack.push(temp);
        }

        else if (display_command == "rotate")
        {
            double angle, ax, ay, az;
            input >> angle >> ax >> ay >> az;

            Transformation rotationTransformation;
            rotationTransformation.RotationMatrixGeneration(angle, ax, ay, az);

            Transformation top_matrix = MatrixStack.top();
            Transformation temp = two_4_by_4_matrix_multiplication(top_matrix, rotationTransformation);
            MatrixStack.pop();
            MatrixStack.push(temp);
        }

        else if (display_command == "end")
        {
            break;
        }
    }
    input.close();
    output.close();

    /* ############################# stage2 ############################################*/

    // input.open("./test-cases/"+testCaseDir+"/stage1.txt");
    input.open(currentPath + testCaseDir + "/stage1.txt");
    if (!input.is_open())
    {
        exit(EXIT_FAILURE);
         cout<<"input<<stage1.txt file open failed."<<endl;
    }

    // output.open("./test-cases/"+testCaseDir+"/stage2.txt");
    output.open(currentPath + testCaseDir + "/stage2.txt");
    if (!output.is_open())
    {
        exit(EXIT_FAILURE);
        cout<<"output >>stage2.txt file open failed."<<endl;
    }

    Transformation viewTransformation;
    Point eye(eyeX, eyeY, eyeZ);
    Point look(lookX, lookY, lookZ);
    Point up(upX, upY, upZ);
    viewTransformation.ViewMatrixGeneration(eye, look, up);

    for (int i = 0; i < triangle_counter; i++)
    {
        Point p1, p2, p3;

        input >> p1.x >> p1.y >> p1.z;
        input >> p2.x >> p2.y >> p2.z;
        input >> p3.x >> p3.y >> p3.z;

        p1 = one_4_by_4_and_one_4_by_1_matrix_multiplication(viewTransformation, p1);
        p1.scale();

        p2 = one_4_by_4_and_one_4_by_1_matrix_multiplication(viewTransformation, p2);
        p2.scale();

        p3 = one_4_by_4_and_one_4_by_1_matrix_multiplication(viewTransformation, p3);
        p3.scale();

        output << fixed << setprecision(7) << p1.x << ' ' << p1.y << ' ' << p1.z << endl;
        output << fixed << setprecision(7) << p2.x << ' ' << p2.y << ' ' << p2.z << endl;
        output << fixed << setprecision(7) << p3.x << ' ' << p3.y << ' ' << p3.z << endl;
        output << endl;
    }
    input.close();
    output.close();

    /* ############################# stage3 ############################################*/
    input.open(currentPath + testCaseDir + "/stage2.txt");
    if (!input.is_open())
    {
        exit(EXIT_FAILURE);
        cout<<"stage2.txt file open failed."<<endl;
    }

    output.open(currentPath + testCaseDir + "/stage3.txt");
    if (!output.is_open())
    {
        exit(EXIT_FAILURE);
        cout<<"stage3.txt file open failed."<<endl;
    }

    Transformation projectionTransformation;
    projectionTransformation.ProjectionMatrixGeneration(fovY, aspectRatio, near, far);

    for (int i = 0; i < triangle_counter; i++)
    {
        Point p1, p2, p3;

        input >> p1.x >> p1.y >> p1.z;
        input >> p2.x >> p2.y >> p2.z;
        input >> p3.x >> p3.y >> p3.z;

        p1 = one_4_by_4_and_one_4_by_1_matrix_multiplication(projectionTransformation, p1);
        p1.scale();
        p2 = one_4_by_4_and_one_4_by_1_matrix_multiplication(projectionTransformation, p2);
        p2.scale();
        p3 = one_4_by_4_and_one_4_by_1_matrix_multiplication(projectionTransformation, p3);
        p3.scale();

        output << fixed << setprecision(7) << p1.x << ' ' << p1.y << ' ' << p1.z << endl;
        output << fixed << setprecision(7) << p2.x << ' ' << p2.y << ' ' << p2.z << endl;
        output << fixed << setprecision(7) << p3.x << ' ' << p3.y << ' ' << p3.z << endl;
        output << endl;
    }
    input.close();
    output.close();

      /* ############################# stage4 ############################################*/

    input.open(currentPath + testCaseDir + "/config.txt");
    if (!input.is_open())
    {
        exit(EXIT_FAILURE);
        cout<<"config.txt file open failed."<<endl;
    }

    int screenWidth, screenHeight;
    double leftLimitX, rightLimitX, bottomLimitY, topLimitY, frontLimitZ, rearLimitZ;

    input >> screenWidth >> screenHeight;
    input >> leftLimitX;
    input >> bottomLimitY;
    input >> frontLimitZ >> rearLimitZ;

    input.close();

    rightLimitX = -leftLimitX;
    topLimitY = -bottomLimitY;

    input.open(currentPath + testCaseDir + "/stage3.txt");
    if (!input.is_open())
    {
        exit(EXIT_FAILURE);
    }

    Triangle triangles[triangle_counter];
    srand(time(0));

    for (int i = 0; i < triangle_counter; i++)
    {
        input >> triangles[i].corners[0].x >> triangles[i].corners[0].y >> triangles[i].corners[0].z;
        input >> triangles[i].corners[1].x >> triangles[i].corners[1].y >> triangles[i].corners[1].z;
        input >> triangles[i].corners[2].x >> triangles[i].corners[2].y >> triangles[i].corners[2].z;
        triangles[i].rgb.redValue = rand() % 256;
        triangles[i].rgb.greenValue = rand() % 256;
        triangles[i].rgb.blueValue = rand() % 256;
    }
    input.close();

    double dx, dy, topY, bottomY, leftX, rightX;

    dx = get_dx(rightLimitX,leftLimitX,screenWidth);
    dy = get_dy(topLimitY,bottomLimitY,screenHeight);
    topY = get_topY(topLimitY,dy);
    bottomY = get_bottomY(bottomLimitY,dy);
    leftX = get_leftX(leftLimitX,dx);
    rightX = get_rightX(rightLimitX,dx);

    double **zBuffer = new double *[screenHeight];
    Color **frameBuffer = new Color *[screenHeight];

    for (int i = 0; i < screenHeight; i++)
    {
        zBuffer[i] = new double[screenWidth];
    }

    for (int row = 0; row < screenHeight; row++)
    {
        for (int column = 0; column < screenWidth; column++)
        {
            zBuffer[row][column] = rearLimitZ;
        }
    }

    for (int i = 0; i < screenHeight; i++)
    {
        frameBuffer[i] = new Color[screenWidth];
    }

    for (int row = 0; row < screenHeight; row++)
    {
        for (int column = 0; column < screenWidth; column++)
        {
            frameBuffer[row][column].redValue = 0;
            frameBuffer[row][column].blueValue = 0;
            frameBuffer[row][column].greenValue = 0;
        }
    }


    int topScanline, bottomScanline;

    for (int i = 0; i < triangle_counter; i++)
    {

        double corner_1_y = triangles[i].corners[1].y;
        double corner_2_y = triangles[i].corners[2].y;
        double corner_0_y = triangles[i].corners[0].y;
        double max_between_1_2 = max(corner_1_y, corner_2_y);
        double min_between_1_2 = min(corner_1_y, corner_2_y);

        double maxY = max(corner_0_y, max_between_1_2);

        double minY = min(corner_0_y, min_between_1_2);

        topScanline = get_top_scanline(maxY, topY,dy);
        bottomScanline = get_bottom_scanline(minY, bottomY,dy,screenHeight);

        int leftIntersectingColumn, rightIntersectingColumn;
        for (int row = topScanline; row <= bottomScanline; row++)
        {

            double y_side = topY - row * dy;

            Point intersectingPoints_list[3];
            Point p1(INF, y_side, 0, 1);
            Point p2(INF, y_side, 1, 2);
            Point p3(INF, y_side, 2, 0);

            intersectingPoints_list[0] = p1;
            intersectingPoints_list[1] = p2;
            intersectingPoints_list[2] = p3;

            for (int j = 0; j < 3; j++)
            {
                Point p1;
                Point p2;

                int corner_index_from_z = (int)intersectingPoints_list[j].z;
                int corner_index_from_w = (int)intersectingPoints_list[j].w;
                p1 = triangles[i].corners[corner_index_from_z];
                p2 = triangles[i].corners[corner_index_from_w];

                if (p1.y == p2.y)
                {
                }
                else
                {

                    double factor1 = (y_side - p1.y) * (p1.x - p2.x);
                    double factor2 = (p1.y - p2.y);
                    intersectingPoints_list[j].x = p1.x + factor1 / factor2;
                }
            }
            for (int j = 0; j < 3; j++)
            {


                Point p1;
                Point p2;

                int corner_index_from_z = (int)intersectingPoints_list[j].z;
                int corner_index_from_w = (int)intersectingPoints_list[j].w;
                p1 = triangles[i].corners[corner_index_from_z];
                p2 = triangles[i].corners[corner_index_from_w];

                if (intersectingPoints_list[j].x != INF)
                {
                    double max_val_x = max(p1.x, p2.x);
                    double min_val_x = min(p1.x, p2.x);
                    double max_val_y = max(p1.y, p2.y);
                    double min_val_y = min(p1.y, p2.y);

                    double inter_sect_point_xvalue = intersectingPoints_list[j].x;
                    double inter_sect_point_yvalue = intersectingPoints_list[j].y;

                    if (inter_sect_point_xvalue > max_val_x || inter_sect_point_xvalue < min_val_x || inter_sect_point_yvalue > max_val_y || inter_sect_point_yvalue < min_val_y)
                    {

                        intersectingPoints_list[j].x = INF;
                    }
                }
            }


            int maxIndex = -1;
            int minIndex = -1;
            double maxX, minX;

            for (int j = 0; j < 3; j++)
            {
                if (maxIndex == -1 && minIndex == -1)
                {
                    if (intersectingPoints_list[j].x != INF)
                    {
                        maxIndex = j;
                        minIndex = j;
                        maxX = intersectingPoints_list[j].x;
                        minX = intersectingPoints_list[j].x;
                    }
                }
                else
                {
                    if (intersectingPoints_list[j].x != INF)
                    {

                        if (intersectingPoints_list[j].x > maxX)
                        {
                            maxX = intersectingPoints_list[j].x;
                            maxIndex = j;
                        }
                        if (intersectingPoints_list[j].x < minX)
                        {
                            minX = intersectingPoints_list[j].x;
                            minIndex = j;
                        }
                    }
                }
            }


            double inter_sect_point_x_min = intersectingPoints_list[minIndex].x;
            leftIntersectingColumn = get_left_inter_sec_col(inter_sect_point_x_min,leftX,dx);

            double inter_sect_point_x_max = intersectingPoints_list[maxIndex].x;
            rightIntersectingColumn = get_right_inter_sec_col(inter_sect_point_x_max,rightX,dx,screenHeight);

            int corner_index_from_z = (int)intersectingPoints_list[minIndex].z;
            int corner_index_from_w = (int)intersectingPoints_list[minIndex].w;
            p1 = triangles[i].corners[corner_index_from_z];
            p2 = triangles[i].corners[corner_index_from_w];
            double factor1 = (intersectingPoints_list[minIndex].y - p1.y);
            double factor2 = (p2.z - p1.z);
            double factor3 = (p2.y - p1.y);
            double za = p1.z + factor1*factor2 /factor3;

            corner_index_from_z = (int)intersectingPoints_list[maxIndex].z;
            corner_index_from_w = (int)intersectingPoints_list[maxIndex].w;
            p1 = triangles[i].corners[corner_index_from_z];
            p2 = triangles[i].corners[corner_index_from_w];
            double factor4 = (intersectingPoints_list[maxIndex].y - p1.y);
            double factor5 = (p2.z - p1.z);
            double factor6 = (p2.y - p1.y);
            double zb = p1.z + factor4*factor5 /factor6;


            double zp;
            double term = (intersectingPoints_list[maxIndex].x - intersectingPoints_list[minIndex].x);
            double constantTerm = dx * (zb - za) / term;

            for (int column = leftIntersectingColumn; column <= rightIntersectingColumn; column++)
            {
                if (column == leftIntersectingColumn)
                {
                    double f1 = ((leftX + leftIntersectingColumn * dx) - intersectingPoints_list[minIndex].x);
                    double f2 =  (intersectingPoints_list[maxIndex].x - intersectingPoints_list[minIndex].x);
                    zp = za + f1*(zb - za) / f2;
                }
                else
                {
                    zp = constantTerm + zp;
                }

                if (zp > frontLimitZ)
                {
                    if(zp < zBuffer[row][column])
                    {

                    zBuffer[row][column] = zp;
                    frameBuffer[row][column].blueValue = triangles[i].rgb.blueValue;
                    frameBuffer[row][column].greenValue = triangles[i].rgb.greenValue;
                    frameBuffer[row][column].redValue = triangles[i].rgb.redValue;


                    }

                }
            }
        }
    }

    bitmap_image bitmapImage(screenWidth, screenHeight);

    for (int row = 0; row < screenHeight; row++)
    {
        for (int column = 0; column < screenWidth; column++)
        {
            bitmapImage.set_pixel(column, row, frameBuffer[row][column].redValue, frameBuffer[row][column].greenValue, frameBuffer[row][column].blueValue);
        }
    }
    bitmapImage.save_image(currentPath + testCaseDir + "/out.bmp");


    output.open(currentPath + testCaseDir + "/zbuffer.txt");
    if (!output.is_open())
    {
        exit(EXIT_FAILURE);
        cout<<"zbuffer.txt file open failed"<<endl;
    }

    for (int row = 0; row < screenHeight; row++)
    {
        for (int column = 0; column < screenWidth; column++)
        {
            if (zBuffer[row][column] < rearLimitZ)
            {

                //cout << zBuffer[row][column] << '\t';
                output << fixed << setprecision(7) <<zBuffer[row][column] << '\t';
            }
        }
        output << endl;
        //cout<<endl;
    }
    output.close();

    for (int i = 0; i < screenHeight; i++)
    {
        delete[] frameBuffer[i];
    }
    delete[] frameBuffer;
    cout<<"frameBuffer memory freeing successfully"<<endl;

    for (int i = 0; i < screenHeight; i++)
    {

        delete[] zBuffer[i];
    }
    delete[] zBuffer;
    cout<<"zBuffer memory freeing successfully"<<endl;

    return 0;
}
