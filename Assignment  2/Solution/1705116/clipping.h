




class clipping
{

public:
    int screenWidth, screenHeight;
    double leftLimitX, rightLimitX, bottomLimitY, topLimitY, frontLimitZ, rearLimitZ;
    clipping()
    {
         screenHeight=0;
         screenWidth =0;
         leftLimitX=0;
         rightLimitX =0;
         bottomLimitY =0;
         topLimitY =0;
         frontLimitZ =0;
         rearLimitZ =0;
    }


int
get_top_scanline(double maxY, double topY, double dy)
{
    int topScanline;
    if (maxY >= topY)
    {
        topScanline = 0;
    }
    else
    {
        double factor = (topY - maxY) / dy;
        topScanline = (int)round(factor);
    }

    return topScanline;
}


     void show_variable_status()
    {
                cout <<screenHeight<<endl;
                cout<<screenWidth<<endl;
                cout<<leftLimitX<<endl;
                cout<<rightLimitX<<endl;
                cout<<bottomLimitY<<endl;
                cout<<topLimitY<<endl;
                cout<<frontLimitZ<<endl;
                cout<<rearLimitZ<<endl;

        }


     void showMatrix()
    {
        cout << "Transformation class matrix : " << endl;
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                //cout << matrix[i][j] << "--";
            }
            cout << endl;
        }
    }


     void showzbuffer()
    {

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                //cout << matrix[i][j] << "--";
            }
            cout << endl;
        }
    }

 void showframebuffer()
    {

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                //cout << matrix[i][j] << "--";
            }
            cout << endl;
        }
    }


};

int
get_top_scanline(double maxY, double topY, double dy)
{
    int topScanline;
    if (maxY >= topY)
    {
        topScanline = 0;
    }
    else
    {
        double factor = (topY - maxY) / dy;
        topScanline = (int)round(factor);
    }

    return topScanline;
}

int get_bottom_scanline(double minY, double bottomY, double dy, double screenHeight)
{
    int bottomScanline;
    if (minY <= bottomY)
    {
        bottomScanline = screenHeight - 1;
    }
    else
    {
        double factor1 = (minY - bottomY) / dy;
        int factor2 = (1 + ((int)round(factor1)));
        bottomScanline = screenHeight - factor2;
    }

    return bottomScanline;
}

int get_left_inter_sec_col(double inter_sect_point_x, double leftX, double dx)
{
    int leftIntersectingColumn;
    if (inter_sect_point_x <= leftX)
    {
        leftIntersectingColumn = 0;
    }
    else
    {
        double factor = (inter_sect_point_x - leftX) / dx;
        leftIntersectingColumn = (int)round(factor);
    }

    return leftIntersectingColumn;
}

int get_right_inter_sec_col(double inter_sect_point_x_max, double rightX, double dx,double screenWidth)
{
            int rightIntersectingColumn;
            if ( inter_sect_point_x_max >= rightX)
            {
                rightIntersectingColumn = screenWidth - 1;
            }
            else
            {
                rightIntersectingColumn = screenWidth - (1 + ((int)round((rightX - inter_sect_point_x_max) / dx)));
            }

            return rightIntersectingColumn;
}



  double get_dx(double rightLimitX,double leftLimitX,double screenWidth)
    {
         return (rightLimitX - leftLimitX) / screenWidth;
    }


    double get_dy(double topLimitY,double bottomLimitY,double screenHeight)
    {
        return (topLimitY - bottomLimitY) / screenHeight;
    }


    double get_topY(double topLimitY,double dy)
    {
        return topLimitY - dy / 2.0;
    }


    double get_bottomY(double bottomLimitY,double dy)
    {
        return bottomLimitY + dy / 2.0;
    }


    double get_leftX(double leftLimitX,double dx)
    {
        return leftLimitX + dx / 2.0;
    }


    double get_rightX(double rightLimitX,double dx)
    {
        return rightLimitX - dx / 2.0;
    }

