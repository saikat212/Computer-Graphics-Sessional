

class Point
{
public:
    double x, y, z, w;

    Point()
    {
        x = 0.0;
        y = 0.0;
        z = 0.0;
        w = 1.0;
    }

    Point(double x1, double y1, double z1)
    {
        x = x1;
        y = y1;
        z = z1;
        w = 1.0;
    }

    Point(double x1, double y1, double z1, double w1)
    {
        x = x1;
        y = y1;
        z = z1;
        w = w1;
    }

    void normalize()
    {
        double val = pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0);
        double factor = sqrt(val);
        x = x / factor;
        y = y / factor;
        z = z / factor;
    }

    void scale()
    {
        x = x / w;
        y = y / w;
        z = z / w;
        w = w / w;
    }

    ~Point()
    {
        x = 0.0;
        y = 0.0;
        z = 0.0;
        w = 1.0;
    }
};

Point AdditionOfPoint(Point p1, Point p2)
{
    Point temp;
    temp.x = p1.x + p2.x;
    temp.y = p1.y + p2.y;
    temp.z = p1.z + p2.z;

    return temp;
}

Point SubtractionOfPoint(Point p1, Point p2)
{
    Point temp;
    temp.x = p1.x - p2.x;
    temp.y = p1.y - p2.y;
    temp.z = p1.z - p2.z;

    return temp;
}

Point vector_scalar_multiplication(Point p, double factor)
{
    Point temp;
    temp.x = factor * p.x;
    temp.y = factor * p.y;
    temp.z = factor * p.z;

    return temp;
}

double vector_dot_multiplication(Point p1, Point p2)
{

    double x_val = p1.x * p2.x;
    double y_val = p1.y * p2.y;
    double z_val = p1.z * p2.z;

    return x_val + y_val + z_val;
}

Point vector_cross_multiplication(Point p1, Point p2)
{
    Point temp;
    temp.x = p1.y * p2.z - p1.z * p2.y;
    temp.y = p1.z * p2.x - p1.x * p2.z;
    temp.z = p1.x * p2.y - p1.y * p2.x;

    return temp;
}
