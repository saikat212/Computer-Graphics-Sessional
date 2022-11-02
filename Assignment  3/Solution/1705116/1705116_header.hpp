#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <limits>

using namespace std;

#define PI 2 * acos(0.0)
#define INF numeric_limits<double>::infinity()

class Vector
{
public:
    double x, y, z;

    Vector()
    {
        x = 0.0;
        y = 0.0;
        z = 0.0;
    }

    Vector(double x1, double y1, double z1)
    {
        x = x1;
        y = y1;
        z = z1;
    }

    void normalize()
    {
        double val = pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0);
        double factor = sqrt(val);
        x = x / factor;
        y = y / factor;
        z = z / factor;
    }

    double computeDistanceBetween(Vector _vector)
    {
        double t1 = pow(x - _vector.x, 2.0);
        double t2 = pow(y - _vector.y, 2.0);
        double t3 = pow(z - _vector.z, 2.0);
        return sqrt(t1 + t2 + t3);
    }

    ~Vector()
    {
        x = 0.0;
        y = 0.0;
        z = 0.0;
    }


};


Vector AdditionOfVector(Vector p1, Vector p2)
{
    Vector temp;
    temp.x = p1.x + p2.x;
    temp.y = p1.y + p2.y;
    temp.z = p1.z + p2.z;

    return temp;
}

Vector SubtractionOfVector(Vector p1, Vector p2)
{
    Vector temp;
    temp.x = p1.x - p2.x;
    temp.y = p1.y - p2.y;
    temp.z = p1.z - p2.z;

    return temp;
}

Vector vector_scalar_multiplication(Vector p, double factor)
{
    Vector temp;
    temp.x = factor * p.x;
    temp.y = factor * p.y;
    temp.z = factor * p.z;

    return temp;
}

double vector_dot_multiplication(Vector p1, Vector p2)
{

    double x_val = p1.x * p2.x;
    double y_val = p1.y * p2.y;
    double z_val = p1.z * p2.z;

    return x_val + y_val + z_val;
}

Vector vector_cross_multiplication(Vector p1, Vector p2)
{
    Vector temp;
    temp.x = p1.y * p2.z - p1.z * p2.y;
    temp.y = p1.z * p2.x - p1.x * p2.z;
    temp.z = p1.x * p2.y - p1.y * p2.x;

    return temp;
}


class Ray
{

public:
    Vector rO; 
    Vector rD; 
    Ray()
    {
       
    }

    Ray(Vector rO1, Vector rD1)
    {
        rO = rO1;
        rD = rD1;
        rD.normalize();
    }

    Vector getRO() const
    {
        return rO;
    }

    Vector getRD() const
    {
        return rD;
    }

    ~Ray()
    {
       
    }
};


struct Color
{
    double red;
    double green;
    double blue;

    Color()
    {
        red = green = blue = 0.0;
    }

    Color(double red1, double green1, double blue1)
    {
        red = red1;
        green = green1;
        blue = blue1;
    }

   

    ~Color()
    {
       
    }
};


class Light
{

public:
    Vector position;
    Color color;
   
    double radius;
    int segments;
    int stacks;

    Light()
    {
        radius = 0.0;
        segments = stacks = 0;
    }

    Light(Vector position1, Color color1, double radius1, int segments1, int stacks1)
    {
        position = position1;
        color = color1;
        radius = radius1;
        segments = segments1;
        stacks = stacks1;
    }

    Vector getPosition() const
    {
        return position;
    }

    Color getColor() const
    {
        return color;
    }

    void PointGeneration_And_Connecting()
    {
        Vector points[stacks + 1][segments + 1];
        double height, _radius;

        for (int i = 0; i <= stacks; i++)
        {
            double factor = ((double)i / (double)stacks);
            height = radius * sin(factor * (PI / 2));
            _radius = radius * cos(factor * (PI / 2));

            for (int j = 0; j <= segments; j++)
            {
                double fac = ((double)j / (double)segments);
                double x = _radius * cos(fac * 2 * PI);
                double y = _radius * sin(fac * 2 * PI);
                points[i][j] = Vector(x, y, height);
            }
        }

        glColor3f(color.red, color.green, color.blue);

        for (int i = 0; i < stacks; i++)
        {
            for (int j = 0; j < segments; j++)
            {
                glBegin(GL_QUADS);
                {

                    /* upper hemisphere */
                    glVertex3f(AdditionOfVector(position, points[i][j]).x, AdditionOfVector(position, points[i][j]).y, AdditionOfVector(position, points[i][j]).z);
                    glVertex3f(AdditionOfVector(position, points[i][j + 1]).x, AdditionOfVector(position, points[i][j + 1]).y, AdditionOfVector(position, points[i][j + 1]).z);
                    glVertex3f(AdditionOfVector(position, points[i + 1][j + 1]).x, AdditionOfVector(position, points[i + 1][j + 1]).y, AdditionOfVector(position, points[i + 1][j + 1]).z);
                    glVertex3f(AdditionOfVector(position, points[i + 1][j]).x, AdditionOfVector(position, points[i + 1][j]).y, AdditionOfVector(position, points[i + 1][j]).z);

                    /* lower hemisphere */

                    glVertex3f(AdditionOfVector(position, points[i][j]).x, AdditionOfVector(position, points[i][j]).y, SubtractionOfVector(position, points[i][j]).z);
                    glVertex3f(AdditionOfVector(position, points[i][j + 1]).x, AdditionOfVector(position, points[i][j + 1]).y, SubtractionOfVector(position, points[i][j + 1]).z);
                    glVertex3f(AdditionOfVector(position, points[i + 1][j + 1]).x, AdditionOfVector(position, points[i + 1][j + 1]).y, SubtractionOfVector(position, points[i + 1][j + 1]).z);
                    glVertex3f(AdditionOfVector(position, points[i + 1][j]).x, AdditionOfVector(position, points[i + 1][j]).y, SubtractionOfVector(position, points[i + 1][j]).z);
                }
                glEnd();
            }
        }
    }

    void draw()
    {
        PointGeneration_And_Connecting();
    }

    ~Light()
    {
        
    }
};

struct ReflectionCoefficient
{
    double ambient;
    double diffuse;
    double specular;
    double recursive;

    ReflectionCoefficient()
    {
        ambient= diffuse = specular = recursive = 0.0;
    }

    ReflectionCoefficient(double ambient1, double diffuse1, double specular1, double recursive1)
    {
        ambient = ambient1;
        diffuse = diffuse1;
        specular = specular1;
        recursive = recursive1;
    }

    ~ReflectionCoefficient()
    {
       
    }
};



class Object
{

public:
    Color color;
    ReflectionCoefficient reflectionCoefficient;
    int shininess;

    Object()
    {
        shininess = 0;
    }

    Color getColor() const
    {
        return color;
    }

    void setColor(Color color1)
    {
        color = color1;
    }

    ReflectionCoefficient getReflectionCoefficient() const
    {
        return reflectionCoefficient;
    }

    void setReflectionCoefficient(ReflectionCoefficient reflectionCoefficient1)
    {
        reflectionCoefficient = reflectionCoefficient1;
    }

    int getShininess() const
    {
        return shininess;
    }

    void setShininess(int shininess1)
    {
        shininess = shininess1;
    }

    void Ambient_LightComponent_Computation(Color &color, Color intersectionPointColor)
    {

        double factor = reflectionCoefficient.ambient;

        color.red = intersectionPointColor.red * factor ;
        color.green = intersectionPointColor.green * factor;
        color.blue = intersectionPointColor.blue * factor;
    }

    double GetDiffuseFactor(double lambertValue)
    {
        double r_d = reflectionCoefficient.diffuse;
        double positive_val = max(lambertValue, 0.0);
        return r_d * positive_val;
    }
    double GetSpecularFactor(double phongValue)
    {
        double r_s = reflectionCoefficient.specular;
        double positive_phongValue = max(phongValue, 0.0);
        double power_value = pow(positive_phongValue, shininess);
        return r_s * power_value;
    }

    void Reflection_Components_Computation(Ray ray, Color &color, Vector intersectionPoint, Color intersectionPointColor, Vector normal, Light light, Ray incidentRay)
    {
       
        double lambertValue = vector_dot_multiplication(vector_scalar_multiplication(incidentRay.getRD(), -1.0), normal);

        
        double result = vector_dot_multiplication(incidentRay.getRD(), normal) * 2.0;
        Vector v1 = vector_scalar_multiplication(normal, result);
        Vector v2 = SubtractionOfVector(incidentRay.getRD(), v1);
        Ray reflectedRay(intersectionPoint, v2);
       
        Vector v3 = vector_scalar_multiplication(ray.getRD(), -1.0);
        double phongValue = vector_dot_multiplication(v3, reflectedRay.getRD());

        double diffuse_factor = GetDiffuseFactor(lambertValue);

        color.red += light.color.red * intersectionPointColor.red * diffuse_factor;
        color.green += light.color.green * intersectionPointColor.green * diffuse_factor;
        color.blue += light.color.blue * intersectionPointColor.blue * diffuse_factor;

        double specular_factor = GetSpecularFactor(phongValue);

        color.red += light.color.red * intersectionPointColor.red * specular_factor;
        color.green += light.color.green * intersectionPointColor.green * specular_factor;
        color.blue += light.color.blue * intersectionPointColor.blue * specular_factor;

    }

    void Recursive_Reflection_Component_Computation(Color &color, Color reflectedColor)
    {
        double factor = reflectionCoefficient.recursive;
        color.red += reflectedColor.red * factor;
        color.green += reflectedColor.green * factor;
        color.blue += reflectedColor.blue * factor;
    }

    virtual void draw() = 0;
    virtual double intersect(Ray, Color &, int) = 0;

    ~Object()
    {
    }
};


Vector position;
int recursionLevel = 0;
vector<Object *> objects;
vector<Light> lights;

class Sphere : public Object
{
    Vector center;
    double radius;

    int segments;
    int stacks;

public:
    Sphere()
    {
        radius = 0.0;
        segments = stacks = 0;
    }

    Sphere(Vector center1, double radius1, int segments1, int stacks1)
    {
        center = center1;
        radius = radius1;

        segments = segments1;
        stacks = stacks1;
    }

    void Sphere_PointGeneration_And_Connecting()
    {

        Vector points[stacks + 1][segments + 1];
        double height, _radius;

        for (int i = 0; i <= stacks; i++)
        {
            double angle = ((double)i / (double)stacks);
            height = radius * sin(angle * (PI / 2));
            _radius = radius * cos(angle * (PI / 2));

            for (int j = 0; j <= segments; j++)
            {
                double theta = ((double)j / (double)segments);
                double x = _radius * cos(theta * 2 * PI);
                double y = _radius * sin(theta * 2 * PI);
                points[i][j] = Vector(x, y, height);
            }
        }


        glColor3f(getColor().red, getColor().green, getColor().blue);

        for (int i = 0; i < stacks; i++)
        {
            for (int j = 0; j < segments; j++)
            {
                glBegin(GL_QUADS);
                {

                    /* upper hemisphere */
                    glVertex3f(AdditionOfVector(center, points[i][j]).x, AdditionOfVector(center, points[i][j]).y, AdditionOfVector(center, points[i][j]).z);
                    glVertex3f(AdditionOfVector(center, points[i][j + 1]).x, AdditionOfVector(center, points[i][j + 1]).y, AdditionOfVector(center, points[i][j + 1]).z);
                    glVertex3f(AdditionOfVector(center, points[i + 1][j + 1]).x, AdditionOfVector(center, points[i + 1][j + 1]).y, AdditionOfVector(center, points[i + 1][j + 1]).z);
                    glVertex3f(AdditionOfVector(center, points[i + 1][j]).x, AdditionOfVector(center, points[i + 1][j]).y, AdditionOfVector(center, points[i + 1][j]).z);

                    /* lower hemisphere */

                    glVertex3f(AdditionOfVector(center, points[i][j]).x, AdditionOfVector(center, points[i][j]).y, SubtractionOfVector(center, points[i][j]).z);
                    glVertex3f(AdditionOfVector(center, points[i][j + 1]).x, AdditionOfVector(center, points[i][j + 1]).y, SubtractionOfVector(center, points[i][j + 1]).z);
                    glVertex3f(AdditionOfVector(center, points[i + 1][j + 1]).x, AdditionOfVector(center, points[i + 1][j + 1]).y, SubtractionOfVector(center, points[i + 1][j + 1]).z);
                    glVertex3f(AdditionOfVector(center, points[i + 1][j]).x, AdditionOfVector(center, points[i + 1][j]).y, SubtractionOfVector(center, points[i + 1][j]).z);
                }
                glEnd();
            }
        }
    }

    void draw()
    {
        Sphere_PointGeneration_And_Connecting();
    }

    double intersect(Ray ray, Color &color, int level)
    {
       
        double a, b, c, tMin;

       
        a = vector_dot_multiplication(ray.getRD(), ray.getRD());

        
        double m1 = vector_dot_multiplication(ray.getRO(), ray.getRD());
        double m2 = vector_dot_multiplication(ray.getRD(), center);
        b = (m1 - m2) * 2.0;

        double m3 = vector_dot_multiplication(ray.getRO(), ray.getRO());
        double m4 = vector_dot_multiplication(center, center);
        double m5 = vector_dot_multiplication(ray.getRO(), center);
        c = m3 + m4 - m5 * 2.0 - radius * radius;

        double discriminant = b * b - 4.0 * a * c;

        if (discriminant < 0.0)
        {
            tMin = INF;
        }
        else if (discriminant > 0.0)
        {

            double val = sqrt(discriminant) / (2.0 * a);
            double tMax = -b / (2.0 * a) + val;
            tMin = -b / (2.0 * a) - val;
            tMin = (tMin > 0.0) ? tMin : tMax;
        }
        else
        {
            tMin = -b / (2.0 * a);
        }

        if (level == 0)
        {
            return tMin;
        }

       
        Vector v1 = vector_scalar_multiplication(ray.getRD(), tMin);
        Vector intersectionPoint = AdditionOfVector(ray.getRO(), v1);

        Color intersectionPointColor = getColor();


        Vector normal = SubtractionOfVector(intersectionPoint, center);

        normal.normalize();

        
        Vector v2 = vector_scalar_multiplication(normal, -1.0);
        normal = (position.computeDistanceBetween(center) > radius) ? normal : v2;

        Ambient_LightComponent_Computation(color, intersectionPointColor);


        for (int i = 0; i < lights.size(); i++)
        {

            Ray incidentRay(lights[i].getPosition(), SubtractionOfVector(intersectionPoint, lights[i].getPosition()));

            
            double t, tMinimum = INF;

            for (int j = 0; j < objects.size(); j++)
            {
                Color dummyColor; 
                t = objects[j]->intersect(incidentRay, dummyColor, 0);


                if (t > 0.0)
                {
                    if (t < tMinimum)
                    {
                        tMinimum = t;
                    }
                }
            }

           
            Vector v6 = vector_scalar_multiplication(incidentRay.getRD(), tMinimum);
            Vector shadowIntersectionPoint = AdditionOfVector(incidentRay.getRO(), v6);

            double epsilon = 0.0000001; 

            if (intersectionPoint.computeDistanceBetween(incidentRay.getRO()) - epsilon > shadowIntersectionPoint.computeDistanceBetween(incidentRay.getRO()))
            {
                
                continue;
            }

           
            Reflection_Components_Computation(ray, color, intersectionPoint, intersectionPointColor, normal, lights[i], incidentRay);
        }

       
        if (level >= recursionLevel)
        {
            return tMin;
        }

      
        double m7 = vector_dot_multiplication(ray.getRD(), normal) * 2.0;
        Vector v7 = vector_scalar_multiplication(normal, m7);
        Vector reflectionDirection = SubtractionOfVector(ray.getRD(), v7);

        reflectionDirection.normalize();
        
        Ray reflectedRay(AdditionOfVector(intersectionPoint, reflectionDirection), reflectionDirection);

        int nearest = INT_MAX;
        double t, tMinimum = INF;

        for (int i = 0; i < objects.size(); i++)
        {
            Color dummyColor; 
            t = objects[i]->intersect(reflectedRay, dummyColor, 0);


            if (t > 0.0)
            {
                if (t < tMinimum)
                {
                    tMinimum = t;
                    nearest = i;
                }
            }
        }

      
        Color reflectedColor; 

        if (nearest != INT_MAX)
        {
            tMinimum = objects[nearest]->intersect(reflectedRay, reflectedColor, level + 1);
        }

       
        Recursive_Reflection_Component_Computation(color, reflectedColor);

      
        if (color.red > 1.0)
        {
            color.red = 1.0;
        }
        else if (color.red < 0.0)
        {
            color.red = 0.0;
        }
        else
        {
            color.red = color.red;
        }

        if (color.green > 1.0)
        {
            color.green = 1.0;
        }
        else if (color.green < 0.0)
        {
            color.green = 0.0;
        }
        else
        {
            color.green = color.green;
        }

        if (color.blue > 1.0)
        {
            color.blue = 1.0;
        }
        else if (color.blue < 0.0)
        {
            color.blue = 0.0;
        }
        else
        {
            color.blue = color.blue;
        }

        return tMin;
    }

    ~Sphere()
    {
       
    }
};


class Triangle : public Object
{
    Vector a, b, c;

public:
    Triangle()
    {
        
    }

    Triangle(Vector a1, Vector b1, Vector c1)
    {
        a = a1;
        b = b1;
        c = c1;
    }

    void draw()
    {
       
        glColor3f(color.red, color.green, color.blue);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(a.x, a.y, a.z);
            glVertex3f(b.x, b.y, b.z);
            glVertex3f(c.x, c.y, c.z);
        }
        glEnd();
    }
    double get_determinantBeta(Ray ray)
    {

        double term1 = (a.x - ray.getRO().x) * ((a.y - c.y) * ray.getRD().z - (a.z - c.z) * ray.getRD().y);
        double term2 = (a.x - c.x) * ((a.z - ray.getRO().z) * ray.getRD().y - (a.y - ray.getRO().y) * ray.getRD().z);
        double term3 =  ray.getRD().x * ((a.y - ray.getRO().y) * (a.z - c.z) - (a.z - ray.getRO().z) * (a.y - c.y));
     
        return  term1 + term2 + term3;
        

    }
    double get_determinantBase (Ray ray)
    {
        double term1 = (a.x - b.x) * ((a.y - c.y) * ray.getRD().z - (a.z - c.z) * ray.getRD().y);
        double term2 = (a.x - c.x) * ((a.z - b.z) * ray.getRD().y - (a.y - b.y) * ray.getRD().z);
        double term3 = ray.getRD().x * ((a.y - b.y) * (a.z - c.z) - (a.z - b.z) * (a.y - c.y));
     
        return  term1 + term2 + term3;
        

    }
    double get_determinantGamma(Ray ray)
    {
        double term1 = (a.x - b.x) * ((a.y - ray.getRO().y) * ray.getRD().z - (a.z - ray.getRO().z) * ray.getRD().y);
        double term2 = (a.x - ray.getRO().x) * ((a.z - b.z) * ray.getRD().y - (a.y - b.y) * ray.getRD().z);
        double term3 = ray.getRD().x * ((a.y - b.y) * (a.z - ray.getRO().z) - (a.z - b.z) * (a.y - ray.getRO().y));
     
        return  term1 + term2 + term3;
        

    }
    double get_determinantT(Ray ray)
    {
        double term1 = (a.x - b.x) * ((a.y - c.y) * (a.z - ray.getRO().z) - (a.z - c.z) * (a.y - ray.getRO().y));
        double term2 = (a.x - c.x) * ((a.z - b.z) * (a.y - ray.getRO().y) - (a.y - b.y) * (a.z - ray.getRO().z));
        double term3 = (a.x - ray.getRO().x) * ((a.y - b.y) * (a.z - c.z) - (a.z - b.z) * (a.y - c.y));
     
        return  term1 + term2 + term3;
    

    }


    double intersect(Ray ray, Color &color, int level)
    {
     
        double determinantBase, determinantBeta, determinantGamma, determinantT, tMin;

        determinantBase = get_determinantBase(ray);
        determinantBeta = get_determinantBeta(ray); 
        determinantGamma = get_determinantGamma(ray);
        determinantT = get_determinantT(ray);

        if (determinantBase == 0.0)
        {
         
            tMin = INF;
        }
        else
        {
           
            if (determinantBeta / determinantBase > 0.0 && determinantGamma / determinantBase > 0.0 && determinantBeta / determinantBase + determinantGamma / determinantBase < 1.0)
            {
                
                tMin = determinantT / determinantBase;
            }
            else
            {
             
                tMin = INF;
            }
        }

        if (level == 0)
        {
            return tMin;
        }

      
        Vector intersectionPoint = AdditionOfVector(ray.getRO(), vector_scalar_multiplication(ray.getRD(), tMin));

        Color intersectionPointColor = getColor();

      
        Vector normal = vector_cross_multiplication(SubtractionOfVector(b, a), SubtractionOfVector(c, a));

        normal.normalize();
    
        normal = (vector_dot_multiplication(vector_scalar_multiplication(ray.getRD(), -1.0), normal) > 0.0) ? normal : vector_scalar_multiplication(normal, -1.0);

       
        Ambient_LightComponent_Computation(color, intersectionPointColor);

       
        for (int i = 0; i < lights.size(); i++)
        {

            Ray incidentRay(lights[i].getPosition(), SubtractionOfVector(intersectionPoint, lights[i].getPosition()));

            double t, tMinimum = INF;

            for (int j = 0; j < objects.size(); j++)
            {
                Color dummyColor; 
                t = objects[j]->intersect(incidentRay, dummyColor, 0);

                if (t > 0.0)
                {
                    if (t < tMinimum)
                    {
                        tMinimum = t;
                    }
                }
            }

          
            Vector shadowIntersectionPoint = AdditionOfVector(incidentRay.getRO(), vector_scalar_multiplication(incidentRay.getRD(), tMinimum));

            double epsilon = 0.0000001; 

            if (intersectionPoint.computeDistanceBetween(incidentRay.getRO()) - epsilon > shadowIntersectionPoint.computeDistanceBetween(incidentRay.getRO()))
            {
                
                continue;
            }

           
            Reflection_Components_Computation(ray, color, intersectionPoint, intersectionPointColor, normal, lights[i], incidentRay);
        }

  
        if (level >= recursionLevel)
        {
            return tMin;
        }

        double m1 = vector_dot_multiplication(ray.getRD(), normal) * 2.0;
        Vector vec1 = vector_scalar_multiplication(normal, m1);
        Vector reflectionDirection = SubtractionOfVector(ray.getRD(), vec1);

        reflectionDirection.normalize();
       
        Ray reflectedRay(AdditionOfVector(intersectionPoint, reflectionDirection), reflectionDirection);

        int nearest = INT_MAX;
        double t, tMinimum = INF;

        for (int i = 0; i < objects.size(); i++)
        {
            Color dummyColor; 
            t = objects[i]->intersect(reflectedRay, dummyColor, 0);

            if (t > 0.0)
            {

                if (t < tMinimum)
                {
                    tMinimum = t;
                    nearest = i;
                }
            }
        }

      
        Color reflectedColor; 

        if (nearest != INT_MAX)
        {
            tMinimum = objects[nearest]->intersect(reflectedRay, reflectedColor, level + 1);
        }

       
        Recursive_Reflection_Component_Computation(color, reflectedColor);

      

        if (color.red > 1.0)
        {
            color.red = 1.0;
        }
        else if (color.red < 0.0)
        {
            color.red = 0.0;
        }
        else
        {
            color.red = color.red;
        }

        if (color.green > 1.0)
        {
            color.green = 1.0;
        }
        else if (color.green < 0.0)
        {
            color.green = 0.0;
        }
        else
        {
            color.green = color.green;
        }

        if (color.blue > 1.0)
        {
            color.blue = 1.0;
        }
        else if (color.blue < 0.0)
        {
            color.blue = 0.0;
        }
        else
        {
            color.blue = color.blue;
        }

        return tMin;
    }

    ~Triangle()
    {
       
    }
};





struct GeneralQuadricSurfaceCoefficient
{
    double a, b, c, d, e, f, g, h, i, j;

  
};

class GeneralQuadricSurface : public Object
{
   

public:
   GeneralQuadricSurfaceCoefficient coefficient;
    Vector cubeReferencePoint;
    double length;
    double width;
    double height;
    GeneralQuadricSurface()
    {
        coefficient.a = 0.0;
        coefficient.b = 0.0;
        coefficient.c = 0.0;
        coefficient.d = 0.0;
        coefficient.e = 0.0;
        coefficient.f = 0.0;
        coefficient.g = 0.0;
        coefficient.h = 0.0;
        coefficient.i = 0.0;
        coefficient.j = 0.0;
        length =0.0;
        width =0.0;
        height = 0.0;
    }

    GeneralQuadricSurface(GeneralQuadricSurfaceCoefficient coefficient1, Vector cubeReferencePoint1, double length1, double width1, double height1)
    {
        coefficient = coefficient1;
        cubeReferencePoint = cubeReferencePoint1;
        length = length1;
        width = width1;
        height = height1;
    }

    void draw()
    {
        
    }

    double get_a(Ray ray)
    {
        double a;
        double term1 = coefficient.a * ray.getRD().x * ray.getRD().x;
        double term2 = coefficient.b * ray.getRD().y * ray.getRD().y;
        double term3 = coefficient.c * ray.getRD().z * ray.getRD().z;
        double term4 = coefficient.d * ray.getRD().x * ray.getRD().y;
        double term5 = coefficient.e * ray.getRD().x * ray.getRD().z;
        double term6 = coefficient.f * ray.getRD().y * ray.getRD().z;
        a = term1 + term2 + term3 + term4 + term5 + term6;
        return a;
    }
    double get_b(Ray ray)
    {
        double b;
        double term1 = 2.0 * coefficient.a * ray.getRO().x * ray.getRD().x;
        double term2 = 2.0 * coefficient.b * ray.getRO().y * ray.getRD().y;
        double term3 = 2.0 * coefficient.c * ray.getRO().z * ray.getRD().z;
        double term4 = coefficient.d * (ray.getRO().x * ray.getRD().y + ray.getRD().x * ray.getRO().y);
        double term5 = coefficient.e * (ray.getRO().x * ray.getRD().z + ray.getRD().x * ray.getRO().z);
        double term6 = coefficient.f * (ray.getRO().y * ray.getRD().z + ray.getRD().y * ray.getRO().z);
        double term7 = coefficient.g * ray.getRD().x + coefficient.h * ray.getRD().y + coefficient.i * ray.getRD().z;

        b = term1 + term2 + term3 + term4 + term5 + term6 + term7;

        return b;
    }
    double get_c(Ray ray)
    {
        double c;

        double term1 = coefficient.a * ray.getRO().x * ray.getRO().x + coefficient.b * ray.getRO().y * ray.getRO().y + coefficient.c * ray.getRO().z * ray.getRO().z;
        double term2 = coefficient.d * ray.getRO().x * ray.getRO().y + coefficient.e * ray.getRO().x * ray.getRO().z + coefficient.f * ray.getRO().y * ray.getRO().z;
        double term3 = coefficient.g * ray.getRO().x + coefficient.h * ray.getRO().y + coefficient.i * ray.getRO().z + coefficient.j;

        c = term1 + term2 + term3;

        return c;
    }

    double intersect(Ray ray, Color &color, int level)
{
   
    double a, b, c, tMin, tMax;

    a = get_a(ray);
    b = get_b(ray);
    c = get_c(ray);

   
    if (a == 0.0)
    {
        tMin = (b == 0.0) ? INF : -c / b;
        tMax = INF;
    }
    else
    {
        double discriminant = b * b - 4.0 * a * c;

        if (discriminant < 0.0)
        {
            tMin = INF;
            tMax = INF;
        }
        else if (discriminant > 0.0)
        {
            double factor = sqrt(discriminant) / (2.0 * a);
            tMax = -b / (2.0 * a) + factor;
            tMin = -b / (2.0 * a) - factor;
        }
        else
        {
            tMax = INF;
            tMin = -b / (2.0 * a);
          
        }
    }

    if (tMin < INF)
    {
        if (tMax < INF)
        {
            if (tMin > 0.0)
            {
               
                Vector intersectionPoint = AdditionOfVector(ray.getRO(), vector_scalar_multiplication(ray.getRD(), tMin));

                if ((length != 0.0 && (intersectionPoint.x < cubeReferencePoint.x || intersectionPoint.x > cubeReferencePoint.x + length)) || (width != 0.0 && (intersectionPoint.y < cubeReferencePoint.y || intersectionPoint.y > cubeReferencePoint.y + width)) || (height != 0.0 && (intersectionPoint.z < cubeReferencePoint.z || intersectionPoint.z > cubeReferencePoint.z + height)))
                {
                    tMin = INF;
                }
            }
            if (tMax > 0.0)
            {
               
                Vector intersectionPoint = AdditionOfVector(ray.getRO(), vector_scalar_multiplication(ray.getRD(), tMax));

                if ((length != 0.0 && (intersectionPoint.x < cubeReferencePoint.x || intersectionPoint.x > cubeReferencePoint.x + length)) || (width != 0.0 && (intersectionPoint.y < cubeReferencePoint.y || intersectionPoint.y > cubeReferencePoint.y + width)) || (height != 0.0 && (intersectionPoint.z < cubeReferencePoint.z || intersectionPoint.z > cubeReferencePoint.z + height)))
                {
                    tMax = INF;
                }
            }
            tMin = (tMin > 0.0 && tMin < tMax) ? tMin : tMax;
        }
        else
        {
            if (tMin > 0.0)
            {
                
                Vector intersectionPoint = AdditionOfVector(ray.getRO(), vector_scalar_multiplication(ray.getRD(), tMin));

                if ((length != 0.0 && (intersectionPoint.x < cubeReferencePoint.x || intersectionPoint.x > cubeReferencePoint.x + length)) || (width != 0.0 && (intersectionPoint.y < cubeReferencePoint.y || intersectionPoint.y > cubeReferencePoint.y + width)) || (height != 0.0 && (intersectionPoint.z < cubeReferencePoint.z || intersectionPoint.z > cubeReferencePoint.z + height)))
                {
                    tMin = INF;
                }
            }
        }
    }

    if (level == 0)
    {
        return tMin;
    }

    Vector intersectionPoint = AdditionOfVector(ray.getRO(), vector_scalar_multiplication(ray.getRD(), tMin));
    Color intersectionPointColor = getColor();

   
    double xNormal, yNormal, zNormal;


    double term1 = 2.0 * coefficient.a * intersectionPoint.x;
    double term2 = coefficient.d * intersectionPoint.y;
    double term3 = coefficient.e * intersectionPoint.z + coefficient.g;
    xNormal = term1 + term2 + term3 ;
   
    
    double term4 = 2.0 * coefficient.b * intersectionPoint.y;
    double term5 = coefficient.d * intersectionPoint.x;
    double term6 = coefficient.f * intersectionPoint.z + coefficient.h;
    yNormal = term4 + term5 + term6;
   
    double term7 =  2.0 * coefficient.c * intersectionPoint.z;
    double term8 = coefficient.e * intersectionPoint.x;
    double term9 =coefficient.f * intersectionPoint.y + coefficient.i;
    zNormal = term7 + term8 + term9 ;
   

    Vector normal(xNormal, yNormal, zNormal);
    normal.normalize();
    normal = (vector_dot_multiplication(vector_scalar_multiplication(ray.getRD(), -1.0), normal) > 0.0) ? normal : vector_scalar_multiplication(normal, -1.0);

   
    Ambient_LightComponent_Computation(color, intersectionPointColor);

   
    for (int i = 0; i < lights.size(); i++)
    {
        Ray incidentRay(lights[i].getPosition(), SubtractionOfVector(intersectionPoint, lights[i].getPosition()));

        double t, tMinimum = INF;

        for (int j = 0; j < objects.size(); j++)
        {
            Color dummyColor; 
            t = objects[j]->intersect(incidentRay, dummyColor, 0);

            if (t > 0.0 )
            {
                 if (t < tMinimum)
            {
                tMinimum = t;
            }
            }
        }

        
        Vector shadowIntersectionPoint = AdditionOfVector(incidentRay.getRO(), vector_scalar_multiplication(incidentRay.getRD(), tMinimum));
        double epsilon = 0.0000001; 

        if (intersectionPoint.computeDistanceBetween(incidentRay.getRO()) - epsilon > shadowIntersectionPoint.computeDistanceBetween(incidentRay.getRO()))
        {
            
            continue;
        }

       
        Reflection_Components_Computation(ray, color, intersectionPoint, intersectionPointColor, normal, lights[i], incidentRay);
    }

    
    if (level >= recursionLevel)
    {
        return tMin;
    }

   
    double m1 = vector_dot_multiplication(ray.getRD(), normal) * 2.0;
    Vector vec1 = vector_scalar_multiplication(normal, m1);
    Vector reflectionDirection = SubtractionOfVector(ray.getRD(), vec1);

    reflectionDirection.normalize();
   
    Ray reflectedRay(AdditionOfVector(intersectionPoint, reflectionDirection), reflectionDirection);
    int nearest = INT_MAX;
    double t, tMinimum = INF;

    for (int i = 0; i < objects.size(); i++)
    {
        Color dummyColor; 
        t = objects[i]->intersect(reflectedRay, dummyColor, 0);

        if (t > 0.0)
        {
            if (t < tMinimum)
            {
                tMinimum = t;
                nearest = i;
            }
        }
    }

    
    Color reflectedColor;

    if (nearest != INT_MAX)
    {
        tMinimum = objects[nearest]->intersect(reflectedRay, reflectedColor, level + 1);
    }

   
    Recursive_Reflection_Component_Computation(color, reflectedColor);


     if (color.red > 1.0)
        {
            color.red = 1.0;
        }
        else if (color.red < 0.0)
        {
            color.red = 0.0;
        }
        else
        {
            color.red = color.red;
        }

        if (color.green > 1.0)
        {
            color.green = 1.0;
        }
        else if (color.green < 0.0)
        {
            color.green = 0.0;
        }
        else
        {
            color.green = color.green;
        }

        if (color.blue > 1.0)
        {
            color.blue = 1.0;
        }
        else if (color.blue < 0.0)
        {
            color.blue = 0.0;
        }
        else
        {
            color.blue = color.blue;
        }

    return tMin;
}


    ~GeneralQuadricSurface()
    {
        
    }
};


/* Floor class */
class Floor : public Object
{
   

public:
    double floorWidth;
    double tileWidth;

    Color foregroundColor;

    Floor()
    {
        floorWidth = tileWidth = 0.0;
    }

    Floor(double floorWidth1, double tileWidth1, Color foregroundColor1)
    {
        floorWidth = floorWidth1;
        tileWidth = tileWidth1;
        foregroundColor = foregroundColor1;
    }

    void PointConnect(Vector leftBottomCorner)
    {
                glBegin(GL_QUADS);
                {
                    glVertex3f(leftBottomCorner.x, leftBottomCorner.y, leftBottomCorner.z);
                    glVertex3f(leftBottomCorner.x + tileWidth, leftBottomCorner.y, leftBottomCorner.z);
                    glVertex3f(leftBottomCorner.x + tileWidth, leftBottomCorner.y + tileWidth, leftBottomCorner.z);
                    glVertex3f(leftBottomCorner.x, leftBottomCorner.y + tileWidth, leftBottomCorner.z);
                }
                glEnd();

    }

    void draw()
    {
        for (int i = 0, row = (int)floorWidth / tileWidth, column = (int)floorWidth / tileWidth; i < row; i++)
        {
            for (int j = 0; j < column; j++)
            {
              
                glColor3f(((i + j) % 2 == 0) ? getColor().red : foregroundColor.red, ((i + j) % 2 == 0) ? getColor().green : foregroundColor.green, ((i + j) % 2 == 0) ? getColor().blue : foregroundColor.blue);
                Vector leftBottomCorner(-floorWidth / 2.0 + tileWidth * j, -floorWidth / 2.0 + tileWidth * i, 0.0);

                PointConnect(leftBottomCorner);

            }
        }
    }

  
double intersect(Ray ray, Color &color, int level)
{

    Vector normal(0.0, 0.0, 1.0);
    normal = (vector_dot_multiplication(position, normal) > 0.0) ? normal : vector_scalar_multiplication(normal, -1.0);

  
    double tMin = INF;

    if (vector_dot_multiplication(normal, ray.getRD()) != 0.0)
    {

       

        double m1 = vector_dot_multiplication(normal, ray.getRO());
        double m2 = vector_dot_multiplication(normal, ray.getRD());
        tMin = (-1.0) * m1 / m2;
    }

    if (tMin > 0.0 && tMin < INF)
    {
      
        Vector intersectionPoint = AdditionOfVector(ray.getRO(), vector_scalar_multiplication(ray.getRD(), tMin));

        if (!((intersectionPoint.x > -floorWidth / 2.0 && intersectionPoint.x < floorWidth / 2.0) && (intersectionPoint.y > -floorWidth / 2.0 && intersectionPoint.y < floorWidth / 2.0)))
        {
          
            tMin = INF;
        }
    }

    if (level == 0)
    {
        return tMin;
    }

    Vector intersectionPoint = AdditionOfVector(ray.getRO(), vector_scalar_multiplication(ray.getRD(), tMin));

   
    Vector referencePosition = SubtractionOfVector(intersectionPoint, Vector(-floorWidth / 2.0, -floorWidth / 2.0, 0.0));

    Color intersectionPointColor = (((int)(floor(referencePosition.x / tileWidth) + floor(referencePosition.y / tileWidth))) % 2 == 0) ? getColor() : foregroundColor;

  
    Ambient_LightComponent_Computation(color, intersectionPointColor);

   
    for (int i = 0; i < lights.size(); i++)
    {
        Ray incidentRay(lights[i].getPosition(), SubtractionOfVector(intersectionPoint, lights[i].getPosition()));

       
        double t, tMinimum = INF;

        for (int j = 0; j < objects.size(); j++)
        {
            Color dummyColor; 
            t = objects[j]->intersect(incidentRay, dummyColor, 0);

            if (t > 0.0)
            {
                if (t < tMinimum)
                {
                    tMinimum = t;
                }
            }
        }

       

        Vector shadowIntersectionPoint = AdditionOfVector(incidentRay.getRO(), vector_scalar_multiplication(incidentRay.getRD(), tMinimum));

        double epsilon = 0.0000001; 
        if (intersectionPoint.computeDistanceBetween(incidentRay.getRO()) - epsilon > shadowIntersectionPoint.computeDistanceBetween(incidentRay.getRO()))
        {
          
            continue;
        }

        
        Reflection_Components_Computation(ray, color, intersectionPoint, intersectionPointColor, normal, lights[i], incidentRay);
    }

    
    if (level >= recursionLevel)
    {
        return tMin;
    }

    
    double m1 = vector_dot_multiplication(ray.getRD(), normal) * 2.0;
    Vector v1 = vector_scalar_multiplication(normal, m1);
    Vector reflectionDirection = SubtractionOfVector(ray.getRD(), v1);

    reflectionDirection.normalize();
    Ray reflectedRay(AdditionOfVector(intersectionPoint, reflectionDirection), reflectionDirection);

    
    int nearest = INT_MAX;
    double t, tMinimum = INF;

    for (int i = 0; i < objects.size(); i++)
    {
        Color dummyColor; 
        t = objects[i]->intersect(reflectedRay, dummyColor, 0);

        if (t > 0.0)
        {
            if (t < tMinimum)
            {
                tMinimum = t;
                nearest = i;
            }
        }
    }

    
    Color reflectedColor; 

    if (nearest != INT_MAX)
    {
        tMinimum = objects[nearest]->intersect(reflectedRay, reflectedColor, level + 1);
    }

    
    Recursive_Reflection_Component_Computation(color, reflectedColor);



     if (color.red > 1.0)
        {
            color.red = 1.0;
        }
        else if (color.red < 0.0)
        {
            color.red = 0.0;
        }
        else
        {
            color.red = color.red;
        }

        if (color.green > 1.0)
        {
            color.green = 1.0;
        }
        else if (color.green < 0.0)
        {
            color.green = 0.0;
        }
        else
        {
            color.green = color.green;
        }

        if (color.blue > 1.0)
        {
            color.blue = 1.0;
        }
        else if (color.blue < 0.0)
        {
            color.blue = 0.0;
        }
        else
        {
            color.blue = color.blue;
        }

    return tMin;
}


    ~Floor()
    {
       
    }
};


