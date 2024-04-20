#include <list>
using namespace std;

class Color;
class Light;
class Figure;
class Plane;
class Sphere;

class Vec
{
    friend Vec operator*(double num, const Vec& v);

private:
    double x;
    double y;
    double z;

public:
    void normalize();
    double norm() const;
    double dot(const Vec& other) const;
    Vec operator+(const Vec& other) const;
    Vec operator-(const Vec& other) const;
    const double& operator[](int i) const;
    double& operator[](int i);
    Vec();
    Vec(ifstream& ifs);
    Vec(double xx, double yy, double zz);
    Vec(const Vec& other);
 //   explicit Vec(const Vect& v);
};


class Color
{
  friend Color operator*(double num, const Color& c);

  protected:
   double red, green, blue;

  public:
	Color();
    Color(ifstream& ifs);
    Color(double r, double g, double b);

};

class Light
{
  public:
    Light(ifstream& ifs);
  private:
    Vec position;
    Color shading;
    double c0, c1, c2;
};


class Figure
{
 protected:
   Color ambient;
   Color diffuse;
   Color specular;
   Color reflectivity;
   Color transmissivity;
   double shininess;
   double indexOfRefraction;
   int rFlag, tFlag;

 public:
   Figure();
   void initFigure(ifstream& ifs);
};



class Plane : public Figure
{
  private:
    Vec abcVector;
    double dScalar;
    Vec direction1;
    Vec direction2;
  public:
    Plane(ifstream& ifs);
};

class Sphere : public Figure
{
  private:
    Vec center;
    double radius;
  public:
    Sphere(ifstream& ifs);
};

