#include <list>
using namespace std;

class Color;
class Light;
class Figure;
class Plane;
class Sphere;

class Vec
{
    //friend Vec operator*(double num, const Vec& v);

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
    Vec operator*(double t) const;
    Vec operator*(Vec abc) const;
    const double& operator[](int i) const;
    double& operator[](int i);
    Vec();
    Vec(ifstream& ifs);
    Vec(double xx, double yy, double zz);
    Vec(const Vec& other);
 //   explicit Vec(const Vect& v);
};

class Ray
{
    private:
        Vec v0;
        Vec v1;
    public:
        Ray(Vec v1);
        Ray(Vec v0, Vec v1);
        Vec getOrigin() const;
        Vec getDirection() const;
        Vec rayPoint(double t) const;
};


class Color
{
  //friend Color operator*(double num, const Color& c);

  protected:
   double red, green, blue;

  public:
	Color();
    Color(ifstream& ifs);
    Color(double r, double g, double b);
    Color operator+(const Color& c) const;
    Color operator*(double c) const;
    Color operator*(const Color& c) const;
    void writeOut(ofstream& ofs);

    void printOut();

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
   Color getColor(double c);
   virtual double intersection(const Ray& r, double minT, double maxT) const = 0;
   virtual Vec getNormal(Vec direction) const = 0;
   Color getColorAmbient();
   Color getColorDiffuse();
   Color getColorSpecular();
   Color getColorReflectivity();
   Color getColorTransmissivity();
   double getShininess();
   double getIndexOfRefraction();
   int getRFlag();
   int getTFlag();
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
    virtual double intersection(const Ray& r, double minT, double maxT) const;
    virtual Vec getNormal(Vec direction) const;
    Vec getABC();
    double getD();
    Vec getD1();
    Vec getD2();
    
};

class Sphere : public Figure
{
  private:
    Vec center;
    double radius;
  public:
    Sphere(ifstream& ifs);
    virtual double intersection(const Ray& r, double minT, double maxT) const;
    virtual Vec getNormal(Vec direction) const;
};

