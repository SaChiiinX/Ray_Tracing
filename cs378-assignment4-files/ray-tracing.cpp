#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <fstream>
#include <cassert>
#include <list>
#include <iostream>

using namespace std;

#include "ray-tracing.h"

list<Figure*> shapeList;
list<Light*> lightList;



double cameraX, cameraY, cameraZ;
int horizontalResolution, verticalResolution;
double zCoor;
double minX, maxX;
double minY, maxY;

int maxDepth;
Color backgroundColor;
Color ambient;

Color** pixelColors;

void Vec::normalize()
{
    double distance = norm();
    x = x / distance;
    y = y / distance;
    z = z / distance;
}

double Vec::norm() const
{
    return sqrt(x * x + y * y + z * z);
}

double Vec::dot(const Vec& other) const
{
    double sumX = x * (other.x);
    double sumY = y * (other.y);
    double sumZ = z * (other.z);
    return sumX + sumY + sumZ;
}

Vec Vec::operator+(const Vec& other) const
{
    double x2 = x + other.x;
    double y2 = y + other.y;
    double z2 = z + other.z;
    return Vec(x2, y2, z2);
}

Vec Vec::operator-(const Vec& other) const
{
    double x3 = x - other.x;
    double y3 = y - other.y;
    double z3 = z - other.z;
    return Vec(x3, y3, z3);
}

Vec::Vec() : x(0.0), y(0.0), z(0.0) {}

Vec::Vec(ifstream& ifs)
{
    ifs >> x >> y >> z;
}

Vec::Vec(double xx, double yy, double zz) : x(xx), y(yy), z(zz) {}

Vec::Vec(const Vec& other) : x(other.x), y(other.y), z(other.z) {}

// Vec::Vec(const Vect& v) : x(v[0]), y(v[1]), z(v[2]) {}


Vec operator*(double num, const Vec& v)
{
    Vec newVec(num * v.x, num * v.y, num * v.z);
    return newVec;
}

const double& Vec::operator[](int i) const
{
    switch (i)
    {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    default: return x;
    }
}


double& Vec::operator[](int i)
{
    switch (i)
    {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    default: return x;
    }
}


Color::Color(): red(0.0), green(0.0), blue(0.0) {};

Color::Color(ifstream& ifs)
{
  ifs >> red >> green >> blue;
}

Color::Color(double r, double g, double b) : red(r), green(g), blue(b)
{}

Color Color::add(Color color)
{
    return Color(this->red + color.red, this->blue + color.blue, this->green + color.green);
}

Color Color::scale(double c)
{
    return Color(this->red * c, this->blue * c, this->green * c);
}

void Color::writeOut(ofstream& ofs)
{
    ofs << this->red << "\t";
    ofs << this->blue << "\t";
    ofs << this->green << "\t";
}

Figure::Figure(){}

void Figure::initFigure(ifstream& ifs)
{
 ambient = Color(ifs);
 diffuse = Color(ifs);
 specular = Color(ifs);
 reflectivity = Color(ifs);
 transmissivity = Color(ifs);
 ifs >> shininess >> indexOfRefraction >> rFlag >> tFlag;
}

Color Figure::getColor(double c)
{
    Color combined = diffuse.add(specular);
    combined = combined.scale(c);
    combined = ambient.add(combined);
    return combined;
}

Light::Light(ifstream& ifs) : position(ifs), shading(ifs)
{
  ifs >> c0 >> c1 >> c2;
}


Sphere::Sphere(ifstream& ifs) : center(ifs)
{
  ifs >> radius;
  initFigure(ifs);
}

double Sphere::intersection(const Ray& r, double minT, double maxT) const
{
}


Plane::Plane(ifstream& ifs) : abcVector(ifs)
{
  ifs >> dScalar;
  initFigure(ifs);
  direction1 = Vec(ifs);
  direction2 = Vec(ifs);
}

double Plane::intersection(const Ray& r, double minT, double maxT) const
{
}


void parseSceneFile(char* sceneName)
{
  double bgr, bgg, bgb;
  double ar, ag, ab;
  ifstream ifs;
  assert (sceneName != 0);
  ifs.open(sceneName);
  ifs >> cameraX;
  ifs >> cameraY;
  ifs >> cameraZ;
  ifs >> zCoor;
  ifs >> minX >> maxX;
  ifs >> minY >> maxY;
  ifs >> bgr >> bgg >> bgb;
  backgroundColor = Color(bgr,bgg,bgb);
  ifs >> ar >> ag >> ab;
  ambient = Color(ar,ag,ab);
  ifs >> maxDepth;
  ifs >> horizontalResolution >> verticalResolution;
  int numLights, numSpheres, numPlanes;
  ifs >> numLights;
  for (int i=0; i<numLights; ++i) lightList.push_front(new Light(ifs));
  ifs >> numSpheres;
  for (int i=0; i<numSpheres; ++i) shapeList.push_front(new Sphere(ifs));
  ifs >> numPlanes;
  for (int i=0; i<numPlanes; ++i) shapeList.push_front(new Plane(ifs));
  ifs.close();
}

Ray::Ray(Vec* v1)
{
    this->v0 = new Vec(0, 0, zCoor);
    // This may be a problem later
    this->v1 = v1;
}

Ray::Ray(Vec* v0, Vec* v1)
{
    // May run into pointer issues later
    this->v0 = v0;
    this->v1 = v1;
}




pair<double, Figure*> nearestIntersection(const Ray& r, double minT, double maxT, bool mayBeTransparent = true)
{
    for (Figure* figure : shapeList)
    {
        figure->intersection(r, minT, maxT);
        double c = 0.0;

    }
}


Color RT_trace(const Ray& r, double depth)
{
    double epsilon = 0.0;
    double maxT = 0.0;
    pair<double, Figure*> intersection =  nearestIntersection(r, epsilon, maxT);
    ... RT_shade(figure, r, i, normal, entering, depth);
}


void RT_algorithm() 
{
    for (int u = 0; u < horizontalResolution; u++) {
        for (int v = 0; v < verticalResolution; v++) {

            double wp = (maxX - minX) / u;
            double hp = (maxY - minY) / v;


            double x = minX + wp/2 + u*wp;
            double y = minY + hp/2 + v*hp;

            Vec* point = new Vec(x,y,zCoor);
            Ray* ray = new Ray(point);
            
            pixelColors[u][v] = RT_trace(*ray, 1);
        }
    }
}

void initializeImage()
{
    pixelColors = new Color * [horizontalResolution];

    for (int u = 0; u < horizontalResolution; u++) {
        pixelColors[u] = new Color[verticalResolution];
    }
}

void writeImageFile()
{
    ofstream ofs;
    ofs.open("test.ppm");
    ofs << "P3\n";
    ofs << horizontalResolution << "\t" << verticalResolution << "\n";

    // Color of each pixel
    for (int u = 0; u < horizontalResolution; u++) {
        for (int v = 0; v < verticalResolution; v++) {
            if (v % 4 == 0) {
                ofs << "\n";
            }
            pixelColors[u][v].writeOut(ofs);
        }
    }
    ofs.close();
}


int main(int, char *argv[])
{
    parseSceneFile(argv[1]);
    initializeImage();
    RT_algorithm();
    writeImageFile();
    return 0;
}


