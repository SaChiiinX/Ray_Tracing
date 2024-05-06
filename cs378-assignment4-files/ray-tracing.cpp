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

Vec Vec::operator*(double num) const
{
    Vec newVec(num * this->x, num * this->y, num * this->z);
    return newVec;
}

Vec Vec::operator*(Vec abc) const
{
    return Vec();
}

Vec::Vec() : x(0.0), y(0.0), z(0.0) {}

Vec::Vec(ifstream& ifs)
{
    ifs >> x >> y >> z;
}

Vec::Vec(double xx, double yy, double zz) : x(xx), y(yy), z(zz) {}

Vec::Vec(const Vec& other) : x(other.x), y(other.y), z(other.z) {}

// Vec::Vec(const Vect& v) : x(v[0]), y(v[1]), z(v[2]) {}


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

Color Color::operator+(const Color & c) const
{
    return Color(this->red + c.red, this->green + c.green, this->blue + c.blue);
}

Color Color::operator*(double c) const
{
    return Color(this->red * c, this->green * c, this->blue * c);
}

Color Color::operator*(const Color& c) const
{
    return Color(this->red * c.red, this->green * c.green, this->blue * c.blue);
}

void Color::writeOut(ofstream& ofs)
{

    ofs << (int)this->red << " ";
    ofs << (int)this->green << " ";
    ofs << (int)this->blue << " ";
}

void Color::printOut() {
    cout << "r: " << this->red << " g: " << this->green << " b: " << this->blue << "\n";
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
    Color combined = diffuse+specular;
    combined = combined*c;
    combined = ambient+combined;
    return combined;
}

Color Figure::getColorAmbient()
{
    return this->ambient;
}

Color Figure::getColorDiffuse()
{
    return this->diffuse;
}

Color Figure::getColorSpecular()
{
    return this->specular;
}

Color Figure::getColorReflectivity()
{
    return this->reflectivity;
}

Color Figure::getColorTransmissivity()
{
    return this->transmissivity;
}


double Figure::getShininess()
{
    return this->shininess;
}


double Figure::getIndexOfRefraction()
{
    return this->indexOfRefraction;
}

int Figure::getRFlag()
{
    return this->rFlag;
}

int Figure::getTFlag()
{
    return this->tFlag;
}

Light::Light(ifstream& ifs) : position(ifs), shading(ifs)
{
  ifs >> c0 >> c1 >> c2;
}

Vec Light::getPosition()
{
    return this->position;
}

Color Light::getShading()
{
    return this->shading;
}

double Light::getFatt(double d)
{
    return 1/(c0 + c1*d+c2*d*d);
}


Sphere::Sphere(ifstream& ifs) : center(ifs)
{
  ifs >> radius;
  initFigure(ifs);
}

double Sphere::intersection(const Ray& r, double minT, double maxT) const
{
    Vec p0 = r.getFirst();
    Vec p1 = r.getSecond();
    Vec h = p1-p0;
    Vec k = p0-center;
    double a = h.dot(h);
    double b = 2.0 * (h.dot(k));
    double c = k.dot(k) - (radius * radius);


    double underRoot = b * b - (4 * a * c);

    if (underRoot < 0 || a == 0.0) {
        return -1;
    }

    double sqrtroot = sqrt(underRoot);
    double denominator = 2.0 * a;
    double posQuad = (double)((-b + sqrtroot) / denominator);
    double negQuad = (double)((-b - sqrtroot) / denominator);

    if ((posQuad >= 1) && (negQuad >= 1)) {
        return min(posQuad, negQuad);
    }
    else if (posQuad >= 1) {
        return posQuad;
    }
    else if (negQuad >= 1) {
        return negQuad;
    }
    else {
        return minT;
    }
}

Vec Sphere::getNormal(Vec direction) const
{
    Vec normal = (direction - center) * (1.0 / radius);
    normal.normalize();
    return normal;
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
    Vec p0 = r.getFirst();
    Vec p1 = r.getSecond();
    double denominator = (p1 - (p0)).dot(abcVector);

    if (denominator == 0.0) {
        return minT;
    }
    else {
        double numerator = dScalar - (p1.dot(abcVector));
        return numerator / denominator;
    }
}

Vec Plane::getNormal(Vec direction) const
{
    double magn = abcVector.norm();
    Vec normal = abcVector * (1.0 / magn);
    return normal;
}

Vec Plane::getABC()
{
    return this->abcVector;
}

double Plane::getD()
{
    return this->dScalar;
}

Vec Plane::getD1()
{
    return this->direction1;
}

Vec Plane::getD2()
{
    return this->direction2;
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

Ray::Ray(Vec v1)
{
    this->v0 = Vec(cameraX, cameraY, cameraZ);
    this->v1 = v1;
}

Ray::Ray(Vec v0, Vec v1)
{
    this->v0 = v0;
    this->v1 = v1;
}

Vec Ray::getFirst() const
{
    return this->v0;
}

Vec Ray::getSecond() const
{
    return this->v1;
}

Vec Ray::getDirection() const
{
    Vec direction = (v1-v0);
    direction.normalize();
    return direction;
}

Vec Ray::rayPoint(double t) const
{
    return v0+(v1-v0)*t;
}

Color getDiffuse(Color iColor, Color kColor, double fatt){
    if (fatt > 0.0){
        return iColor * kColor * fatt;
    }
    return Color();
}

Color getSpecular(Color iColor, Color kColor, double shininess, double fatt){
    if(fatt > 0.0){
        return iColor * kColor * pow(fatt, shininess);
    }
    return Color();
}

pair<double, Figure*> nearestIntersection(const Ray& r, double minT, double maxT, bool mayBeTransparent = true)
{
    pair<double, Figure*> ret (maxT+minT, NULL);
    for (Figure* figure : shapeList)
    {
        double val = figure->intersection(r, minT, maxT);
        if (val >= minT && val < ret.first) {
            ret.first = val;
            ret.second = figure;
        }   
    }

    return ret;
}

Color RT_trace(const Ray& r, double depth)
{
    double epsilon = 0.00001;
    double maxT = 5000.0;
    pair<double, Figure*> intersection =  nearestIntersection(r, epsilon, maxT);
    Figure* obj = intersection.second;

    if (!obj) {
        return backgroundColor;
    }

    Vec i = r.rayPoint(intersection.first);
    Vec normal = intersection.second->getNormal(r.getDirection());
    bool entering = true;

    return RT_shade(obj, r, i, normal, entering, depth);
}

Color RT_lights(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal)
{
    Color cumLight = Color();
    Color diffuseColor = obj->getColorDiffuse();
    Color specularColor = obj->getColorSpecular();
    Vec inverseDirection = ray.getDirection()*-1.0;

    double shininess = obj->getShininess();
    for (Light* light : lightList) {
        Vec dir = light->getPosition();
        Ray LRay = Ray(i, dir);
        double distance = (i-dir).norm();
        double fattD = light->getFatt(distance);
        Vec LRayDirection = LRay.getDirection();
        double LN = LRayDirection.dot(normal);
        Color iColor = light->getShading();

        bool notBlocking = !(nearestIntersection(LRay, 0.00001, 1.0, false).second);
        if(notBlocking){
            Color diffuse = getDiffuse(iColor, diffuseColor, fattD);
            double fattS = inverseDirection.dot((normal*2.0*fattD)-LRayDirection); 
            Color specular = getSpecular(iColor, specularColor, shininess, fattS);
            cumLight = cumLight + diffuse + specular;
        } 
    }

    return cumLight;
}

Color RT_reflect(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal, double depth)
{
    if(/*(depth <= maxDepth) && */(obj->getRFlag() == 1)){
        Vec inverseDirection = ray.getDirection()*-1.0;
        double s = inverseDirection.dot(normal);
        Ray reflection = Ray(i, i+(normal * 2.0*s)-inverseDirection);
        return obj->getColorReflectivity()*RT_trace(reflection,depth+1);
    }
    return Color();
}

Color RT_transmit(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal, bool entering, double depth)
{
    if((obj->getTFlag() != 1)){
        return Color();
    }

    Vec inverseDirection = ray.getDirection()*-1.0;
    double s =  inverseDirection.dot(normal);
    double indexOfRefraction = obj->getIndexOfRefraction();
    double temp1 = indexOfRefraction * s;
    double temp2 = 1.0 - indexOfRefraction + indexOfRefraction*indexOfRefraction*s*s;

    if(temp2 < 0.0){
        return Color();
    }

    double temp3 = temp1 - sqrt(temp2);
    Ray TRay(i,i+(normal*temp3 - inverseDirection*indexOfRefraction));
    return obj->getColorTransmissivity() * RT_trace(TRay,depth+1);
}

Color RT_shade(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal, bool entering, double depth)
{
    Color newColor = ambient * obj->getColorAmbient();
    Color l = RT_lights(obj, ray, i, normal);
    if (depth < maxDepth){
        Color r = RT_reflect(obj, ray, i, normal, depth);
        Color t = RT_transmit(obj, ray, i, normal, entering, depth);
        return newColor + l + r + t;
    }
    
    return newColor + l;
}

void RT_algorithm() 
{
    for (int u = 0; u < horizontalResolution; u++) {
        double dx = (maxX - minX)/(double) horizontalResolution;
        double dy = (maxY - minY)/(double) verticalResolution;
        for (int v = 0; v < verticalResolution; v++) {

            double x = minX + dx/2.0 + (double)u*dx;
            double y = minY + dy / 2.0 + (double)v * dy;

            Vec point = Vec(x,y,zCoor);
            Ray ray = Ray(point);
            
            pixelColors[u][v] = RT_trace(ray, 1);
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
    ofs << "P3 \n";
    ofs << horizontalResolution << " " << verticalResolution << " \n";
    ofs << 255 << " \n";
    for (int v = verticalResolution-1; v >= 0; v--) {
        for (int u = 0; u < horizontalResolution; u++) {
            (pixelColors[u][v]* 255.0).writeOut(ofs);
            ofs << "\n";
        }
    }
    
    ofs.close();
}


int main(int, char *argv[])
{
    parseSceneFile("scenes\\scene.03");//argv[1]);
    initializeImage();
    RT_algorithm();
    writeImageFile();
    return 0;
}


