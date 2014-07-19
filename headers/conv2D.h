#ifndef CONV2D_H
#define CONV2D_H

#include<iostream>
#include<cstdlib>
#include<fstream>
#include<cmath>
#include<sstream>

class conv2D
{
  double *x,*y;
  int size;
  std::string zName;

  typedef struct Polyring
  {
    double a,b,c;
    int id;
    struct Polyring *next;
    Polyring(){next=NULL;}
    bool contains(double xp,double yp,double xc,double yc);
  }polyr;

  typedef struct Polygon
  {
    polyr *ring;
    double xcen,ycen;
    bool contains(double xp,double yp);
  }polyg;

  typedef struct Polyzones
  {
    polyg zone;
    struct Polyzones *next;
  }polyz;

  polyz *zo;
  polyr *ringo;

  int zcon;

  void makeRing();
  void decomposeRing();
  bool angleCheck(const double x0,const double y0,const double x1,const double y1);
  void computeLine(const double x0,const double y0,const double x1,const double y1,double& a,double& b,double& c);
  bool checkInclusion(const polyr *lseg,const polyr *pnt,const double xc,const double yc);
  bool checkInclusion(const double a,const double b,const double c,const polyr *pnt,const double xc,const double yc);
  void addPolyz(polyr *ring,const double xcen,const double ycen);
  void priZ(polyr *rng);
public:
  conv2D(int n,double *xv=NULL,double *yv=NULL,std::string zoname="region");
  void copyData(double *xv,double *yv);
  bool contains(const double xp,const double yp);
  void printZones();
};

#endif
