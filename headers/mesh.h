#ifndef MESH_H
#define MESH_H

#include<iostream>
#include<cmath>
#include<fstream>
#include<sstream>
#include<cstdlib>

//#define code_B printf("Hi\n");
#define code_B								\
  circumcentre(x,y,cur->id,cur->next->id,j,xc,yc,rad);			\
  jStat=true;								\
  for(k=0;k<nPoints;k++)						\
    {									\
      if(k==cur->id || k==cur->next->id || k==j)continue;		\
      if(distance(x[cur->id],y[cur->id],x[k],y[k])>maxDis)continue;	\
      dis=distance(xc,yc,x[k],y[k]);					\
      if(dis<rad)							\
	{								\
	  jStat=false;							\
	  break;							\
	}								\
      if(cur->id<nBPoints && cur->next->id<nBPoints && j<nBPoints && k<nBPoints) \
	continue;							\
      if(dis==rad)							\
	{								\
	  std::cout<<"\nDegenrate condition faced.\n";			\
	  std::cout<<"Centroid: X = "<<xc<<" Y = "<<yc<<" rad = "<<rad<<"\n"; \
	  std::cout<<"Point: X = "<<x[k]<<" Y = "<<y[k]<<"\n";		\
	  std::cout<<"Distance point-centroid: "<<dis<<"\n";		\
	  std::cout<<"Point Ids: J = "<<j<<" K = "<<k<<"\n";		\
	  std::cout<<"Terminating.\n\n";				\
	  std::fstream G;						\
	  G.open("deg",std::ios::out);					\
	  G<<x[cur->id]<<" "<<y[cur->id]<<"\n";				\
	  G<<x[j]<<" "<<y[j]<<"\n";					\
	  G<<x[cur->next->id]<<" "<<y[cur->next->id]<<"\n";		\
	  G<<x[cur->id]<<" "<<y[cur->id]<<"\n\n";			\
	  G<<x[k]<<" "<<y[k]<<"\n";					\
	  G.close();							\
	  exit(0);							\
	}								\
    }									\
  ;									\


class mesh
{
  typedef struct Ring
  {
    int id;
    struct Ring *next;
    struct Ring *prev;
  }ring;

  double *x,*y;
  int nPoints,nBPoints;
  int nCells,nFaces;
  int numDupPoints;
  int faceCtr,numIter;
  double maxDis;
  double xmin,ymin,xmax,ymax;
  ring *rng;
  int **faceId;
  bool *touched;

  bool isValidCorner(double x0,double y0,double x1,double y1);
  void getIntersectionPoints(double a0,double b0,double c0,double a1,double b1,double c1,double& xp,double &yp);
  double distance(double x0,double y0,double x1,double y1);
  void circumcentre(double *xv,double *yv,int i,int j,int k,double& xc,double &yc,double &r);
  void collapse();
  void writeRing(ring *curRng,int id);
  void writeFaces(int id);
  std::string fileName(std::string name,int id);

 public:
  mesh(int nP,int nBP,double *xv=NULL,double *yv=NULL);
  void copyPoints(double *xv,double *yv);
  void setMaxDist(double d);
  void makeAnim(std::string name="ani",std::string tim="0.1");
  void setDomainBounds(double xmn,double xmx,double ymn,double ymx);
};

#endif
