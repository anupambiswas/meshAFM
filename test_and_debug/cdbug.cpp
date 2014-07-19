#include<iostream>
#include<cstdlib>
#include<cmath>
#include<fstream>
#include<sstream>

using namespace std;

void getIntersectionPoints(double a0,double b0,double c0,double a1,double b1,double c1,double& xp,double &yp)
{
  double den=a0*b1-a1*b0;
  xp=(b0*c1-b1*c0)/den;
  yp=(a1*c0-a0*c1)/den;
}

void circumcentre(double *xv,double *yv,int i,int j,int k,double& xc,double &yc,double &r)
{
  double a0,b0,c0,a1,b1,c1;
  double xp,yp;

  a0=xv[j]-xv[i];
  b0=yv[j]-yv[i];
  c0=-0.5*(xv[j]*xv[j]+yv[j]*yv[j]-xv[i]*xv[i]-yv[i]*yv[i]);
  a1=xv[k]-xv[i];
  b1=yv[k]-yv[i];
  c1=-0.5*(xv[k]*xv[k]+yv[k]*yv[k]-xv[i]*xv[i]-yv[i]*yv[i]);
  getIntersectionPoints(a0,b0,c0,a1,b1,c1,xc,yc);
  r=sqrt(pow(xv[i]-xc,2.0)+pow(yv[i]-yc,2.0));
}

void makeCircle(double xc,double yc,double r,string file)
{
  fstream F;
  int i;
  double tht,pi=4.0*atan(1.0),xp,yp;
  F.open(file.c_str(),ios::out);
  for(i=0;i<100;i++)
    {
      tht=i*2*pi/99;
      xp=r*cos(tht)+xc;
      yp=r*sin(tht)+yc;
      F<<xp<<" "<<yp<<endl;
    }
  F<<(r+xc)<<" "<<yc<<endl;
  F.close();
}

string fileName(string name,int id)
{
  std::ostringstream os;
  os<<id;
  name+="_";
  name+=os.str();
  return(name);
}

int main(int argc,char *argv[])
{
  double *x,*y;
  int i,j,k;

  int n=100;
  int nb=40;
  if(argc>3)
    n=atoi(argv[3]);
  if(argc>4)
    nb=atoi(argv[4]);

  x=new double [n];
  y=new double [n];

  fstream F("pdata",ios::in);
  for(i=0;i<n;i++)
    F>>x[i]>>y[i];
  F.close();

  double xc,yc,r;

  i=atoi(argv[1]);
  j=atoi(argv[2]);

  string file;
  system("rm -rf circl");
  system("mkdir circl");
  for(k=0;k<n;k++)
    {
      if(k==i || k==j)continue;
      file=fileName("circl/circle",k);
      circumcentre(x,y,i,j,k,xc,yc,r);
      makeCircle(xc,yc,r,file);
    }

  string tim;
  cout<<"Enter time gap: ";
  cin>>tim;
  cout<<"\n";
  F.open("cirplo",ios::out);
  F<<"set xrange [-1:1]; set yrange [-1:1]; set size ratio 1;\n";
  for(k=0;k<n;k++)
    {
      if(k==i || k==j)continue;
      F<<"p 'bdata' w l,'pdata','"<<fileName("circl/circle",k)<<"' w l;pause "<<tim<<"\n";
    }
  F.close();

  return(0);
}
