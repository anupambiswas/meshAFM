#include<iostream>
#include<cstdlib>
#include"../headers/mesh.h"

using namespace std;

void test_0(double *xv,double *yv,int n,int nb)
{
  system("reset");
  double r=1.0,tht;
  int i,j;
  double pi=4.0*atan(1.0);
  double dx=2*pi*r/nb;
  double dtht=2*pi/nb;

  for(i=0;i<nb;i++)
    {
      tht=i*dtht;
      xv[i]=r*cos(tht);
      yv[i]=r*sin(tht);
    }

  double minbd=0.8*dx,mind=0.8*minbd,dist;
  double xp,yp,rad;
  double rm=r-minbd;
  int dec;
  i=nb;
  cout<<"Creating internal points starting from id: "<<nb<<"\n";
  while(i<n)
    {
      cout<<i<<" ";//endl;
      rad=(rand()/(RAND_MAX+1.0))*rm;
      tht=(rand()/(RAND_MAX+1.0))*2*pi;
      xp=rad*cos(tht);
      yp=rad*sin(tht);
      dec=0;
      for(j=nb;j<i;j++)
	{
	  dist=sqrt(pow(xp-xv[j],2.0)+pow(yp-yv[j],2.0));
	  if(dist<mind)
	    {
	      dec=1;
	      break;
	    }
	}
      if(dec)
	continue;
      xv[i]=xp;
      yv[i]=yp;
      i++;
    }
  cout<<"\n\n";

  fstream F("bdata",ios::out);
  for(i=0;i<nb;i++)
    F<<xv[i]<<" "<<yv[i]<<endl;
  F<<xv[0]<<" "<<yv[0]<<endl;
  F.close();

  F.open("pdata",ios::out);
  for(i=0;i<n;i++)
    F<<xv[i]<<" "<<yv[i]<<endl;
  F.close();
  cout<<"Test_0 completed\n\n";
}

int main(int argc,char *argv[])
{
  double *xp,*yp;
  int n=100,nb=40;

  if(argc>1)
    n=atoi(argv[1]);
  if(argc>2)
    nb=atoi(argv[2]);

  xp=new double [n];
  yp=new double [n];
  test_0(xp,yp,n,nb);
  mesh msh(n,nb,xp,yp);
  msh.makeAnim("dyn","0.01");
  return(0);
}
