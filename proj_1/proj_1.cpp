#include<iostream>
#include<cstdlib>
#include"../headers/mesh.h"

using namespace std;

void test_0(double *xv,double *yv,int n,int ns)
{
  system("reset");
  double s=1.0;
  int i,j;
  double dx=s/ns;
  int cnt=0;

  for(i=0;i<ns;i++)
    {
      xv[cnt]=(i+1)*dx;
      yv[cnt]=0;
      cnt++;
    }

  for(i=0;i<ns;i++)
    {
      xv[cnt]=s;
      yv[cnt]=(i+1)*dx;
      cnt++;
    }

  for(i=ns-1;i>=0;i--)
    {
      xv[cnt]=i*dx;
      yv[cnt]=s;
      cnt++;
    }

  for(i=ns-1;i>=0;i--)
    {
      xv[cnt]=0;
      yv[cnt]=i*dx;
      cnt++;
    }

  double minbd=0.8*dx,mind=0.8*minbd,dist;
  double xp,yp,rad;
  double sm=s-2.0*minbd;
  int dec;
  int nb=4*ns;
  i=nb;
  cout<<"Creating internal points starting from id: "<<nb<<"\n";
  while(i<n)
    {
      cout<<i<<" ";
      xp=minbd+(rand()/(RAND_MAX+1.0))*sm;
      yp=minbd+(rand()/(RAND_MAX+1.0))*sm;
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
  int n=100,ns=10;

  if(argc>1)
    n=atoi(argv[1]);
  if(argc>2)
    ns=atoi(argv[2]);

  int nb=4*ns;

  xp=new double [n];
  yp=new double [n];
  test_0(xp,yp,n,ns);//return(0);
  mesh msh(n,nb,xp,yp);
  msh.setDomainBounds(-0.5,1.5,-0.5,1.5);
  msh.makeAnim("dyn","0.01");
  return(0);
}
