#include"../headers/mesh.h"

mesh::mesh(int nP,int nBP,double *xv,double *yv)
{
  nPoints=nP;
  nBPoints=nBP;
  nFaces=10000; // will be corrected
  nCells=1; // will be corrected
  numIter=0;
  xmin=ymin=-1;
  xmax=ymax=1;
  std::cout<<"\nSetting Domain bounds: xmin = "<<xmin<<" xmax = "<<xmax<<" ymin = "<<ymin<<" ymax = "<<ymax<<"\n\n";
  x=new double [nPoints];
  y=new double [nPoints];
  faceId=new int* [nFaces];
  for(int i=0;i<nFaces;i++)faceId[i]=new int [2];
  rng=new ring;
  faceCtr=0;
  maxDis=1.0;
  numDupPoints=0;
  if(!(xv && yv))return;
  copyPoints(xv,yv);
}

void mesh::copyPoints(double *xv,double *yv)
{
  std::cout<<"\nCopying points ... \n\n";
  int i;
  for(i=0;i<nPoints;i++)
    {
      x[i]=xv[i];
      y[i]=yv[i];
    }

  std::cout<<"Creating ring ...\n";
  ring *cur=rng;
  int nxt;
  for(i=0;i<nBPoints;i++)
    {
      cur->id=i;
      if(i<nBPoints-1)
	{
	  cur->next=new ring;
	  nxt=i+1;
	}
      else
	{
	  cur->next=rng;
	  nxt=0;
	}
      cur->next->prev=cur;
      cur=cur->next;
      faceId[i][0]=i;
      faceId[i][1]=nxt;
      faceCtr++;
    }

  touched=new bool [nPoints];
  for(i=0;i<nPoints;i++)
    if(i<nBPoints)
      touched[i]=true;
    else
      touched[i]=false;
  faceCtr=nBPoints;

  system("rm -rf results");
  system("mkdir results");

  collapse();
}

void mesh::setMaxDist(double d)
{
  maxDis=d;
}

bool mesh::isValidCorner(double x0,double y0,double x1,double y1)
{
  return(x0*y1>x1*y0);
}

void mesh::getIntersectionPoints(double a0,double b0,double c0,double a1,double b1,double c1,double& xp,double &yp)
{
  double den=a0*b1-a1*b0;
  xp=(b0*c1-b1*c0)/den;
  yp=(a1*c0-a0*c1)/den;
}

double mesh::distance(double x0,double y0,double x1,double y1)
{
  double disq=pow(x0-x1,2.0)+pow(y0-y1,2.0);
  return(sqrt(disq));
}

void mesh::circumcentre(double *xv,double *yv,int i,int j,int k,double& xc,double &yc,double &r)
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

void mesh::collapse()
{
  std::cout<<"\nStarting ring collapse ...\n";
  int i,j,k;

  int bef,aft;
  bool jStat,curStat;
  ring *cur=rng;
  double xc,yc,rad;
  double dis;
  int rngCnt=nBPoints;
  ring *toDel;
  int iter=0;

  while(rngCnt>3)
    {
      writeRing(cur,iter);
      writeFaces(iter);
      //std::cout<<"here!\n";
      std::cout<<"ring count = "<<rngCnt<<" face count = "<<faceCtr<<" iter = "<<iter<<"\n";
      bef=cur->prev->id;
      aft=cur->next->next->id;
      std::cout<<"bef = "<<bef<<" cur = "<<cur->id<<" next = "<<cur->next->id<<" aft = "<<aft<<"\n";
      curStat=false;

      for(j=nBPoints;j<nPoints;j++)
	{
	  if(touched[j])continue;
	  if(j==cur->id || j==cur->next->id)continue;
	  //if(distance(x,y,cur->id,j)>maxd)continue;//replaced by following line
	  if(distance(x[cur->id],y[cur->id],x[j],y[j])>maxDis)continue;
	  code_B;
	  if(jStat)
	    {
	      curStat=true;
	      if(j>=nBPoints)
		{
		  std::cout<<"kokokoko\n";
		  break;
		}
	    }
	}
      std::cout<<"After J loop\n";
      if(!curStat)
	{
	  std::cout<<"In bef-if\n";
	  int fir,sec,thi;
	  fir=bef;
	  sec=cur->id;
	  thi=cur->next->id;
	  if(!isValidCorner(x[sec]-x[fir],y[sec]-y[fir],x[thi]-x[sec],y[thi]-y[sec]))
	    {
	      std::cout<<"Invalid corner in bef-if\n";
	    }
	  else
	    {
	      j=bef;
	      code_B;
	      curStat=jStat;
	    }
	}
      if(!curStat)
	{
	  std::cout<<"In aft-if\n";
	  int fir,sec,thi;
	  fir=cur->id;
	  sec=cur->next->id;
	  thi=aft;
	  if(!isValidCorner(x[sec]-x[fir],y[sec]-y[fir],x[thi]-x[sec],y[thi]-y[sec]))
	    {
	      std::cout<<"Invalid corner in aft-if\n";
	    }
	  else
	    {
	      j=aft;
	      code_B;
	      curStat=jStat;
	    }
	}
      if(curStat)
	{
	  std::cout<<"J found = "<<j<<"\n";
	  if(!touched[j])
	    {
	      faceId[faceCtr][0]=cur->id;
	      faceId[faceCtr][1]=j;
	      faceId[faceCtr+1][0]=cur->next->id;
	      faceId[faceCtr+1][1]=j;
	      faceCtr+=2;
	      //addFace(&f,cur->id,j);
	      //addFace(&f,cur->next->id,j);
	      touched[j]=true;
	      ring *rinEl=new ring;
	      rinEl->next=cur->next;
	      rinEl->prev=cur;
	      rinEl->id=j;
	      cur->next=rinEl;
	      rinEl->next->prev=rinEl;
	      rngCnt++;
	      std::cout<<"section internal\n";
	    }
	  else
	    {
	      if(j==bef)
		{
		  faceId[faceCtr][0]=cur->next->id;
		  faceId[faceCtr][1]=j;
		  faceCtr++;
		  //addFace(&f,cur->next->id,j);
		  cur->prev->next=cur->next;
		  cur->next->prev=cur->prev;
		  toDel=cur;
		  cur=cur->next;
		  delete toDel;
		  std::cout<<"section boundary bef\n";
		}
	      else
		{
		  faceId[faceCtr][0]=cur->id;
		  faceId[faceCtr][1]=j;
		  faceCtr++;
		  //addFace(&f,cur->id,j);
		  toDel=cur->next;
		  cur->next=cur->next->next;
		  toDel->next->prev=cur;
		  delete toDel;
		  std::cout<<"section boundary aft\n";
		}
	      rngCnt--;
	    }
	}

      cur=cur->next;
      iter++;
      std::cout<<"cur->id = "<<cur->id<<" ring count at end = "<<rngCnt<<"\n\n";
    }
  writeRing(cur,iter);
  writeFaces(iter);
  numIter=iter;
}

void mesh::writeRing(ring *curRng,int id)
{
  std::string fName=fileName("results/RN",id);
  std::fstream F(fName.c_str(),std::ios::out);
  ring *cur=curRng;
  do
    {
      F<<x[cur->id]<<" "<<y[cur->id]<<"\n";
      cur=cur->next;
    }
  while(cur!=curRng);
  F<<x[cur->id]<<" "<<y[cur->id]<<"\n";
  F.close();
}

void mesh::writeFaces(int id)
{
  std::string fName=fileName("results/FL",id);
  std::fstream F(fName.c_str(),std::ios::out);
  for(int i=0;i<faceCtr;i++)
    {
      F<<x[faceId[i][0]]<<" "<<y[faceId[i][0]]<<"\n";
      F<<x[faceId[i][1]]<<" "<<y[faceId[i][1]]<<"\n\n";
    }
  F.close();
}

std::string mesh::fileName(std::string name,int id)
{
  std::ostringstream os;
  os<<id;
  name+="_";
  name+=os.str();
  return(name);
}

void mesh::makeAnim(std::string name,std::string tim)
{
  std::fstream F;
  name="results/"+name;
  F.open(name.c_str(),std::ios::out);
  F<<"set size ratio 1\nset xrange ["<<xmin<<":"<<xmax<<"]\nset yrange ["<<ymin<<":"<<ymax<<"]\n";
  for(int i=0;i<=numIter;i++)
    F<<"p '"<<fileName("FL",i)<<"' w l,'"<<fileName("RN",i)<<"' w l;pause "<<tim<<"\n";
  F.close();
}

void mesh::setDomainBounds(double xmn,double xmx,double ymn,double ymx)
{
  xmin=xmn;
  xmax=xmx;
  ymin=ymn;
  ymax=ymx;
  std::cout<<"\nSetting Domain bounds: xmin = "<<xmin<<" xmax = "<<xmax<<" ymin = "<<ymin<<" ymax = "<<ymax<<"\n\n";
}
