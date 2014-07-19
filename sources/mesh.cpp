#include"../headers/mesh.h"

mesh::mesh(int nP,int nBP,double *xv,double *yv)
{
  pi=4.0*atan(1.0);
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
  //printAngles();return;
  //test();
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
	  if(!checkCorner(cur->next->id,bef,cur->id,j))continue;
	  if(!checkCorner(aft,cur->id,cur->next->id,j))continue;
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
    F<<"p 'results/"<<fileName("FL",i)<<"' w l,'results/"<<fileName("RN",i)<<"' w l;pause "<<tim<<"\n";
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

double mesh::getAngle(double xp,double yp)
{
  if(xp==0 && yp==0)
    {
      std::cout<<"\nInderminate point (0. 0). Terminating from getAngle.\n\n";
      exit(0);
    }
  double tht=atan(fabs(yp)/fabs(xp));
  if(xp>0.0)
    {
      if(yp>=0.0)
	return(tht);
      else
	return(2.0*pi-tht);
    }
  else
    {
      if(yp>0.0)
	return(pi-tht);
      else
	return(pi+tht);
    }
}

bool mesh::isWithinCorner(double xf,double yf,double xb,double yb,double xo,double yo,double xp,double yp)
{
  double xdf=xf-xo;
  double ydf=yf-yo;
  double xdb=xb-xo;
  double ydb=yb-yo;
  double xdp=xp-xo;
  double ydp=yp-yo;
  double thtf=getAngle(xdf,ydf);
  double thtb=getAngle(xdb,ydb);
  double thtp=getAngle(xdp,ydp);

  bool cond;
  cond=(thtf<thtb && thtp>thtf && thtp<thtb) || (thtf>thtb && (thtp>thtf || thtp<thtb));

  return(cond);
}

bool mesh::checkCorner(int nxt,int prv,int cur,int pid)
{
  //std::cout<<"\ncheckCorner("<<nxt<<","<<prv<<","<<cur<<","<<pid<<")\n\n";
  return(isWithinCorner(x[nxt],y[nxt],x[prv],y[prv],x[cur],y[cur],x[pid],y[pid]));
}

void mesh::printAngles()
{
  srand(time(NULL));
  double rad,tht;
  bool cond;
  std::fstream F;

  double tht0,tht1;
  tht0=rand()/(RAND_MAX+1.0)*2*pi;
  tht1=rand()/(RAND_MAX+1.0)*2*pi;
  std::cout<<"\nFirst Angle = "<<tht0*180.0/pi<<"\n";
  std::cout<<"Second angle = "<<tht1*180.0/pi<<"\n\n";
  F.open("bounds",std::ios::out);
  F<<cos(tht0)<<" "<<sin(tht0)<<"\n";
  F<<"0 0\n";
  F<<cos(tht1)<<" "<<sin(tht1)<<"\n";
  F.close();

  double xb,yb,xf,yf;
  xf=cos(tht0);
  yf=sin(tht0);
  xb=cos(tht1);
  yb=sin(tht1);

  F.open("points",std::ios::out);
  for(int i=0;i<1000;i++)
    {
      rad=rand()/(RAND_MAX+1.0);
      tht=rand()/(RAND_MAX+1.0)*2*pi;
      if(isWithinCorner(xf,yf,xb,yb,0,0,rad*cos(tht),rad*sin(tht)))
	F<<rad*cos(tht)<<" "<<rad*sin(tht)<<"\n";
    }
  F.close();

  return;

  // std::cout<<pi<<"\n";
  // std::fstream F("angf",std::ios::out);
  // std::fstream G("angb",std::ios::out);
  // std::fstream H("ang",std::ios::out);
  // int cnt=0;
  // ring *cur=rng;
  // double angf,angb;
  // do
  //   {
  //     angf=getAngle(x[cur->next->id]-x[cur->id],y[cur->next->id]-y[cur->id])*180.0/pi;
  //     angb=getAngle(x[cur->prev->id]-x[cur->id],y[cur->prev->id]-y[cur->id])*180.0/pi;
  //     F<<cnt<<" "<<angf<<"\n";
  //     G<<cnt<<" "<<angb<<"\n";
  //     H<<cnt<<"     "<<angf<<"     "<<angb<<"\n";
  //     cur=cur->next;
  //     cnt++;
  //   }
  // while(cur!=rng);
  // F.close();
  // G.close();
  // H.close();
}

bool mesh::isWithinDomain(int i,int j)
{
  if(i<nBPoints && j<nBPoints)
    {
      //std::cout<<"i = "<<i<<" and j = "<<j<<std::endl;
      int nxt,prv;
      bool cond0,cond1;
      nxt=i+1;
      prv=i-1;
      if(i==0)
	prv=nBPoints-1;
      if(i==nBPoints-1)
	nxt=0;
      //std::cout<<"one\n";
      if(j==prv || j==nxt)return(true);
      cond0=checkCorner(nxt,prv,i,j);
      nxt=j+1;
      prv=j-1;
      if(j==0)
	prv=nBPoints-1;
      if(j==nBPoints-1)
	nxt=0;
      //if(i==prv || i==nxt)return(true);//not needed actually
      cond1=checkCorner(nxt,prv,j,i);
      //std::cout<<"two "<<cond0<<" "<<cond1<<"\n";
      return(cond0 && cond1);
    }
  else if(i>=nBPoints && j>=nBPoints)
    {
      return(!checkIntersectionWithBoundary(i,j));
    }
  else
    {
      return(!checkIntersectionWithBoundary(i,j,true));
    }
}

// void mesh::convertToLine(double x0,double y0,double x1,double y1,double& a,double& b,double& c)
// {
//   a=y1-y0;
//   b=x0-x1;
//   c=x1*y0-x0*y1;
// }

bool mesh::checkIntersectionWithBoundary(int i,int j,bool cond)
{
  int k,nxt;
  for(k=0;k<nBPoints;k++)
    {
      nxt=k+1;
      if(k==nBPoints-1)
	nxt=0;
      if(cond)
	if(nxt==i || nxt==j || k==i || k==j)
	  continue;
      if(checkIntersectionOfLines(i,j,k,nxt))
	return(true);
    }
  return(false);
}

bool mesh::checkIntersectionOfLines(int i00,int i01,int i10,int i11)
{
  double m[2];
  double a[2][2];
  double b[2];

  a[0][0]=x[i01]-x[i00];
  a[0][1]=x[i10]-x[i11];
  a[1][0]=y[i01]-y[i00];
  a[1][1]=y[i10]-y[i11];
  b[0]=x[i10]-x[i00];
  b[1]=y[i10]-y[i00];
  double det=a[0][0]*a[1][1]-a[1][0]*a[0][1];
  m[0]=(a[1][1]*b[0]-a[0][1]*b[1])/det;
  m[1]=(-a[1][0]*b[0]+a[0][0]*b[1])/det;

  if(m[0]>=0 && m[0]<=1 && m[1]>=0 && m[1]<=1)
    return(true);
  return(false);
}

// void mesh::test()
// {
//   if(isWithinDomain(36,38))
//     std::cout<<"36 and 38 within domain\n";
//   if(isWithinDomain(38,36))
//     std::cout<<"38 and 36 within domain\n";
//   exit(0);
// }
