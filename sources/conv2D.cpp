#include "../headers/conv2D.h"

void conv2D::makeRing()
{
  std::string s = zName + "_vid";
  remove(s.c_str());
  int i, nxt;
  ringo = new polyr;
  polyr *cur = ringo;
  for(i = 0; i < size; i++)
  {
    nxt = i + 1;
    if(i == size - 1)
    {
      nxt = 0;
    }
    computeLine(x[i], y[i], x[nxt], y[nxt], cur->a, cur->b, cur->c);
    cur->id = i;
    if(i < size - 1)
    {
      cur->next=new polyr;
    }
    else
    {
      cur->next=ringo;
    }
    cur=cur->next;
  }
}

void conv2D::decomposeRing()
{
  polyr *befStart, *start, *end, *ins, *befEnd, *startNew, *endNew;
  int numSE;
  double xcenSE, ycenSE, xcenSI, ycenSI;
  double a, b, c;
  bool insStat;
  double x0, y0, x1, y1;
  polyr *cur, *cur2;
  int dec;

  befStart = ringo;
  start = befStart->next;
  befEnd = start;
  end = start->next;
  numSE = 2;
  xcenSE = 0.5 * (x[start->id] + x[end->id]);
  ycenSE = 0.5 * (y[start->id] + y[end->id]);
  ins = NULL;

  while(1)
    {
      if(ins==start)
	{
	  addPolyz(start,xcenSE,ycenSE);
	  priZ(start);
	  zcon++;
	  break;
	}
      ins=end->next;
      xcenSI=(xcenSE*numSE+x[ins->id])/(numSE+1);
      ycenSI=(ycenSE*numSE+y[ins->id])/(numSE+1);
      computeLine(x[ins->id],y[ins->id],x[start->id],y[start->id],a,b,c);

      insStat=true;
      while(1)
	{
	  // Check 0
	  x0=x[end->id]-x[befEnd->id];
	  y0=y[end->id]-y[befEnd->id];
	  x1=x[ins->id]-x[end->id];
	  y1=y[ins->id]-y[end->id];
	  if(!angleCheck(x0,y0,x1,y1))
	    {
	      insStat=false;
	      break;
	    }
	  // Check 1
	  x0=x[start->id]-x[ins->id];
	  y0=y[start->id]-y[ins->id];
	  x1=x[start->next->id]-x[start->id];
	  y1=y[start->next->id]-y[start->id];
	  if(!angleCheck(x0,y0,x1,y1))
	    {
	      insStat=false;
	      break;
	    }
	  // Check 2
	  cur=ins->next;
	  while(cur!=start)
	    {
	      cur2=start;
	      dec=0;
	      while(cur2!=ins)
		{
		  if(!checkInclusion(cur2,cur,xcenSI,ycenSI))
		    {
		      dec=1;
		      break;
		    }
		  cur2=cur2->next;
		}
	      if(dec){cur=cur->next;continue;}
	      if(checkInclusion(a,b,c,cur,xcenSI,ycenSI))
		{
		  insStat=false;
		  break;
		}
	      cur=cur->next;
	    }
	  // Breaking out of while
	  break;
	}

      if(insStat)
	{
	  xcenSE=xcenSI;
	  ycenSE=ycenSI;
	  numSE++;
	  befEnd=end;
	  end=ins;
	  ins=ins->next;
	}
      else
	{
	  if(numSE==2)
	    {
	      befStart=start;
	      start=end;
	      befEnd=end;
	      end=ins;
	      ins=ins->next;
	      xcenSE=0.5*(x[start->id]+x[end->id]);
	      ycenSE=0.5*(y[start->id]+y[end->id]);
	    }
	  else
	    {
	      startNew=new polyr;
	      endNew=new polyr;
	      startNew->id=start->id;
	      startNew->next=endNew;
	      endNew->id=end->id;
	      endNew->next=ins;
	      befStart->next=startNew;

	      end->next=start;
	      computeLine(x[end->id],y[end->id],x[start->id],y[start->id],a,b,c);
	      end->a=a;
	      end->b=b;
	      end->c=c;
	      addPolyz(start,xcenSE,ycenSE);
	      priZ(start);
	      zcon++;

	      start=startNew;
	      end=endNew;
	      start->a=a;
	      start->b=b;
	      start->c=c;
	      computeLine(x[end->id],y[end->id],x[ins->id],y[ins->id],a,b,c);
	      end->a=a;
	      end->b=b;
	      end->c=c;
	      befEnd=start;
	      numSE=2;
	      xcenSE=0.5*(x[start->id]+x[end->id]);
	      ycenSE=0.5*(y[start->id]+y[end->id]);
	    }
	}
    }
}

bool conv2D::angleCheck(const double x0,const double y0,const double x1,const double y1)
{
  return(x0*y1>=x1*y0);
}

void conv2D::computeLine(const double x0,const double y0,const double x1,const double y1,double& a,double& b,double& c)
{
  a=y0-y1;
  b=x1-x0;
  c=x0*y1-x1*y0;
}

bool conv2D::checkInclusion(const polyr *lseg,const polyr *pnt,const double xc,const double yc)
{
  return(checkInclusion(lseg->a,lseg->b,lseg->c,pnt,xc,yc));
}

bool conv2D::checkInclusion(const double a,const double b,const double c,const polyr *pnt,const double xc,const double yc)
{
  double valp=a*x[pnt->id]+b*y[pnt->id]+c;
  double valc=a*xc+b*yc+c;
  return(valp*valc>0);
}

void conv2D::addPolyz(polyr *ring,const double xcen,const double ycen)
{
  polyz *cur=new polyz;
  cur->zone.ring=ring;
  cur->zone.xcen=xcen;
  cur->zone.ycen=ycen;

  cur->next=NULL;
  if(zo==NULL)
    zo=cur;
  else
    {
      cur->next=zo;
      zo=cur;
    }
}

bool conv2D::polyr::contains(double xp,double yp,double xc,double yc)
{
  double valp=a*xp+b*yp+c;
  double valc=a*xc+b*yc+c;
  return(valp*valc>=0.0);
}

bool conv2D::polyg::contains(double xp,double yp)
{
  polyr *cur=ring;
  do
    {
      if(!cur->contains(xp,yp,xcen,ycen))
	return(false);
      cur=cur->next;
    }
  while(cur!=ring);
  return(true);
}

conv2D::conv2D(int n,double *xv,double *yv,std::string zoname)
{
  zName=zoname;
  size=n;
  zo=NULL;
  zcon=0;
  x=new double [n];
  y=new double [n];
  if(xv==NULL || yv==NULL)
    return;
  copyData(xv,yv);
}

void conv2D::copyData(double *xv,double *yv)
{
  for(int i=0;i<size;i++)
    {
      x[i]=xv[i];
      y[i]=yv[i];
    }
  makeRing();
  decomposeRing();
}

bool conv2D::contains(const double xp,const double yp)
{
  polyz *cur=zo;
  while(cur)
  {
    if(cur->zone.contains(xp,yp))
      return(true);
    cur=cur->next;
  }
  return(false);
}

void conv2D::printZones()
{
  std::string file,cent;
  polyz *cur=zo;
  polyr *rng;
  int numz=0;
  if(!cur)
    {
      std::cout<<"No zones present.\n";
      return;
    }
  file=zName+"_zones";
  cent=zName+"_centroids";
  std::fstream F(file.c_str(),std::ios::out);
  std::fstream G(cent.c_str(),std::ios::out);
  while(cur)
    {
      rng=cur->zone.ring;
      do
	{
	  F<<x[rng->id]<<" "<<y[rng->id]<<"\n";
	  rng=rng->next;
	}
      while(rng!=cur->zone.ring);
      F<<x[rng->id]<<" "<<y[rng->id]<<"\n\n";
      G<<cur->zone.xcen<<" "<<cur->zone.ycen<<"\n";
      cur=cur->next;
      numz++;
    }
  F.close();
  G.close();
  std::cout<<"NUMZ = "<<numz<<"\n";
}

void conv2D::priZ(polyr *rng)
{
  std::ostringstream os;
  os<<zName<<"Zone"<<zcon;
  std::fstream F(os.str().c_str(),std::ios::out);
  polyr *cur=rng;
  do
    {
      F<<x[cur->id]<<" "<<y[cur->id]<<"\n";
      cur=cur->next;
    }
  while(cur!=rng);
  F<<x[cur->id]<<" "<<y[cur->id]<<"\n";
  F.close();

  std::string s=zName+"_vid";
  F.open(s.c_str(),std::ios::out|std::ios::app);
  if(zcon==0)
    F<<"set xrange [-1:1]\nset yrange [-1:1]\nset size ratio 1\nunset key\n";
  F<<"p ";
  for(int i=0;i<=zcon;i++)
    {
      std::ostringstream ostr;
      ostr<<zName<<"Zone"<<i;
      F<<"'"<<ostr.str()<<"' w l lw 3";
      if(i<zcon)
	F<<",";
    }
  F<<"\npause 0.3\n";
  F.close();
}
