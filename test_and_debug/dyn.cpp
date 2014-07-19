#include<iostream>
#include<fstream>
#include<cstdlib>
#include<sstream>

using namespace std;

string fileName(string name,int id)
{
  ostringstream os;
  os<<id;
  name+="_";
  name+=os.str();
  return(name);
}

int main(int argc,char *argv[])
{
  string tim="0.1\n";
  if(argc>2)
    tim=argv[2];
  fstream F;
  F.open("ani",ios::out);
  F<<"set size ratio 1\nset xrange [-1:1]\nset yrange [-1:1]\n";
  for(int i=0;i<atoi(argv[1]);i++)
    F<<"p 'results/"<<fileName("FL",i)<<"' w l,'results/"<<fileName("RN",i)<<"' w l,'bdata' w l,'pdata';pause "<<tim<<"\n";
  F.close();
  return(0);
}
