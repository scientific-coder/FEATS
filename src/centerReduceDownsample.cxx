#include <boost/lexical_cast.hpp>
#include <iterator>
#include <iostream>
#include <sstream>
#include <vector>

#include<algorithm>
#include "utility.hxx"
#include "feats/utility/na_reader.hxx"

using namespace feats::utility;

//g++ centerReduceDownsample.cxx -std=c++0x -o centerReduceDownsample -I../include -I/home/bernard/Code/repositories/boost-trunk -lstdc++ -Wall -Wno-sign-compare -O4

int main(int argc, char* argv[]){
  typedef double data_type;
  bool badArgs(argc!=4);

  if(!badArgs){
    bool const shouldCenter(boost::lexical_cast<bool>(argv[1]));
    bool const shouldReduce(boost::lexical_cast<bool>(argv[2]));
    int downsampleStep(boost::lexical_cast<int>(argv[3]));
    typedef std::vector<data_type> container_type;
    std::string lineBuffer;

    while(std::getline(std::cin, lineBuffer)){
      std::istringstream tmp(lineBuffer);
      std::istream_iterator<na_reader<data_type> >b(tmp),e ;
      container_type s(b,e);
      if(downsampleStep>1){
	s.resize(std::distance(s.begin(),downsample(downsampleStep,s.begin(),s.end(),s.begin())));
      }
      if(shouldCenter){
	center(s.begin(), s.end(), s.begin());
      }
      if(shouldReduce){
	reduce(s.begin(), s.end(), s.begin());
      }
      std::copy(s.begin(), s.end(), std::ostream_iterator<data_type>(std::cout,"\t"));
    std::cout<<'\n';
    }
  }
    return 0;
  }

