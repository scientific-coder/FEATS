#include <iterator>
#include <iostream>
#include<vector>
#include <functional>

#include <numeric>
#include <string>
#include <sstream>
#include <boost/lexical_cast.hpp>

// g++ splitter.cxx -o splitter -lstdc++
int main(int argc, char* argv[]){
  int nbElts(boost::lexical_cast<int>(argv[1]));
  typedef double data_type;
  typedef std::vector<data_type> container_type;
  
  std::string lineBuffer;

  container_type c;
  while(std::getline(std::cin, lineBuffer)){
    std::istringstream tmp(lineBuffer);
    std::istream_iterator<data_type> b(tmp),e ;
    std::copy(b,e,std::back_inserter(c));
    unsigned int lostLasts(c.size()%nbElts);
    lostLasts=lostLasts?nbElts-lostLasts:lostLasts;
    std::cerr<<"losing "<<lostLasts<<"last elements\n";
    for (unsigned int i(0); i!=(c.size()-lostLasts);++i){
      //      std::cerr<<"((i+1)%nbElts)="<<((i+1)%nbElts)<<'\n';
      std::cout<<c[i]<<(((i+1)%nbElts)?'\t':'\n');
    }
    { container_type eraser; c.swap(eraser);}
  }
  
  return 0;
}
