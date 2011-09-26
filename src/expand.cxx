// expand segmentation of levels
// first column is the length of segments
// other columns are dupplicated on output
#include <vector>
#include <iterator>
#include <string>
#include <iostream>
#include <sstream>
// g++ expand.cxx -o expand -lstdc++
int main(){
  int total(0);
  std::string lineBuffer;
  while(std::getline(std::cin, lineBuffer)){
    std::istringstream tmp(lineBuffer);
    int length;
    tmp>>length;
    std::string restOfLine;
    std::getline(tmp,restOfLine);
    //    std::cerr<<lineBuffer<<'|'<<length<<'\t'<<(total+=length)<<'\n';
    while(length >0 && length--){ std::cout<<restOfLine<<'\n';}
  }
  return 0;
}
