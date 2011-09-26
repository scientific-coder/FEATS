#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>

// g++ tabsym.cxx -o tabsym -lstdc++
int main(){
  /*
  std::vector<std::vector<int> >v(4,std::vector<int>(4));
  for (unsigned int i(0); i!=4; ++i){
      for (unsigned int j(i); j!=4; ++j){
	v[i][j]=i+j;
      }
  }
  for (unsigned int i(0); i!=4; ++i){
      for (unsigned int j(i); j!=4; ++j){
	v[j][i]=-v[i][j];
      }
  }
  for (unsigned int i(0); i!=4; ++i){
      for (unsigned int j(0); j!=4; ++j){
	std::cout<<"v["<<i<<"]["<<j<<"]="<< v[i][j]<<'\n';
      }
  }
  */
  std::vector<int> v;
  for (unsigned int i(0); i!=4; ++i){ v.push_back(i);}
  std::copy(v.begin(), v.end(), std::ostream_iterator<int>(std::cout, "\t"));
  v.erase(v.begin());
  std::cout<<std::endl;
  std::copy(v.begin(), v.end(), std::ostream_iterator<int>(std::cout, "\t"));
  std::cout<<std::endl;

  return 0;
}
