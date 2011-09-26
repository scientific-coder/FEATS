#include <iostream>

#include <vector>

int main(){
  typedef std::vector<int> c_type;
  std::vector<c_type> v1(3,c_type(2,1)), v2(2,c_type(3,555));
  v1.swap(v2);
  std::cout<<v1.size()<<' '<<v1.front().size();
}
