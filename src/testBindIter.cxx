#include <vector>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/bind.hpp>
struct toto{
  int operator()(int t)const{return t+5;}
};
int main(){

  std::vector<toto> tmp(10);
  std::vector<int> 
    res(boost::make_transform_iterator(tmp.begin(), boost::bind(&toto::operator(),_1,2))
	,boost::make_transform_iterator(tmp.end(), boost::bind(&toto::operator(),_1,2)));
  return 0;
}
