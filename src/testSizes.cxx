#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <vector>

//struct titi{
struct toto{
  int test()const{return 5;}
};
//  static int value_from_toto(toto t){return t.test();}
//};

int main(int argc, char* argv[]){
  typedef std::vector<toto> container_type;
  typedef container_type::const_iterator c_iter;
  container_type c(5);
  //  typedef int (*letype)(const toto&);
  typedef boost::function<int (toto const *const)> b_f;
  typedef boost::transform_iterator<b_f, c_iter> trans_it;

  trans_it it(c.begin(),b_f(&toto::test));

  /*
  bool t(boost::lexical_cast<bool>(argv[1]));
  bool f(boost::lexical_cast<bool>(argv[2]));
  if(t)std::cout<<"t is "<<t<<std::endl;
  if(!f)std::cout<<"f is "<<f<<std::endl;
  */
  /*
  int nbSegs(boost::lexical_cast<int>(argv[1]));
  int seqSize(boost::lexical_cast<int>(argv[2]));

  int currentSize(0);
  for(int remain=seqSize; remain !=0; --nbSegs, remain-=currentSize){
    currentSize=remain/nbSegs;
    std::cout<<currentSize<<' ';
  }
  */

  return 0;
}
