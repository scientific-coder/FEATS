#include <iterator>
#include <limits>
#include <algorithm>
#include <ext/algorithm> // random_sample is a SGI extension
#include <iostream>

#include <vector>
#include <functional>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/iterator/transform_iterator.hpp>

//g++ -std=c++0x test-interlace.cxx -O4 -I../include -I/home/bernard/Code/repositories/boost-trunk  -lstdc++ -o test-interlace

struct Test{
  typedef int size_type;
  typedef double data_type;
  typedef std::vector<data_type> container_type;
  typedef container_type::iterator values_iterator;
  
  Test(size_type nbProtos, size_type nbPoints):protosVals_(nbProtos*nbPoints){
    typedef std::vector<container_type> series_set;
    series_set initPop(nbProtos*2, container_type(nbPoints));
    typedef boost::tuple<series_set::const_iterator,series_set::const_iterator> serie_iterators_pair;
    typedef boost::tuple<values_iterator,values_iterator> values_iterators_pair;
    for(size_type i(0); i!= initPop.size(); ++i)
      { std::fill_n(initPop[i].begin(),nbPoints, i);}
    
    std::vector<values_iterators_pair> allIts(initPop.size()),protosIt(nbProtos);
    for(size_type i(0); i!= allIts.size(); ++i)
      { allIts[i]=boost::make_tuple(initPop[i].begin(), initPop[i].end());}
    
    __gnu_cxx::random_sample(allIts.begin(), allIts.end(), protosIt.begin(), protosIt.end());
    set_proto_values(boost::make_tuple(protosIt.begin(), protosIt.end()));
  }

  template<typename InitItTuple>
  void set_proto_values(const InitItTuple& initP){ 
    typedef typename boost::tuples::element<0, InitItTuple>::type proto_it_iterator;
    typedef typename std::iterator_traits<proto_it_iterator>::value_type proto_it_tuple;
    typedef typename boost::tuples::element<0, proto_it_tuple>::type proto_it;
    size_type nbProtos(std::distance(boost::get<0>(initP),boost::get<1>(initP))); 
    values_iterator protosValsIt(protosVals_.begin()), protoValsInitIt(protosVals_.begin());

    for(proto_it_iterator it(boost::get<0>(initP)); it != boost::get<1>(initP)
	  ; ++it, ++protoValsInitIt, protosValsIt=protoValsInitIt){
      for(proto_it sit(boost::get<0>(*it)); sit!= boost::get<1>(*it); ++sit, protosValsIt+=nbProtos){
	std::cerr<<*sit<<'|';
	*protosValsIt=*sit;
      }
      std::cerr<<'\n';
    }
    std::cerr<<"\n\n";
  }

  void print()const{ std::copy(protosVals_.begin(), protosVals_.end(),std::ostream_iterator<data_type>(std::cout,"\t")); std::cout<<std::endl;}
    
private:
  container_type protosVals_;
};

int main(){

  Test t(10,3);
  t.print();
  return 0;
}
