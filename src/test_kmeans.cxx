#include <iterator>
#include <iostream>
#include <vector>
#include <functional>

#include <numeric>
#include <string>
#include <sstream>
#include <boost/utility.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/progress.hpp>

#include "utility.hxx"

#include "feats/clustering/kmeans.hxx"

//g++ -std=c++0x test_kmeans.cxx -O4 -I../include -I/home/bernard/Code/repositories/boost-trunk  -lstdc++ -o test_kmeans

int main(int argc, char* argv[]){
  int nbClusters(boost::lexical_cast<int>(argv[1]));
  int nbRetries(boost::lexical_cast<int>(argv[2]));
  //  bool onlyCost(boost::lexical_cast<true>(argv[3])); no use yet

  typedef double data_type;
  typedef std::vector<data_type> container_type;
  typedef container_type::const_iterator c_iter_type;
  typedef std::vector<container_type> series_set_type;
  
  std::string lineBuffer;
  series_set_type series;

  while(std::getline(std::cin, lineBuffer)){
    series.push_back(container_type());

    std::istringstream tmp(lineBuffer);
    std::istream_iterator<data_type> b(tmp),e ;
    std::copy(b,e,std::back_inserter(series.back()));

  }
  std::cerr<<"reading "<<series.size()<<" series\n";
  
  typedef std::vector<boost::tuple<c_iter_type, c_iter_type> > sit_type;
  typedef sit_type::iterator set_it;
  typedef boost::tuple<set_it, set_it> sit_it_pair;

  sit_type seriesIts;
  for( series_set_type::const_iterator it(series.begin()); it != series.end(); ++it){
    //    sit_it_pair tmp(it->begin(), it->end());
    seriesIts.push_back(boost::make_tuple(it->begin(), it->end()));
  }
  
  typedef kmeans<sit_it_pair> kmean_type;
  kmean_type km(boost::make_tuple(boost::make_tuple(seriesIts.begin(), seriesIts.end()), nbClusters));
  km.compute();
  km.output_sizes(std::ostream_iterator<int>(std::cout," "));
  for(int i(0); i!=nbClusters; ++i){km.output_prototype(i,std::ostream_iterator<data_type>(std::cout," "));std::cout<<std::endl;}

  boost::shared_ptr<kmean_type> kmp(best_kmean<kmean_type>(boost::make_tuple(boost::make_tuple(seriesIts.begin(), seriesIts.end()), nbClusters), nbRetries));
  kmp->output_sizes(std::ostream_iterator<int>(std::cout," "));std::cout<<"\n";
  for(int i(0); i!=nbClusters; ++i){kmp->output_prototype(i,std::ostream_iterator<data_type>(std::cout," "));std::cout<<std::endl;}


  typedef std::vector<int> w_seq;
  typedef weighted_kmeans<sit_it_pair,w_seq::const_iterator > w_kmean_type;

  w_seq wSeq(series.back().size(),1);wSeq[2]=100;
  w_kmean_type wkm(boost::make_tuple(boost::make_tuple(seriesIts.begin(), seriesIts.end()), boost::make_tuple(wSeq.begin(), wSeq.end()),nbClusters));
  wkm.compute();
  wkm.output_sizes(std::ostream_iterator<int>(std::cout," "));
  for(int i(0); i!=nbClusters; ++i){wkm.output_prototype(i,std::ostream_iterator<data_type>(std::cout," "));std::cout<<std::endl;}

  
  boost::shared_ptr<w_kmean_type> wkmp(best_kmean<w_kmean_type>(boost::make_tuple(boost::make_tuple(seriesIts.begin(), seriesIts.end()), boost::make_tuple(wSeq.begin(), wSeq.end()), nbClusters), nbRetries));
  wkmp->output_sizes(std::ostream_iterator<int>(std::cout," "));std::cout<<"\n";
  for(int i(0); i!=nbClusters; ++i){wkmp->output_prototype(i,std::ostream_iterator<data_type>(std::cout," "));std::cout<<std::endl;}
 
  
  return 0;
}


      

