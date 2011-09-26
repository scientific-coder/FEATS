#include <iterator>
#include <iostream>
#include <fstream>
#include<vector>
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

#include "feats/misc/croeuc.hxx"


/*
for NBC in `seq 2 12`; do for NBS in `seq 2  12`; do cat courbes144.txt |~/Feats/src/test_croeucG $NBC $NBS 5  ${NBC}_${NBS}_initSeg.txt ${NBC}_${NBS}_protos.txt ${NBC}_${NBS}_segProtos.txt ${NBC}_${NBS}_segSeries.txt 2>/dev/null >${NBC}_${NBS}_initProtos.txt; done; done




variantes pour essayer d'avoir des classes de variables plus intéressantes:
- une segmentation par classe, puis Inf des partitions (union des points de soupure)
- pb si bcp de classes : on risque de ne plus avoir du tout de segmentation globale: il faut choisir un sous ensemble des segmentations qui repréente au mieux les individus.

- faire des classes de partitions et des prototypes de partitions


1) nb de segments variable suivant les classes
*/

template<typename InSeg, typename InSize, typename Out>
Out expand(InSeg in,InSeg ie, InSize is, Out o){
  while(in!=ie){
    for(typename std::iterator_traits<InSize>::value_type i(*is); i!=0; --i, ++o){
      *o=*in;
    }
    ++in; ++is;
  }
  return o;
}

int main(int argc, char* argv[]){
  if(argc==8){
  int nbClusters(boost::lexical_cast<int>(argv[1]));
  int nbSegments(boost::lexical_cast<int>(argv[2]));

  int nbRetries(boost::lexical_cast<int>(argv[3]));


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
  for( series_set_type::const_iterator it(series.begin()); it != series.end(); ++it)
    { seriesIts.push_back(boost::make_tuple(it->begin(), it->end())); }

  //  std::cout<<"clusters sizes original series:"; 
  croeuc<sit_it_pair> cr(boost::make_tuple(seriesIts.begin(), seriesIts.end()), nbClusters,nbSegments, nbRetries, true);
  std::vector<int> sizes(nbSegments);
  std::vector<data_type> tmp(nbSegments);
  cr.optimize();
  //std::cout<<"segments sizes individual series:";
  cr.output_individual_segments_size(sizes.begin());
  {
    std::ofstream f(argv[4]);
      for(int i(0); i!=series.size(); ++i){
	means(series[i].begin(), series[i].end(), sizes.begin(),tmp.begin());
	expand(tmp.begin(), tmp.end(),sizes.begin(),std::ostream_iterator<data_type>(f,"\t"));f<<"\n";}
  }
      //  std::cout<<"segments sizes clustered series:";
      //  cr.output_segments_size(sizes.begin());
      //  std::cout<<std::endl;

  //  cr.next_segmentation();
  // cr.next_clustering();
      //  std::cout<<"prototypes"<<std::endl;
  {
    std::ofstream f(argv[5]);
  for(int i(0); i!=nbClusters; ++i){

cr.output_prototype(i,std::ostream_iterator<data_type>(f,"\t"));f<<"\n";}
  }
  //  std::cout<<"segmented prototypes"<<std::endl;
  {
    std::ofstream f(argv[6]);
  for(int i(0); i!=nbClusters; ++i){
    cr.output_segmented_prototype(i,tmp.begin());
    expand(tmp.begin(), tmp.end(),sizes.begin(),std::ostream_iterator<data_type>(f,"\t"));f<<"\n";}
  }
  //  std::cout<<"segmented_series"<<std::endl;
  {
    std::ofstream f(argv[7]);
    for(int i(0); i!=series.size(); ++i){
      cr.output_segmented_series(i,tmp.begin());
 expand(tmp.begin(), tmp.end(),sizes.begin(),std::ostream_iterator<data_type>(f,"\t"));f<<"\n";}
  }
  
  }
  return 0;
}


      

