#include <iterator>
#include <iostream>
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

variantes pour essayer d'avoir des classes de variables plus intéressantes:
- une segmentation par classe, puis Inf des partitions (union des points de soupure)
- pb si bcp de classes : on risque de ne plus avoir du tout de segmentation globale: il faut choisir un sous ensemble des segmentations qui repréente au mieux les individus.

- faire des classes de partitions et des prototypes de partitions


1) nb de segments variable suivant les classes
 */
// g++ -std=c++0x test_croeuc.cxx -I../include -I/home/bernard/Code/repositories/boost-trunk -lstdc++
int main(int argc, char* argv[]){
  int nbClusters(boost::lexical_cast<int>(argv[1]));
  int nbSegments(boost::lexical_cast<int>(argv[2]));

  int nbRetries(boost::lexical_cast<int>(argv[3]));

  int printSeries(boost::lexical_cast<int>(argv[4]));

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

   std::cout<<"clusters sizes original series:"; 
  croeuc<sit_it_pair> cr(boost::make_tuple(seriesIts.begin(), seriesIts.end()), nbClusters,nbSegments, nbRetries, false);
  std::cout<<"clusters sizes segmented series:";
  cr.output_clusters_size(std::ostream_iterator<data_type>(std::cout,"\t"));
  std::cout<<std::endl;
  std::cout<<"segments sizes global series:";
  cr.output_global_segments_size(std::ostream_iterator<data_type>(std::cout,"\t"));
  std::cout<<std::endl;
  std::cout<<"segments sizes individual series:";
  cr.output_individual_segments_size(std::ostream_iterator<data_type>(std::cout,"\t"));
  std::cout<<std::endl;
  std::cout<<"segments sizes clustered series:";
  cr.output_segments_size(std::ostream_iterator<data_type>(std::cout,"\t"));
  std::cout<<std::endl;

  //  cr.next_segmentation();
  // cr.next_clustering();
  /*
  for(int i(0); i!=nbClusters; ++i){cr.output_prototype(i,std::ostream_iterator<data_type>(std::cout,"\t"));std::cout<<"@";}
  std::cout<<std::endl;
  std::cout<<cr.segmentation_cost()<<'\t'<<cr.clustering_cost()<<std::endl;
  for(int i(0); i!=nbClusters; ++i){cr.output_segmented_prototype(i,std::ostream_iterator<data_type>(std::cout,"\t"));std::cout<<"@";}
   std::cout<<std::endl;
   if(printSeries){
   for(int i(0); i!=series.size(); ++i){cr.output_segmented_series(i,std::ostream_iterator<data_type>(std::cout,"\t"));std::cout<<"@";}
   std::cout<<std::endl;
   }
  */
  return 0;
}


      

