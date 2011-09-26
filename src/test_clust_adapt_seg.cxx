#include <iterator>
#include <iostream>
#include <iomanip>

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

#include "feats/segmentations/segmentations_for_total_number.hxx"
#include "feats/segmentations/models.hxx"
#include "feats/segmentations/segment.hxx"
#include "feats/segmentations/optimal_splits.hxx"

#include "feats/misc/clusters_of_adaptive_segmentations.hxx"
#include "feats/utility/na_reader.hxx"

using namespace feats::segmentations;
using namespace feats::utility;

// g++ -std=c++0x test_clust_adapt_seg.cxx -I../include/ -o  test_clust_adapt_seg -lstdc++ -O4

/*

variantes pour essayer d'avoir des classes de variables plus intéressantes:
1) nb de segments variable suivant les classes
 */

int main(int argc, char* argv[]){
  int nbClustersMin(boost::lexical_cast<int>(argv[1]));
  int nbClustersMax(boost::lexical_cast<int>(argv[2]));
  int nbSegmentsMin(boost::lexical_cast<int>(argv[3]));
  int nbSegmentsMax(boost::lexical_cast<int>(argv[4]));

  int nbRetries(boost::lexical_cast<int>(argv[5]));

  int printSeries(boost::lexical_cast<int>(argv[6]));

  typedef double data_type;
  typedef std::vector<data_type> container_type;
  typedef container_type::const_iterator c_iter_type;
  typedef std::vector<container_type> series_set_type;
  
  std::string lineBuffer;
  series_set_type series;

  while(std::getline(std::cin, lineBuffer)){
    series.push_back(container_type());
    std::istringstream tmp(lineBuffer);
    std::istream_iterator<na_reader<data_type> >b(tmp),e ;
    std::copy(b,e,std::back_inserter(series.back()));
  }
  std::cerr<<"reading "<<series.size()<<" series\n";
  
  typedef std::vector<boost::tuple<c_iter_type, c_iter_type> > sit_type;
  typedef sit_type::iterator set_it;
  typedef boost::tuple<set_it, set_it> sit_it_pair;

  sit_type seriesIts;
  for( series_set_type::const_iterator it(series.begin()); it != series.end(); ++it)
    { seriesIts.push_back(boost::make_tuple(it->begin(), it->end())); }

  std::cerr.setf(std::ios::fixed);
  std::cerr<<std::setprecision(4);

   std::cout.setf(std::ios::fixed);
   std::cout<<std::setprecision(4);
   for(int i(nbClustersMin); i!=nbClustersMax; ++i){
     for(int j(std::max(nbClustersMin,nbSegmentsMin)); j!=nbSegmentsMax; ++j){
       feats::misc::clusters_of_adaptive_segmentations<sit_it_pair> clustAdapSegs(boost::make_tuple(seriesIts.begin(), seriesIts.end()), i, j, nbRetries, false);
       std::cout<<clustAdapSegs.nb_iter()<<'\t'<<clustAdapSegs.cost()<<'\t';
       if(printSeries){
	 for(int k(0); k!=i; ++k, std::cout<<'@'<<k<<'|')
	   { clustAdapSegs.output_prototype(k,std::ostream_iterator<data_type>(std::cout,"\t")); }
	 std::cout<<std::endl;
       }
     }
   }
   return 0;
}


      

