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

#include "feats/segmentations/segmentations_for_total_number.hxx"
#include "feats/segmentations/models.hxx"
#include "feats/segmentations/segment.hxx"
#include "feats/segmentations/optimal_splits.hxx"

#include "feats/misc/clusters_of_adaptive_segmentations.hxx"
#include "feats/utility/na_reader.hxx"

using namespace feats::segmentations;
using namespace feats::utility;

// g++ -std=c++0x test_seg_with_total.cxx -I../include/ -o test_seg_with_total -lstdc++

/*

variantes pour essayer d'avoir des classes de variables plus intéressantes:
1) nb de segments variable suivant les classes
 */

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

   std::cout<<"nb d'épisodes par serie"; 
   typedef segment<linear0_model<data_type>,c_iter_type > segment_type;
   segmentations_for_total_number<optimal_split<segment_type > > segm(boost::make_tuple(seriesIts.begin(), seriesIts.end()));

   while(segm.nb_segments()!=nbSegments) { segm.inc_segments(); }

   segm.nb_episodes(std::ostream_iterator<int>(std::cout,"\t"));
   std::cout<<std::endl;

   typedef  feats::segmentations::seg_to_tuple<segment_type> convert_type;
   typedef std::ostream_iterator< segment_type::data_type> os_it_type;
   for(unsigned int i(0); i!=series.size(); ++i){
     std::vector<segment_type> res;
     segm.segments(i,std::back_inserter(res));
     
     std::copy(boost::make_transform_iterator(res.begin(),convert_type(res.front().begin()))
	       , boost::make_transform_iterator(res.end(), convert_type(res.front().begin()))
	       , std::ostream_iterator<convert_type::result_type>(std::cout,"\t"));
     std::cout<<std::endl; 
   }

   feats::misc::clusters_of_adaptive_segmentations<sit_it_pair> clustAdapSegs(boost::make_tuple(seriesIts.begin(), seriesIts.end()), nbClusters,nbSegments, nbRetries, false);

   return 0;
}


      

