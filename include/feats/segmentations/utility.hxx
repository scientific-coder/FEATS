#ifndef FEATS_SEGMENTATIONS_UTILITY_HXX
#define FEATS_SEGMENTATIONS_UTILITY_HXX

#include <vector>
#include <boost/utility.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_container_iterator.hpp>

#include "feats/segmentations/models.hxx"

namespace feats{
  namespace segmentations{


/***
    template functors to dispatch initial strivial segmentation for algorithms
    first one (default for topdown alog (binary_splits and optimal_splits) computes an initial segmentation of 1 segment.
    second one (sepcialized for bottomup (merges) computes an initial segmentation of N segments where N is the number of points
 ***/
template<typename Segmentation> struct init_segments{
  typedef typename Segmentation::segment_type segment_type;
  typedef typename segment_type::model_type model_type;
  typedef typename segment_type::iterator data_iterator;
  typedef boost::shared_container_iterator<std::vector<segment_type> >seg_iterator;
  boost::tuple<seg_iterator, seg_iterator> operator()( boost::tuple<data_iterator, data_iterator> begEnd, const model_type& initM=model_type()){
    typedef std::vector<segment_type> c_type;
    boost::shared_ptr< c_type > ptr(new c_type());
    ptr->push_back(segment_type(boost::get<0>(begEnd), boost::get<1>(begEnd), initM));
    return boost::make_tuple(seg_iterator(ptr->begin(),ptr), seg_iterator(ptr->end(),ptr));
  }
};
/// could be specialized for linear1_model to have init segs of size 2
template<typename Seg, typename Nb> struct init_segments<merges<Seg, Nb> >{
  typedef Seg segment_type;
  typedef typename segment_type::model_type model_type;
  typedef typename segment_type::iterator data_iterator;
  typedef boost::shared_container_iterator<std::vector<segment_type> >seg_iterator;
  boost::tuple<seg_iterator, seg_iterator> operator()( boost::tuple<data_iterator, data_iterator> begEnd, const model_type& initM=model_type()){
    typedef std::vector<segment_type> c_type;
    boost::shared_ptr< c_type > ptr(new c_type());
    // segs of 1 elts should assert distance(beg,end)>1
    for(data_iterator segEnd(boost::next(boost::get<0>(begEnd))); boost::get<0>(begEnd) !=boost::get<1>(begEnd)
	  ;boost::get<0>(begEnd)=segEnd, segEnd=(segEnd==boost::get<1>(begEnd))?segEnd:++segEnd){
      ptr->push_back(segment_type(boost::get<0>(begEnd),segEnd, initM));
    }
    return boost::make_tuple(seg_iterator(ptr->begin(),ptr), seg_iterator(ptr->end(),ptr));
  }
};


  }
}

#endif
