#include <limits>
#include <vector>
#include <tuple>
#include <memory>
#include <iterator>
#include <iostream>
#include <algorithm>

#include <boost/lexical_cast.hpp>
// TODO find best abstractions
// g++-snapshot -DPROTO -std=c++11 simple-proto-seg.cxx -o simple-proto-seg.S -g3 -O4 -Wall -S -march=native -fverbose-asm
/*
We used an idiomatic basic C++ iterator based API.
(could also use a C++11 container absed API with move)
The data to segment in K segments is in [beg, end).
The segments for segmentations in K-1 segments of [i,end) with i \in [beg,end)
are in [computed_beg, computed_end)
segments are recursive tuples
 */

/*
  The segment data structure contains the cost (sum of square errors), the episod size and a pointer to the next segment.
 */

template<typename CostType, typename SizeType=std::size_t>
struct segment {
  typedef CostType cost_type;
  typedef SizeType size_type;
  explicit segment(cost_type c= std::numeric_limits<cost_type>::max(), size_type s=0, segment const* n=0)
    :cost_(c), size_(s), next_(n){}
  cost_type cost() const{return cost_;}
  size_type size() const{return size_;}
  segment const* next()const{return next_;}
  cost_type cost_;
  size_type size_;
  segment const* next_;
};
// fake segment iterator used by models for initial segmentation
template<typename SegmentType> struct fake_segment_iterator {
  typedef SegmentType value_type;
  fake_segment_iterator& operator++() { return *this;}
  struct fake_segment {
    typename SegmentType::cost_type cost() const {return 0.;}
    typename SegmentType::size_type size() const {return 0;}
    SegmentType* operator&() const { return 0;}
  };
  fake_segment const& operator*()const{ return instance;}
  static fake_segment instance;
};
template<typename SegmentType> 
typename fake_segment_iterator<SegmentType>::fake_segment fake_segment_iterator<SegmentType>::instance= 
  typename fake_segment_iterator<SegmentType>::fake_segment();

namespace std
{ template<typename S> struct iterator_traits<fake_segment_iterator<S> > : iterator_traits<S*> {}; }

// needed because std::back_insert_iterator<> do not define ::value_type
template <typename It> struct generic_traits : std::iterator_traits<It> {};
template <typename C> struct generic_traits<std::back_insert_iterator<C> >
: std::iterator_traits<typename C::iterator> {};

template<typename DataIt, typename SegIt, typename ModelData=int> struct regular_model {
  typedef typename std::iterator_traits<SegIt>::value_type seg_type;
  typedef typename seg_type::cost_type cost_type;
  typedef typename seg_type::size_type size_type;
  // dummy arg is here for API compatibility with models that need them
  template<typename Dummy=int>
  regular_model(DataIt data_beg, SegIt seg_beg, Dummy const& unused = Dummy())
    : data_it(data_beg), seg_it(seg_beg), sum(0.), sum_sq(0.), cost(0.), size(0.) {}
  regular_model& operator++(){
    //    std::cerr<<"adding "<<*data_it<<std::endl;
    sum+= *data_it;
    sum_sq+= (*data_it) * (*data_it);
    ++size;
    ++data_it;
    ++seg_it;
    cost= sum_sq-(sum*sum) / size;
    return *this;
  }
  bool operator<(regular_model const& o) const
  { return (cost + seg_it->cost()) < (o.cost + o.seg_it->cost());}
  explicit operator seg_type ()const {return seg_type(cost+ (*seg_it).cost(), size, &(*seg_it));}
  explicit operator typename std::iterator_traits<DataIt>::value_type () const {return sum/size;}// size==0 -> na ?
  SegIt next() const {return seg_it;} 
  DataIt data_it;
  SegIt seg_it;
  cost_type sum, sum_sq, cost;
  size_type size;
};

template<typename DataIt, typename SegIt, typename Prototypes> struct prototypes_model {
  typedef typename std::iterator_traits<SegIt>::value_type seg_type;
  typedef typename seg_type::cost_type cost_type;
  typedef typename seg_type::size_type size_type;
  prototypes_model(DataIt data_beg, SegIt seg_beg, Prototypes const& p)
    : data_it(data_beg), seg_it(seg_beg), protos(p), size(0.), costs(p.size(), 0.)
    , best_cost_it(costs.begin()){}
  prototypes_model(prototypes_model const& o): data_it(o.data_it), seg_it(o.seg_it), protos(o.protos), size(o.size), costs(o.costs), best_cost_it(costs.begin()+ (o.best_cost_it-o.costs.begin() )){}

  prototypes_model& operator=(prototypes_model const& o){
    // assert same prototypes
    data_it= o.data_it;
    seg_it= o.seg_it;
    size=  o.size;
    costs= o.costs;
    best_cost_it= costs.begin()+(o.best_cost_it-o.costs.begin());
    return *this;
  }

  prototypes_model& operator++(){
    //    std::cerr<<"adding "<<*data_it<<std::endl;
    for(std::size_t i(0); i != protos.size(); ++i){
      cost_type const tmp(protos[i]- (*data_it));
      costs[i]+= tmp * tmp;
    }
    ++size;
    ++data_it;
    ++seg_it;
    best_cost_it= std::min_element(costs.begin(), costs.end());
    return *this;
  }
  bool operator<(prototypes_model const& o) const
  { //std::cerr<<"comparing "<< *best_cost_it <<" + "<<seg_it->cost() << " = "<<(*best_cost_it + seg_it->cost())<<" to "<<(*(o.best_cost_it))<<" + "<<o.seg_it->cost()<<" = "<<(*(o.best_cost_it) + o.seg_it->cost())<<std::endl;
    return (*best_cost_it + seg_it->cost()) < (*(o.best_cost_it) + o.seg_it->cost());}
  explicit operator seg_type () const { return seg_type(*best_cost_it + (*seg_it).cost(), size, &(*seg_it));}
  explicit operator typename std::iterator_traits<DataIt>::value_type () const 
  {return protos[best_cost_it-costs.begin()];}// size==0 -> na ?
  SegIt next() const {return seg_it;} 
  
  DataIt data_it;
  SegIt seg_it;
  Prototypes const& protos;
  size_type size;
  std::vector<cost_type> costs;
  typename std::vector<cost_type>::iterator best_cost_it;
};

/* The generic optimal segmentation algorithm for one more segment. Dynamic programing makes it O(n^2).
 */
template<template<class DI, class SI, class M> class ModelType
         , typename ItData, typename ItSeg, typename ModelData, typename Out>
Out compute_seg(ItData beg, ItData end, ItSeg computed_beg, ItSeg computed_end
                , ModelData const& model_data, Out result_out) {
  typedef ModelType<ItData, ItSeg, ModelData> model_type;
  typedef typename std::iterator_traits<ItSeg>::value_type seg_type;
  typedef typename seg_type::cost_type cost_type;
  typedef typename std::iterator_traits<ItData>::difference_type size_type;

  ItData start_data(beg); // we won't compute N=(end-beg) segmentations because there is no seg in K segments for the last K-1 points, so we iterate in computed_ instead
   // we compute a seg starting from start_data
  for(ItSeg start_seg(computed_beg); start_seg != computed_end; ++start_seg, ++start_data, ++result_out){
    model_type best_model(start_data, start_seg, model_data);
    // for each seg, we test starting from start_data 
    for(model_type current_model(best_model); (current_model).next() != computed_end; ++current_model){
      if( current_model < best_model){
        //        std::cerr<<"optim found for size "<<current_model.size<<" and cost "<< (seg_type(current_model).cost())<<"next size is "<<(current_model.next()->size())<<std::endl;
        best_model= current_model;
      }
    }
    //    std::cerr<<"seg done\n";
    *result_out= seg_type(best_model);
  }
  return result_out;
}


/* the initial step computing "segmentations" in one segment.
 */
template<template<class ID, class IS, class MD> class ModelType, typename ItData, typename ModelData, typename ItSeg>
ItSeg compute_initial_seg(ItData beg, ItData end, ModelData const& model_data, ItSeg out) {
  typedef typename generic_traits<ItSeg>::value_type seg_type;
  typedef typename seg_type::cost_type cost_type;
  typedef typename std::iterator_traits<ItData>::difference_type size_type;
  typedef std::reverse_iterator<ItData> rev_it_type;
  ModelType<rev_it_type, fake_segment_iterator<seg_type>, ModelData > model(rev_it_type(end), fake_segment_iterator<seg_type>(), model_data);
  // for efficiency reasons, we compute the segments in reverse order
  // so we store them in a temporary buffer and output the buffer in reverse order at the end.
  std::vector<seg_type> tmp;
  for(rev_it_type it(end), rend(beg); it != rend; ++it, ++model, tmp.push_back(seg_type(model)))
    {}
  return std::copy(tmp.rbegin(), tmp.rend(), out);
}
/* Out put segmented series values "unpacking" the list of segments.
As the mean is not stored into the segment data structure, we have to compute it.
 */
template<template<class ID, class IS, class MD> class ModelType
         , typename Seg, typename ItData, typename ModelData, typename Out>
Out print_segmentation(Seg seg, ItData it_data, ModelData const& model_data, Out out) {
  for(Seg const* s(&seg); s != 0; s=s->next()){
    ModelType<ItData, fake_segment_iterator<Seg>, ModelData > model(it_data, fake_segment_iterator<Seg>(), model_data);
    for(typename Seg::size_type i(0); i != s->size(); ++i, ++model, ++it_data)
    {}
    std::cerr<<"output seg of size "<<s->size()<<std::endl;
    for(typename Seg::size_type i(0); i != s->size(); ++i, ++out){
      *out= typename std::iterator_traits<ItData>::value_type(model);
    }
  }
  return out;
}

template<template<class ID, class IS, class MD> class ModelType
         , typename ItData, typename ModelData=int>
void do_it(ItData beg, ItData end, std::size_t nb_segs, ModelData const& model_data= ModelData()){
  typedef typename std::iterator_traits<ItData>::value_type data_type;
  typedef segment<data_type> segment_type;
  typedef std::vector<segment_type> segmentation_type;
  typedef std::vector<segmentation_type> segs_type;
  segs_type all_segs(nb_segs);
  for(std::size_t i(0); i != nb_segs; ++i){
    if(i==0){
      compute_initial_seg<ModelType>(beg, end, model_data, std::back_inserter(all_segs.front()));
    }else {
      compute_seg<ModelType>(beg, end, all_segs[i-1].begin(), all_segs[i-1].end(), model_data
                             , std::back_inserter(all_segs[i]));
    }
  }
  print_segmentation<ModelType>(all_segs.back().front(), beg, model_data, std::ostream_iterator<data_type>(std::cout, " "));
}  

// main program, reading times series from std::cin and outputing the optimal segmented time series to std::cout in the requested nb of segments (given as the arg, defaults to 2)
int main(int argc, char* argv[]){
  typedef double data_type;
  typedef std::vector<data_type> container_t;
  typedef std::istream_iterator<data_type> input_t;
  container_t data(input_t(std::cin), (input_t()));

  // nb segs is given as argument (default to 2) but must not exceed the nb of points
  std::size_t const nb_segs(std::min(argc >1 ? boost::lexical_cast<std::size_t>(argv[1]):2
                                     , data.size())); 

  container_t protos(2);
  protos[0]= 0; protos[1]=10;
  if(true){ do_it<regular_model>(data.begin(), data.end(), nb_segs); }
  else { do_it<prototypes_model>(data.begin(), data.end(), nb_segs, protos); }
  return 0;
}
