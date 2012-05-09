#include <limits>
#include <vector>
#include <tuple>
#include <iterator>
#include <iostream>
#include <algorithm>

#include <boost/lexical_cast.hpp>
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
  segment():cost(std::numeric_limits<cost_type>::max()), size(0), next(0){}
  cost_type cost;
  size_type size;
  segment* next;
};

/* The generic optimal segmentation algorithm for one more segment. Dynamic programing makes it O(n^2).
 */
template<typename ItData, typename ItSeg, typename Out>
Out compute_seg(ItData beg, ItData end, ItSeg computed_beg, ItSeg computed_end
                   , Out result_out) {
  typedef typename std::iterator_traits<ItSeg>::value_type seg_type;
  typedef typename seg_type::cost_type cost_type;
  typedef typename std::iterator_traits<ItData>::difference_type size_type;

  ItData start_data(beg); // we won't compute N=(end-beg) segmentations because there is no seg in K segments for the last K-1 points, so we iterate in computed_ instead
   // we compute a seg starting from start_data
  for(ItSeg start_seg(computed_beg); start_seg != computed_end; ++start_seg, ++start_data, ++result_out){
    seg_type current_seg;
    cost_type sum(0.), sum_sq(0.);
    size_type nb_elts(1);
    ItData it_data= start_data;
    // for each seg, we test strating from start_data 
    for(ItSeg it_seg(start_seg); it_seg != computed_end;++it_seg, ++it_data, ++nb_elts){
      sum+= *it_data; sum_sq+= (*it_data)* (*it_data);
      cost_type const current_cost((sum_sq-(sum*sum)/nb_elts) + it_seg->cost);
      //      std::cerr<<" cost is "<<(sum_sq-(sum*sum)/nb_elts)<< " + "<<it_seg->cost <<" = "<<current_cost<<" ";
      if(current_cost < current_seg.cost) {
        //std::cerr<<"optim found for size "<<nb_elts<<" next size is "<<it_seg->size<<std::endl;
        current_seg.cost= current_cost;
        current_seg.next= &(*it_seg);
        current_seg.size= nb_elts;
      }
    }
    *result_out= current_seg;
  }
  return result_out;
}
template <typename It> struct generic_traits : std::iterator_traits<It> {};
template <typename C> struct generic_traits<std::back_insert_iterator<C> >
: std::iterator_traits<typename C::iterator> {};

/* the initial step computing "segmentations" in one segment.
 */
template<typename ItData, typename ItSeg>
ItSeg compute_initial_seg(ItData beg, ItData end, ItSeg out) {
  typedef typename generic_traits<ItSeg>::value_type seg_type;
  typedef typename seg_type::cost_type cost_type;
  typedef typename std::iterator_traits<ItData>::difference_type size_type;
  typedef std::reverse_iterator<ItData> rev_it_type;
  seg_type res;
  cost_type sum(0.), sum_sq(0.);
  // for efficiency reasons, we compute the segments in reverse order
  // so we store them in a temporary buffer and output the buffer in reverse order at the end.
  std::vector<seg_type> tmp;
  for(rev_it_type it(end), rend(beg); it != rend; ++it, tmp.push_back(res)){
    sum+= *it; sum_sq+= (*it)* (*it);
    ++res.size;
    res.cost= sum_sq-(sum*sum)/ res.size;
  }
  return std::copy(tmp.rbegin(), tmp.rend(), out);
}
/* Out put segmented series values "unpacking" the list of segments.
As the mean is not stored into the segment data structure, we have to compute it.
 */
template< typename Seg, typename ItData, typename Out>
Out print_segmentation(Seg seg, ItData it_data, Out out) {
  for(Seg* s(&seg); s != 0; s=s->next){
    typename std::iterator_traits<ItData>::value_type current_mean=0.;
    for(typename Seg::size_type i(0); i != s->size; ++i, ++it_data)
    { current_mean+= *it_data;}// TODO Kahan sum
    current_mean/= s->size;
    for(typename Seg::size_type i(0); i != s->size; ++i, ++out){
      *out= current_mean;
    }
  }
  return out;
}

// main program, reading times series from std::cin and outputing the optimal segmented time series to std::cout in the requested nb of segments (given as the arg, defaults to 2)
int main(int argc, char* argv[]){
  typedef double data_type;
  typedef std::istream_iterator<data_type> input_t;
  std::vector<data_type> data(input_t(std::cin), (input_t()));
  typedef std::vector<segment<data_type> > seg_type;

  // nb segs is given as argument (default to 2) but must not exceed the nb of points
  std::size_t const nb_segs(std::min(argc >1 ? boost::lexical_cast<std::size_t>(argv[1]):2
                                     , data.size())); 

  typedef std::vector<seg_type> segs_type;
  segs_type all_segs(nb_segs);
  for(std::size_t i(0); i != nb_segs; ++i){
    if(i==0){
      compute_initial_seg(data.begin(), data.end(), std::back_inserter(all_segs.front()));
      //      std::reverse(all_segs.front().begin(), all_segs.front().end());
    }else {
      compute_seg(data.begin(), data.end(), all_segs[i-1].begin()+1, all_segs[i-1].end()
                  , std::back_inserter(all_segs[i]));
    }
  }
  print_segmentation(all_segs.back().front(), data.begin(), std::ostream_iterator<data_type>(std::cout, " "));
  return 0;
}
