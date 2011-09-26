#ifndef FEATS_LOCAL_OPTIMIZATION_HXX
#define FEATS_LOCAL_OPTIMIZATION_HXX

template< typename Algo> void next_segmentation(const Algo& a){ a.inc_segments();}
template< typename Seg, typename Nb> void next_segmentation(const merges<Seg,Nb>& a){ a.dec_segments();}


/*
  this struct is an Algo that performs a segmentation with Algo and then does a local optimisation on a sliding window of segments (until convergence or not) with OptimAlgo.

Parameters are
 */
template <typename Algo, typename OptimAlgo>
struct local_optimization {
  typedef Algo seg_algo;
  typedef OptimSegAlgo optim_algo;
  typedef typedef init_segments<optim_algo> init_seg_optim;
  typedef typename seg_algo::segment_type segment_type;
  typedef std::vector<segment_type> segments_cont; 


  template<typename In>
  local_optimization( boost::tuple<In,In> segsBeginEnd)
    :segments_(boost::get<0>(segsBeginEnd),boost::get<1>(segsBeginEnd)), cost_(std::accumulate())
  {}
  
  cost_type compute(size_type winSize, size_type step=1, bool converge=true){
    cost_type odlCost(cost_);
    
    

private:
  template<typename It>
  void optimize_window(It begin, It end){
    seg_algo seg(initial_segmentation(boost::make_tuple(beg->begin(), boost::prior(end)->end()), beg->model()));
    while(seg.nb_segs() != std::distance(begin, end))
      { next_segmentation(seg);}
    seg.segments(begin);
  }
  
  segments_cont segments_;
  
};

#endif
