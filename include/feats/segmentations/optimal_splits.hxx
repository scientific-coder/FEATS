#ifndef FEATS_SEGMENTATION_OPTIMAL_SPLITS_HXX
#define FEATS_SEGMENTATION_OPTIMAL_SPLITS_HXX

#include <vector>
#include <queue>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <boost/bind.hpp>
#include <boost/checked_delete.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include "feats/utility/static_caster.hxx"
namespace feats{
  namespace segmentations{


    /***
	Optimal segmentation. One problem is that I cannot know delta_cost before computing the next segmentation. I could use another class to provide the cost, using optimal for computing n+1 segmentation but it would have to cache the n segmentation. optimal can provide the segmentation for previous segmentation so I will do just that: always compute n+1, but return seg and cost of n segments. Unless the segmentation is full: nb seg = nb points. I cannot then compute further segmentation so I stop and decide next_cost is max possible value.

Optimal segmentation with P segments of p points is computed with complexity O(PN-A²)times the complexity to update the segment models (with is O(1) for linear models order 0 and 1). The algorithm uses dynamic programming. We have to compute (and store) optimal segmentation on P-1 segments-b
	
This is problematic with models having an initial state (such as prototypes), because I need to preserve this information even while not storing segments, hence the need for emptyModel_
    ***/
    template<typename Seg, typename Nb=int> struct optimal_split {
      typedef Seg segment_type;
      typedef Nb nb_segments_type;
      typedef typename segment_type::cost_type cost_type;
      typedef std::vector<cost_type> costs_container;
      typedef typename costs_container::iterator costs_iterator;
      typedef typename costs_container::const_iterator costs_const_iterator;
      template<typename S, typename N=nb_segments_type> struct rebind{typedef optimal_split<S,N> other;};

      typedef typename segment_type::iterator iterator;
      typedef std::vector<iterator> ends_container;
      typedef std::vector<ends_container> ends_containers_container;

      explicit optimal_split(Seg initSeg): begin_(initSeg.begin()), end_(initSeg.end()), nbPoints_(initSeg.size())
					 , costs_(nbPoints_), endPoints_(1,ends_container(nbPoints_?nbPoints_:1,end_)), nextCost_(initSeg.model().cost()), full_(false), emptyModel_(empty_model(initSeg.model())){
	// compute costs with one segment, first cutPoints_ have all be set to end_
	for(costs_iterator it=costs_.begin();initSeg.size(); ++it, initSeg.shrink_begin()){*it=initSeg.model().cost(); }
	totalCosts_.push_back(nextCost_);
	inc_segments(); // compute n+1
      }
      
      template<typename Out>
      Out segments(Out o)const {
	iterator segBeg(begin_), segEnd;
	//	std::cerr<<"optim segments:"<<end_-begin_<<"points et "<<endPoints_.size()<<"segments ";
	for(typename ends_containers_container::const_reverse_iterator endsIt=boost::next(endPoints_.rbegin()) // boost::next because we computed n+1 segs and only want n segs, so we don't take the last endPoints unless we are full
	      ;  endsIt != endPoints_.rend() ; segBeg = segEnd, ++o, ++endsIt){ 
	  //	  std::cerr<<segBeg-begin_<<'|';
	  segEnd=(*endsIt)[segBeg-begin_];
	  *o=segment_type(segBeg, segEnd, emptyModel_);
	}
	return o;
      }

      template<typename Out>
      Out segments(nb_segments_type nbS, Out o){
	while (nbS>endPoints_.size()){ inc_segments(); }
	typename ends_containers_container::reverse_iterator endsIt(endPoints_.rbegin());
	std::advance(endsIt, endPoints_.size()-nbS);
	iterator segBeg(begin_), segEnd;
	for(      ;  endsIt != endPoints_.rend() ; segBeg = segEnd, ++o, ++endsIt){ 
	  //	  std::cerr<<segBeg-begin_<<'|';
	  segEnd=(*endsIt)[segBeg-begin_];
	  *o=segment_type(segBeg, segEnd, emptyModel_);
	}
	return o;
      }

    
      void inc_segments(){
	if(endPoints_.size()<nbPoints_){
	  cost_=nextCost_;
	  nb_segments_type nbS(endPoints_.size());
	  costs_container nextCosts(nbPoints_ - nbS);
	  endPoints_.push_back(ends_container(nbPoints_-nbS));//maybe compiler cannot optimize away the copy, should check with :
	  // endPoints_.push_back(ends_container(0)); endPoints_.back().resize(nbPoints_-nbS)
	
	  costs_const_iterator beginCosts(costs_.begin());
	  costs_const_iterator const endCosts(costs_.end());
	  costs_iterator currentCostIt(nextCosts.begin());
	  typename ends_container::iterator currentEndIt(endPoints_.back().begin());

	  for(iterator beginSegmentation(begin_), endBeginSegmentation(end_-nbS)
		;beginSegmentation!=endBeginSegmentation
		; ++beginSegmentation, ++beginCosts, ++currentEndIt, ++currentCostIt){
	    segment_type tmp(beginSegmentation,emptyModel_); 
	    *currentCostIt= tmp.optimize_with_costs(beginCosts, endCosts);
	    *currentEndIt=tmp.end();
	  }
	  std::swap(nextCosts, costs_);
	  nextCost_=costs_.front();
	}else{

	  // full_=true; cost_=nextCost_; nextCost_=std::numeric_limits<cost_type>::max();
	  // add an empty segment
	  cost_=nextCost_;
	  endPoints_.push_back(ends_container(1,begin_));

	}
	totalCosts_.push_back(nextCost_);
      }

      cost_type cost()const{return cost_;}
      cost_type cost(nb_segments_type nbS){
	while (nbS>endPoints_.size()) { inc_segments(); }
	return totalCosts_[nbS-1];
      }
      cost_type delta_cost()const{ return nextCost_-cost_;}

      nb_segments_type nb_segs()const{ return endPoints_.size()-(full_?0:1);} //-1 because I compute n+1 only for delta_cost
    private:
      iterator begin_, end_;
      typename std::iterator_traits<iterator>::difference_type nbPoints_;
      costs_container costs_, totalCosts_;
      ends_containers_container endPoints_;
      cost_type cost_, nextCost_;
      bool full_;
      const typename segment_type::model_type emptyModel_;
    };




    /***
	This one can split more than on segments, using one optimal_split per inital segment.
	I need to access the optimal_split with the best delta_cost (smallest because delta_cost is <0)
	and update this when the best optimal split is done. For maximum flexibility (with either lots of optimal splits or big optimal splits, I want to avoid a list where updating position would be slow (O(N)) and a priority_queue of optimal_split that would be expensive to copy. I can either use my letfish_tree, or a priority_queue of pointers. I do the later. I also need to access the set of optimal_split in chronological order to ouput segments. TO that effect, I store them in a vector (preservedOrder_)

    ***/
    template<typename Seg, typename Nb=int> struct optimal_splits {
      typedef Seg segment_type;
      typedef typename segment_type::model_type model_type;
      typedef Nb nb_segments_type;
      typedef typename segment_type::cost_type cost_type;
      typedef optimal_split<segment_type> splitter_type;
      typedef std::vector<splitter_type*> seq_type;
      template<typename S, typename N=nb_segments_type> struct rebind{typedef optimal_splits<S,N> other;};
      struct SplittersComparator:std::binary_function<splitter_type*, splitter_type*, bool>{
	bool operator()(splitter_type const* const s1, splitter_type const* const s2)const{return s1->delta_cost()<s2->delta_cost();}
      };

      typedef std::priority_queue<splitter_type*, seq_type, SplittersComparator> queue_type; 

      template <typename In>
      optimal_splits(boost::tuple<In,In> initSegsBegEnd)
	:preservedOrder_(boost::make_transform_iterator(boost::get<0>(initSegsBegEnd), feats::utility::static_caster<segment_type,splitter_type*>())
			 ,boost::make_transform_iterator(boost::get<1>(initSegsBegEnd), feats::utility::static_caster<segment_type,splitter_type*>()))
	,splitters_(preservedOrder_.begin(), preservedOrder_.end())

	/*	,cost_(std::accumulate(boost::make_indirect_iterator(preservedOrder_.begin())
			       ,boost::make_indirect_iterator(preservedOrder_.end())
			       ,0.
			       ,boost::bind(std::plus<typename segment_type::cost_type>(),_1
					    , boost::bind(std::mem_fun_ref(&splitter_type::cost), _2))))
	*/
	, nbSegs_(splitters_.size()) {
	for (typename seq_type::const_iterator it(preservedOrder_.begin()); it!= preservedOrder_.end(); ++it)
	  { cost_+=(*it)->cost(); }
      }
      // could be faster: non need to sort again
      optimal_splits(const optimal_splits& other)
	:preservedOrder_(boost::make_transform_iterator(boost::make_indirect_iterator(other.preservedOrder_.begin())
							,feats::utility::static_caster<splitter_type, splitter_type*>())
			 ,boost::make_transform_iterator(make_indirect_iterator(other.preservedOrder_.end())))
	,splitters_(preservedOrder_.begin(), preservedOrder_.end()), cost_(other.cost_), nbSegs_(other.nbSegs_){}

      // not exception safe !
      optimal_splits& operator=(const optimal_splits& other){
	if(this!=&other){
	  delete_splitters();
	  {
	    seq_type tmp(boost::make_transform_iterator(boost::make_indirect_iterator(other.preservedOrder_.begin())
							,feats::utility::static_caster<splitter_type, splitter_type*>())
			 ,boost::make_transform_iterator(make_indirect_iterator(other.preservedOrder_.end())));
	    preservedOrder_.swap(tmp);
	  }
	  // priority_queue has no swap() :-(
	  while(!splitters_.empty()){ splitters_.pop();}
	  std::for_each(preservedOrder_.begin(), preservedOrder_.end(), boost::bind(&queue_type::push,boost::ref(splitters_),_1));
	  cost_=other.cost_;
	  nbSegs_=other.nbSegs_;
	}
	return *this;
      }
 
      void inc_segments(){
	splitter_type* bestSplit(splitters_.top());
	splitters_.pop();
	cost_+=bestSplit->delta_cost(); // should monitor roundoff errors
	bestSplit->inc_segments();
	splitters_.push(bestSplit);
	++nbSegs_;
      }
      cost_type cost()const{return cost_;}
      cost_type delta_cost()const{return splitters_.front()->delta_cost();}
      nb_segments_type nb_segs()const{ return nbSegs_;}

      template<typename Out>
      Out segments(Out out)const{
	for(typename seq_type::const_iterator it(preservedOrder_.begin()); it!= preservedOrder_.end(); ++it)
	  { *out=(*it)->segments(out);}
	return out;
      }

      ~optimal_splits(){delete_splitters();}

    private:
      void delete_splitters(){std::for_each(preservedOrder_.begin(), preservedOrder_.end(), boost::checked_deleter<splitter_type>());}

      seq_type preservedOrder_;
      queue_type splitters_;
      cost_type cost_;
      nb_segments_type nbSegs_;
    };
  }
}


#endif
