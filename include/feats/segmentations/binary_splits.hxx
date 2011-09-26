#ifndef FEATS_SEGMENTATIONS_BINARY_SPLITS_HXX
#define  FEATS_SEGMENTATIONS_BINARY_SPLITS_HXX

#include<list>
#include<queue>
#include <algorithm>
#include<boost/operators.hpp>
#include<boost/tuple/tuple.hpp>
#include <boost/type_traits/add_reference.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include "feats/utility/static_caster.hxx"
#include "assert.h"
namespace feats{
  namespace segmentations{


    template<typename Seg, typename Nb=int> struct binary_splits {
      typedef Seg segment_type;
      typedef typename segment_type::model_type model_type;
      typedef Nb nb_segments_type;
      typedef typename Seg::cost_type cost_type;
      template<typename S, typename N=nb_segments_type> struct rebind{typedef binary_splits<S,N> other;};
    private:
      typedef std::list<Seg> list_type;
      typedef typename list_type::iterator list_iterator;
      //      typedef typename boost::add_reference<Seg>::type ref_seg;
      typedef boost::tuple<Seg,Seg> segs_type;

      struct Split: boost::less_than_comparable<Split>{
	explicit Split( list_iterator it) : it_(it),segs_(split_in_two(*it_))
	{
#if 0
	  std::cerr<<"Constructing Split from split in two cost and size ptr"<<boost::get<0>(segs_).model().cost()
		   <<" "<<boost::get<0>(segs_).model().size()
		   <<" "<<&(boost::get<0>(segs_).model())
		   <<" "<<" and cost and size ptr"
		   <<boost::get<1>(segs_).model().cost()
		   <<" "<<boost::get<1>(segs_).model().size()
		   <<" "<<&(boost::get<1>(segs_).model())<<'\n';
#endif
 boost::get<0>(segs_).optimize_with(boost::get<1>(segs_)); cost_=  (boost::get<0>(segs_).model().cost()+boost::get<1>(segs_).model().cost()) - it_->model().cost();
 assert(cost_==cost_);
	}
	//check
	bool operator<(const Split& other) const{ return other.cost_ < cost_;}
	cost_type cost()const{return cost_;}
	list_iterator l_iterator()const{return it_;}
	const segs_type& segs()const{return segs_;}

      private:
	list_iterator it_;
	segs_type  segs_;
	cost_type cost_;
      };

    public:
      binary_splits( Seg initSeg):segments_(1,initSeg),splits_(), cost_(initSeg.model().cost()), nbSegs_(1)
      { splits_.push(Split(segments_.begin()));}

      template<typename In>
      binary_splits( boost::tuple<In,In> segsBeginEnd)
	:segments_(boost::get<0>(segsBeginEnd),boost::get<1>(segsBeginEnd))
	,splits_(boost::make_transform_iterator(boost::make_counting_iterator(segments_.begin())
						, feats::utility::static_caster<list_iterator, Split>())
		 ,boost::make_transform_iterator(boost::make_counting_iterator(segments_.end())
						 , feats::utility::static_caster<list_iterator, Split>()))
	/*	,cost_(std::accumulate(segments_.begin(),segments_.end(), 0.
			       ,boost::bind(std::plus<typename segment_type::cost_type>(),_1
					    , boost::bind(std::mem_fun_ref(&model_type::cost)
							  ,boost::bind(std::mem_fun_ref(&segment_type::model)
							       , _2)))))
	*/, nbSegs_(segments_.size())
      { 
	for (list_iterator it(segments_.begin()); it!= segments_.end(); ++it)
	  { cost_+=it->model().cost(); }
	
      }

      void inc_segments() {
	const Split& spl(splits_.top());
	list_iterator it0(segments_.insert(segments_.erase(spl.l_iterator()),boost::get<1>(spl.segs())));
	list_iterator it1(segments_.insert(it0, boost::get<0>(spl.segs())));
	assert(cost_==cost_);
	cost_+= delta_cost(); // if roundoff errors, I should use real costs
	assert(cost_==cost_);
	splits_.pop();
	// I cannot push before pop
	splits_.push(Split(it0));
	splits_.push(Split(it1));
	++nbSegs_;
      }
      nb_segments_type nb_segs()const{ return nbSegs_;}
      cost_type cost()const{return cost_;};
      cost_type delta_cost()const{return splits_.top().cost();}
      
      template<typename Out>
      Out segments(Out out)const{return   std::copy(segments_.begin(), segments_.end(), out);}
    private:
      list_type segments_;
      std::priority_queue<Split> splits_;
      cost_type cost_;
      nb_segments_type nbSegs_;
    };
  }
}
#endif
