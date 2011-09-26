#ifndef FEATS_SEGMENTATIONS_MERGES_HXX
#define FEATS_SEGMENTATIONS_MERGES_HXX
#include <functional>
#include <algorithm>

#define DEBUG 1

#if DEBUG
#include <iostream>
#endif
#include <boost/utility.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "feats/segmentations/bottom_up_segmentation.hxx"

#include "feats/utility/leftish_tree.hxx"



namespace feats{
  namespace segmentations{

    /*
      This is more complex than binary_split because when merging two segments, I have to remove them and change the Merge of the previous segment
      The problem is that I cannot remove elements from priority queues, so i made a leftish_tree as described in The Art of Computer Programming vol.3
    */

    template<typename Seg, typename Nb=int> struct merges: bottom_up_segmentation<Seg, Nb> {
      typedef Seg segment_type;
      typedef Nb nb_segments_type;
      typedef typename segment_type::cost_type cost_type;

      template<typename S, typename N=nb_segments_type> struct rebind{typedef merges<S,N> other;};
    private:
      /***
	  Generic utility dereference functor, should be in some library.
      ***/

      template <typename F>
      class dereference : std::binary_function<typename F::first_argument_type,typename F::second_argument_type,typename F::result_type>
      {
	typedef std::binary_function<typename F::first_argument_type,typename F::second_argument_type,typename F::result_type> parent;
	F m_op;
      public:
	typename parent::result_type operator()(typename parent::first_argument_type* a1, typename parent::second_argument_type* a2)
 	{ return m_op( *a1, *a2);}
      };

      /// forward declaration needed because the Merge need to handle nodes of leftish_tree of Merge pointers
      struct Merge;

      typedef feats::utility::leftish_tree<Merge*, dereference<std::less<Merge> > >heap_type;
      typedef typename heap_type::node_type node_type;

      /***
	  A Merge computes the merging of a segment (original)with the next one (merged).
	  I should really try to avoid storing both merged segments, but it really helps to have a clean implmentation...
	  I store the cost because it is both small and often used. Maybe I should also store result_segment...
      ***/
      struct Merge {
	typedef typename Seg::model_type::value_type value_type;
	typedef typename heap_type::node_type node_type;
	// should use boost fastest argument passing type call_traits
	explicit Merge( const Seg& seg0, const Seg& seg1, node_type prev  ) : seg0_(seg0),seg1_(seg1)
									    , prev_node_(prev), next_node_(0)
	{ cost_=result_segment().model().cost()-seg0_.model().cost()-seg1_.model().cost(); }
	Seg result_segment()const{Seg result_seg(seg1_); result_seg.merge_begin(seg0_); return result_seg;}
	Seg get_original_segment()const {return seg0_;}

	void update_from_next(Merge const*const m)
	{seg1_=m->result_segment(); update_cost(); next_node_=m->next_node_;}
	void update_from_prev(Merge const * const m)
	{seg0_=m->result_segment(); update_cost(); prev_node_=m->prev_node_;}

	void set_next_node(node_type n){next_node_=n;}
	node_type get_next_node()const{return next_node_;}
	node_type get_prev_node()const{return prev_node_;}

	value_type cost()const{return cost_;}

	bool operator<(const Merge& other)const{ return cost()>other.cost();}

      private:
	void update_cost(){ Seg result_seg(seg1_); result_seg.merge_begin(seg0_); cost_=result_seg.model().cost()-seg0_.model().cost()-seg1_.model().cost(); }
#if DEBUG
      public:
#endif
	Seg seg0_,seg1_;
	value_type cost_;
	node_type  prev_node_,next_node_;
      };


    public:
      template <typename It> merges( boost::tuple<It,It> initSegBeginEnd):lastSeg_(*boost::get<0>(initSegBeginEnd)), firstMerge_(0), cost_(0.), nbSegs_(0)
#if DEBUG
								      ,dbg_(lastSeg_.begin())
#endif
      {

	--boost::get<1>(initSegBeginEnd); // no merging for the last seg
	node_type prevNode(0);
	for(It current(boost::get<0>(initSegBeginEnd)); boost::get<0>(initSegBeginEnd)++ != boost::get<1>(initSegBeginEnd)
	      ; current=boost::get<0>(initSegBeginEnd), ++nbSegs_){
	  node_type currentNode(merges_.push(new Merge(*current, *boost::get<0>(initSegBeginEnd), prevNode)));
	  if(prevNode)(**prevNode)->set_next_node(currentNode);
	  if(!firstMerge_)firstMerge_=**currentNode;
	  prevNode=currentNode;
	  cost_+=current->model().cost();
	}
	lastSeg_=*boost::prior(boost::get<0>(initSegBeginEnd)); cost_+= lastSeg_.model().cost();
	++nbSegs_;
      }

#if DEBUG
      void dbg(const Merge&m)const{	std::cerr<<"merge of"<<dbg_(m.get_original_segment())<<" with "<<dbg_(m.seg1_)<<" result is "<<dbg_(m.result_segment())<<'\n';}
#endif

      void dec_segments() {
	if(!merges_.empty()){
	Merge* current(merges_.top());
	cost_+= delta_cost(); // maybe should use the segments costs to avoid roundoff errors
	merges_.pop();

	node_type toUpdate(current->get_prev_node());
	if(toUpdate){ merges_.update(toUpdate, boost::bind(&Merge::update_from_next,_1, current));}
	else{firstMerge_=(current->get_next_node())?**(current->get_next_node()):0;}

	toUpdate=current->get_next_node();
	if(toUpdate){ merges_.update(toUpdate, boost::bind(&Merge::update_from_prev, _1, current)); }
	else{ lastSeg_=current->result_segment(); }
	delete current;
	}
	--nbSegs_;
      }
  
      template<typename Out>
      Out segments(Out out)const{
	if(!merges_.empty()){
	  Merge* current(firstMerge_);
	  while(current){
	    *out=current->get_original_segment(); ++out;
	    current=(current->get_next_node())?**(current->get_next_node()):0;
	  }
	}
	*out=lastSeg_;
	return ++out;
      }
      nb_segments_type nb_segs()const{ return nbSegs_;}
      cost_type cost()const{ return cost_;}
      cost_type delta_cost()const{ return merges_.top()->cost();}

      ~merges(){ while(!merges_.empty()){ delete(merges_.top()); merges_.pop();}}

    private:
      heap_type merges_;
      Seg lastSeg_;
      Merge* firstMerge_;
      cost_type cost_;
      nb_segments_type nbSegs_;
#if DEBUG
      seg_to_tuple<Seg> dbg_;
#endif
    };

  }
}

#endif
