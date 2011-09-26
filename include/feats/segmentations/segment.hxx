#ifndef FEATS_SEGMENTATIONS_SEGMENT_HXX
#define FEATS_SEGMENTATIONS_SEGMENT_HXX

#include <algorithm>
#include <functional>
#include <boost/tuple/tuple.hpp>
#include <boost/bind.hpp>

namespace feats{
  namespace segmentations{


    template<typename Model, typename Iter> struct segment {
      // should static assert for Model::value_type and Iter::value_type, Model::size_type and Iter::difference_type It could be better to have a more precise type for Model::value_type than Iter::value_type
      typedef Iter iterator;
      typedef Model model_type;
      typedef typename model_type::cost_type cost_type;
      typedef typename model_type::data_type data_type;
      typedef boost::tuple<model_type, boost::tuple<iterator, iterator> > tuple_type;

      segment(iterator begEnd, const model_type& initM=model_type()): model_(initM), begin_(begEnd), end_(begEnd){}


      // althought not stricly necessary a constructor for seg and model taking (beging end) could be (marginally?) faster by not updating size and end_for each point
      segment(iterator beg, iterator end, const model_type&initM=model_type()): model_(initM), begin_(beg), end_(end){
	while(beg!=end){ model_+=*beg; ++beg;}
      }

      segment& operator=(const segment& other){
	if(this != &other){
	  model_=other.model_;
	  begin_=other.begin_;
	  end_=other.end_;
	}
	return *this;
      }
      // even better compiler might be better at unrolling
      void shrink_begin(){ model_-=*begin_; ++begin_;}
      void grow_begin(){ ++begin_; model_ += *begin_;}
      void grow_end() {model_+=*end_; ++end_;}
      void shrink_end() {--end_; model_-=*end_;}

      void merge_begin( const segment& other){
	// should assert other.end_==this->begin_
	model_+=other.model_;
	begin_=other.begin_;
      }
      void merge_end( const segment& other){
	// should assert other.begin_==this->end_
	model_+=other.model_;
	end_=other.end_;
      }
  
      void optimize_with( segment& other){
	typename model_type::cost_type bestCost(model_.cost()+other.model_.cost());
	iterator testingEnd(end_);
	typename boost::remove_reference<model_type>::type thisTestingModel(model_), otherTestingModel(other.model_);
	if(thisTestingModel.size()!=0) {otherTestingModel+=thisTestingModel; thisTestingModel=empty_model(thisTestingModel); testingEnd=begin_;}
	while( otherTestingModel.size()!=0){
	  thisTestingModel.transfer_from_model(*testingEnd, otherTestingModel);
	  ++testingEnd;
	  if(thisTestingModel.cost()>bestCost)break;
	  if(thisTestingModel.cost()+otherTestingModel.cost()<bestCost)
	    {bestCost=thisTestingModel.cost()+otherTestingModel.cost(); model_=thisTestingModel; other.model_=otherTestingModel; end_=testingEnd;}
	}
	other.begin_=end_;
	//	std::cerr<<"end of optimizeWith bestCost="<<bestCost<<'\n';
      }
      
      template<typename CostsIt> typename model_type::cost_type optimize_with_costs(CostsIt beg, CostsIt end) {
	typename model_type::cost_type bestCost(model_.cost()+*beg);
	model_type testingModel(model_);
	//	std::cerr<<"entering optimize with costs"<<testingModel.cost()<<"or "<<model_.cost()<<"+"<<*beg<<"="<<bestCost<<"\n";
	if(beg!=end)++beg;// will be adding one point 
	for( iterator testingEnd(end_); beg != end; ++beg){
	  testingModel+=*testingEnd;
	  ++testingEnd;
	  cost_type tmpCost(testingModel.cost());
	  if(tmpCost>bestCost)break;
	  tmpCost+=*beg;
	  if(tmpCost <bestCost)
	    { bestCost=tmpCost; model_=testingModel; end_=testingEnd; }
	}
	//	std::cerr<<"out of optimize with cost best:"<<bestCost<<'\n';
	return bestCost;
      }
	
      iterator begin() const { return begin_;}
      iterator end() const { return end_;}

      // this is really inefficient (esp for Linear1Model) If I want to be fast, I will have to add a Model method to return a Modeller or something that would cahe the parameters. Maybe my Models should really be ModelBuilders
      template<typename Out> Out values(Out out)const
      {	return std::transform(begin_, end_, out, std::bind1st(std::mem_fun_ref(&model_type::value), model_));}

      // must export model instead to allow segmentations of model series
      // size (end-begin)will then be different from model.size (nb of points)

      const model_type& model()const{return model_;}
      const typename std::iterator_traits<iterator>::difference_type size()const{return std::distance(begin_,end_);}
      tuple_type tuple()const{return boost::make_tuple(model_, boost::make_tuple(begin_,end_));}

    private:
      model_type model_;
      iterator begin_, end_;
    };

    template<typename Model, typename Iter> struct models_segment {
      // should static assert for Model::value_type and Iter::value_type, Model::size_type and Iter::difference_type It could be better to have a more precise type for Model::value_type than Iter::value_type
      typedef Iter iterator;
      typedef typename Model::value_type model_type;
      typedef typename model_type::cost_type cost_type;
      typedef typename model_type::data_type data_type;
      typedef boost::tuple<model_type, boost::tuple<iterator, iterator> > tuple_type;

      template<typename M, typename I> struct rebind{ typedef models_segment<M,I> other; };
      // trouble with model not beeing model for moedls_segment : band aid solution :
      template<typename I> struct rebind_iter{ typedef models_segment<Model,I> other; };

      models_segment(iterator begEnd, const model_type& initM=model_type()): model_(initM), begin_(begEnd), end_(begEnd){}


      // althought not stricly necessary a constructor for seg and model taking (beging end) could be (marginally?) faster by not updating size and end_for each point
      models_segment(iterator beg, iterator end, const model_type&initM=model_type()): model_(initM), begin_(beg), end_(end){
	while(beg!=end){ model_+=*beg; ++beg;}
      }

      models_segment& operator=(const models_segment& other){
	if(this != &other){
	  model_=other.model_;
	  begin_=other.begin_;
	  end_=other.end_;
	}
	return *this;
      }
      // even better compiler might be better at unrolling
      void shrink_begin(){ model_-=*begin_; ++begin_;}
      void grow_begin(){ ++begin_; model_ += *begin_;}
      void grow_end() {model_+=*end_; ++end_;}
      void shrink_end() {--end_; model_-=*end_;}

      void merge_begin( const models_segment& other){
	// should assert other.end_==this->begin_
	model_+=other.model_;
	begin_=other.begin_;
      }
  
      void optimize_with( models_segment& other){
	typename model_type::cost_type bestCost(model_.cost()+other.model_.cost());
	iterator testingEnd(end_);
	typename boost::remove_reference<model_type>::type thisTestingModel(model_), otherTestingModel(other.model_);
	if(thisTestingModel.size()!=0) {otherTestingModel+=thisTestingModel; thisTestingModel=empty_model(thisTestingModel); testingEnd=begin_;}
	while( otherTestingModel.size()!=0){
	  thisTestingModel.transfer_from_model(*testingEnd, otherTestingModel);
	  ++testingEnd;
	  if(thisTestingModel.cost()>bestCost)break;
	  if(thisTestingModel.cost()+otherTestingModel.cost()<bestCost)
	    {bestCost=thisTestingModel.cost()+otherTestingModel.cost(); model_=thisTestingModel; other.model_=otherTestingModel; end_=testingEnd;}
	}
	other.begin_=end_;
	//	std::cerr<<"end of optimizeWith bestCost="<<bestCost<<'\n';
      }
      
      template<typename CostsIt> typename model_type::cost_type optimize_with_costs(CostsIt beg, CostsIt end) {
	typename model_type::cost_type bestCost(model_.cost()+*beg);
	model_type testingModel(model_);
	if(beg!=end)++beg;// will be adding one point 
	for( iterator testingEnd(end_); beg != end; ++beg){
	  testingModel+=*testingEnd;
	  ++testingEnd;
	  cost_type tmpCost(testingModel.cost());
	  if(tmpCost>bestCost)break;
	  tmpCost+=*beg;
	  if(tmpCost <bestCost)
	    { bestCost=tmpCost; model_=testingModel; end_=testingEnd; }
	}
	return bestCost;
      }
	
      iterator begin() const { return begin_;}
      iterator end() const { return end_;}

      // this is really inefficient (esp for Linear1Model) If I want to be fast, I will have to add a Model method to return a Modeller or something that would cahe the parameters. Maybe my Models should really be ModelBuilders
      // A DEBUGGER NE COMPILE PAS !!!!!
      template<typename Out> Out values(Out out)const
      {	return std::transform(begin_, end_, out, std::bind1st(std::mem_fun_ref(&model_type::value), model_));}

      // must export model instead to allow segmentations of model series
      // size (end-begin)will then be different from model.size (nb of points)

      const model_type& model()const{return model_;}
      const typename std::iterator_traits<iterator>::difference_type size()const{return std::distance(begin_,end_);}
      tuple_type tuple()const{return boost::make_tuple(model_, boost::make_tuple(begin_,end_));}

    private:
      model_type model_;
      iterator begin_, end_;
    };



    /***
	Segmentations algorithms use two basic block : basic segment creation and optimisation.
	for basic segment creation and segmentation, we have to decide weither nb of segments is a static or dynamic property
	better to allow both of course.
	basic segment creation has to allow statioc nb of segs and segmentation has to allow dynamic nb of segs.
	I intend to use something like Boost Fusion to allow both. The segs sequence will be an object : either a tuple of values for static size of a tuple of iterators for dynamic one. Maybe I should revert to pair for range to ease sstatic dispatch (overload resolution), or always use a tuple of (even one) tuple.
	To allow static computing, I also need to have a template meta function ..., will do this later.
    
    ***/

    template<typename Seg> boost::tuple<Seg,Seg> split_in_two(const Seg& seg){ 
      //      typename Seg::model_type m=1;
      return boost::make_tuple(Seg(seg.begin(),empty_model(seg.model())), seg);}

    /* I do not use std::accumulate because I want to use += instead of + for efficiency reasons
     */
    /*
    template<typename It> typename std::iterator_traits<typename boost::tuples::element<0,It>::type>::value_type
    merge_segments(boost::tuple<It,It> beginEnd){
      typedef typename std::iterator_traits<typename boost::tuples::element<0,It>::type>::value_type segment_type;
      segment_type  res(*boost::get<0>(beginEnd));
      std::for_each(boost::next(boost::get<0>(beginEnd)), boost::get<1>(beginEnd)
		    ,boost::bind(segment_type::merge_end, boost::ref(res), _1));
      return res;
    }
    */
    
    template<typename Seg> struct seg_to_tuple {
      typedef const Seg& argument_type;
      typedef typename Seg::tuple_type seg_tuple_type;
      typedef typename Seg::iterator iterator;
      typedef typename Seg::model_type model_type;
  
      typedef typename std::iterator_traits<iterator>::difference_type index_type;
      typedef boost::tuple<typename model_type::tuple_type,typename boost::tuple<index_type, index_type> > result_type;

      seg_to_tuple(iterator base):base_(base){}

      result_type operator()(const Seg& seg)const{
	const seg_tuple_type & tmp(seg.tuple());

	return boost::make_tuple(boost::get<0>(tmp).tuple(),boost::make_tuple( boost::get<0>(boost::get<1>(tmp))-base_
									       ,boost::get<1>(boost::get<1>(tmp))-base_));
      }

    private:
      iterator base_;
    };

  }
}
#endif
