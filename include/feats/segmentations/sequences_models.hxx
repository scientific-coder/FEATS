#ifndef FEATS_SEGMENTATIONS_SEQUENCES_MODELS_HXX
#define FEATS_SEGMENTATIONS_SEQUENCES_MODELS_HXX

#include "feats/segmentations/models.hxx"

#include <utility>
#include <limits>
/***

 ***/


namespace feats{
  namespace segmentations{
    template<typename WeightIt, typename Val, typename Size=int> struct linear0_weighted_sequence_model:model_without_time
    {
      typedef Val value_type;
      typedef value_type data_type;
      typedef value_type cost_type;
      typedef Size size_type;

      typedef boost::tuple<value_type> parameters_type;


      typedef std::vector<std::pair<value_type,value_type> > container_type;
      typedef typename container_type::iterator iterator;
      typedef typename container_type::const_iterator const_iterator;


      typedef boost::tuple<parameters_type, cost_type> tuple_type; // first tuple for params, last element is cost


      linear0_weighted_sequence_model(WeightIt weightBeg, WeightIt weightEnd)
	: weightBegin_(weightBeg),sums_(std::distance(weightBeg, weightEnd),std::make_pair(0.,0.)),size_(0),cost_(0.){}

      linear0_weighted_sequence_model(const linear0_weighted_sequence_model& other)
	:weightBegin_(other.weightBegin_),sums_(other.sums_), size_(other.size_), cost_(other.cost_) {}
  

      linear0_weighted_sequence_model& operator=(const linear0_weighted_sequence_model& other){
	if(this != &other){
	  weightBegin_=other.weightBegin_;sums_=other.sums_;
	  size_=other.size_; cost_=other.cost_;
	}
	return *this;
      }

      template<typename ValIt>
      const linear0_weighted_sequence_model& operator+=( ValIt vi){
	WeightIt wIt(weightBegin_);
	cost_=0.; ++size_;
	for(iterator it(sums_.begin()); it!= sums_.end(); ++it, ++vi,++wIt){
	  it->first+=*vi;
	  it->second+=*vi * *vi;
	  cost_+=(it->second - (it->first * it->first)/size_) * *wIt;
	  //	  std::cerr<<"adding "<<*vi<<" cost is now "<<cost_<<'\n';  
	}
	return *this;
      }

      template<typename ValIt>
      const linear0_weighted_sequence_model& operator-=( ValIt vi){
	WeightIt wIt(weightBegin_);
	cost_=0.; --size_;
	for(iterator it(sums_.begin()); it!= sums_.end(); ++it, ++vi,++wIt){
	  it->first -= *vi;
	  it->second -= *vi * *vi;
	  cost_+=(it->second - (it->first * it->first)/size_) * *wIt;
	}
	return *this;
      }

      const linear0_weighted_sequence_model& operator+=( const linear0_weighted_sequence_model& lm){
	const_iterator otherIt(lm.sums_.begin());
	WeightIt wIt(weightBegin_);
	cost_=0.; size_+=lm.size_;
	for(iterator selfIt(sums_.begin()); selfIt!=sums_.end(); ++selfIt, ++otherIt, ++wIt){
	  selfIt->first += otherIt->first; selfIt->second += otherIt->second;
	  cost_+=(selfIt->second - (selfIt->first * selfIt->first)/size_) * *wIt;
	}
	return *this;
      }
      const linear0_weighted_sequence_model& operator-=( const linear0_weighted_sequence_model& lm){
	const_iterator otherIt(lm.sums_.begin());
	WeightIt wIt(weightBegin_);
	cost_=0.; size_-=lm.size_;
	for(iterator selfIt(sums_.begin()); selfIt!=sums_.end(); ++selfIt, ++otherIt, ++wIt){
	  selfIt->first -= otherIt->first; selfIt->second -= otherIt->second;
	  cost_+=(selfIt->second - (selfIt->first * selfIt->first)/size_) * *wIt;
	}
	return *this;
      }

      template<typename ValIt>
      const void transfer_from_model(ValIt vIt, linear0_weighted_sequence_model& lm){

	WeightIt wIt(weightBegin_);
	++size_; --lm.size_;
	cost_= lm.cost_= 0. ;
	for(iterator selfIt(sums_.begin()), otherIt(lm.sums_.begin()); selfIt!=sums_.end(); ++selfIt, ++otherIt, ++vIt, ++wIt){
	  const value_type& v(*vIt);
	  selfIt->first += v; otherIt->first -= v; selfIt->second += v*v; otherIt->second -= v*v;
	  cost_+=(selfIt->second - (selfIt->first * selfIt->first)/size_) * *wIt;
	  lm.cost_+=(otherIt->second - (otherIt->first * otherIt->first)/size_) * *wIt;
	}
      }
      const void transfer_from_model(const linear0_weighted_sequence_model& t, linear0_weighted_sequence_model& lm){(*this)+=t; lm-=t;}
      
      parameters_type operator()()const {return  boost::make_tuple(-555.555);}// FIXME could be an iterator computing the means (boost iterator adaptor divider by size_ over sums_ .first but it's not value type then (it should not anyway)

      value_type value( value_type unused) const{return boost::get<0>(operator()());}
      cost_type cost() const { 	return size_?cost_:0.;}
      
      size_type size()const {return size_;}
      tuple_type tuple()const {return boost::make_tuple(operator()(), cost());}
    private:
      WeightIt weightBegin_;
      container_type sums_;
      size_type size_;
      value_type cost_;
    };

    template< typename Val, typename Size=int> struct linear0_sequence_model:model_without_time
    {
      typedef Val value_type;
      typedef value_type data_type;
      typedef value_type cost_type;
      typedef Size size_type;


      typedef boost::tuple<value_type, value_type, Size> params_type;
      typedef std::vector<params_type > container_type;
      typedef typename container_type::iterator iterator;
      typedef typename container_type::const_iterator const_iterator;


      struct value_from_params: std::unary_function<params_type,value_type>{
	value_type operator()(const params_type&p)const{
	  return boost::get<2>(p)
	    ?boost::get<0>(p)/boost::get<2>(p)
	    :-std::numeric_limits<value_type>::max();//;std::numeric_limits<value_type>::quietNaN();
	}
      };
      typedef boost::transform_iterator<value_from_params, const_iterator> values_iterator;

      typedef boost::tuple<values_iterator> parameters_type;



      typedef boost::tuple<parameters_type, cost_type> tuple_type; // first tuple for params, last element is cost


      linear0_sequence_model(Size seqSize)
	:params_(seqSize,params_type(0.,0.,0)),size_(0),cost_(0.){}

      linear0_sequence_model(const linear0_sequence_model& other)
	:params_(other.params_),size_(other.size_),cost_(other.cost_) {}
  

      linear0_sequence_model& operator=(const linear0_sequence_model& other){
	if(this != &other){ params_=other.params_; cost_=other.cost_; size_= other.size_; }
	return *this;
      }

      template<typename ValIt>
      const linear0_sequence_model& operator+=( ValIt vi){
	cost_=0.; ++size_;
	for(iterator it(params_.begin()); it!= params_.end(); ++it, ++vi){
	  const value_type& v(*vi);
	  params_type& p(*it);
	  if(v==v){
	    boost::get<0>(p)+=v;
	    boost::get<1>(p)+=v*v;
	    ++boost::get<2>(p);
	  }
	  if(boost::get<2>(p)){
	    cost_+=boost::get<1>(p) 
	      - (boost::get<0>(p)* boost::get<0>(p))/boost::get<2>(p);
	  }
	  //	  std::cerr<<"adding "<<*vi<<" cost is now "<<cost_<<'\n';  
	}
	return *this;
      }

      template<typename ValIt>
      const linear0_sequence_model& operator-=( ValIt vi){
	cost_=0.; --size_;
	for(iterator it(params_.begin()); it!= params_.end(); ++it, ++vi){
	  const value_type& v(*vi);
	  params_type& p(*it);
	  if(v==v){
	    boost::get<0>(p)-=v;
	    boost::get<1>(p)-=v*v;
	    --boost::get<2>(p);
	  }
	  if(boost::get<2>(p)){
	    cost_+=boost::get<1>(p) 
	      - (boost::get<0>(p)* boost::get<0>(p))/boost::get<2>(p);
	  }
	}
	return *this;
      }

      const linear0_sequence_model& operator+=( const linear0_sequence_model& lm){
	const_iterator otherIt(lm.params_.begin());
	cost_=0.; size_+=lm.size_;
	for(iterator selfIt(params_.begin()); selfIt!=params_.end(); ++selfIt, ++otherIt){
	  params_type& p(*selfIt);
	  boost::get<0>(p) += boost::get<0>(*otherIt);
	  boost::get<1>(p) += boost::get<1>(*otherIt);
	  boost::get<2>(p) += boost::get<2>(*otherIt);
	  if(boost::get<2>(p)){
	    cost_+=boost::get<1>(p) 
	      - (boost::get<0>(p)* boost::get<0>(p))/boost::get<2>(p);
	  }
	}
	return *this;
      }
      const linear0_sequence_model& operator-=( const linear0_sequence_model& lm){
	const_iterator otherIt(lm.params_.begin());
	cost_=0.; size_-=lm.size_;
	for(iterator selfIt(params_.begin()); selfIt!=params_.end(); ++selfIt, ++otherIt){
	  params_type& p(*selfIt);
	  boost::get<0>(p) -= boost::get<0>(*otherIt);
	  boost::get<1>(p) -= boost::get<1>(*otherIt);
	  boost::get<2>(p) -= boost::get<2>(*otherIt);
	  if(boost::get<2>(p)){
	    cost_+=boost::get<1>(p) 
	      - (boost::get<0>(p)* boost::get<0>(p))/boost::get<2>(p);
	  }
	}
	return *this;
      }

      template<typename ValIt>
      const void transfer_from_model(ValIt vIt, linear0_sequence_model& lm){
	cost_= lm.cost_= 0. ; ++size_; --lm.size_;
	for(iterator selfIt(params_.begin()), otherIt(lm.params_.begin()); selfIt!=params_.end(); ++selfIt, ++otherIt, ++vIt){
	  const value_type& v(*vIt);
	  params_type& sp(*selfIt);
	  params_type& op(*otherIt);
	  if(v==v){
	    boost::get<0>(sp)+=v;
	    boost::get<1>(sp)+=v*v;
	    ++boost::get<2>(sp);
	    boost::get<0>(op)-=v;
	    boost::get<1>(op)-=v*v;
	    --boost::get<2>(op);
	  }
	  if(boost::get<2>(sp)){
	    cost_+=boost::get<1>(sp) 
	      - (boost::get<0>(sp)* boost::get<0>(sp))/boost::get<2>(sp);
	  }
	  if(boost::get<2>(op)){
	    lm.cost_+=boost::get<1>(op) 
	      - (boost::get<0>(op)* boost::get<0>(op))/boost::get<2>(op);
	  }
	}
      }
      const void transfer_from_model(const linear0_sequence_model& t, linear0_sequence_model& lm){(*this)+=t; lm-=t;}
      
      parameters_type operator()()const 
      {return  boost::make_tuple(values_iterator(params_.begin()));}

      values_iterator value( value_type unused) const
	  {return values_iterator(params_.begin());}
      cost_type cost() const { 	return cost_;}

      template<typename Out> Out means(Out o)const{return std::copy(values_iterator(params_.begin()),values_iterator(params_.end()), o);} 
      
      size_type size()const {return size_;}
      size_type sequence_size()const{return params_.size();}
      tuple_type tuple()const {return boost::make_tuple(operator()(), cost());}
    private:
      container_type params_;
      size_type size_;
      value_type cost_;
    };
  }
}

#endif
