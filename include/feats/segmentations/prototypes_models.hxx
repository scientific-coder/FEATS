#ifndef FEATS_SEGMENTATIONS_PROTOTYPES_MODELS_HXX
#define FEATS_SEGMENTATIONS_PROTOTYPES_MODELS_HXX

#include "feats/segmentations/models.hxx"

/***
    For symbolic representations based on segmentation, I use prototype parameters. That means that some segmentations will be done with restricted prototype values for some paramters. I need:

1) Models to allows those segmentations for a given set of prototype parameters
This is done with linear_prototype_level_model and linear1_prototype_slope_model

2) Quantizers to compute a set of prototype parameters given a segmentation
This is done by segmenting models stored according to the paramter to be quantized. Note that I use models to set the boundaries, but paramters are computed "as if" from original data, and not on model parameters. This, however, prevetns me from quantizing directly *prototype*_model becaus they do not store the required partial sums. I use an adapter the same way I do to create initial N-segmentations for bottom_up (with a pair of shared_container_iterator )

3) Optimizers to optimize both prototypes values and segmentations
The idea is as follows:
1) do segmentation on free (unrestricted parameters) model
2) Quantize a parameters of the previous segmentation
3) do a restricted segmentation with prototype values for parameter
4) loop to 2 until convergence
 ***/


namespace feats{
  namespace segmentations{

    /***
	Prototypes models : one (for now) parameter is restricted to a limited set of prototypes values.
	For linear0_prototypes_model it is the level (of course)
	For linear1_prototypes_model it is the slope
	I do not copy the prototypes into each model, to spare memory,  but changing prototype values invalidates the models
	ProtoIt *must* be RandomAccessIterator or the models would be slow
    ***/
    /***
	I am using traits to stor relationship between prototype_model and other related types:
    ***/
    struct undefined_prototypes_model_error{};
    template<template<typename,typename,typename> class ProtoModelT, typename Data, typename Size=int> struct prototypes_model_traits{
      typedef undefined_prototypes_model_error free_model;
      typedef undefined_prototypes_model_error meta_model;
    };
      



    template<typename ProtoIt, typename ValYX, typename Size=int> struct linear1_prototype_slope_model;
    template<typename ValYX, typename Size=int> struct linear1_prototype_slope_meta_model;


    template<typename ProtoIt, typename Val=typename std::iterator_traits<ProtoIt>::value_type, typename Size=int> struct linear0_prototype_level_model:model_without_time
    {
      typedef Val value_type;
      typedef value_type data_type;
      typedef value_type cost_type;
      typedef Size size_type;
      typedef ProtoIt prototypes_iterator;
      typedef boost::tuple<value_type> parameters_type;

      typedef std::vector<cost_type> container_type;
      typedef typename container_type::iterator iterator;
      typedef typename container_type::const_iterator const_iterator;

      typedef boost::tuple<parameters_type, cost_type> tuple_type; // first tuple for params, last element is cost

      // L0Model can be either linear0_model or linear1_prototype_level_model
      // this template functor is used to get the parameter (here slope) to quantize for prototypes
      template<typename L0Model> struct parameter_selector:std::unary_function<L0Model, value_type>  {	value_type operator()(const L0Model& lm)const{return boost::get<0>(lm());} };





      linear0_prototype_level_model(ProtoIt protoLvlBeg, ProtoIt protoLvlEnd)
	: protoBegin_(protoLvlBeg), costs_(std::distance(protoBegin_, protoLvlEnd),0), bestCostIt_(costs_.begin()),size_(0) {}

      linear0_prototype_level_model(const linear0_prototype_level_model& other)
	:protoBegin_(other.protoBegin_),costs_(other.costs_)
	, bestCostIt_(costs_.begin()+(other.bestCostIt_-other.costs_.begin())), size_(other.size_) {}
  

      linear0_prototype_level_model& operator=(const linear0_prototype_level_model& other){
	if(this != &other){
	  protoBegin_=other.protoBegin_; costs_=other.costs_;
	  bestCostIt_=costs_.begin()+(other.bestCostIt_-other.costs_.begin()); size_=other.size_;  
	}
	return *this;
      }

      const linear0_prototype_level_model& operator+=( const Val& v){
	if(v==v){
	ProtoIt pIt(protoBegin_);
	bestCostIt_=costs_.begin();
	for(iterator it(costs_.begin()); it!= costs_.end(); ++it, ++pIt){
	  const Val tmp(*pIt-v); 
	  if((*it+= tmp*tmp)< *bestCostIt_){ bestCostIt_=it;}
	}
	++size_; 
	}
	return *this;
      }
      const linear0_prototype_level_model& operator-=( const Val& v){
	if(v==v){
	  ProtoIt pIt(protoBegin_);
	bestCostIt_=costs_.begin();
	for(iterator it(costs_.begin()); it!= costs_.end(); ++it, ++pIt){
	  const Val tmp(*pIt-v); 
	  if((*it-= tmp*tmp)< *bestCostIt_){ bestCostIt_=it; }
	}
	--size_; 
	}
	return *this;
      }

      const linear0_prototype_level_model& operator+=( const linear0_prototype_level_model& lm){
	// could check that we have same prototypes in same order
	const_iterator otherIt(lm.costs_.begin());	bestCostIt_=costs_.begin();
	for(iterator selfIt(costs_.begin()); selfIt!=costs_.end(); ++selfIt, ++otherIt){
	  if((*selfIt+= *otherIt)<*bestCostIt_){ bestCostIt_=selfIt;}
	}
	size_+=lm.size_; return *this;
      }
      const linear0_prototype_level_model& operator-=( const linear0_prototype_level_model& lm){
	// could check that we have same prototypes in same order
	const_iterator otherIt(lm.costs_.begin());	bestCostIt_=costs_.begin();
	for(iterator selfIt(costs_.begin()); selfIt!=costs_.end(); ++selfIt, ++otherIt){
	  if((*selfIt-= *otherIt)<*bestCostIt_){ bestCostIt_=selfIt;} 
	}
	size_-=lm.size_; return *this;
      }

      const void transfer_from_model(Val v, linear0_prototype_level_model& lm){
	// could check that we have same prototypes in same order
	if(v==v){
	ProtoIt pIt(protoBegin_);
	bestCostIt_=costs_.begin();
	for(iterator selfIt(costs_.begin()), otherIt(lm.costs_.begin()); selfIt!=costs_.end(); ++selfIt, ++otherIt,++pIt){
	  Val tmp(*pIt-v); tmp *=tmp;
	  if((*selfIt+= tmp)<*bestCostIt_){ bestCostIt_=selfIt;}
	  if((*otherIt-= tmp)<*(lm.bestCostIt_)){ lm.bestCostIt_=otherIt;}
	}
	++size_; --lm.size_;
	}
      }
      const void transfer_from_model(const linear0_prototype_level_model& t, linear0_prototype_level_model& lm){(*this)+=t; lm-=t;}
      
      parameters_type operator()()const {return  boost::make_tuple(*(protoBegin_+(bestCostIt_-costs_.begin())));}

      size_type get_proto_nb()const{ return bestCostIt_-costs_.begin();}

      prototypes_iterator get_prototypes_begin()const{return protoBegin_;};

      prototypes_iterator get_prototypes_end()const{return protoBegin_+costs_.size();};

      value_type value( value_type unused) const{return boost::get<0>(operator()());}
      cost_type cost() const { 	return *bestCostIt_;}

      template<typename OutCosts>
      OutCosts costs( OutCosts o)const{return std::copy(costs_.begin(),costs_.end(),o);}
      
      size_type size()const {return size_;}
      tuple_type tuple()const {return boost::make_tuple(operator()(), cost());}
    private:
      prototypes_iterator protoBegin_;
      container_type costs_;
      iterator bestCostIt_;
      size_type size_;
    };

      // linear0_model could be used as the meta modelbecause linear0 only has one parameter 
    template<typename Data, typename Size> struct linear0_prototype_level_meta_model : model_without_time{
      typedef linear0_model<Data, Size> data_type;
      typedef typename data_type::value_type value_type;
      typedef typename data_type::cost_type cost_type;
      typedef typename data_type::data_type free_model_data_type;
      typedef Size size_type;
      typedef boost::tuple<value_type> parameters_type; // the slope is the only parameter
      typedef boost::tuple<parameters_type,value_type> tuple_type; // first tuple for params, last element is cost

      template<typename L0Model> struct parameter_selector:std::unary_function<L0Model, Data>
      {	Data operator()(const L0Model& lm)const{return boost::get<0>(lm());} };


      linear0_prototype_level_meta_model():sumY_(),sumYY_(),size_(0){}

      const linear0_prototype_level_meta_model& operator+=( const data_type& v){
	sumY_+=v.sum_y(); sumYY_+=v.sum_yy(); size_+=v.size();
	return *this;
      }
      const linear0_prototype_level_meta_model& operator-=( const data_type& v){
	sumY_-=v.sum_y(); sumYY_-=v.sum_yy(); size_-=v.size();
	return *this;
      }
      const linear0_prototype_level_meta_model& operator+=( const linear0_prototype_level_meta_model& lm)
      { sumY_+=lm.sumY_; sumYY_+=lm.sumYY_; size_+=lm.size_; return *this;}
      const linear0_prototype_level_meta_model& operator-=( const linear0_prototype_level_meta_model& lm)
      { sumY_-=lm.sumY_; sumYY_-=lm.sumYY_; size_-=lm.size_; return *this;}

      const void transfer_from_model(const data_type& v, linear0_prototype_level_meta_model& lm){sumY_+=v.sum_y(); lm.sumY_-=v.sum_y(); sumYY_+=v.sum_yy(); lm.sumYY_-=v.sum_yy(); size_+=v.size(); lm.size_-=v.size();}
      const void transfer_from_model(const linear0_prototype_level_meta_model& t, linear0_prototype_level_meta_model& lm){(*this)+=t; lm-=t;}

      parameters_type operator()()const {return  boost::make_tuple((size_==0)? value_type() : sumY_/size_ );}
      free_model_data_type value( value_type unused) const{return boost::get<0>(operator()());}
      cost_type cost() const { return (size_==0)? value_type() :  sumYY_ - (sumY_*sumY_)/size_ ;}

      size_type size()const {return size_;}
      tuple_type tuple()const{return boost::make_tuple(operator()(),cost());}
    private:
      value_type sumY_, sumYY_;
      size_type size_;



    };

    template< typename Data,typename Size> struct prototypes_model_traits<linear0_prototype_level_model, Data,Size>{
      typedef linear0_model<Data, Size> free_model;
      typedef linear0_prototype_level_meta_model<Data, Size> meta_model;
    };
    /***
	for this prototype model there are trade-offs because I can precompute partial costs for he known slopes. I do not do it because I prefer to save memory: memory bus should be the bottleneck for huge time-series (those for which speed is of essence).
	This model is very much like linear1_model, but slopes are known, I only have to compute offset but I have to do it for each slope
	Computing cost is expensive, so I cache the best cost.
     ***/

    template<typename ProtoIt, typename ValYX, typename Size> struct linear1_prototype_slope_model : model_with_time
    {
      // should differencite data_type taken and returnd for constness ...
      typedef ValYX data_type_ref;
      typedef typename boost::remove_const<typename boost::remove_reference<typename boost::tuples::element<0, data_type_ref>::type>::type>::type value_type;
      typedef  typename boost::remove_const<typename boost::remove_reference<typename boost::tuples::element<1, data_type_ref>::type>::type>::type time_type;
      typedef boost::tuple<value_type, time_type> data_type;
      typedef value_type cost_type;

      typedef ProtoIt prototypes_iterator;

      typedef Size size_type;
      typedef boost::tuple<value_type, value_type> parameters_type; // tricky to define :-(
      typedef boost::tuple<parameters_type,value_type> tuple_type; // first tuple for params, last element is cost

      // L1Model can be either linear1_model or linear1_prototype_slope_model
      // this template functor is used to get the parameter (here slope) to quantize for prototypes
      template<typename L1Model> struct parameter_selector:std::unary_function<L1Model, value_type>
      {	value_type operator()(const L1Model& lm)const{return boost::get<1>(lm());} };


      linear1_prototype_slope_model(ProtoIt protoSlopesBeg, ProtoIt protoSlopesEnd)
	: sumX_(0.),sumXX_(0.),sumXY_(0.),sumYY_(0.),size_(0)
	, protoSlopesBegin_(protoSlopesBeg), protoSlopesEnd_(protoSlopesEnd), bestSlopeIt_(protoSlopesBeg), bestCost_(0.){}

      /*
      linear1_prototype_slope_model(const linear1_prototype_slope& other):sumX_(other.sumX_), sumXX_(other.sumXX_), sumXY_(other.sumXY_), sumY_(other.sumY_), sumYY_(other.sumYY_), size_(other.size_), protoSlopesBegin_(other.protoSlopesBegin_), protoSlopesEnd_(other.protoSlopesEnd_), bestSlopeIt(
      */
      const linear1_prototype_slope_model& operator+=( const data_type& t){
	const value_type& y(boost::get<0>(t)); 
	if(y==y){
	  const time_type& x(boost::get<1>(t)); sumX_+=x; sumXX_+=x*x; sumY_+=y; sumYY_+=y*y; sumXY_+=x*y; ++size_; 
	update();
	}
	return *this;
      }

      const linear1_prototype_slope_model& operator-=( const data_type& t){ 
	const value_type& y(boost::get<0>(t)); 
	if(y==y){
	  const time_type& x(boost::get<1>(t)); sumX_-=x; sumXX_-=x*x; sumY_-=y; sumYY_-=y*y; sumXY_-=x*y; --size_;
	  update();
	}
	return *this;
      }
      const void transfer_from_model(const data_type& t, linear1_prototype_slope_model& lm){
	value_type y(boost::get<0>(t)); 
	if(y==y){
	time_type x(boost::get<1>(t));
	sumX_+=x; lm.sumX_-=x; sumY_+=y; lm.sumY_-=y; sumXY_+=x*y; lm.sumXY_-=x*y;
	y*=y; x*=x;
	sumXX_+=x; lm.sumXX_-=x; sumYY_+=y; lm.sumYY_-=y;
	++size_; --lm.size_;
	update(); lm.update();
	}
      }

      const linear1_prototype_slope_model& operator+=( const linear1_prototype_slope_model& lm){ 
	sumX_+=lm.sumX_; sumXX_+=lm.sumXX_;sumY_+=lm.sumY_; sumYY_+=lm.sumYY_; sumXY_+=lm.sumXY_; size_+=lm.size_; 
	update();
	return *this;
      }
      const linear1_prototype_slope_model& operator-=( const linear1_prototype_slope_model& lm){
	sumX_-=lm.sumX_; sumXX_-=lm.sumXX_;sumY_-=lm.sumY_; sumYY_-=lm.sumYY_; sumXY_-=lm.sumXY_;size_-=lm.size_; 
	update();
	return *this;
      }
      const void transfer_from_model(const linear1_prototype_slope_model& t, linear1_prototype_slope_model& lm){(*this)+=(t);lm-=t; update(); lm.update();}

      parameters_type operator()()const { return boost::make_tuple(offset(*bestSlopeIt_),*bestSlopeIt_); }
      value_type cost() const { return bestCost_ ; }
      size_type size()const {return size_;}

      size_type get_proto_nb()const{ return bestSlopeIt_-protoSlopesBegin_;}

      prototypes_iterator get_prototypes_begin()const{return protoSlopesBegin_;};

      prototypes_iterator get_prototypes_end()const{return protoSlopesEnd_;};

      data_type value(data_type d)const {
	const parameters_type& p(operator()()); 
	return boost::make_tuple(boost::get<0>(p)+boost::get<1>(p)*boost::get<1>(d), boost::get<1>(d));
      }

      tuple_type tuple()const{return boost::make_tuple(operator()(),cost());}

      value_type sum_x()const{return sumX_;}
      value_type sum_y()const{return sumY_;}
      value_type sum_xx()const{return sumXX_;}
      value_type sum_yy()const{return sumYY_;}
      value_type sum_xy()const{return sumXY_;}
    private:
      void update(){
	if(size_<2){bestCost_=0; bestSlopeIt_=protoSlopesBegin_; return;}
	bestCost_=std::numeric_limits<cost_type>::max();
	for(ProtoIt it(protoSlopesBegin_); it!=protoSlopesEnd_; ++it){ 
	  const value_type c(compute_cost(*it));
	  if(c<bestCost_){ bestCost_=c; bestSlopeIt_=it;}
	} 
      }
      cost_type compute_cost(const value_type& a1)const 
      {value_type a0(offset(a1));  return sumYY_-2*(a1*(sumXY_-a0*sumX_)+a0*sumY_)+a1*a1*sumXX_+size_*a0*a0; }
      value_type offset(const value_type& a1)const{return (sumY_-a1*sumX_)/size_;}

      value_type sumX_, sumXX_, sumXY_, sumY_, sumYY_;
      size_type size_;
      prototypes_iterator protoSlopesBegin_, protoSlopesEnd_, bestSlopeIt_;
      cost_type bestCost_;
    };

    /***
	to compute best prototypes, I have to cluster models of segments. As I only have prototypes of one parameter, I can cluster as a segmentation over an ordered set of models (order on the parameter used for the clustering).
	For linear0_prototype_level_model, I can use a linear0_model( and reuse directly the models of segments)
	For linear1_prototype_slope_model, I have to use a new meta-model of models:
	the cost of the meta-model is the sum of costs over each model, for the best slope over the whoel set of models but the best offset for each model
	After some tedious maths, I came up with :
	slope=(SwXY-SXY)/(SwXX-SXX) , where 
	SwXY=sum((1/n_i)*sum(y,startofi,endofi)*sum(x,startofi,endofi),firsti,lasti)
with firsti lasti first and last models of the meta model and startofi (resp.endofi) first(resp.last) point of model i and n_i size of model i.
 SXY=sum(x*y,startoffirsti, endoflasti)
 SwXX=sum((1/n_i)*(sum(x,startofi, endofi)).AŽ²,firsti, lasti)
 SXX=sum(xŽ²,,startoffirsti, endoflasti)

and then cost=(SYY-SwYY)-2*slope*(SXY-SwXY)+slopeŽ²*(SXX-SwXX) with 
SYY=sum(yŽ²,startoffirsti, endoflasti)
 SwYY=sum((1/n_i)*(sum(y,startofi, endofi))Ž²,firsti, lasti)
 SwXY=sum((1/n_i)*sum(x,startofi, endofi)*sum(y,startofi, endofi),firsti, lasti)


     ***/

    template<typename ValYX, typename Size> struct linear1_prototype_slope_meta_model :model_with_time{
      typedef linear1_model<ValYX,Size> data_type;
      typedef typename data_type::value_type value_type;
      typedef typename data_type::data_type free_model_data_type;
      typedef typename data_type::time_type time_type;
      typedef typename data_type::cost_type cost_type;
      typedef Size size_type;
      typedef boost::tuple<value_type> parameters_type; // the slope is the only parameter
      typedef boost::tuple<parameters_type,value_type> tuple_type; // first tuple for params, last element is cost

      template<typename L1Model> struct parameter_selector:std::unary_function<L1Model, value_type>
      {	value_type operator()(const L1Model& lm)const{return boost::get<1>(lm());} };


      linear1_prototype_slope_meta_model():sumYY_(),sumXX_(), sumXY_(), sumWYY_(), sumWXX_(), sumWXY_(), size_(0), slope_(0), cost_(0){}

      // L1Model will either be  linear1_prototype_slope_model or linear1_model

      template<typename L1Model>
      const linear1_prototype_slope_meta_model& operator+=( const L1Model& m){
	if(m.size()){
	sumYY_+=m.sum_yy(); sumXX_+=m.sum_xx(); sumXY_+=m.sum_xy();
	// should check that optimiser is good enough
	sumWYY_+=(m.sum_y()*m.sum_y())/m.size(); sumWXX_+=(m.sum_x()*m.sum_x())/m.size(); sumWXY_+=(m.sum_x()*m.sum_y())/m.size();
	++size_; update();
	}
	return *this;
      }
      template<typename L1Model>
      const linear1_prototype_slope_meta_model& operator-=( const L1Model& m){
	if(m.size()){
	sumYY_-=m.sum_yy(); sumXX_-=m.sum_xx(); sumXY_-=m.sum_xy();
	// should check that optimiser is good enough
	sumWYY_-=(m.sum_y()*m.sum_y())/m.size(); sumWXX_-=(m.sum_x()*m.sum_x())/m.size(); sumWXY_-=(m.sum_x()*m.sum_y())/m.size();
	--size_; update(); 
	}
	return *this;
      }
      template<typename L1Model>      
      const void transfer_from_model(const L1Model& m, linear1_prototype_slope_meta_model& other){
	if(m.size()){
	sumYY_+=m.sum_yy(); sumXX_+=m.sum_xx(); sumXY_+=m.sum_xy();
	other.sumYY_-=m.sum_yy(); other.sumXX_-=m.sum_xx(); other.sumXY_-=m.sum_xy();
	{ const value_type tmp((m.sum_y()*m.sum_y())/m.size()); sumWYY_+=tmp; other.sumWYY_-=tmp; }
	{ const value_type tmp((m.sum_x()*m.sum_x())/m.size()); sumWXX_+=tmp; other.sumWXX_-=tmp; }
	{ const value_type tmp((m.sum_x()*m.sum_y())/m.size()); sumWXY_+=tmp; other.sumWXY_-=tmp; }
	++size_; update(); --(other.size_); other.update();
	}
      }

      const linear1_prototype_slope_meta_model& operator+=( const linear1_prototype_slope_meta_model& lpsm){
	if(lpsm.size()){
	sumYY_+=lpsm.sumYY_; sumXX_+=lpsm.sumXX_; sumXY_+=lpsm.sumXY_;
	sumWYY_+=lpsm.sumWYY_;sumWXX_+=lpsm.sumWXX_; sumWXY_+=lpsm.sumWXY_;
	size_+=lpsm.size_; update(); 
	}
	return *this;
      }
      const linear1_prototype_slope_meta_model& operator-=( const  linear1_prototype_slope_meta_model& lpsm){
	if(lpsm.size()){
	sumYY_-=lpsm.sumYY_; sumXX_-=lpsm.sumXX_; sumXY_-=lpsm.sumXY_;
	sumWYY_-=lpsm.sumWYY_;sumWXX_-=lpsm.sumWXX_; sumWXY_-=lpsm.sumWXY_;
	size_-=lpsm.size_; update();
	}
	return *this;
      }

      const void transfer_from_model(const linear1_prototype_slope_meta_model& t, linear1_prototype_slope_meta_model& lm)
      {(*this)+=(t);lm-=t; }

      parameters_type operator()()const {return boost::make_tuple(slope_);}
      value_type cost() const { return cost_;}
      // is size useful ?
      size_type size()const {return size_;}
      // is value() meaningful ?
      free_model_data_type value(free_model_data_type d)const { boost::get<0>(d)= slope_ *boost::get<1>(d); return d;}
      tuple_type tuple()const{return boost::make_tuple(operator()(),cost());}

    private:
      void update(){ 
	slope_=(sumWXX_!=sumXX_)?(sumWXY_-sumXY_)/(sumWXX_-sumXX_):999; 
	cost_= (sumYY_-sumWYY_) - 2*slope_*(sumXY_-sumWXY_)+slope_*slope_*(sumXX_-sumWXX_);
#if 0
	std::cerr<<"sumWXX_="<<sumWXX_<<" sumXX_="<<sumWXX_<< "size_= "<<size_;
	std::cerr<<" slope_"<<slope_<<" cost_"<<cost_<<'\n';
#endif
      }

      value_type sumYY_, sumXX_, sumXY_, sumWYY_, sumWXX_, sumWXY_ ;
      size_type size_;
      value_type slope_, cost_;
    };


    template< typename Data,typename Size> struct prototypes_model_traits<linear1_prototype_slope_model, Data,Size>{
      typedef linear1_model<Data, Size> free_model;
      typedef linear1_prototype_slope_meta_model<Data, Size> meta_model;
    };


  }
}






















#if 0
    template<typename ProtoIt, typename Val=std::iterator_traits<ProtoIt>::value_type, typename Size=int> struct linear0_prototype_level_model
    {
      typedef Val value_type;
      typedef value_type data_type;
      typedef value_type cost_type;
      typedef Size size_type;
      typedef linear0_model<Val, Size> free_model;

      typedef boost::tuple<value_type,cost_type> lvl_cost_type;
      typedef std::vector<lvl_cost_type> container_type;
      typedef typename container_type::iterator iterator;

      struct initCost0{lvl_cost_type operator()(const value_type& protoLvl)const{return lvl_cost_type(protoLvl,0);} };

      typedef 
      typedef boost::tuple<boost::tuple<value_type>,value_type> tuple_type; // first tuple for params, last element is cost

      linear0_prototype_level_model(In protoLvlBeg, In protoLvlEnd)
	: size_(0), protoLvlsCosts_(boost::make_transform_iterator(protoLvlBeg, initCost0())
				    , boost::make_transform_iterator(protoLvlEnd, initCost0()))
	, bestLvlCost_(protoLvlsCost_.begin()) {}


      const linear0_prototype_level_model& operator+=( const Val& v){
	for(iterator it(protoLvlsCosts_.begin()); it!= protoLvlsCosts_.end(); ++it){
	  const Val tmp(boost::get<0>(*it)-v); 
	  if((boost::get<1>(*it)+= tmp*tmp)<boost::get<1>(*bestLvlCost_)){ bestLvlCost_=it; }
	}
	++size_; return *this;
      }
      const linear0_prototype_level_model& operator-=( const Val& v){
	for(typename constainer_type::iterator it(protoLvlsCosts_.begin()); it!= protoLvlsCosts_.end(); ++it){
	  const Val tmp(boost::get<0>(*it)-v); 
	  if((boost::get<1>(*it)-= tmp*tmp)<boost::get<1>(*bestLvlCost_)){ bestLvlCost_=it; }
	}
	--size_; return *this;
      }
      const linear0_prototype_level_model& operator+=( const linear0_prototype_level_model& lm){
	// could check that we have same prototypes in same order
	for(iterator selfIt(protoLvlsCosts_.begin()), otherIt(lm.protoLvlsCosts_.begin()); selfIt!=protoLvlsCosts_.end(); ++selfIt, ++otherIt){
	  if((boost::get<1>(*selfIt)+= boost::get<1>(*otherIt))<boost::get<1>(*bestLvlCost_)){ bestLvlCost_=selfIt;} 
	}
	size_+=lm.size_; return *this;
      }
      const linear0_prototype_level_model& operator-=( const linear0_prototype_level_model& lm){
	// could check that we have same prototypes in same order
	for(iterator selfIt(protoLvlsCosts_.begin()), otherIt(lm.protoLvlsCosts_.begin()); selfIt!=protoLvlsCosts_.end(); ++selfIt, ++otherIt){
	  if((boost::get<1>(*selfIt)-= boost::get<1>(*otherIt))<boost::get<1>(*bestLvlCost_)){ bestLvlCost_=selfIt;} 
	}
	size_-=lm.size_; return *this;
      }

      const void transfer_from_model(Val v, linear0_prototype_level_model& lm){
	// could check that we have same prototypes in same order
	for(iterator selfIt(protoLvlsCosts_.begin()), otherIt(lm.protoLvlsCosts_.begin()); selfIt!=protoLvlsCosts_.end(); ++selfIt, ++otherIt){
	    Val tmp(boost::get<0>(*selfIt)-v);
	    tmp*=tmp;
	    if((boost::get<1>(*selfIt)+= tmp)<boost::get<1>(*bestLvlCost_)){ bestLvlCost_=selfIt;}
	    if((boost::get<1>(*otherIt)-= tmp)<boost::get<1>(*(lm.bestLvlCost_))){ lm.bestLvlCost_=otherIt;}
	  }
	++size_; --lm.size_;

      }
      const void transfer_from_model(const linear0_prototype_level_model& t, linear0_prototype_level_model& lm){(*this)+=t; lm-=t;}

      value_type operator()()const {return  boost::get<0>(*bestLvlCost_);}
      free_model_data_type value( value_type unused) const{return operator()();}
      cost_type cost() const { return boost::get<1>(*bestLvlCost_);}

      size_type size()const {return size_;}
      tuple_type tuple()const{return boost::make_tuple(boost::make_tuple(operator()()),cost());}
    private:
      container_type protoLvlsCosts_;
      iterator bestLvlCost_;
      size_type size_;
    };

#endif




#endif
