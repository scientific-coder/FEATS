#ifndef FEATS_SEQUENCES_SEGMENTATIONS_PROTOTYPES_MODELS_HXX
#define FEATS_SEQUENCES_SEGMENTATIONS_PROTOTYPES_MODELS_HXX

#include "feats/segmentations/models.hxx"
#include "feats/segmentations/sequences_models.hxx"
#include "feats/segmentations/prototypes_models.hxx"
#include <assert.h>

namespace feats{
  namespace segmentations{

    template<typename ProtoIt, typename Val=typename std::iterator_traits<ProtoIt>::value_type, typename Size=int> struct linear0_prototype_level_sequence_model:model_without_time
    {
      typedef Val value_type;
      typedef value_type data_type;
      typedef value_type cost_type;
      typedef Size size_type;

      typedef boost::tuple<value_type, value_type, Size, ProtoIt> params_type;
      typedef std::vector<params_type > container_type;

      typedef typename container_type::iterator iterator;
      typedef typename container_type::const_iterator const_iterator;

      struct value_from_params: std::unary_function<params_type,value_type>{
	value_type operator()(const params_type&p)const{
	  return boost::get<2>(p)
	    ? *boost::get<3>(p) : std::numeric_limits<value_type>::quiet_NaN();
	}
      };

      typedef boost::transform_iterator<value_from_params, const_iterator> values_iterator;
      typedef boost::tuple<value_type> model_parameters;
      struct parameters_from_params: std::unary_function<params_type,model_parameters>{
	model_parameters operator()(const params_type&p)const{
	  value_from_params v;  return boost::make_tuple(v(p));
	}
      };

      typedef boost::transform_iterator<parameters_from_params, const_iterator> parameters_iterator;

      typedef boost::tuple<model_parameters,cost_type> model_tuple_type;
      struct tuple_from_params: std::unary_function<params_type,model_tuple_type>{

	model_tuple_type operator()(const params_type&p)const{
	  value_type value(std::numeric_limits<value_type>::quiet_NaN());
	  cost_type cost(0.);
	  if(boost::get<2>(p)){
	    value=*boost::get<3>(p);
	    cost=boost::get<1>(p) - 2*value*(boost::get<0>(p))+boost::get<2>(p)*value*value;
	  }
	  return boost::make_tuple(boost::make_tuple(value),cost);
	}
      };

      typedef boost::transform_iterator<tuple_from_params, const_iterator> tuples_iterator;




      typedef boost::tuple<value_type,cost_type> value_cost_tuple;
      struct value_cost_from_params: std::unary_function<params_type, value_cost_tuple>{
	value_cost_tuple operator()(const params_type&p)const{
	  value_type value(std::numeric_limits<value_type>::quiet_NaN());
	  cost_type cost(0.);
	  if(boost::get<2>(p)){
	    value=*boost::get<3>(p);
	    cost=boost::get<1>(p) - 2*value*(boost::get<0>(p))+boost::get<2>(p)*value*value;
	  }
	  return value_cost_tuple(value,cost);
	}
      };
      typedef boost::transform_iterator<value_cost_from_params, const_iterator> value_cost_iterator;




      linear0_prototype_level_sequence_model(ProtoIt protoLvlBeg, ProtoIt protoLvlEnd, Size nbP)
	: protoBegin_(protoLvlBeg), protoEnd_(protoLvlEnd),params_(nbP,params_type(0.,0.,0,protoEnd_)),cost_(0.),size_(0) {
	std::cerr<<"in proto seq mod with"<<std::distance(protoBegin_,protoEnd_)<<" protos and sequences of size"<<nbP<<'\n'; 
}

      linear0_prototype_level_sequence_model( const linear0_prototype_level_sequence_model& other):
	protoBegin_(other.protoBegin_),protoEnd_(other.protoEnd_), params_(other.params_), cost_(other.cost_),size_(other.size_){}
      linear0_prototype_level_sequence_model& operator=( const linear0_prototype_level_sequence_model& other){
	if(this != &other){
	  protoBegin_=other.protoBegin_; protoEnd_=other.protoEnd_;
	  params_=other.params_; cost_=other.cost_; size_=other.size_;
	}
	return *this;
      }
    
      template<typename ValIt>
      const linear0_prototype_level_sequence_model& operator+=( ValIt vi){
	cost_=0.;++size_;
	for(iterator it(params_.begin()); it!= params_.end(); ++it, ++vi){
	  const value_type& v(*vi);
	  params_type& p(*it);
	  if(v==v){
	    boost::get<0>(p)+=v;
	    boost::get<1>(p)+=v*v;
	    ++boost::get<2>(p);
	  }
	  updateCost(p);
	}
	return *this;
      }

      template<typename ValIt>
      const linear0_prototype_level_sequence_model& operator-=( ValIt vi){
	cost_=0.; --size_;
	for(iterator it(params_.begin()); it!= params_.end(); ++it, ++vi){
	  const value_type& v(*vi);
	  params_type& p(*it);
	  if(v==v){
	    boost::get<0>(p)-=v;
	    boost::get<1>(p)-=v*v;
	    --boost::get<2>(p);
	  }
	  updateCost(p);
	}
	return *this;
      }

      const linear0_prototype_level_sequence_model& operator+=( const linear0_prototype_level_sequence_model& lm){
	const_iterator otherIt(lm.params_.begin());
	cost_=0.; size_+=lm.size_;
	for(iterator selfIt(params_.begin()); selfIt!=params_.end(); ++selfIt, ++otherIt){
	  params_type& p(*selfIt);
	  boost::get<0>(p) += boost::get<0>(*otherIt);
	  boost::get<1>(p) += boost::get<1>(*otherIt);
	  boost::get<2>(p) += boost::get<2>(*otherIt);
	  updateCost(p);
	}
	return *this;
      }
      const linear0_prototype_level_sequence_model& operator-=( const linear0_prototype_level_sequence_model& lm){
	const_iterator otherIt(lm.params_.begin());
	cost_=0.; size_-=lm.size_;
	for(iterator selfIt(params_.begin()); selfIt!=params_.end(); ++selfIt, ++otherIt){
	  params_type& p(*selfIt);
	  boost::get<0>(p) -= boost::get<0>(*otherIt);
	  boost::get<1>(p) -= boost::get<1>(*otherIt);
	  boost::get<2>(p) -= boost::get<2>(*otherIt);
	  updateCost(p);
	}
	return *this;
      }

      template<typename ValIt>
      const void transfer_from_model(ValIt vIt, linear0_prototype_level_sequence_model& lm){

	cost_= lm.cost_= 0. ;++size_, --lm.size_;
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
	  updateCost(sp);
	  lm.updateCost(op);
	}
      }
      const void transfer_from_model(const linear0_prototype_level_sequence_model& t, linear0_prototype_level_sequence_model& lm){(*this)+=t; lm-=t;}
      
      parameters_iterator operator()()const 
      {return  parameters_iterator(params_.begin());}

      value_cost_iterator value_cost_begin()const{return value_cost_iterator(params_.begin());}
      value_cost_iterator value_cost_end()const{return value_cost_iterator(params_.end());}
      values_iterator value( value_type unused) const{values_iterator(params_.begin());}
      cost_type cost() const {return cost_;}
      // ! sans doutes supprimer size de l'interface des models, utiliser celle de segment
      // pas possible car utilise dans optim de segments
      size_type size()const {return size_;}

      typedef tuples_iterator tuple_type;
      tuples_iterator tuple()const {return tuples_iterator(params_.begin());}

      size_type sequence_size()const{return params_.size();}

      typedef boost::tuple<boost::tuple<std::string>,std::string> tuple_names_type;
      tuple_names_type tuple_names()const
      {return tuple_names_type(boost::make_tuple("level"),"cost");}

      template<typename Out> Out mean_costs(size_type n, Out o){
	for(ProtoIt i(protoBegin_); i!=protoEnd_; ++i, ++o){ *o=costForValue(params_[n], *i)/size_; }
	return o;
      }
      cost_type cost(size_type n) {return costForValue(params_[n], *boost::get<3>(params_[n]));}

      template<typename Out> Out protos_number(Out o) const{ 
	for(const_iterator i(params_.begin()); i!=params_.end(); ++i, ++o){ *o=std::distance(protoBegin_,boost::get<3>(*i)); }
	return o;
      }

    private:

      value_type costForValue( params_type& p, value_type v) const
	{ return boost::get<1>(p) - 2 * v *boost::get<0>(p) + boost::get<2>(p)*v*v; }

      value_type costForValue( params_type& p, value_type v) 
	{ return boost::get<1>(p) - 2 * v *boost::get<0>(p) + boost::get<2>(p)*v*v; }
	
      void updateCost(params_type& p){
	if(boost::get<2>(p)){
	cost_type best(std::numeric_limits<cost_type>::max()), tmp;
	for(ProtoIt pIt(protoBegin_); pIt != protoEnd_; ++pIt){
	  if((tmp =costForValue(p,*pIt))<best)
	    { best=tmp; boost::get<3>(p)=pIt;}
	}
	//	std::cout<<"best:"<<best<<"proto"<<*boost::get<3>(p)<<std::endl;
	//	std::cerr<<"updating cost"<<best<<"for mean"<< it->first/size_<<'\n';
	cost_+=best;
	assert(!(cost_!=cost_));// detect nan cost_
	}
      }
      ProtoIt protoBegin_, protoEnd_;
      container_type params_;
      value_type cost_;
      size_type size_;
    };

    
template< typename Data,typename Size> struct prototypes_model_traits<linear0_prototype_level_sequence_model, Data,Size>{
  typedef typename std::iterator_traits<Data>::value_type value_type;
  typedef linear0_model<value_type, Size> free_points_model;
      typedef linear0_sequence_model<value_type, Size> free_model;
      typedef linear0_prototype_level_meta_model<value_type, Size> meta_model;
    };



#ifdef BOF
    //////////////////////////////////////////////
    template<typename ProtoIt, typename value_type,typename time_type, typename Size=int> struct linear1_prototype_slope_sequence_model:model_with_time
    {
      typedef ProtoIt prototypes_iterator;
      typedef boost::tuple<value_type, time_type> data_type;
      typedef value_type cost_type;
      typedef boost::tuple<value_type, value_type> model_parameters_type;
      typedef boost::tuple<parameters_type, cost_type> model_tuple_type;
      typedef Size size_type;
      /*
	I cannot factor sumX and sumXX because an x should not be taken into account when the y is NA
       */
      struct model_data{
	explicit model_data(prototypes_iterator pit)
	  :sumX_(0.),sumXX_(0.),sumXY_(0.), sumY_(0.),sumYY_(0.),size_(0),bestSlopeIt_(pit){}
	value_type sumX_;
	value_type sumXX_;
	value_type sumXY_;
	value_type sumY_;
	value_type sumYY_;
	size_type size_;
	cost_type bestCost_;
	prototypes_iterator bestSlopeIt_;

	cost_type add(value_type const y, time_type const x){
	  if(y==y){
	    sumX_+=x; sumXX_+=x*x; sumY_+=y; sumYY_+=y*y; sumXY_+=x*y; ++size_; 
	    update();
	  }
	  return bestCost_;
	}

	cost_type add( const model_data& other){
	    sumX_+=other.sumX_; sumXX_+=other.sumXX_; sumY_+=other.sumY_; sumYY_+=other.sumYY_; sumXY_+=other_.sumXY_; size_+=other_.size_; 
	    update();
	  }
	  return bestCost_;
	}


	cost_type remove(value_type const y, time_type const x){
	  if(y==y){
	    sumX_-=x; sumXX_-=x*x; sumY_-=y; sumYY_-=y*y; sumXY_-=x*y; --size_; 
	    update();
	  }
	  return bestCost_;
	}

	cost_type remove( const model_data& other){
	    sumX_-=other.sumX_; sumXX_-=other.sumXX_; sumY_-=other.sumY_; sumYY_-=other.sumYY_; sumXY_-=other_.sumXY_; size_-=other_.size_; 
	    update();
	  }
	  return bestCost_;
	}

	value_type offset(const value_type a1)const{ return (sumY_-a1*sumX_)/size_;}
	value_type offset()const{ return offset(*bestSlopeIt_);}
	cost_type compute_cost( const value_type& a1) const {
	  value_type a0(offset(a1));
	  return sumYY_-2*(a1*(sumXY_-a0*sumX_)+a0*sumY_)+a1*a1*sumXX_+size_*a0*a0; 
	}

	void update(prototypes_iterator beg, prototypes_iterator en){
	  if(size_>=2){
	    bestCost_=std::numeric_limits<cost_type>::max();
	    for(prototypes_iterator it(beg); it!=en; ++it){ 
	      const value_type c(compute_cost(*it));
	      if(c<bestCost_){ bestCost_=c; bestSlopeIt_=it;}
	    } 
	  }else{ bestCost_=0; bestSlope=beg;}
	}
	value_type value(time_type const t)const{
	  return size_
	    ? *(bestSlopeIt_)*t+offset()
	    : std::numeric_limits<value_type>::quiet_NaN();
	}

	model_parameters_type operator()()const
	{return model_parameters_type(offset(),*bestSlopeIt);}

    model_tuple_type tuple()const{return model_tuple_type(operator()(),cost_);}

	prototypes_iterator bestSlopeIt_;
	cost_type bestCost_;
      };
      typedef std::vector<model_data> container_type;

      typedef typename container_type::iterator iterator;
      typedef typename container_type::const_iterator const_iterator;

      struct value_from_model: std::unary_function<model_data,value_type>{
	explicit value_from_model(const time_type& t):t_(t){}
	value_type operator()(const model_data& m)const{ return m.value(t_);}
	time_type t_;
      };
      typedef boost::transform_iterator<value_from_model, const_iterator> values_iterator;
  typedef boost::tuple<values_iterator, time_type> values_it_time;
      typedef boost::function<model_parameters_type, (model_data*)> parameters_from_model;
      typedef boost::transform_iterator<parameters_from_model, const_iterator> parameters_iterator;



      typedef boost::function<model_tuple_type, (model_data*)> tuples_from_model;
      typedef boost::transform_iterator<tuples_from_model, const_iterator> tuples_iterator;


      linear1_prototype_slope_sequence_model(ProtoIt protoSlopeBeg, ProtoIt protoSlopeEnd, Size nbP)
	: protoBegin_(protoLvlBeg), protoEnd_(protoLvlEnd),cost_(0.),size_(0),models_(nbP,model_data(protoEnd_)) {
	std::cerr<<"in proto seq mod with"<<std::distance(protoBegin_,protoEnd_)<<" protos and sequences of size"<<nbP<<'\n'; 
}

      linear1_prototype_slope_sequence_model( const linear1_prototype_slope_sequence_model& other):
	protoBegin_(other.protoBegin_),protoEnd_(other.protoEnd_), cost_(other.cost_),size_(other.size_),models_(other.models_){}
      linear1_prototype_slope_sequence_model& operator=( const linear1_prototype_slope_sequence_model& other){
	if(this != &other){
	  protoBegin_=other.protoBegin_; protoEnd_=other.protoEnd_; 
	  cost_=other.cost_; size_=other.size_; models_=other.models_; 
	}
	return *this;
      }
      // DataIt is a boost::tuple<begin iterator over seqSize value_type, time_type>
      template<typename DataIt>
      const linear1_prototype_slope_sequence_model& operator+=( const DataIt& di){
	typename boost::remove_const<typename boost::remove_reference<typename boost::tuples::element<0, DataIt>::type>::type>::type vi(boost::get<0>(di));
	const time_type& x(boost::get<1>(di));
	cost_=0.;++size_;
	for(iterator it(models_.begin()); it!= models_.end(); ++it, ++vi){
	  cost_+=it->add(*vi,x);
	}
	return *this;
      }

      template<typename DataIt>
      const linear1_prototype_slope_sequence_model& operator-=( const DataIt& di){
	typename boost::remove_const<typename boost::remove_reference<typename boost::tuples::element<0, DataIt>::type>::type>::type vi(boost::get<0>(di));
	conat time_type& x(boost::get<1>(di));
	cost_=0.;--size_;
	for(iterator it(models_.begin()); it!= models_.end(); ++it, ++vi){
	  cost_+=it->remove(*vi,x);
	}
	return *this;
      }

      const linear1_prototype_slope_sequence_model& operator+=( const linear1_prototype_slope_sequence_model& lm){
	const_iterator otherIt(lm.models_.begin());
	cost_=0.; size_+=lm.size_;
	for(iterator selfIt(models_.begin()); selfIt!=models_.end(); ++selfIt, ++otherIt){
	  cost_+=selfIt->add(*otherIt);
	}
	return *this;
      }
      const linear1_prototype_slope_sequence_model& operator-=( const linear1_prototype_slope_sequence_model& lm){
	const_iterator otherIt(lm.models_.begin());
	cost_=0.; size_-=lm.size_;
	for(iterator selfIt(models_.begin()); selfIt!=models_.end(); ++selfIt, ++otherIt){
	  cost_+=selfIt->remove(*otherIt);
	}
	return *this;
      }


      template<typename DataIt>
      const linear1_prototype_slope_sequence_model& operator-=( const DataIt& di){
	typename boost::remove_const<typename boost::remove_reference<typename boost::tuples::element<0, DataIt>::type>::type>::type vi(boost::get<0>(di));
	const time_type& x(boost::get<1>(di));
	cost_=0.;--size_;
	for(iterator it(models_.begin()); it!= models_.end(); ++it, ++vi){
	  cost_+=it->remove(*vi,x);
	}
	return *this;
      }

      template<typename DataIt>
      const void transfer_from_model(DataIt di, linear1_prototype_slopesequence_model& lm){
	cost_= lm.cost_= 0. ;++size_, --lm.size_;
	typename boost::remove_const<typename boost::remove_reference<typename boost::tuples::element<0, DataIt>::type>::type>::type vi(boost::get<0>(di));
	const time_type& x(boost::get<1>(di));
	for(iterator selfIt(models_.begin()), otherIt(lm.models_.begin()); selfIt!=models_.end(); ++selfIt, ++otherIt, ++vi){
	  cost_+=selfIt->add(*vi,x);
	  lm.cost_+=otherIt->remove(*vi,x);
	}
      }
  const void transfer_from_model(const linear1_prototype_slope_sequence_model& t, linear1_prototype_slope_sequence_model& lm){(*this)+=t; lm-=t;}
      
  parameters_type operator()()const 
      {return  parameters_iterator(models_.begin(),parameters_from_model(&model_type::operator()));}
  template<typename DataIt>
      values_it_time value( const DataIt& dit) const  {
    const time_type x(boost::get<1>(dit));
    return values_it_time(values_iterator(models_.begin(), value_from_model(x)),x);
  }
      cost_type cost() const { return cost_;}
      // ! sans doutes supprimer size de l'interface des models, utiliser celle de segment
      // pas possible car utilise dans optim de segments
      size_type size()const {return size_;}
  tuples_iterator tuple()const {return tuples_iterator(models_.begin(),tuples_from_model(&model_data::tuple));}

      size_type sequence_size()const{return params_.size();}
    private:

      ProtoIt protoBegin_, protoEnd_;
      value_type cost_;
      size_type size_;
      container_type models_;
    };

    template< typename Data,typename Size> struct prototypes_model_traits<linear1_prototype_level_sequence_model, Data,Size>{
      typename std::iterator_traits<typename boost::remove_const<typename boost::remove_reference<typename boost::tuples::element<0, DataIt>::type>::type>::type >::value_type value_type;
      typename boost::remove_const<typename boost::remove_reference<typename boost::tuples::element<1, DataIt>::type>::type>::type time_type;
      typedef boost::tuple<value_type,time_type> data_type;
      typedef linear1_model<data_type, Size> free_points_model;
      typedef linear1_sequence_model<Data, Size> free_model;
      typedef linear1_prototype_slope_meta_model<data_type, Size> meta_model;
    };

#endif
  }
}

#endif
