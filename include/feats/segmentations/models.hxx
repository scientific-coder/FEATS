#ifndef FEATS_SEGMENTATIONS_MODELS_HXX
#define FEATS_SEGMENTATIONS_MODELS_HXX

#include <boost/tuple/tuple.hpp>
#include <boost/type_traits.hpp>

namespace feats{
  namespace segmentations{

struct model_with_time{  static bool const with_time=true; };
struct model_without_time{ static bool const with_time=false;};

      

    /***
	Models

	rem pour eviter les pb d'arrondi, il faudra maintenir plusieurs sous-sommes.
	Peut-être que deux suffiraient, avec accumulations des valeurs dans la somme intermédiaire et accumulation de la somme intermédiaire dans la plus grande somme.

have tocheck boost::tuple<Val> overhead over Val : shoudl model_without_time take tuple<Val> as data_type ?
    ***/

    template<typename Val, typename Size=int> struct linear0_model:model_without_time
    {
      typedef Val value_type;
      typedef value_type cost_type;
      typedef value_type data_type;
      typedef boost::tuple<value_type> parameters_type;
      typedef Size size_type;

      typedef boost::tuple<parameters_type,value_type> tuple_type; // first tuple for params, last element is cost

      linear0_model():sumY_(),sumYY_(),size_(0){}

      const linear0_model& operator+=( const Val& v){ 
	if(v==v){sumY_+=v; sumYY_+=v*v; ++size_; }
	//	else{std::cerr<<"nan value passed\n";}
	return *this;
      }
      const linear0_model& operator-=( const Val& v){ 
	if(v==v){sumY_-=v; sumYY_-=v*v; --size_;}
	//	else{std::cerr<<"nan value passed\n";}
	return *this;
      }

      const linear0_model& operator+=( const linear0_model& lm){ sumY_+=lm.sumY_; sumYY_+=lm.sumYY_; size_+=lm.size_; return *this;}
      const linear0_model& operator-=( const linear0_model& lm){ sumY_-=lm.sumY_; sumYY_-=lm.sumYY_; size_-=lm.size_; return *this;}

      const void transfer_from_model(Val v, linear0_model& lm){
	if(v==v){
	sumY_+=v; lm.sumY_-=v;
	v *= v;
	sumYY_+=v; lm.sumYY_-=v;
	++size_; --lm.size_;
	//	}else{std::cerr<<"nan value passed\n";}
	}
      }
      const void transfer_from_model(const linear0_model& t, linear0_model& lm)
      {(*this)+=t; lm-=t;}

      parameters_type operator()()const {return  boost::make_tuple((size_==0)? value_type() : sumY_/size_ );}
      data_type value( value_type unused) const{return boost::get<0>(operator()());}
      cost_type cost() const {return (size_==0)? value_type() :  sumYY_ - (sumY_*sumY_)/size_;}
      value_type sum_y()const{ return sumY_;}
      value_type sum_yy()const{return sumYY_;}
      size_type size()const {return size_;}
      tuple_type tuple()const{return boost::make_tuple(operator()(),cost());}
    private:
      value_type sumY_, sumYY_;
      size_type size_;
      };
    
    // ValXY should be a boost::tuple
    template<typename ValYX, typename Size=int> struct linear1_model :model_with_time
    {
      typedef ValYX data_type;
      typedef typename boost::remove_const<typename boost::remove_reference<typename boost::tuples::element<0, data_type>::type>::type>::type value_type;
      typedef  typename boost::remove_const<typename boost::remove_reference<typename boost::tuples::element<1, data_type>::type>::type>::type time_type;
      typedef value_type cost_type;
      typedef Size size_type;
      typedef boost::tuple<value_type, value_type> parameters_type; // tricky to define :-(

      typedef boost::tuple<parameters_type,value_type> tuple_type; // first tuple for params, last element is cost

      linear1_model():sumX_(),sumXX_(), sumXY_(), sumY_(), sumYY_(),size_(0){}

      const linear1_model& operator+=( const data_type& t){ const value_type& y(boost::get<0>(t));
	if(y==y){
	  const time_type& x(boost::get<1>(t)); sumX_+=x; sumXX_+=x*x; sumY_+=y; sumYY_+=y*y; sumXY_+=x*y; ++size_;
	}
	return *this;
      }
	const linear1_model& operator-=( const data_type& t){ const value_type& y(boost::get<0>(t)); 
	  if(y==y){
	    const time_type& x(boost::get<1>(t)); sumX_-=x; sumXX_-=x*x; sumY_-=y; sumYY_-=y*y; sumXY_-=x*y; --size_;
	  }
	  return *this;
	}

      const void transfer_from_model(const data_type& t, linear1_model& lm){
	value_type y(boost::get<0>(t));
	if(y==y){
	  time_type x(boost::get<1>(t));
	  sumX_+=x; lm.sumX_-=x; sumY_+=y; lm.sumY_-=y; sumXY_+=x*y; lm.sumXY_-=x*y;
	  y*=y; x*=x;
	  sumXX_+=x; lm.sumXX_-=x; sumYY_+=y; lm.sumYY_-=y;
	  ++size_; --lm.size_;
	}
      }

      const linear1_model& operator+=( const linear1_model& lm){ sumX_+=lm.sumX_; sumXX_+=lm.sumXX_;sumY_+=lm.sumY_; sumYY_+=lm.sumYY_; sumXY_+=lm.sumXY_; size_+=lm.size_; return *this;}
      const linear1_model& operator-=( const linear1_model& lm){ sumX_-=lm.sumX_; sumXX_-=lm.sumXX_;sumY_-=lm.sumY_; sumYY_-=lm.sumYY_; sumXY_-=lm.sumXY_;size_-=lm.size_; return *this;}
      const void transfer_from_model(const linear1_model& t, linear1_model& lm){(*this)+=(t);lm-=t; }
      // tuple (offset, slope)
      parameters_type operator()()const {
	if(size_==0)return boost::make_tuple(value_type(),value_type());
	value_type c1( (size_==1) ? 
		       0
		       : (sumXY_-(sumX_*sumY_)/size_)
		       /(sumXX_-(sumX_*sumX_)/size_));
	return boost::make_tuple((sumY_-c1*sumX_)/size_,c1);
      }
      cost_type cost() const { 
	if(size_<2)return  value_type() ;
	value_type sXY((sumXY_-(sumX_*sumY_)/size_));
	return (sumYY_-(sumY_*sumY_)/size_)-((sXY*sXY)/(sumXX_-(sumX_*sumX_)/size_));
      }
      size_type size()const {return size_;}
      data_type value(data_type d)const
      { const parameters_type& p(operator()()); boost::get<0>(d)=boost::get<0>(p)+boost::get<1>(p)*boost::get<1>(d); return d;}
      tuple_type tuple()const{return boost::make_tuple(operator()(),cost());}

      value_type sum_x() const{return sumX_;}
      value_type sum_y()const{return sumY_;}
      value_type sum_xx()const{return sumXX_;}
      value_type sum_yy()const{return sumYY_;}
      value_type sum_xy()const{return sumXY_;}


    private:
      value_type sumX_, sumXX_, sumXY_, sumY_, sumYY_;
      size_type size_;
    };

    // some Models will have a state even when empty, to be used as initM in constructors, hence this helper function
    template<typename Model> Model empty_model(Model m){ return m-=m;}


    /*
models_sequence is a model that really is a sequence of models used over sequences of data
the length of data sequence is only known at run-time and we take iterators at the beginning of the sequences
parameters are vectors of the parameters of the models, cost is the sum of the costs of the models
     */

#ifdef PASFINI
    template<typename ModelType, typename Size=int > struct models_sequence_with_time : model_with_time  {
      typedef ModelType model_type;
      typedef Size models_sequence_size;
      typedef typename model_type::data_type_ref data_type_ref;
      typedef typename model_type::value_type value_type;
      typedef typename model_type::time_type time_type;
      typedef typename model_type::data_type data_type;
      typedef typename model_type::cost_type cost_type;
      typedef typename model_type::prototypes_iterator prototypes_iterator;
      typedef typename model_type::size_type size_type;
      typedef std::vector<model_type> models_container;
      typedef typename models_container::iterator models_iter;
      typedef typename models_container::const_iterator models_citer;
      // sequences of values are linked to one time value
      // hence we need to reconstruc tuples of data_type 
      // it is unnecessarily expensive to do but values should not be used too much
      

      struct value_from_model : std::unary_function<model_type, value_type>{
	explicit value_from_model(data_type d):d_(d){}
	value_type operator( const model_type& m)const
	{ return boost::get<0>(m.value(d_));}
      };

      
      typedef boost::transform_iterator<value_from_model, models_citer> values_iterator;
      typedef boost::tuple<values_iterator,time_type> seq_data_type;
      typedef typename model_type::tuple_type model_tuple_type;
      typedef boost::function<model_tuple_type (model_type*)> tuple_from_model;
      typedef boost::transform_iterator<tuple_from_model,models_citer> tuple_type;

      typedef typename model_type::params_type model_params_type;

      explicit models_sequence(Size seqSize,const model_type initM&=model_type()):size_(0),cost_(0.),models_(seqSize,initM){}
      template<typename ProtoIt>
      explicit models_sequence(ProtoIt protosBegin, ProtoIt protosEnd,Size seqSize):size_(0),cost_(0.),models_(seqSize,model_type(protosBegin,protosEnd)){}

      models_sequence(const models_sequence& other):size_(other.size_),cost_(other.cost_),models_(other.models_){}

      models_sequence& operator=(const models_sequence& other){
	if(this != &other){ size_=other.size_;cost_=other.cost_; models_=other.models_;}
	return *this;
      }

      template<typename SeqData>
      const models_sequence& operator+=(SeqData sd){
	const time_type& t(boost::get<1>(sd));
	typename boost::tuples::element<0,SeqData>::type vit(boost::get<0>(sd));
	cost_=0.;++size_;
	for(models_iter mit(models_.begin()); mit != models_.end(); ++mit, ++vi)
	  { (*mit)+=boost::make_tuple(*vi,t); cost_+=mit->cost();}
	return *this;
      }

      template<typename ValIt>
      const models_sequence& operator-=(ValIt vi){
	cost_=0.;--size_;
	for(models_iter mit(models_.begin()); mit != models_.end(); ++mit, ++vi)
	  { (*mit)-=*vi; cost_+=mit->cost();}
	return *this;
      }

      const models_sequence& operator+=(const models_sequence& other){
	cost_=0.;size_+=other.size_;
	for(models_iter mit(models_.begin()), omit(other.models_.begin())
	      ; mit != models_.end(); ++mit, ++omit)
	  { (*mit)+=(*omit); cost_+=mit->cost();}
	return *this;
      }

      const models_sequence& operator-=(const models_sequence& other){
	cost_=0.;size_-=other.size_;
	for(models_iter mit(models_.begin()), omit(other.models_.begin())
	      ; mit != models_.end(); ++mit, ++omit)
	  { (*mit)-=(*omit); cost_+=mit->cost();}
	return *this;
      }


      template<typename ValIt>
      const models_sequence& transfer_from_model(ValIt vi, models_sequence& other){
	cost_=0.;other.cost_=0.;++size_;--(other.size_);
	for(models_iter mit(models_.begin()), omit(other.models_.begin()); mit != models_.end(); ++mit, ++omit, ++vi)
	  { (*mit)+=*vi; (*omit)-=*vi; cost_+=mit->cost();other.cost_+=omit->cost();}
	return *this;
      }


      void transfer_from_model(const models_sequence& t, models_sequence& lm){(*this)+=t; lm-=t;}
      
      // pourrait être un tuple de sequence avec tple de dimension celle des params des models : est-ce que ca en vaut la peine ?
      parameters_type operator()()const {
	return boost::make_tuple(boost::make_transform_iterator(models_.begin()
								, boost::bind(&model_type::operator()(),_1)));
      }      }
      value_type value( value_type v) const{
	return 	  boost::make_tuple(values_type 
			    (boost::make_transform_iterator(models_.begin()
							    , boost::bind(&model_type::value(),_1,v))
			     ,boost::make_transform_iterator(tmp.end()
							     , boost::bind(&model_type::value(),_1,v))));
      }

      cost_type cost() const { 	return cost_;}
      
      // quelle valeur pertinente pour le taille ? le moyenne, mediane, min, max, somme ?
      // pas encore de choix pertinent, donc une valeur NaN
      size_type size()const {return std::numeric_limits<size_type>::quiet_NaN();}
      tuple_type tuple()const {return boost::make_tuple(operator()(), cost());}
    private:
      container_type sums_;
      size_type size_;
      value_type cost_;

      
    };
#endif
  }
}
#endif
