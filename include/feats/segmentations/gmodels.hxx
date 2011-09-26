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

have tocheck boost::tuple<Val> overhead over Val : should model_without_time take tuple<Val> as data_type ?
    ***/

    template<typename Val, typename Size=int> struct linear0_model:model_without_time
    {
      typedef Val value_type;
      typedef value_type data_type;
      typedef boost::tuple<value_type> parameters_type;
      typedef Size size_type;

      typedef boost::tuple<parameters_type,value_type> tuple_type; // first tuple for params, last element is cost

      linear0_model():sumY_(),sumYY_(),size_(0){}


      const linear0_model& operator+=( const Val& v){ sumY_+=v; sumYY_+=v*v; ++size_; return *this;}
      const linear0_model& operator-=( const Val& v){ sumY_-=v; sumYY_-=v*v; --size_; return *this;}

      const linear0_model& operator+=( const linear0_model& lm){ sumY_+=lm.sumY_; sumYY_+=lm.sumYY_; size_+=lm.size_; return *this;}
      const linear0_model& operator-=( const linear0_model& lm){ sumY_-=lm.sumY_; sumYY_-=lm.sumYY_; size_-=lm.size_; return *this;}

      const void transfer_from_model(Val v, linear0_model& lm){
	sumY_+=v; lm.sumY_-=v;
	v *= v;
	sumYY_+=v; lm.sumYY_-=v;
	++size_; --lm.size_;
      }
      const void transfer_from_model(const linear0_model& t, linear0_model& lm){(*this)+=t; lm-=t;}

      parameters_type operator()()const {return  boost::make_tuple((size_==0)? value_type() : sumY_/size_ );}
      data_type value( value_type unused) const{return boost::get<0>(operator()());}
      value_type cost() const { return (size_==0)? value_type() :  sumYY_ - (sumY_*sumY_)/size_ ;}

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
      typedef Size size_type;
      typedef boost::tuple<value_type, value_type> parameters_type; // tricky to define :-(

      typedef boost::tuple<parameters_type,value_type> tuple_type; // first tuple for params, last element is cost

      linear1_model():sumX_(),sumXX_(), sumXY_(), sumY_(), sumYY_(),size_(0){}

      const linear1_model& operator+=( const data_type& t){ const value_type& y(boost::get<0>(t)); const time_type& x(boost::get<1>(t)); sumX_+=x; sumXX_+=x*x; sumY_+=y; sumYY_+=y*y; sumXY_+=x*y; ++size_; return *this;}
      const linear1_model& operator-=( const data_type& t){ const value_type& y(boost::get<0>(t)); const time_type& x(boost::get<1>(t)); sumX_-=x; sumXX_-=x*x; sumY_-=y; sumYY_-=y*y; sumXY_-=x*y; --size_; return *this;}

      const void transfer_from_model(const data_type& t, linear1_model& lm){
	value_type y(boost::get<0>(t)); time_type x(boost::get<1>(t));
	sumX_+=x; lm.sumX_-=x; sumY_+=y; lm.sumY_-=y; sumXY_+=x*y; lm.sumXY_-=x*y;
	y*=y; x*=x;
	sumXX_+=x; lm.sumXX_-=x; sumYY_+=y; lm.sumYY_-=y;
	++size_; --lm.size_;
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
      value_type cost() const { 
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

  }
}
#endif
