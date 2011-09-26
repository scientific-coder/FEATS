#ifndef SRC_UTILITY_HXX
#define SRC_UTILITY_HXX

#include <iterator>
#include <string>
#include <iostream>
#include <functional>
#include <numeric>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/bind.hpp>

// should be split to put dependancies into proper include directory
#include "feats/segmentations/segment.hxx"
#include "feats/segmentations/models.hxx"


/***
    template functors to dispatch serie for linear0_model and linear1_model
    the first functor (general case, used by linear0_model) returns the iterator unchanged, while the second (specialized for linear1_model) makes an iterator over tuples to add a "time" to each value.
***/
template<bool, typename It> struct serie_for_model_dispatcher{  typedef It iterator; boost::tuple<It,It> operator()(boost::tuple<It,It> begEnd)const{ return begEnd;} };
template<typename It> struct serie_for_model_dispatcher<true, It>{ 
  typedef boost::counting_iterator<long> time_iter_type;
  typedef  boost::zip_iterator<boost::tuple<It, time_iter_type> > iterator;
  
  boost::tuple<iterator, iterator > operator()(boost::tuple<It,It> begEnd)const {
    time_iter_type first(0),last(0+std::distance(boost::get<0>(begEnd),boost::get<1>(begEnd)));
    return boost::make_tuple(boost::make_zip_iterator(boost::make_tuple(boost::get<0>(begEnd), first))
			     ,boost::make_zip_iterator(boost::make_tuple(boost::get<1>(begEnd), last)));
  }
};

template<typename Model, typename It> struct serie_for_model:serie_for_model_dispatcher<Model::with_time, It>{};


/***
    output functor, is call on the segmentation to display computing time, estimated cost etc...
    may also print the resulting segmentation (depending on outputFormat
***/
template<typename Segmentation> struct some_output{
  typedef typename Segmentation::segment_type segment_type;
  typedef typename segment_type::model_type model_type;
  explicit some_output(const std::string& outputFormat): of_(outputFormat){}
  template<typename Out>
  Out& operator()(const Segmentation& seg, Out& out )const{
    out<<seg.nb_segs()<<'\t'<<seg.cost()<<'\t';
    if(of_!="None"){
      std::vector<segment_type>res; res.reserve(seg.nb_segs());
      seg.segments(std::back_inserter(res));
      out<<res.size()<<'\t'<<std::accumulate(res.begin(),res.end(),0.
					     , boost::bind(std::plus<typename segment_type::cost_type>(),_1
							   , boost::bind(std::mem_fun_ref(&model_type::cost)
									 ,boost::bind(std::mem_fun_ref(&segment_type::model)
										      , _2))))<<'\t';
      typedef  feats::segmentations::seg_to_tuple<segment_type> convert_type;
      typedef std::ostream_iterator<typename segment_type::data_type> os_it_type;
      if(of_=="Segments"){
	std::copy(boost::make_transform_iterator(res.begin(),convert_type(res.front().begin()))
		  , boost::make_transform_iterator(res.end(), convert_type(res.front().begin()))
		  , std::ostream_iterator<typename convert_type::result_type>(out,"\t"));
      }else if(of_=="Points"){
	std::for_each(res.begin(), res.end(),std::bind2nd(std::mem_fun_ref(&segment_type::template values<os_it_type>) 
							  ,os_it_type(out," ")));
      }else{std::cerr<<"unknown output format:"<<of_<<std::endl;}
    }
    return out;
  }
  const std::string& of_;
};

// needed for sequence models: I cannot output the tuple because the models parameters are accessed through an iterator
// the tuple from seg_to_tuple is: (((iteratorOverModelsParams),cost),size)
template <typename TupleType,typename Size> struct output_sequence{
  output_sequence(const TupleType&t , const Size& s):t_(t),s_(s){}
  const TupleType& t_;
  const Size& s_;
};

template<typename Os, typename TupleType, typename Size>
Os& operator<<(Os& os, const output_sequence<TupleType,Size>& tmp){
  typename boost::tuples::element<0,TupleType>::type it(boost::get<0>(tmp.t_));
  os<<"((";
  for(Size i(0); i!=tmp.s_;++i, ++it){
    os<<*it<<'\t';
  }
  os<<')'<<boost::get<1>(tmp.t_)<<')';
  return os;
}

template<typename Segmentation> struct some_sequences_output{
  typedef typename Segmentation::segment_type segment_type;
  typedef typename segment_type::model_type model_type;
  typedef typename model_type::size_type size_type;
  explicit some_sequences_output(const std::string& outputFormat, size_type size): of_(outputFormat),size_(size){}
  
  template<typename Out>
  Out& operator()(const Segmentation& seg, Out& out )const{
    out<<seg.nb_segs()<<'\t'<<seg.cost()<<'\t';
    if(of_!="None"){
      typedef std::vector<segment_type> res_cont_type;
      typedef typename res_cont_type::const_iterator citer_type;
      res_cont_type res; res.reserve(seg.nb_segs());
      seg.segments(std::back_inserter(res));
      out<<res.size()<<'\t'<<std::accumulate(res.begin(),res.end(),0.
					     // , boost::bind(std::plus<typename segment_type::cost_type>(),_1
					     //    	   , boost::bind(std::mem_fun_ref(&model_type::cost)
					     //    			 ,boost::bind(std::mem_fun_ref(&segment_type::model)
					     //    				      , _2)))
                                             ,[](typename segment_type::cost_type c,  segment_type const& s)
                                             {return c+ s.model().cost();}
                                             )<<'\t';
      if(of_=="Segments"){
      typedef  feats::segmentations::seg_to_tuple<segment_type> convert_type;
      typedef  typename convert_type::result_type tuple_type;

      typedef output_sequence<tuple_type, size_type> seq_out_type;
      convert_type c(res.front().begin());
      for(citer_type it(res.begin()); it!=res.end(); ++it){
	seq_out_type tmp(c(*it),size_);
	out<<tmp<<'\t';
      }


      }else if(of_=="Points"){
	/*	std::for_each(res.begin(), res.end(),std::bind2nd(std::mem_fun_ref(&segment_type::template values<os_it_type>) 
							  ,os_it_type(out," ")));
	*/
      }else{std::cerr<<"unknown output format:"<<of_<<std::endl;}
    }
    return out;
  }
  const std::string& of_;
  const size_type size_;
};

template<typename Segmentation, typename Dates> struct some_sequences_output_with_dates{
  typedef typename Segmentation::segment_type segment_type;
  typedef typename segment_type::model_type model_type;
  typedef typename model_type::size_type size_type;
  explicit some_sequences_output_with_dates(const std::string& outputFormat, size_type size, const Dates& d):of_(outputFormat),size_(size), dates_(d){}
  
  template<typename Out>
  Out& operator()(const Segmentation& seg, Out& out )const{
      typedef std::vector<segment_type> res_cont_type;
      typedef typename res_cont_type::const_iterator citer_type;
      res_cont_type res; res.reserve(seg.nb_segs());
      seg.segments(std::back_inserter(res));
      typedef typename segment_type::iterator seg_iter;
      seg_iter base(res.front().begin());
      size_type count(0);
      for(citer_type it(res.begin()); it!=res.end(); ++it,out<<std::endl){
	const model_type& m(it->model());
	
	const size_type nextCount(count+it->size());
	out<<dates_[count]<<'\t'<<dates_[nextCount]<<'\t'<<it->size()<<'\t';
	std::copy(m.value_cost_begin(),m.value_cost_end(),std::ostream_iterator<typename model_type::value_cost_tuple>(out<<boost::tuples::set_open(' ')<<boost::tuples::set_delimiter('\t')<<boost::tuples::set_close(' '),"\t"));
	count=nextCount;
      }
      return out;
  }
  const std::string& of_;
  const size_type size_;
  const Dates& dates_;
};


template<typename In, typename Out>
Out center(In begIn, In begEnd, Out begOut){
  typename std::iterator_traits<In>::value_type mean(0.);
  typename std::iterator_traits<In>::difference_type count(0);
  for(In it(begIn); it!=begEnd; ++it){
    if(*it==*it){
      mean+=*it;
      ++count;
    }
  }
  mean/=(count? count:1);
  //  std::cerr<<"mean:"<<mean<<" count:"<<count<<'\n';
  while(begIn != begEnd)
    { *begOut=(*begIn)-mean; ++begOut; ++begIn; }
  return begOut;
}

template<typename In, typename Out>
Out reduce(In begIn, In begEnd, Out begOut){
  typedef typename std::iterator_traits<In>::value_type value_type;
  feats::segmentations::segment<feats::segmentations::linear0_model<value_type>, In > seg(begIn, begEnd);
  value_type sd(sqrt(seg.model().cost()/seg.model().size()));
  if(sd==0.)sd=1.;
  while(begIn != begEnd)
    { *begOut=(*begIn)/sd; ++begOut; ++begIn;}
  return begOut;
}

// downsampleing using mean, could be median, could discard NaN
template<typename Size, typename In, typename Out>
Out downsample(Size step, In begIn, In begEnd, Out begOut){
  typedef typename std::iterator_traits<In>::value_type value_type;
  In currentEnd(begIn);
  for(;begIn!=begEnd; begIn=currentEnd, ++begOut){
    value_type currentSum(0.);
    Size currentStep(0);
    for(; currentEnd!=begEnd && currentStep!=step; ++currentEnd,++currentStep){
      currentSum+=*currentEnd;
    }
    //last window might not be complete so I cannot use step
    *begOut=currentSum/currentStep;
  }
  return begOut;
}



#endif
