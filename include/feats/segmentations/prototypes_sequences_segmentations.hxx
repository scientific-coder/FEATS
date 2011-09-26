#ifndef FEATS_SEQUENCES_SEGMENTATIONS_PROTOTYPES_SEGMENTATION_HXX
#define FEATS_SEQUENCES_SEGMENTATIONS_PROTOTYPES_SEGMENTATION_HXX

#include <functional>
#include <limits>
#include <boost/bind.hpp>
#include <boost/type_traits.hpp>
#include <boost/shared_ptr.hpp>
#include "feats/segmentations/prototypes_models.hxx"
#include "feats/segmentations/utility.hxx"
#include "feats/utility/nearly_equals.hxx"

//#define DBG 1

namespace feats{
  namespace segmentations{
    // this one is used to compute optimal prototype parameters, given a sequence of models
    // for now I only accept "top down" segmentations that is, with inc_segments()
    // not really sure about the API should I create my own copy of models ?
    // pros: I have to sort it, and usually I would have a Colleciton of Segments not models: copying models makes a dense representation
    // cons: I usually do not copy (for exemple time series)
    // If goal is execution speed I think copying from sparse models of segments is can be more efficient as long as extra memory usage does not cause swapping.

    template<typename meta_model, template<typename, typename> class Segmentation > struct sequence_parameter_quantizer{
      typedef typename boost::tuples::element<0,typename meta_model::parameters_type>::type parameter_type;
      typedef typename meta_model::data_type model_type;
      typedef std::vector<model_type> container_type;
      typedef typename container_type::const_iterator const_iterator;
      typedef segment<meta_model, const_iterator > segment_type;
      typedef Segmentation<segment_type, int> segmentation_type;
      typedef typename segmentation_type::cost_type cost_type;
      typedef typename segmentation_type::nb_segments_type nb_prototypes_type;
      template<typename In>
      sequence_parameter_quantizer(const boost::tuple<In,In>& modelsBeginEnd):models_(boost::get<0>(modelsBeginEnd)
									    ,boost::get<1>(modelsBeginEnd)){
	typename meta_model::template parameter_selector<model_type> selectParam;
	std::sort(models_.begin(), models_.end(), boost::bind(std::less<parameter_type>()
							      ,boost::bind(selectParam,_1)
							      ,boost::bind(selectParam,_2)));
#ifdef DBG
	std::cout<<"in quantizer :";
	for(const_iterator it(models_.begin()); it!=models_.end(); ++it){
	  std::cout<<it->tuple()<<','<<it->size()<<'\t';
	}
	std::cout<<std::endl;
#endif
	init_segments<segmentation_type> initSeg;
	seg_=new segmentation_type(initSeg(boost::make_tuple(models_.begin(), models_.end())));
      }
      // boost bind should replace this
      struct get_parameter_from_segment: std::unary_function<segment_type,parameter_type>
      { parameter_type operator()(const segment_type& s)const{ return boost::get<0>(s.model()());} };
      template<typename Out>
      Out prototypes(Out out)const{
	// There should be a better way !
	std::vector<segment_type> tmp; tmp.reserve(seg_->nb_segs());
#if DBG
	std::cerr<<"calling prototypes with"<<nb_prototypes()<<" protos over"<<models_.size()<<"models\n";
#endif
	seg_->segments(std::back_inserter(tmp));
	return std::transform(tmp.begin(), tmp.end(), out, get_parameter_from_segment());
      }
      nb_prototypes_type nb_prototypes()const{return seg_->nb_segs();}
      void inc_prototypes(){ seg_->inc_segments();}
      cost_type cost()const{return seg_->cost();}
      cost_type delta_cost()const{ return seg_->delta_cost();}
      ~sequence_parameter_quantizer(){delete seg_;}
    private:
      container_type  models_;
      segmentation_type* seg_;

    };

    /***
	API not frozen , I will have to implement more algorithms in order to find right abstractions
	For now, I prefer to give direct access to underlying segmentations, as I do not know if it has inc_segments or dec_segments
	I try to keep quantizer and prototype segmentations around, so I use mutable pointers with the following dependancies: proto_segmentation depends on parameters quantizer which depends on free segmentation
	meta_model is used by quantizer over free_model s to compute prototype parameters to be used by proto_model

I have a problem because prototype models templates depend on the type of the iterator over prototype parameters the later being decided by prototypes_segmentations, so prototypes_segmentation cannot receive the prototype_model class as parameter. 

     ***/

    template<template<typename,typename, typename > class ProtoModelT, typename Iter, template<typename, typename >class Segmentation=optimal_splits, template<typename, typename >class QuantizerSegmentation=optimal_splits> struct prototypes_sequence_segmentations {
      typedef Iter serie_iterator;
      typedef typename std::iterator_traits<serie_iterator>::value_type sequence_type;
      typedef typename std::iterator_traits<sequence_type>::value_type data_type;
      typedef typename std::iterator_traits<sequence_type>::difference_type size_type;
      typedef typename prototypes_model_traits<ProtoModelT,sequence_type,size_type>::free_points_model free_model;
      typedef typename prototypes_model_traits<ProtoModelT,sequence_type,size_type>::free_model free_sequence_model;
      typedef typename prototypes_model_traits<ProtoModelT,sequence_type,size_type>::meta_model meta_model;
      /*
      typedef linear0_model<data_type> free_model;
      typedef linear0_model<data_type> free_model;
      */
      typedef segment<free_sequence_model, serie_iterator> free_segment;
      typedef Segmentation<free_segment, int > free_segmentation_type;
      typedef sequence_parameter_quantizer<meta_model, QuantizerSegmentation> quantizer;
      typedef typename quantizer::parameter_type parameter_type; 
      typedef std::vector<typename boost::remove_reference<typename boost::remove_const<parameter_type>::type>::type> protos_container;
      typedef typename protos_container::const_iterator protos_it;
      typedef ProtoModelT<protos_it, data_type, int> proto_model;
      typedef segment<proto_model, serie_iterator> proto_segment;
      typedef Segmentation<proto_segment, int> prototypes_segmentation;
      typedef typename prototypes_segmentation::cost_type cost_type;
      typedef const boost::tuple<serie_iterator,serie_iterator> serie_boundaries_type;
      typedef int sequences_size;
      prototypes_sequence_segmentations(const serie_boundaries_type& begEnd, sequences_size seqSize)
	:beginEnd_(begEnd),freeSeg_(init_segments<free_segmentation_type>()(begEnd,free_sequence_model(seqSize)))
	,seqSize_(seqSize)
			      //	 , ptrQuantizer_(0)
	 ,protoParams_(0)
			      //	 , ptrProtoSeg_(0)
      {
#ifdef DBG
std::cout<<" in prototypes_sequences_segmentation multiv nbp="<<std::distance(boost::get<0>(begEnd),boost::get<1>(begEnd))<<"seq size"<<seqSize<<std::endl;
#endif
}

      const free_segmentation_type& get_free_segmentation()const{return freeSeg_;}
      const free_segmentation_type& free_segmentation()const{return freeSeg_;}
      free_segmentation_type& free_segmentation(){invalidate_protos(); return freeSeg_;}

      const quantizer& get_parameter_quantizer()const{ return *ensure_quantizer();}
      const quantizer& parameter_quantizer()const{ return *ensure_quantizer();}
      quantizer& parameter_quantizer(){ invalidate_segmentation(); return *ensure_quantizer();}
      
      const prototypes_segmentation& get_segmentation()const{return *ensure_segmentation();}
      const prototypes_segmentation& segmentation()const{return *ensure_segmentation();}
      prototypes_segmentation& segmentation(){return *ensure_segmentation();}
      
      // returns the number of optimization steps
      int optimize(){
#if DBG
	std::cerr<<"in optimize nb Protos"<<get_parameter_quantizer().nb_prototypes()<<" nb segments"<<get_free_segmentation().nb_segs()<<'\n';
#endif
	ensure_segmentation();
	feats::utility::nearly_equals<cost_type> cmp;
	int nbSteps(0);

	for(cost_type lastCost(std::numeric_limits<cost_type>::max()); !cmp(lastCost,get_segmentation().cost());  ++nbSteps ){
#if DBG
	  std::cerr<<"lastCost"<<lastCost<<" currentCost"<<get_segmentation().cost()<<'\n';
#endif
	  lastCost=get_segmentation().cost();
	  optimization_step_subopt();
	}
	return nbSteps;
      }

    private:

      struct free_model_from_proto_segment:std::unary_function<proto_segment, free_model>
      { free_model operator()(const proto_segment& ps)const{ free_segment tmp(ps.begin(), ps.end());return tmp.model();} };

      void optimization_step_opt(){
      // recompute optimal prototypes, but this time over the segments found using prototypes models 
	// note that we have to compute free_model s over thoses segments because the prototype models do not have the partial sums taht we need
	{
	  // first I have to save nb of prototypes because I will invalidate prototypes
	  typename quantizer::nb_prototypes_type nbP(ptrQuantizer_->nb_prototypes());
#if DBG
	  std::cerr<<"in optimization_step nbP "<<nbP;
#endif
	  // must get segments before because segmentation is invalidated by invalidate_protos
	  std::vector<proto_segment>tmp;tmp.reserve(get_segmentation().nb_segs()); get_segmentation().segments(std::back_inserter(tmp));
	  invalidate_protos();
	  size_type sequencesSize(tmp.front().model().sequence_size());
	  std::vector<free_model> modelsToQuantize(sequencesSize*tmp.size());
	  serie_iterator sIt(boost::get<0>(beginEnd_));
	  size_type firstIndex(0);
	  for(size_type iSegment(0); iSegment != tmp.size(); ++iSegment, firstIndex+=sequencesSize){
	    for(size_type iSequence(0); iSequence!=tmp[iSegment].size(); ++iSequence, ++sIt){
	      sequence_type currentSequence(*sIt);
	      for(size_type iPoint(0); iPoint!=sequencesSize; ++iPoint, ++currentSequence){
		modelsToQuantize[firstIndex+iPoint]+=*currentSequence;
	      }
	    }
	  }

	  /*

	   */
	  ptrQuantizer_= boost::shared_ptr<quantizer>(new quantizer(boost::make_tuple(modelsToQuantize.begin(), modelsToQuantize.end())));
	  // pb si un seul segment et 2 protos avec seg optimale!
#if DBG
	  std::cerr<<"nb Protos"<<nbP<<" nb segs"<<get_free_segmentation().nb_segs()<<'\n';
#endif
	  while(--nbP) ptrQuantizer_->inc_prototypes();
	}
	protoParams_.resize(ptrQuantizer_->nb_prototypes());// useful if optimize_step was called before any valid PtrQuantizer_ was ensured
	ptrQuantizer_->prototypes(protoParams_.begin());
#if DBG
	std::cerr<<"protos:";
	std::copy(protoParams_.begin(), protoParams_.end(), std::ostream_iterator<double>(std::cerr,"\t"));
#endif
	// invalidate_segmentation() was called by invalidate_protos()
	prototypes_segmentation* ptrS=ensure_segmentation();
	while(ptrS->nb_segs()<freeSeg_.nb_segs())ptrS->inc_segments();
      }

      // compile time template should allow to call directly optimize_step_opt if Seg is optimal but it is not possible to specialize a member template function without specializing the template class :-(
      // this method is necessary because sub optimal segmentation can worsen the solution, so I make backup of the current solution and restore it is cost in worse after optimization_step_opt()

      void optimization_step_subopt(){
	cost_type oldCost=get_segmentation().cost();
	boost::shared_ptr<quantizer> backQt(ptrQuantizer_);
	boost::shared_ptr<prototypes_segmentation> backSeg(ptrProtoSeg_);
	optimization_step_opt();
	if(get_segmentation().cost()>oldCost){ std::cerr<<"Worse cost:rolling back optimisation!\n";ptrQuantizer_=backQt; ptrProtoSeg_=backSeg;}
      }
      quantizer* ensure_quantizer()const{
	if(!ptrQuantizer_){
	  std::vector<free_segment> tmp;tmp.reserve(freeSeg_.nb_segs()); freeSeg_.segments(std::back_inserter(tmp));

	  size_type sequencesSize(seqSize_);
	  std::vector<free_model> modelsToQuantize(sequencesSize*tmp.size());
	  serie_iterator sIt(boost::get<0>(beginEnd_));
	  size_type firstIndex(0);
	  for(size_type iSegment(0); iSegment != tmp.size(); ++iSegment, firstIndex+=sequencesSize){
	    for(size_type iSequence(0); iSequence!=tmp[iSegment].size(); ++iSequence, ++sIt){
	      sequence_type currentSequence(*sIt);
	      for(size_type iPoint(0); iPoint!=sequencesSize; ++iPoint, ++currentSequence){
		modelsToQuantize[firstIndex+iPoint]+=*currentSequence;
	      }
	    }
	  }

	  /*

	   */
	  ptrQuantizer_= boost::shared_ptr<quantizer>(new quantizer(boost::make_tuple(modelsToQuantize.begin(), modelsToQuantize.end())));
	}
	return ptrQuantizer_.get();
      }
      prototypes_segmentation* ensure_segmentation()const{
	if(!ptrProtoSeg_){
	  if(protoParams_.size()!=get_parameter_quantizer().nb_prototypes()){
	    protos_container tmp; 
	    tmp.reserve(ptrQuantizer_->nb_prototypes());
	    ptrQuantizer_->prototypes(std::back_inserter(tmp));
	    protoParams_.swap(tmp);
	  }
	  init_segments<prototypes_segmentation> initS;
	  ptrProtoSeg_=boost::shared_ptr<prototypes_segmentation>(new prototypes_segmentation(initS(beginEnd_, proto_model(protoParams_.begin(), protoParams_.end(), seqSize_))));
	  while(ptrProtoSeg_->nb_segs()!=freeSeg_.nb_segs()) ptrProtoSeg_->inc_segments();
	}
	return ptrProtoSeg_.get();
      }
      void invalidate_protos(){	ptrQuantizer_.reset();invalidate_segmentation(); }
      void invalidate_segmentation(){  ptrProtoSeg_.reset(); }
      const serie_boundaries_type beginEnd_;
      free_segmentation_type freeSeg_;
      int seqSize_;
      mutable boost::shared_ptr<quantizer> ptrQuantizer_;
      mutable protos_container protoParams_;
      mutable boost::shared_ptr<prototypes_segmentation> ptrProtoSeg_;
    
    };



  }
}
#endif
