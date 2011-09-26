#ifndef FEATS_MISC_CROEUC_HXX
#define FEATS_MISC_CROEUC_HXX

#include <vector>
#include <iterator>
#include <algorithm> // std::copy
#include <boost/tuple/tuple.hpp>
#include <boost/shared_ptr.hpp>

#include "feats/segmentations/models.hxx"
#include "feats/segmentations/sequences_models.hxx"
#include "feats/segmentations/segment.hxx"
#include "feats/segmentations/optimal_splits.hxx"
#include "feats/clustering/kmeans.hxx"

/**
   zoonekind zoonek site d'un matheu sous linux
   rusell.vcharite.univ-mrs.fr
   implementing croeuc algorithm to cluster time series and segment them at the same time

   parameters are
   1) number of time series clusters
   2) number of segments

   result of the algorithm should be twofold
   1) clusters defined by their values (values are thoses of the segments)
   2) segments sizes of the segmentation

   The croeuc algorithm is iterative because it optimizes both the clusters and the segmentation according to the other, one at a time in each iteration until convergence.

   The CROEUC algorithm is as follows:
   1) clustering time series with kmeans, to compute prototypes
   2) linear0 segmentation of prototypes (weighted according to the size of their respectives clusters), computing only one segmentation for all prototypes
   3) clustering time series segmented with segmentation found in 2) to compute prototypes and loop to 2)

   repeat last steps until convergence.

   It should be noted that kmeans except the first one are not with random init but with prototypes of the previous pass.

   result will be retrieved with 
   output_weights(Out o)
   output_prototype(Size clusterNumber, Out o)

   output_classification(Out o) using the same order for the time series given at construction time, croeuc will output clusternumber of each time series.

   maybe last one does not have to be part of the croeuc if I can factor the functionnality out. prototypes and weights will always be dwarfed by time series for memory size (in reasonable uses) so requiring the user to copy them is no bad. User can dispose of the croeuc after copying to avoid duplicating info during program liftetime.

   kmeans in step 1) is done with kmeans template struct
   To solve segmentation in step 2), we will compute linear0 segmentation of linear0_models. One model for each instant. Each 
   model containing the weighted values of the prototypes.

   kmeans in step 3) is done with weighted_kmeans template struct

   One trouble is that kmeans in step 3) is to be started with prototypes from the last iteration. However, the last iteration was done with prototypes found with an other segmentation, so I have to resegment the prototypes. TO that effect, i use the struct described below. However, I also need to store segments lenght from the last iteration. In order to save space, I reuse clusersSize_ as previousSegmentsLength_.

   After kmeans in step 3), I have to "unsegment" the prototypes. Prototypes are "segmented" because they were computed
   from segmented time series. Hoever, I will use the prototypes to compute the next segmentation, so I need
   "unsegmented" prototypes.

   for illustration purposes, I have to be able to output intermediate results of
   each step of the algorithm


*/




template<typename SerieIt, typename LengthsIt, typename Out>
boost::tuple<Out,typename std::iterator_traits<SerieIt>::value_type>
means(SerieIt sB, SerieIt sE, LengthsIt lB, Out o){
  typedef typename std::iterator_traits<SerieIt>::value_type value_type;
  value_type res;
  for(res=0.; sB != sE; ++o,++lB){
    feats::segmentations::linear0_model<value_type , typename std::iterator_traits<SerieIt>::difference_type> m;
    for(   ;	m.size()!=*lB;  ++sB)
      {m+=*sB;}
    *o=boost::get<0>(m());res+=m.cost();
  }
  return boost::make_tuple(o,res);
}


template<typename ValItPairsItPair>
struct croeuc {
  typedef kmeans<ValItPairsItPair> kmeans_type;
  typedef typename kmeans_type::it_iterator it_iterator;
  typedef typename kmeans_type::size_type size_type;
  typedef typename kmeans_type::serie_iterator serie_iterator;
  typedef typename kmeans_type::value_type value_type;
  typedef typename kmeans_type::cost_type cost_type;

  typedef std::vector<value_type> val_seq;
  typedef std::vector< val_seq > val_seq_seq;
  typedef typename val_seq_seq::iterator val_seq_seq_it;
  typedef feats::segmentations::linear0_model<value_type, size_type> model_type;
  typedef std::vector<model_type> protos_seq;
    typedef feats::segmentations::models_segment<feats::segmentations::linear0_model<model_type>, typename protos_seq::iterator> seg_type;
  typedef feats::segmentations::optimal_split<seg_type> algo_type;

  typedef std::vector<size_type> sizes_seq;

  typedef feats::segmentations::linear0_weighted_sequence_model<typename sizes_seq::const_iterator, value_type, size_type> ws_model_type;
  typedef std::vector<typename val_seq::const_iterator> transposed_it_seq;
    typedef feats::segmentations::segment<ws_model_type, typename transposed_it_seq::const_iterator> ws_seg_type;
  typedef feats::segmentations::optimal_split<ws_seg_type> ws_algo_type;


  // should take args for kmean compute and for best_kmeans also
  explicit croeuc(const ValItPairsItPair& vipip, size_type nbClusts, size_type nbSegs, size_type nbRetries=20, bool illustration=false)
    :seriesBegin_(boost::get<0>(vipip)), seriesEnd_(boost::get<1>(vipip))
    , nbSeries_(std::distance(seriesBegin_, seriesEnd_))
    , nbPrototypes_(nbClusts), nbSegments_(nbSegs)
    , seriesLength_(std::distance(boost::get<0>(*seriesBegin_), boost::get<1>(*seriesBegin_)))
    ,segmentedSeries_(nbSeries_, val_seq(nbSegments_,0.))
    ,prototypes_(nbPrototypes_, val_seq(seriesLength_,0.))
    , segmentedPrototypes_(nbPrototypes_,val_seq(nbSegments_,0.))
    , segmentedPrototypesItPairs_(nbPrototypes_),segmentedSeriesItPairs_(nbSeries_)
    ,segmentsSize_(nbSegments_), clustersSize_(nbPrototypes_)
    ,segmentationCost_(std::numeric_limits<cost_type>::max()),clusteringCost_(std::numeric_limits<cost_type>::max())
  {
    // segmentedSeries_ and segmentedPrototypes_ will not be resized, so iterators are stable
    // I store them once now
    for(size_type i(0); i!=nbSeries_; ++i){
      segmentedSeriesItPairs_[i]=boost::make_tuple(segmentedSeries_[i].begin(), segmentedSeries_[i].end());
    }
    for(size_type i(0); i!=nbPrototypes_; ++i){
      segmentedPrototypesItPairs_[i]=boost::make_tuple(segmentedPrototypes_[i].begin(),segmentedPrototypes_[i].end());
    }
    // compute initial clustering
    boost::shared_ptr<kmeans_type> kmp(best_kmean<kmeans_type>(boost::make_tuple(vipip, nbPrototypes_), nbRetries));
    // store prototypes and clusters size
    kmp->output_sizes(clustersSize_.begin());
    std::cerr<<"initial clusters sizes\n";
    kmp->output_sizes(std::ostream_iterator<size_type>(std::cerr," "));
    std::cerr<<'\n';
    std::cerr<<"\ninitial prototypes:\n";
    for(size_type i(0); i!=nbPrototypes_; ++i){kmp->output_prototype(i,prototypes_[i].begin());
      if(illustration){
      kmp->output_prototype(i,std::ostream_iterator<value_type>(std::cout,"\t"));std::cout<<std::endl; 
      }
    }
    
    // if illustration is true, I do not compute final result at once to let caller get access to intermediate results
    // otherwise, do the whole computation
    if(!illustration){
      optimize();
    }
  }  

  void optimize(){
      next_segmentation(),next_clustering();
      size_type i(0);
      for(cost_type oldCost(std::numeric_limits<cost_type>::max()); (oldCost-cost())/cost()>1e-6
	    ;next_segmentation(),next_clustering(),++i)
	{      std::cerr<<"segmentation cost :"<<segmentation_cost()<<"clustering cost:"<<clustering_cost()<<'\n';oldCost=cost();}
      std::cerr<<"nb iter:"<<i<<'\n';
     std::cerr<<"segmentation cost :"<<segmentation_cost()<<"clustering cost:"<<clustering_cost()<<'\n';
  }    
  cost_type cost() const{ return segmentationCost_ + clusteringCost_;}
  cost_type segmentation_cost()const{return segmentationCost_;}
  cost_type clustering_cost()const{return clusteringCost_;}


  template<typename Out>
  Out output_clusters_size(Out o)const
  { return std::copy(clustersSize_.begin(), clustersSize_.end(),o);}

  template<typename Out>
  Out output_segments_size(Out o)const
  { return std::copy(segmentsSize_.begin(), segmentsSize_.end(),o);}

  template<typename Out>
  Out output_individual_segments_size(Out o)const{
    val_seq_seq transposedSeries(seriesLength_,val_seq(nbSeries_,0.));
    size_type j(0);
    for(it_iterator sIt(seriesBegin_); sIt !=seriesEnd_; ++sIt, ++j){
      size_type i(0);
      for(serie_iterator it(boost::get<0>(*sIt)); it!= boost::get<1>(*sIt);++it,++i){
	transposedSeries[i][j]=*it;
      }
    }
    // do segmentation
    {
  typedef feats::segmentations::linear0_sequence_model<value_type, size_type> seq_model_type;
    typedef feats::segmentations::segment<seq_model_type, typename transposed_it_seq::const_iterator> seq_seg_type;
  typedef feats::segmentations::optimal_split<seq_seg_type> seq_algo_type;

      transposed_it_seq transposedSeriesIterSeq;
      for(val_seq_seq_it it(transposedSeries.begin()); it!=transposedSeries.end(); ++it)
	{ transposedSeriesIterSeq.push_back(it->begin()); }
      seq_seg_type initSeg(transposedSeriesIterSeq.begin(), transposedSeriesIterSeq.end(),seq_model_type(nbSeries_));
      seq_algo_type seg(initSeg);
      while(seg.nb_segs()!=nbSegments_) seg.inc_segments();
      std::cerr<<"done with cost "<<seg.cost()<<"\n";
      std::vector<seq_seg_type>res;
      seg.segments(std::back_inserter(res));

      // get segments size
      for(size_type i(0); i!=nbSegments_; ++i,++o)
	{ *o=std::distance(res[i].begin(), res[i].end()); }
    }
    return o;
  }



  template<typename Out>
  Out output_global_segments_size(Out o)const
  {
    protos_seq prototypesModelsToSeg(seriesLength_);
    std::cerr<<"computing models for GLOBAL segmentation\n";
    for(it_iterator sIt(seriesBegin_); sIt !=seriesEnd_; ++sIt){
      size_type i(0);
      for(serie_iterator it(boost::get<0>(*sIt)); it!= boost::get<1>(*sIt);++it,++i){
	prototypesModelsToSeg[i]+=*it;
      }
    }
    // do segmentation
    {
      seg_type initSeg(prototypesModelsToSeg.begin(), prototypesModelsToSeg.end());
      algo_type seg(initSeg);
      std::cerr<<"\ncomputing segmentation of "<<prototypesModelsToSeg.size()<<" models...";
      while(seg.nb_segs()!=nbSegments_) seg.inc_segments();
      std::cerr<<"done with cost "<<seg.cost()<<"\n";
      std::vector<seg_type>res;
      seg.segments(std::back_inserter(res));

      // get segments size
      for(size_type i(0); i!=nbSegments_; ++i,++o)
	{ *o=std::distance(res[i].begin(), res[i].end()); }
    }
    return o;
  }
  
  template<typename Out>
  Out output_segmented_prototype(size_type n, Out o)const
  { return std::copy(segmentedPrototypes_[n].begin(),segmentedPrototypes_[n].end(),o);}

  template<typename Out>
  Out output_segmented_series(size_type n, Out o)const
  { return std::copy(segmentedSeries_[n].begin(),segmentedSeries_[n].end(),o);}
  template<typename Out>
  Out output_prototype(size_type n, Out o)const
  { return std::copy(prototypes_[n].begin(),prototypes_[n].end(),o);}

  //private: those methods should be private, but for illustration purpuses (illustration=true in constructor)
  // maybe I should make them protected and make a derived class illustration_croeuc exposing them.
  // trouble is that the ancestor's constructor would be doing too much, so I would have to make another cstor doing less
  cost_type next_segmentation() {
    // preconditions: prototypes_ and clustersSize_ are valid
    // postcondition: segmentsSize_, segmentedPrototypes_, segmentedSeries_ and  segmentationCost_ are updated
    std::cerr<<"entering next_segmentation...";
    // compute models used for segmentation
    val_seq_seq transposedPrototypes(seriesLength_,val_seq(nbPrototypes_,0.));
    std::cerr<<"transposing prototypes for segmentation\n";
    for(size_type i(0); i!= seriesLength_; ++i){
      for(size_type j(0); j!=nbPrototypes_; ++j){
	if(clustersSize_[j]){
	  transposedPrototypes[i][j]=prototypes_[j][i];
	}
      }
    }
    // do segmentation
    {
      transposed_it_seq transposedProtoIterSeq;
      for(val_seq_seq_it it(transposedPrototypes.begin()); it!=transposedPrototypes.end(); ++it)
	{ transposedProtoIterSeq.push_back(it->begin()); }
      ws_seg_type initSeg(transposedProtoIterSeq.begin(), transposedProtoIterSeq.end(),ws_model_type(clustersSize_.begin(),clustersSize_.end()));
      ws_algo_type seg(initSeg);
      while(seg.nb_segs()!=nbSegments_) seg.inc_segments();
      std::cerr<<"done with cost "<<seg.cost()<<"\n";
      std::vector<ws_seg_type>res;
      seg.segments(std::back_inserter(res));

      // get segments size
      for(size_type i(0); i!=nbSegments_; ++i)
	{ segmentsSize_[i]=std::distance(res[i].begin(), res[i].end()); }
    }
    std::cerr<<"segments size\n";
    output_segments_size(std::ostream_iterator<value_type>(std::cerr,"\t"));
std::cerr<<"\n";
 
    // compute segmented prototypes

//std::cerr<<"segmented prototypes:\n";
    for(size_type i(0); i!=nbPrototypes_; ++i)
      { means(prototypes_[i].begin(), prototypes_[i].end(), segmentsSize_.begin(),segmentedPrototypes_[i].begin());
	//	std::copy(segmentedPrototypes_[i].begin(),segmentedPrototypes_[i].end() , std::ostream_iterator<value_type>(std::cerr," "));  std::cerr<<std::endl;

}

    val_seq_seq_it segIt(segmentedSeries_.begin());
    // compute segmented series
    //std::cerr<<"segmented series:\n";
    cost_type res(0.);
    for(it_iterator it(seriesBegin_); it!=seriesEnd_; ++it, ++segIt)
      { res+=boost::get<1>(means(boost::get<0>(*it), boost::get<1>(*it),segmentsSize_.begin(),segIt->begin()));
	//	std::copy(segIt->begin(),segIt->end() , std::ostream_iterator<value_type>(std::cerr," "));  std::cerr<<std::endl;
}
    std::cerr<<"segmentation done with cost "<<res<<"\n";
    return segmentationCost_=res;
  }

  cost_type next_clustering() {
    // preconditions segmentedSeries_ segmentsSize_ and segmentedPrototypes_ are valid
    // postconditions segmentedPrototypes_, prototypes_,  are updated
    std::cerr<<"entering next clustering...";
    value_type res(0.);
    
    //         std::cerr<<"segmented prototypes before clustering:\n";
    for(size_type i(0); i!=nbPrototypes_;++i){
      //      std::copy(boost::get<0>(segmentedPrototypesItPairs_[i]),boost::get<1>(segmentedPrototypesItPairs_[i]), std::ostream_iterator<value_type>(std::cerr," ")); std::cerr<<std::endl;
	    }
    /*
    weighted_kmeans<ValItPairsItPair, typename sizes_seq::const_iterator> 
      wkm(boost::make_tuple(boost::make_tuple(segmentedSeriesItPairs_.begin(), segmentedSeriesItPairs_.end())
			    , boost::make_tuple(segmentsSize_.begin(), segmentsSize_.end())
			    ,boost::make_tuple(segmentedPrototypesItPairs_.begin(),segmentedPrototypesItPairs_.end())));
    */
//    std::cerr<<"wkm";

    weighted_kmeans<ValItPairsItPair, typename sizes_seq::const_iterator> 
      wkm(boost::make_tuple(boost::make_tuple(segmentedSeriesItPairs_.begin(), segmentedSeriesItPairs_.end())
			    ,boost::make_tuple(segmentedPrototypesItPairs_.begin(),segmentedPrototypesItPairs_.end()))
	  , segmentsSize_.begin(), segmentsSize_.end());
    wkm.compute();

    //    wkm.output_sizes(clustersSize_.begin());
    std::cerr<<"clusters size\n";
wkm.output_sizes(std::ostream_iterator<size_type>(std::cerr," "));
//    std::cerr<<"segmented prototypes resulting from clustering\n";
    for(size_t i(0); i!=nbPrototypes_; ++i)
      { wkm.output_prototype(i, segmentedPrototypes_[i].begin());
	//std::copy(segmentedPrototypes_[i].begin(),segmentedPrototypes_[i].end() , std::ostream_iterator<value_type>(std::cerr," "));  std::cerr<<std::endl;
}
    // now I must "unsegment" prototypes for next resegmentation
    val_seq distances(nbPrototypes_);
    // zero prototypes_
    for(size_type i(0); i!=nbPrototypes_; ++i)
      {std::fill_n(prototypes_[i].begin(),seriesLength_,0.);}
    std::fill_n(clustersSize_.begin(),nbPrototypes_,0);
    size_type segmentedSerieIndex(0);
    for(it_iterator s(seriesBegin_); s!=seriesEnd_; ++s,++segmentedSerieIndex){
      for(size_type i(0); i!=nbPrototypes_; ++i){
	serie_iterator sIt(boost::get<0>(*s));
	value_type d(0.);
	for(size_type j(0); j!=nbSegments_; ++j){
	  value_type tmp(segmentedPrototypes_[i][j]-segmentedSeries_[segmentedSerieIndex][j]);
	  d+=tmp*tmp*segmentsSize_[j];
	}
	distances[i]=d;
      }
      //      std::cerr<<"distances";      std::copy(distances.begin(),distances.end()		, std::ostream_iterator<value_type>(std::cerr," "));      std::cerr<<std::endl;

      typename val_seq::iterator cIt(std::min_element(distances.begin(), distances.end()));
      res+=*cIt;
      const size_type clusterNumber(std::distance(distances.begin(), cIt));
      ++clustersSize_[clusterNumber];
      serie_iterator sIt(boost::get<0>(*s));
      for(size_type i(0); i!=seriesLength_; ++i,++sIt)
	{ prototypes_[clusterNumber][i]+=*sIt; }
    }

    //    std::cerr<<"unsegmented prototypes:\n";
    for(size_type c(0);c!=nbPrototypes_; ++c){
      for(size_type i(0); i!=seriesLength_; ++i)
	{prototypes_[c][i]/=clustersSize_[c];}
      //      std::copy(prototypes_[c].begin(),prototypes_[c].end()		, std::ostream_iterator<value_type>(std::cerr," "));      std::cerr<<std::endl;
    }
    std::cerr<<"clustering done with cost "<<res<<"\n";
    return clusteringCost_=res;
  }


private:
  it_iterator seriesBegin_, seriesEnd_;
  size_type nbSeries_, nbPrototypes_, nbSegments_, seriesLength_;
  val_seq_seq segmentedSeries_, prototypes_, segmentedPrototypes_;
  std::vector<boost::tuple<typename val_seq::const_iterator, typename val_seq::const_iterator> > segmentedSeriesItPairs_, segmentedPrototypesItPairs_;
  sizes_seq segmentsSize_, clustersSize_;
  cost_type segmentationCost_, clusteringCost_;
};


#endif
