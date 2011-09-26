#ifndef CLUSTERS_OF_ADAPTIVE_SEGMENTATIONS_HXX
#define CLUSTERS_OF_ADAPTIVE_SEGMENTATIONS_HXX

#include <vector>
#include <iterator>
#include <algorithm> // std::copy std::fill
#include <boost/tuple/tuple.hpp>
#include <boost/shared_ptr.hpp>

#include "feats/segmentations/models.hxx"
#include "feats/segmentations/sequences_models.hxx"
#include "feats/segmentations/segment.hxx"
#include "feats/segmentations/optimal_splits.hxx"
#include "feats/segmentations/segmentations_for_total_number.hxx"
#include "feats/clustering/kmeans.hxx"
namespace feats{
  namespace misc{

#ifdef TODO
template<typename KmeanType, typename Arg, typename Nb> 
boost::shared_ptr<KmeanType> best_cluster_of_adaptive_segmentations(const Arg& a, Nb retries){
  boost::shared_ptr<KmeanType> bestPtr(new KmeanType(a));
  bestPtr->compute();
  std::cerr<<"best kmeans:\n";
  while(retries-- !=0){
    boost::shared_ptr<KmeanType> current(new KmeanType(a)); current->compute();
    std::cerr<<"testing:"<<current->cost()<<std::endl;
    if(bestPtr->cost()>current->cost())bestPtr=current;
  }
  std::cerr<<"kept best cost:"<<bestPtr->cost()<<std::endl;
  return bestPtr;
}
#endif


/**
   computing clusters of segmented time series and segmented prototypes for thoses clusters. Each prototype can be segmented with a relevent number of segments. Clusters are computed in a k-means way according to the segmented time series.
However, time series can be segmented according to various segmentations.

First version of the algorithm DOES NOT perform the k-means step, we only tune the clustering without recomputing prototypes

   parameters are
   1) number of time series clusters
   2) total number of segments

   result of the algorithm should be twofold
   1) clusters defined by their values (values are thoses of the segments)
   2) episodes sizes of the segmentations

   The algorithm is iterative because it optimizes both the clusters and the segmentation according to the other, one at a time in each iteration until convergence.

   The  algorithm is as follows:
   1) clustering time series with kmeans, to compute prototypes
   2) segmentations_for_total_number with linear0 model of series, optimizing cost over all the segmentations for a given total number of segments.
   3) clustering time series segmented with segmentations found in 2) to compute prototypes and loop to 2)



   It should be noted that kmeans except the first one are not with random init but with prototypes of the previous pass.

   result will be retrieved with 
   output_weights(Size clusterNumber, Out o) weights are differents for each cluster
   output_prototype(Size clusterNumber, Out o)

   output_classification(Out o) using the same order for the time series given at construction time, croeuc will output clusternumber of each time series.

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


/*
given time series (sB,sE) and episods lengths (starting from lB), output means over episods to o
retuns tuple of last output iterator value and total cost of the model

used to segment prototypes and series over a given set of episods lengths
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


template<typename ValItPairsItPair, template<typename>class SegmentationsType=feats::segmentations::segmentations_for_total_number, template<typename,typename> class AlgoType= feats::segmentations::optimal_split>
struct clusters_of_adaptive_segmentations {
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
  typedef AlgoType<feats::segmentations::optimal_split<seg_type> algo_type;

  typedef std::vector<size_type>  sizes_seq;
/*
we need 2 dimensions array to store episodes size */

 typedef std::vector<sizes_seq > sizes_seq_seq;

/* j'utilisais une seg de weighted_sequence_model :
  typedef feats::segmentations::linear0_weighted_sequence_model<typename sizes_seq::const_iterator, value_type, size_type> ws_model_type;
  typedef std::vector<typename val_seq::const_iterator> transposed_it_seq;
    typedef feats::segmentations::segment<ws_model_type, typename transposed_it_seq::const_iterator> ws_seg_type;
  typedef feats::segmentations::optimal_split<ws_seg_type> ws_algo_type;
maintenant c'est segmentations_for_total_number avec des segmentations utilisant des models_segments */

  typedef std::vector<model_type> protos_seq;
  typedef feats::segmentations::models_segment<feats::segmentations::linear0_model<model_type>, typename protos_seq::iterator> seg_type;
  typedef feats::segmentations::optimal_split<seg_type> algo_type;

  typedef feats::segmentations::segmentations_for_total_number<algo_type> segmentations_type;

  // should take args for kmean compute and for best_kmeans also

  explicit clusters_of_adaptive_segmentations(const ValItPairsItPair& vipip, size_type nbClusts, size_type nbSegs, size_type nbRetries=20, bool initRandom=false, bool illustration=false)
    :seriesBegin_(boost::get<0>(vipip)), seriesEnd_(boost::get<1>(vipip))
    , nbSeries_(std::distance(seriesBegin_, seriesEnd_))
    , nbPrototypes_(nbClusts), nbSegments_(nbSegs)
     , seriesLength_(std::distance(boost::get<0>(*seriesBegin_), boost::get<1>(*seriesBegin_)))
    ,prototypes_(nbPrototypes_, val_seq(seriesLength_,0.))
    , prototypesItPairs_(nbPrototypes_)
     , clusters_(nbSeries_),nbIter_(0)
  {
    // prototypes will not be resized, so iterators will stay valid. I store them once and for all
    for(size_type i(0); i!=nbPrototypes_; ++i)
      { prototypesItPairs_[i]=boost::make_tuple(prototypes_[i].begin(), prototypes_[i].end()); }
    
    if(!initRandom){
    // compute initial clustering
    boost::shared_ptr<kmeans_type>bkm(best_kmean<kmeans_type>(boost::make_tuple(vipip, nbPrototypes_), nbRetries));
    // store prototypes and clusters size
        std::cerr<<"initial clusters sizes\n";
    bkm->output_sizes(std::ostream_iterator<size_type>(std::cerr," "));
    std::cerr<<'\n';
    bkm->output_clusters(clusters_.begin());
    for(size_type i(0); i!=nbPrototypes_; ++i){
      bkm->output_prototype(i,prototypes_[i].begin());
    }
    }else{ // random init clusters : we do not call compute on kmeans
      kmeans_type km(boost::make_tuple(vipip, nbPrototypes_));
      km.output_clusters(clusters_.begin(), true); // we update the prototype
      for(size_type i(0); i!=nbPrototypes_; ++i){
	km.output_prototype(i,prototypes_[i].begin());
      }
    }
    // if illustration is true, I do not compute final result at once to let caller get access to intermediate results
    // otherwise, do the whole computation
    if(!illustration){
      optimize();
    }
  }

  // must FIX API wrt kmeans
  void compute(){ optimize();}

  void optimize(){
      next_segmentation(),next_clustering();
      size_type i(0);
      for(cost_type oldCost(std::numeric_limits<cost_type>::max()); (oldCost-cost())/cost()>1e-6
	    ;next_segmentation(),std::cerr<<"seg cost"<<cost(), next_clustering(),std::cerr<<"clust cost"<<cost(),++i)
	{      std::cerr<<"oldcost:"<<oldCost<<'\n';oldCost=cost();}
      std::cerr<<"nb iter:"<<i<<'\n';nbIter_=i;
      std::cerr<<"cost :"<<cost()<<'\n';
  }    
  cost_type cost() const{ return cost_;}
  size_type nb_iter()const{return nbIter_;}

  template<typename Out>
  Out output_prototype(size_type n, Out o)const
  { return std::copy(prototypes_[n].begin(),prototypes_[n].end(),o);}

  template<typename Out>
  Out output_clusters(Out o)const
  { return std::copy(clusters_.begin(),clusters_.end(), o);}



  //private: those methods should be private, but for illustration purpuses (illustration=true in constructor)
  // maybe I should make them protected and make a derived class illustration_croeuc exposing them.
  // trouble is that the ancestor's constructor would be doing too much, so I would have to make another cstor doing less
  cost_type next_segmentation() {
    // compute clusters series of models
    std::vector<protos_seq> protosModels(nbPrototypes_,protos_seq(seriesLength_));
    size_type j(0);
    for(it_iterator sIt(seriesBegin_); sIt !=seriesEnd_; ++sIt, ++j){
      size_type i(0);
      for(serie_iterator it(boost::get<0>(*sIt)); it!= boost::get<1>(*sIt);++it,++i){
	protosModels[clusters_[j]][i]+=*it;
      }
    }
    std::vector<boost::tuple<typename protos_seq::iterator,typename protos_seq::iterator> > itPairsSeq;
    for(size_type i(0); i!=nbPrototypes_; ++i){
      itPairsSeq.push_back(boost::make_tuple(protosModels[i].begin(), protosModels[i].end()));
    }
    segmentations_type segmentations(boost::make_tuple(itPairsSeq.begin(), itPairsSeq.end()));
    while(segmentations.nb_segments()!=nbSegments_)
      { segmentations.inc_segments(); }
    // reconstruct prototypes
    for(size_type i(0); i!=nbPrototypes_; ++i){
      std::cerr<<"reconstruction du proto "<<i;
      std::vector<seg_type> tmp;
      segmentations.segments(i,std::back_inserter(tmp));
      std::cerr<<" composé de "<<tmp.size()<<" segments#";
      std::fill(prototypes_[i].begin(), prototypes_[i].end(), 0.);
      for(size_type j(0),k(0); j!=tmp.size(); ++j){
	//	  std::cerr<<tmp[j].size()<<" ";
	for(size_type m(0); m!=tmp[j].size();++m, ++k){
	  prototypes_[i][k]=boost::get<0>(tmp[j].model()());
	  std::cerr<<prototypes_[i][k]<<'\t';
	}
      }
      std::cerr<<std::endl;
    }
    return cost_=segmentations.cost();
  }

  cost_type next_clustering(){
    kmeans_type km(boost::make_tuple(boost::make_tuple(seriesBegin_, seriesEnd_)
				     , boost::make_tuple(prototypesItPairs_.begin(), prototypesItPairs_.end())),"toto");
    km.output_clusters(clusters_.begin(),true);
    std::cerr<<"clusters :\n";
    km.output_clusters(std::ostream_iterator<size_type>(std::cerr," "),false);
    std::cerr<<"clusters sizes\n";
    km.output_sizes(std::ostream_iterator<size_type>(std::cerr," "));
    std::cerr<<"clusters means\n";
    for(size_type i(0); i!=nbPrototypes_; ++i, std::cerr<<'\n')
      { km.output_prototype(i, std::ostream_iterator<value_type>(std::cerr," "));}
    return cost_=km.cost();
  }


private:
  it_iterator seriesBegin_, seriesEnd_;
  size_type nbSeries_, nbPrototypes_, nbSegments_, seriesLength_;
  val_seq_seq prototypes_;
  std::vector<boost::tuple<typename val_seq::const_iterator, typename val_seq::const_iterator> > prototypesItPairs_;
  sizes_seq clusters_;
  cost_type cost_;
  size_type nbIter_;
};

  }
}
#endif
