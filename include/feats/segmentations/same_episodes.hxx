#ifndef SAME_EPISODES_HXX
#define SAME_EPISODES_HXX

#include<vector>

#include <boost/shared_ptr.hpp>
namespace feats{
  namespace segmentations{

    /*
      this class computes segmentations over a set of time series 
      in fact it compute the same global semgentation amongst all the time series, as a strawman for segmentations_for_total_number

trouble is that we get a segmentation algorithm with a segment_type ofmodels_segment kind to compute segmentations over each set
we cache an aggregation of all the set and perform segmentation on this, so we have to recompute another type segment type because iterators are not the same :-(
hence I created the inner "template typedef" rebind<>::other in segmentation templates

      Static Parameters
      segmentation algorithm i.e. optimal_split<linear0_model<...>
      Note: for clustering of adaptive segmentations, we segment clusters of time series with  linear0_prototype_level_meta_model over the series of the linear0_model over the population

      WeightIt


      constructeur


    */

    template<typename Segmentations>
    struct same_episodes {
      typedef typename Segmentations::iterator sets_iterator;
      typedef boost::tuple<sets_iterator,sets_iterator> sets_it_pairs;
      typedef std::vector<sets_it_pairs> sets_it_pairs_cont;
      typedef typename sets_it_pairs_cont::const_iterator sets_it_pairs_it;
      typedef typename Segmentations::segment_type original_segment_type;
      // non il y a un truc avec le model de models_segments
      typedef typename original_segment_type::model_type model_type;
      typedef std::vector<model_type> inner_models_container;
      typedef typename inner_models_container::const_iterator c_iter;

      typedef typename Segmentations::segment_type::template rebind_iter<c_iter>::other segment_type;
      typedef typename Segmentations::template rebind<segment_type>::other seg_type;
      typedef boost::shared_ptr<seg_type> ptr_seg_type;

      typedef typename seg_type::cost_type cost_type;
      typedef typename seg_type::nb_segments_type nb_segments_type;
      typedef nb_segments_type size_type;

      template<typename SeriesItPairsIt>
      same_episodes(const boost::tuple<SeriesItPairsIt,SeriesItPairsIt>& seriesItPairsItPair)
	: seriesItPairs_(boost::get<0>(seriesItPairsItPair),boost::get<1>(seriesItPairsItPair))
	  ,nbSeries_(seriesItPairs_.size()),nbSegments_(0)
	  ,modelsAggr_(boost::get<0>(*boost::get<0>(seriesItPairsItPair)),boost::get<1>(*boost::get<0>(seriesItPairsItPair)))
      {
	// init segmentations
	// next() because we init with first series
	for(sets_it_pairs_it itIt(boost::next(seriesItPairs_.begin()));itIt != seriesItPairs_.end();++itIt){
	  typename inner_models_container::iterator aggrIt(modelsAggr_.begin());
	  for(sets_iterator it(boost::get<0>(*itIt)); it!= boost::get<1>(*itIt); ++it, ++aggrIt)
	    {(*aggrIt)+=*it;}
	}
	segPtr_=ptr_seg_type(new seg_type(segment_type(modelsAggr_.begin(), modelsAggr_.end())));
	nbSegments_=nbSeries_; 
      }

     
      template<typename Out>
      Out nb_episodes(Out o){
	for(size_type i(0); i!=nbSeries_; ++i,++o)
	  { *o=nbSegments_/nbSeries_;}
	return o;
      }
    
      template<typename Out>
      Out costs(Out o){ // very inefficient, I have to compute segmentations to return costs
	for(size_type i(0); i!=nbSeries_; ++i,++o){
	  std::vector<original_segment_type> tmpSegs;
	  segments(i,std::back_inserter(tmpSegs));
	  cost_type res(0.);
	  for(size_type i(0); i!=tmpSegs.size(); ++i)
	    { res+=tmpSegs[i]->cost();}
	  *o=res;
	}
	return o;
      }

      template<typename Out>
      Out segments(size_type nS, Out o) {
	std::vector<segment_type> segments;
	segPtr_->segments(std::back_inserter(segments));
	original_segment_type tmpSeg(boost::get<0>(seriesItPairs_[nS]));
	for(typename std::vector<segment_type>::const_iterator it(segments.begin()); it!=segments.end(); ++it,++o){
	  while(tmpSeg.size()!=it->size()){tmpSeg.grow_end();}
	  *o=tmpSeg;
	  tmpSeg=original_segment_type(tmpSeg.end());
	}
	return o;
      }

      cost_type cost()const{return segPtr_->cost(); }

      nb_segments_type nb_segments()const{return nbSegments_;}

      void inc_segments(){ nbSegments_+=nbSeries_; segPtr_->inc_segments(); }
    private:
      sets_it_pairs_cont seriesItPairs_;
      size_type nbSeries_, nbSegments_;
      inner_models_container modelsAggr_;
      ptr_seg_type segPtr_;
    };
  }
}
#endif
