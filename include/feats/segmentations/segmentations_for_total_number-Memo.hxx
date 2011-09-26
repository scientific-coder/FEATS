#ifndef SEGMENTATIONS_FOR_TOTAL_NUMBER_HXX
#define SEGMENTATIONS_FOR_TOTAL_NUMBER_HXX

#include<vector>

#include<boost/shared_prt.hpp>
namespace feats{
  namespace segmentations{

    /*
      this class computes segmentations over a set of time series 

      Static Parameters
      segmentation algorithm i.e. optimal_split<linear0_model<...>
      Note: for clustering of adaptive segmentations, we segment clusters of time series with  linear0_prototype_level_meta_model over the series of the linear0_model over the population

      WeightIt


      constructeur


    */

    template<typename Segmentation>
    struct segmentations_for_total_number{
      typedef Segmentation seg_type;
      typedef boost::shared_ptr<seg_type> ptr_seg_type;
      typedef std::vector<ptr_seg_type> segs_container;

      typedef typename seg_type::cost_type cost_type;
      typedef std::vector<cost_type> costs_cont;
      typedef std::vector<costs_cont> costs_cont_cont;
      typedef typename seg_type::nb_segments_type nb_segments_type;
      typedef nb_segments_type size_type;
      typedef std::vector<nb_segments_type> nb_segs_cont;
      typedef std::vector<nb_segs_cont> nb_segs_cont_cont;
      typedef typename seg_type::segment_type segment_type;
      typedef std::vector<segment_type> segments_cont;
      typedef std::vector<segments_cont> segments_cont_cont;
      typedef std::vector<segments_cont_cont> segments_cont_cont_cont;

      template<typename SeriesItPairsItPair>
      segmentations_for_total_number(const boost::tuple<SeriesItPairsIt,SeriesItPairsIt>& seriesItPairsItPair)
	: nbSeries_(std::distance(boost::get<0>(seriesItPairsItPair),boost::get<1>(seriesItPairsItPair))), nbSegments_(0)
	  ,segs_(nbSeries_), costs_(nbSeries_), segCosts_(nbSeries_),nbSegs_(nbSeries_)
      {
	// init segmentations
	SeriesItPairsIt it(boost::get<0>(seriesItPairsItPair));
	for(size_type i(0); i!=nbSeries_; ++i,++it){
	  segs_[i]=ptr_seg_type(new typename seg_type::segment_type(boost::get<0>(*it), boost::get<1>(*it)));
	  // init memoize for costs
	  segCosts_[i].push_back(segs_[i]->cost());
	  // init memoize for segments
	  segments_[i].push_back(segments_cont(1))
	  segs_[i]->segments(std::back_inserter(segments_[i][1]));
	}
	while(nbSegments_ != nbSeries_)
	  { inc_segments(); }
      }

     
      template<typename Out>
      Out nb_episodes(Out o){
	nb_segments_type nbs(nbSegments_);
	for(size_type i(0); i!=nbSeries_; ++i,++o){
	  nb_segments_type tmp(segs_[nbSeries_-i][nbs]);
	  nbs-=tmp;
	  *o=tmp;
	}
	return o;
      }
    
      template<typename Out>
      Out costs(Out o){
	nb_segments_type nbs(nbSegments_);
	for(size_type i(0); i!=nbSeries_; ++i,++o){
	  nb_segments_type tmp(segs_[nbSeries_-i][nbs]);
	  *o=costs_[nbSeries_-i][nbs];
	  nbs-=tmp;
	}
	return o;
      }

      template<typename Out>
      Out segments(size_type nbSeries, Out o){
	nb_segments_type nbs(nbSegments_);
	for(size_type i(nbSeries_); i!=nbSeries; --i)
	  { nbs-=segs_[i][nbs]; }
	return std::copy(segments_[nbSeries][nbs].begin()o;
      }


      cost_type cost()const{return costs_[nbSeries_][nbSegments_];}
      nb_segments_type nb_segments()const{return nbSegments_;}

      void inc_segments(){
	++nbSegments_;
	// we compute repartition of nbS segments over the first i series
	for(size_type i(0); i!=nbSeries_); ++i){
	// last time series will get j segments
	cost_type bestCost(std::numeric_limits<cost_type>::max());
	nb_segments_type bestNbS(0);
	if(i<nbSegment_){ // it is legal to compute 
	  for(nb_segments_type j(1); j!=nbSegments_-i+1; ++j){
	    currentCost= segmentationCost(i, j);
	    if(i) currentCost+=costs_[i-1][nbSegments_-j];
	    if(currentCost < bestCost){
	      bestCost=currentCost; bestNbS=j;
	    }
	  }
	}
	costs_[i].push_back(bestCost);
	nbSegs_[i].push_back(bestNbS);
      }
    private:
      /* not allowed*/
      weighted_segmentations_for_total_number(const weighted_segmentations_for_total_number& other){/* not allowed*/}
      weighted_segmentations_for_total_number& operator=(const weighted_segmentations_for_total_number& other){/* */}

      // memoize costs for segmentations
      // must also memoize segments
      cost_type segmentationCost(size_type nbS, nb_segments_type nbE) {
	while(seg_[nbS].nb_segs()<nbE) { 
	  segs_[nbS].inc_segments(); 
	  segCosts_[nbS].push_back(segs_[nbS].cost());
	  segments_[i].push_back(segments_cont(1))
	  segs_[nbS]->segments(std::back_inserter(segments_[nbS]));
	}
	return segCosts_[nbS][nbE];
      }

      size_type nbSeries_;
      nb_segments_type nbSegments_;
      segs_container segs_;
      costs_cont_cont costs_,segCosts_;
      nb_segs_cont_cont nbSegs_;  
      segments_cont_cont_cont segments_;
    };
  }
}
#endif
