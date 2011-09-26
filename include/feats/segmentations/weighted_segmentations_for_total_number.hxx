#ifndef WEIGHTED_SEGMENTATION_FOR_TOTAL_NUMBER_HXX
#define WEIGHTED_SEGMENTATION_FOR_TOTAL_NUMBER_HXX

#include<vector>

#include<boost/shared_prt.hpp>
namespace feats{
  namespace segmentations{

    /*
      this class computes segmentations over a set of time series each of the time series has a weigth (the number of elements in the cluster it represents in clusters_of_adaptive_segmentations)

      Static Parameters
      segmentation algorithm i.e. optimal_split<linear_0_model<...>

      WeightIt


      constructeur


    */

    template<typename Segmentation, typename WeightIt>
    struct weighted_segmentations_for_total_number{
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

      typedef WeightIt weights_it;

      template<typename SeriesItIt>
      weighted_segmentations_for_total_number(const boost::tuple<weights_it,weights_it>& weightsBegEnd, SeriesItPairsIt SeriesItPairsIt)
	: weightsBegin_(boost::get<0>(weightsBegEnd)), weightEnd_(boost::get<1>(weightsBegEnd))
	,nbSeries_(std::distance(weightsBegin_,weightsEnd_)), nbSegments_(0),segs_(nbSeries_)
	,costs_(nbSeries_), segCosts_(nbSeries_),nbSegs_(nbSeries_)
      {
	// init segmentations
	for(size_type i(0); i!=nbSeries_; ++i,++seriesItPairsIt){
	  segs_[i]=ptr_seg_type(new typename seg_type::segment_type(boost::get<0>(*seriesItPairsIt), boost::get<1>(*seriesItPairsIt)));
	  segCosts_[i].push_back(segs[i]->cost());
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
    cost_type segmentationCost(size_type nbS, nb_segments_type nbE) {
      while(seg_[nbS].nb_segs()<nbE) { segs_[nbS].inc_segments(); segCosts_[nbS].push_back(segs_[nbS].cost());}
      return segCosts_[nbS][nbE];
    }

    weights_it weightsBegin_, weightsEnd_;
    size_type nbSeries_;
    nb_segments_type nbSegments_;
    segs_container segs_;
    costs_cont_cont costs_,segCosts_;
    nb_segs_cont_cont nbSegs_;  
    };
  }
}
#endif
