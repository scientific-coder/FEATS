#ifndef EQUAL_NB_OF_SEGMENTS_HXX
#define EQUAL_NB_OF_SEGMENTS_HXX

#include<vector>

#include <boost/shared_ptr.hpp>
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
    struct equal_nb_of_segments {
      typedef Segmentation seg_type;
      typedef boost::shared_ptr<seg_type> ptr_seg_type;
      typedef std::vector<ptr_seg_type> segs_container;

      typedef typename seg_type::cost_type cost_type;
      typedef typename seg_type::nb_segments_type nb_segments_type;
      typedef nb_segments_type size_type;

      template<typename SeriesItPairsIt>
      equal_nb_of_segments(const boost::tuple<SeriesItPairsIt,SeriesItPairsIt>& seriesItPairsItPair)
	: nbSeries_(std::distance(boost::get<0>(seriesItPairsItPair),boost::get<1>(seriesItPairsItPair))), nbSegments_(0)
	  ,segs_(nbSeries_)
      {
	// init segmentations
	SeriesItPairsIt it(boost::get<0>(seriesItPairsItPair));
	for(size_type i(0); i!=nbSeries_; ++i,++it){
	  segs_[i]=ptr_seg_type(new seg_type(typename seg_type::segment_type(boost::get<0>(*it), boost::get<1>(*it))));
	}
	nbSegments_=nbSeries_; 
      }

     
      template<typename Out>
      Out nb_episodes(Out o){
	for(size_type i(0); i!=nbSeries_; ++i,++o)
	  { *o=segs_[i]->nb_segs();}
	return o;
      }
    
      template<typename Out>
      Out costs(Out o){
	for(size_type i(0); i!=nbSeries_; ++i,++o)
	  { *o=segs_[i]->cost();}
	return o;
      }

      template<typename Out>
      Out segments(size_type nS, Out o)
      { return segs_[nS]->segments(o);}

      cost_type cost()const{
	cost_type res(0.);
	for ( size_type i(0); i!=nbSeries_; ++i)
	  { res+=segs_[i]->cost();}
	return res;
      }

      nb_segments_type nb_segments()const{return nbSegments_;}

      void inc_segments(){
	++nbSegments_;
	// we find the first series having less segments than the previous one
	size_type i(1);
	for(size_type prev=segs_[0]->nb_segs(); i!=nbSeries_ && prev==segs_[i]->nb_segs(); ++i)
	  { prev= segs_[i]->nb_segs(); }
	if(i==nbSeries_)i=0;
	segs_[i]->inc_segments();
      }
    private:
      size_type nbSeries_, nbSegments_;
      segs_container segs_;

    };
  }
}
#endif
