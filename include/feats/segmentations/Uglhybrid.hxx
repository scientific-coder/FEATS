#ifndef UGLHYBRID_HXX
#define UGLHYBRID_HXX


#include "feats/segmentations/binary_splits.hxx"
#include "feats/segmentations/optimal_splits.hxx"
#include "feats/utility/nearly_equals.hxx"

// Quick and very dirty hack to have something to run for the week-end
// comme c'est vite fait mal fait, je ne gère même pas correctement le deltacost, mais c'est pas grave. Gere comme pour optimal splits en ayant un coup d'avance et en stockant le résultat de la segmentation pour n-1 segments.
namespace feats{
  namespace segmentations{


    template<typename Seg, typename Nb=int> struct hybrid_splits {
      typedef Seg segment_type;
      typedef typename segment_type::model_type model_type;
      typedef Nb nb_segments_type;
      typedef typename segment_type::cost_type cost_type;
      typedef optimal_split<Seg, Nb> optimal_seg;
      typedef binary_splits<Seg, Nb> fast_seg;

      template <typename In>
      hybrid_splits(boost::tuple<In,In> initSegsBegEnd):fastSeg_(initSegsBegEnd), optimized_(boost::get<0>(initSegsBegEnd), boost::get<1>(initSegsBegEnd)),upToDate_(false), cost_(0.), winSize_(3){}// cost should be NAN
      void inc_segments()
      { fastSeg_.inc_segments();upToDate_=false;}

      cost_type cost() const{update();return cost_;}
      cost_type delta_cost()const{return fastSeg_.delta_cost();}// NOT CORRECT
      nb_segments_type nb_segs()const{ return fastSeg_.nb_segs();}

      template<typename Out>
      Out segments(Out out)const
      { update(); return std::copy(optimized_.begin(), optimized_.end(),out); }
	
    private:
      typedef std::vector<segment_type> segs_cont;
      typedef typename segs_cont::iterator segs_it;

      // does the local optimisation
      void local_opt(segs_it b, segs_it e)const {
	nb_segments_type winSize(std::distance(b,e));
	if(winSize<2) return;
	segment_type whole(b->begin(), boost::prior(e)->end(),b->model());
	optimal_seg optim(whole);
	while(optim.nb_segs()!=winSize)
	  {optim.inc_segments();}
	optim.segments(b);
      }
      void update() const{
	if(upToDate_){ return ;}
	std::cerr<<"updating";
	segs_cont tmp;tmp.reserve(fastSeg_.nb_segs());
	fastSeg_.segments(std::back_inserter(tmp));
	optimized_.swap(tmp);
	feats::utility::nearly_equals<cost_type> tester;
	for(  cost_type oldCost(std::numeric_limits<cost_type>::max());(oldCost>cost_) &&!tester(oldCost,cost_)
		; oldCost=cost_
		, updateCost() ){
	  tmp=optimized_;
	  segs_it wb(optimized_.begin()), we(optimized_.begin());
	  int ws(std::min(static_cast<nb_segments_type>(std::distance(we,optimized_.end())),winSize_));
	  std::advance(we,ws);
	  for(; we!=optimized_.end(); ++wb, ++we){ // step =1 
	    local_opt(wb,we);
	    std::cerr<<'.';
	  }
	  std::cerr<<cost_;
	  std::cerr<<'!';
	  if(oldCost<cost_){
	    std::cerr<<"worse cost, rolling back";
	    optimized_.swap(tmp);
	  }
	}
	std::cerr<<cost_;
	std::cerr<<"done\n";
	upToDate_=true;
      }
      void updateCost()const{
	cost_=std::accumulate(optimized_.begin(),optimized_.end(),0.
			      // , boost::bind(std::plus<typename segment_type::cost_type>(),_1
			      //   	    , boost::bind(std::mem_fun_ref(&model_type::cost)
			      //   			  ,boost::bind(std::mem_fun_ref(&segment_type::model)
			      //   				       , _2)))
                              ,[](typename segment_type::cost_type c, segment_type const& s){return c + s.model().cost();});
      }
      fast_seg fastSeg_;
      mutable segs_cont optimized_;
      mutable bool upToDate_;
      mutable cost_type cost_;
      nb_segments_type winSize_;
    };
  }
}

#endif
