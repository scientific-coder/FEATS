#ifndef FEATS_CLUSTERING_KMEANS_HXX
#define FEATS_CLUSTERING_KMEANS_HXX

#include <iterator>
#include <limits>
#include <algorithm>
#include <ext/algorithm> // random_sample is a SGI extension

#include <vector>
#include <functional>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/iterator/transform_iterator.hpp>

//#include <boost/functional.hpp> // for boost::binder2nd with default constructor 

/**
   This class template is primary used by CROEUC algorithms
   for efficiency reasons, data is not copied.
   Also for efficiency reasons, prototypes will be interlaced in next version for better memory transversal.

   Prototypes are in a twofold structure
   1) prototypes values used for computing distance
   2) prototypes models used for aggregating values

   I could use only one structure of models but initial case would be more complex and prototype value would be computed over and over again.

   As part of the CROEUC algorithm, this class template is to be used for clustering two different series:
   1)usual time series during the first step of the algorithm
   2) APAC or weighted_value time series during subsequent steps 

   It should be noted that all the weighted_value times series have the same weight sequences (because CROEUC is used to compute a global segmentation). It would be very wasteful to store the weights for every time series, so I use an external data structure storing the weigths.

   The only difference between the two uses of kmeans is during distance and cost computing, of course, the weights must be registered during construction.

   We have a weighted_kmeans inheriting from kmeans taking care of the differences. No dynamic polymorphism needed, I should be able to use a kind of Barton Nackman trick to select nearest a compile time .

*/

template<typename KmeanType, typename Arg, typename Nb> 
boost::shared_ptr<KmeanType> best_kmean(const Arg& a, Nb retries){
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
  


template <typename ValItPairsItPair>
struct kmeans{
  typedef typename boost::tuples::element<0, ValItPairsItPair>::type it_iterator;
  typedef typename std::iterator_traits<it_iterator>::difference_type size_type;
  typedef typename std::iterator_traits<it_iterator>::value_type serie_iterators_pair;
  typedef typename boost::tuples::element<0, serie_iterators_pair>::type serie_iterator;
  typedef typename std::iterator_traits<serie_iterator>::value_type value_type;

  typedef value_type cost_type;

  typedef std::vector<value_type> values_seq;
  typedef typename values_seq::iterator values_iterator;
  typedef boost::tuple<values_seq,size_type> proto; // heavy tuple beware of copies !!!
  typedef std::vector<proto> protos_seq;
  typedef typename protos_seq::iterator protos_seq_iterator;

  typedef boost::tuple<cost_type, size_type> compute_type;
  // rem maybe double/int slower than double/double
  struct divider_type: std::unary_function<value_type, value_type>{
    explicit divider_type(size_type s):s_(s){}
    value_type operator()(value_type v)const{return v/s_;}
    size_type s_;
  };
    
  typedef boost::transform_iterator<divider_type, typename values_seq::const_iterator> divider_iterator;


  // first constructor only gives the number of prototypes
  // they will be chosen at random
  kmeans( const boost::tuple<const ValItPairsItPair&, size_type>& vippInit):seriesBegin_(boost::get<0>(boost::get<0>(vippInit))), seriesEnd_(boost::get<1>(boost::get<0>(vippInit))), seriesLength_(std::distance(boost::get<0>(*seriesBegin_), boost::get<1>(*seriesBegin_))),  protoModels_(boost::get<1>(vippInit), proto(values_seq(seriesLength_),0)),protosVals_(seriesLength_*boost::get<1>(vippInit)), cost_(std::numeric_limits<cost_type>::max())
  { 
    std::vector<serie_iterators_pair> protosIt(boost::get<1>(vippInit));
    __gnu_cxx::random_sample(seriesBegin_, seriesEnd_, protosIt.begin(), protosIt.end());
    //    std::cerr<<"init proto are:\n";
    for(typename std::vector<serie_iterators_pair>::const_iterator it(protosIt.begin()); it!= protosIt.end(); ++it){
      //	  std::copy(boost::get<0>(*it), boost::get<1>(*it), std::ostream_iterator<value_type>(std::cerr," "));	  std::cerr<<std::endl;
    }
    set_proto_values(boost::make_tuple(protosIt.begin(), protosIt.end()));
  }

  template<typename InitIt>
  kmeans( const boost::tuple<const ValItPairsItPair&, boost::tuple<InitIt,InitIt>& >& vippInitP)
    :seriesBegin_(boost::get<0>(boost::get<0>(vippInitP))), seriesEnd_(boost::get<1>(boost::get<0>(vippInitP)))
    , seriesLength_(std::distance(boost::get<0>(*seriesBegin_), boost::get<1>(*seriesBegin_)))
    ,protoModels_(std::distance(boost::get<0>(boost::get<1>(vippInitP)), boost::get<1>(boost::get<1>(vippInitP)))
		  , proto(values_seq(seriesLength_),0))
    ,protosVals_(seriesLength_*protoModels_.size())
    ,cost_(std::numeric_limits<cost_type>::max())  { init(boost::get<1>(vippInitP));  }

  template<typename Shit>
  kmeans( const Shit& vippInitP, char*dummy)
    :seriesBegin_(boost::get<0>(boost::get<0>(vippInitP))), seriesEnd_(boost::get<1>(boost::get<0>(vippInitP)))
    , seriesLength_(std::distance(boost::get<0>(*seriesBegin_), boost::get<1>(*seriesBegin_)))
    ,protoModels_(std::distance(boost::get<0>(boost::get<1>(vippInitP)), boost::get<1>(boost::get<1>(vippInitP)))
		  , proto(values_seq(seriesLength_),0))
    ,protosVals_(seriesLength_*std::distance(boost::get<0>(boost::get<1>(vippInitP)), boost::get<1>(boost::get<1>(vippInitP))))
    ,cost_(std::numeric_limits<cost_type>::max())  { 

    //    std::cerr<<"kmeans from wkm: seriesLEngth_:"<<seriesLength_<<" nbProtos:"<<protoModels_.size()<<" protos size:"<<std::distance(boost::get<0>(*boost::get<0>(boost::get<1>(vippInitP))),boost::get<1>(*boost::get<0>(boost::get<1>(vippInitP))))<<"\n";
    init(boost::get<1>(vippInitP));  }
  
  

  compute_type compute(cost_type epsi=1e-6, size_type nbItMax=500){ 
    size_type curIt(0);
    cost_type curDeltaCost(std::numeric_limits<cost_type>::max());
    for(; curIt != nbItMax && curDeltaCost >epsi; ++curIt, curDeltaCost=one_run())
      {/*std::cerr<<"curDeltaCost:"<<curDeltaCost<<"\n";*/}
    return compute_type(curDeltaCost, curIt);
  }
  
  template<typename Out>
  Out output_prototype(size_type n, Out o)const{
    const proto& protoModel(protoModels_[n]);
    divider_type div( boost::get<1>(protoModel));

    return std::copy(divider_iterator(boost::get<0>(protoModel).begin(), div)
		     ,divider_iterator(boost::get<0>(protoModel).end(),div), o);
  }
  template<typename Out>
  Out output_sizes(Out o)const{
    for(typename protos_seq::const_iterator it(protoModels_.begin()); it!= protoModels_.end(); ++it,++o)
      { *o=boost::get<1>(*it);}
    return o;
  }

  //maybe should make cost_ mutable and have this method const
  template<typename Out>
  Out output_clusters(Out o, bool updateP=false) {
    cost_=0.;
    for(it_iterator it(seriesBegin_); it!= seriesEnd_; ++it,++o)
      { *o=present(*it, updateP);}
    return o;
  }

    
    cost_type cost()const{return cost_;}

    virtual ~kmeans(){}
  private:
    // takes a pair of iterators over pairs of iterators
    // I will interlace
    template<typename InitItTuple>
      void init(const InitItTuple& initP){ 
      set_proto_values(initP);
      clear_proto_models();// maybe not useful if constructor clears things
    }
  
    void clear_proto_models(){
      for(typename protos_seq::iterator it(protoModels_.begin()); it != protoModels_.end(); ++it)
	{ std::fill_n(boost::get<0>(*it).begin(), seriesLength_, 0.); boost::get<1>(*it)=0;}

    }

    template<typename InitItTuple>
      void set_proto_values(const InitItTuple& initP){ 
      typedef typename boost::tuples::element<0, InitItTuple>::type proto_it_iterator;
      typedef typename std::iterator_traits<proto_it_iterator>::value_type proto_it_tuple;
      typedef typename boost::tuples::element<0, proto_it_tuple>::type proto_it;
      size_type nbProtos(std::distance(boost::get<0>(initP),boost::get<1>(initP))); 
      values_iterator protosValsIt(protosVals_.begin()), protoValsInitIt(protosVals_.begin());
      std::cerr<<"protosVals_.size():"<<protosVals_.size()<<"\n";
      std::cerr<<"nb protos:"<<nbProtos<<"\n"; //std::cerr<<"val protos:";
      for(proto_it_iterator it(boost::get<0>(initP)); it != boost::get<1>(initP)
	    ; ++it, ++protoValsInitIt, protosValsIt=protoValsInitIt){
	std::cerr<<"proto size:"<<std::distance(boost::get<0>(*it), boost::get<1>(*it))<<"\n";  //std::copy(boost::get<0>(*it), boost::get<1>(*it) , std::ostream_iterator<value_type>(std::cerr," ")); std::cerr<<std::endl;

	for(proto_it sit(boost::get<0>(*it)); sit!= boost::get<1>(*it); ++sit, protosValsIt+=nbProtos){
	  *protosValsIt=*sit;
	}
      }
    }

    virtual size_t present( const serie_iterators_pair& serie, bool updateProto=true){
      const size_type nbP(protoModels_.size());
      //          std::cerr<<"present km of "<<std::distance(boost::get<0>(serie),boost::get<1>(serie))<<"points\n";
      typedef std::vector<value_type> dist_seq;
      dist_seq distances(nbP,0.);
      values_iterator protosValsIt(protosVals_.begin());
      //            std::cerr<<"presenting "; std::copy(boost::get<0>(serie),boost::get<1>(serie), std::ostream_iterator<value_type>(std::cerr," ")); std::cerr<<"\n";

      for(serie_iterator it(boost::get<0>(serie)); it != boost::get<1>(serie); ++it){
	for(typename std::vector<value_type>::iterator distIt(distances.begin()); distIt!=distances.end(); ++distIt, ++protosValsIt)
	  { const value_type delta(*protosValsIt-*it); *distIt+=delta*delta;}
      }
      //      std::cerr<<"distances "; std::copy(distances.begin(),distances.end(), std::ostream_iterator<value_type>(std::cerr," ")); std::cerr<<"\n";

      typename dist_seq::iterator cIt(std::min_element(distances.begin(),distances.end()));
      cost_+=*cIt;
      const size_t clusterNumber(std::distance(distances.begin(), cIt));
      if(updateProto){
	protos_seq_iterator pIt(protoModels_.begin()+clusterNumber);
	// update cluster size
	++boost::get<1>(*pIt);
	// update sum of cluster
	values_iterator vIt(boost::get<0>(*pIt).begin());
	for(serie_iterator it(boost::get<0>(serie)); it != boost::get<1>(serie); ++it, ++vIt){*vIt+=*it;}
      }
      return clusterNumber;
    }
  
    cost_type one_run(){
      cost_type oldCost=cost_;
      cost_=0.;
      clear_proto_models();
      for(it_iterator it(seriesBegin_); it!= seriesEnd_; ++it)
	{ present(*it);}
      // I should be able to use a custom iterator
      typedef std::vector<boost::tuple<divider_iterator, divider_iterator> > divided_seq;
      divided_seq div;
      //    std::cerr<<"computing proto val iter in wkm";
      for(protos_seq_iterator it(protoModels_.begin()); it !=protoModels_.end(); ++it )
	{
	  //	std::cerr<<"size:"<<boost::get<1>(*it)<<" vals ";std::copy(boost::get<0>(*it).begin(),boost::get<0>(*it).end(), std::ostream_iterator<value_type>(std::cerr," "));      std::cerr<<std::endl;

	  divider_type divider(boost::get<1>(*it));
	  div.push_back(boost::make_tuple(divider_iterator(boost::get<0>(*it).begin(), divider)
					  , divider_iterator(boost::get<0>(*it).end(), divider)));
	  //	std::copy(boost::get<0>(div.back()),boost::get<1>(div.back()), std::ostream_iterator<value_type>(std::cerr," ")); std::cerr<<std::endl;

	}
      set_proto_values(boost::make_tuple(div.begin(), div.end()));
      //    std::cerr<<"cost is"<<cost_<<"\n";
      return (oldCost-cost_)/cost_;
    }
  protected:    
    const it_iterator seriesBegin_, seriesEnd_;
    const size_type seriesLength_;
    values_seq protosVals_;
    protos_seq protoModels_;
    cost_type cost_;
  };


  template<typename ValItPairsItPair, typename WeightIt>
  struct weighted_kmeans: kmeans<ValItPairsItPair>{
    typedef typename  kmeans<ValItPairsItPair>::size_type size_type;
    typedef typename  kmeans<ValItPairsItPair>::serie_iterator serie_iterator;
    typedef typename  kmeans<ValItPairsItPair>::values_iterator values_iterator;
    typedef typename  kmeans<ValItPairsItPair>::value_type value_type;
    typedef typename  kmeans<ValItPairsItPair>::protos_seq_iterator protos_seq_iterator;
    typedef typename  kmeans<ValItPairsItPair>::serie_iterators_pair serie_iterators_pair;
    /*
      weighted_kmeans( const boost::tuple<const ValItPairsItPair&,const boost::tuple<WeightIt, WeightIt>&, size_type>& vippWinit)
      :  kmeans<ValItPairsItPair>(boost::make_tuple(boost::get<0>(vippWinit),boost::get<2>(vippWinit)))
      , weightsBegin_(boost::get<0>(boost::get<1>(vippWinit))),weightsEnd_(boost::get<1>(boost::get<1>(vippWinit))){}
    */ 
    template<typename InitIt>
    explicit   weighted_kmeans( const boost::tuple<const ValItPairsItPair&,const boost::tuple<WeightIt, WeightIt>&
				,boost::tuple<InitIt,InitIt> >& vippWinitP)
      : kmeans<ValItPairsItPair>(boost::make_tuple(boost::get<0>(vippWinitP),boost::get<2>(vippWinitP)))
      ,weightsBegin_(boost::get<0>(boost::get<1>(vippWinitP)))
      , weightsEnd_(boost::get<1>(boost::get<1>(vippWinitP))){}

    typedef typename std::vector<size_type>::const_iterator w_it;
    template<typename Shit>
    explicit   weighted_kmeans( const boost::tuple<const ValItPairsItPair&,const boost::tuple<WeightIt, WeightIt>&
				,Shit >& vippWinitP)
      : kmeans<ValItPairsItPair>(boost::make_tuple(boost::get<0>(vippWinitP),boost::get<2>(vippWinitP)))
      ,weightsBegin_(boost::get<0>(boost::get<1>(vippWinitP)))
      , weightsEnd_(boost::get<1>(boost::get<1>(vippWinitP))){}

    template<typename Shit>
    explicit   weighted_kmeans( const Shit& vippWinitP)
      : kmeans<ValItPairsItPair>(boost::make_tuple(boost::get<0>(vippWinitP),boost::get<2>(vippWinitP)))
      ,weightsBegin_(boost::get<0>(boost::get<1>(vippWinitP)))
      , weightsEnd_(boost::get<1>(boost::get<1>(vippWinitP))){}

    template<typename Arg1,typename Arg2, typename Arg3>
    explicit   weighted_kmeans( const Arg1& a1, const Arg2& a2, const Arg3& a3)
      : kmeans<ValItPairsItPair>(a1,/*"toto"*/0), weightsBegin_(a2), weightsEnd_(a3){}


    virtual size_t present( const serie_iterators_pair& serie, bool updateProto=true){
      const size_type nbP(this->protoModels_.size());
      //    std::cerr<<"present wkm of "<<std::distance(boost::get<0>(serie),boost::get<1>(serie)) <<"points\n";
      typedef std::vector<value_type> dist_seq;
      dist_seq distances(nbP,0.);
      values_iterator protosValsIt(this->protosVals_.begin());
      WeightIt wIt(weightsBegin_);
      for(serie_iterator it(boost::get<0>(serie)); it != boost::get<1>(serie); ++it, ++wIt){
	for(typename std::vector<value_type>::iterator distIt(distances.begin()); distIt!=distances.end(); ++distIt, ++protosValsIt)
	  { const value_type delta(*protosValsIt-*it); *distIt+=delta*delta* (*wIt);}
      }
      typename dist_seq::const_iterator cIt(std::min_element(distances.begin(),distances.end()));
      this->cost_+=*cIt;
      const dist_seq& cref(distances);
      const size_t clusterNumber(std::distance(cref.begin(), cIt));
      if (updateProto){
	//    std::cerr<<"wkm present best cluster number "<<clusterNumber<<"\n";
	protos_seq_iterator pIt(this->protoModels_.begin()+clusterNumber);
	// update cluster size
	++boost::get<1>(*pIt);
	// update sum of cluster
	values_iterator vIt(boost::get<0>(*pIt).begin());
	for(serie_iterator it(boost::get<0>(serie)); it != boost::get<1>(serie); ++it, ++vIt){*vIt+=*it;/*std::cerr<<*vIt<<" ";*/}
      }
      return clusterNumber;
    }



    virtual ~weighted_kmeans(){}
  private:
    weighted_kmeans(const weighted_kmeans& k){}
    const WeightIt weightsBegin_, weightsEnd_;
  };

#endif
