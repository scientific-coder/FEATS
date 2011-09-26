#ifndef FEATS_CLUSTERING_KMEANS_HXX
#define FEATS_CLUSTERING_KMEANS_HXX

#include <iterator
#include <algorithm>

#include "feats/segmentation/models.hxx"

/**
This class template is primary used by CROEUC algorithms
for efficiency reasons, data is not copied.
Also for efficiency reasons, prototypes will be interlaced in next version for better memory transversal.

Prototypes are in a twofold structure
1) prototypes values used for computing distance
2) prototypes models used for aggregating values

I could use only one structure of models but initial case would be more complex and prototype value would be computed over and over again.



 */
template <typename ValItPairsItPair>
struct kmeans{
  typedef typename boost::tuples::element<0, ValItPairsItPair>::type it_iterator;
  typedef typename std::iterator_traits<it_iterator>::difference_type size_type;
  typedef typename std::iterator_traits<it_iterator>::value_type serie_iterators_pair;
  typedef typename boost::tuples::element<0, serie_iterators_pair>::type serie_iterator;
  typedef typename std::iterator_traits<serie_iterator>::value_type value_type;

  typedef linear0_model<value_type> proto_model_type;
  typedef std::vector<value_type> proto_values;
  typedef std::vector<proto_model_type> proto_models; // they are not interlaced
  typedef std::vector<proto_models > proto_models_sequence; // they are not interlaced
  typedef typename proto_models_sequence::iterator proto_models_seq_iterator;

  // first constructor only gives the number of prototypes
  // they will be chosen at random
  kmeans( const ValItPairsPair& vipp, size_type init):seriesBegin_(boost::get<0>(vipp)), seriesEnd_(boost::get<1>(vipp)), seriesLength_(std::distance(boost::get<0>(*seriesBegin_), boost::get<1>(*seriesBegin_))), protosVals(seriesLength_*init) { 
    std::vector<it_iterator> protosIt(init);
    std::random_sample(seriesBegin_, seriesEnd_, protosIt.begin(), protosIt.end());
    init(boost::make_tuple(protosIt.begin(), protosIt.end()));
  }
  template<typename InitIt>
  kmeans( const ValItPairsPair& vipp, boost::tuple<InitIt,InitT> initP):seriesBegin_(boost::get<0>(vipp)), seriesEnd_(boost::get<1>(vipp)), seriesLength_(std::distance(boost::get<0>(*seriesBegin_), boost::get<1>(*seriesBegin_))), protosVals(seriesLength_*init)  { init(initP);  }
  
  

  compute(){}
 private:
  // takes a pair of iterators over pairs of iterators
  // I will interlace
  template<typename InitItTuple>
  void init(const InitItTuple& initP){ 
    set_proto_values();
    // clear protoModels
    proto_model_type emptyModel;
    for(typename proto_models_sequence::iterator it(protoModels.begin()); it != protoModels.end(); ++it){
      for(typename proto_models::iterator innerIt(it->begin()); innerIt!=it->end();++innerIt)
	{ *innerIt=emptyModel;}
    }
  }

  void set_proto_values(const InitItTuple& initP){ 
    typedef typename boost::tuples::element<0, InitItTuple>::type proto_it_iterator;
    // init protoValues must make a fonction of it, to call it later with a focntion iterator over models to update the values.
    size_type nbProtos(std::distance(boost::get<0>(initP),boost::get<1>(initP))); 
    protos_vals_iterator protosValsIt(protosVals_.begin()), protoValsInitIt(protoVals_.begin());
    for(proto_it_iterator it(boost::get<0>(initP)); it != boost::get<1>(initP)
	  ; ++it, ++protoValsInitIt, protoValsIt=protoValsInitIt){
	for(serie_iterator sit(boost::get<0>(*it)); sit!= boost::get<1>(*sit); ++sit, protoValsIt+=nbProtos){
	  *protoValsIt=*sit;
	}
    }
  }

  proto_models_seq_iterator nearest( const serie_iterators_pair& serie)const{
    const size_type nbP(protoModels.size());
    std::vector<value_type> distances(nbP,0.);
    protos_vals_iterator protosValsIt(protosVals_.begin());
    for(serie_iterator it(boost::get<0>(serie)); it != boost::get<1>(serie); ++it){
      for(typename std::vector<value_type>::iterator distIt(distances.begin()); distIt!=distances.end(); ++distIt, ++protoValsIt)
	{ const value_type delta(*protoValsIt-*it); distIt+=delta*delta;}
    }
      return protoModels.begin()+(std::distance(distances.begin()
						,std::max_element(distances.begin(),distances.end())));
  }
  
  void one_run(){
    for(it_iterator it(seriesBegin_); it!= seriesEnd_; ++it)
      { updateProto(*it, nearest(*it));}
    for(

};

#endif
