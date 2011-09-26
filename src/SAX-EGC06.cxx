#include <iterator>
#include <iostream>
#include<vector>
#include <functional>

#include <numeric>
#include <string>
#include <cmath>

#include <boost/utility.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/progress.hpp>


#include "feats/segmentations/models.hxx"
#include "feats/segmentations/segment.hxx"
#include "feats/segmentations/binary_splits.hxx"
#include "feats/segmentations/merges.hxx"
#include "feats/segmentations/optimal_splits.hxx"
#include "feats/segmentations/utility.hxx"
#include "feats/segmentations/prototypes_sequences_models.hxx"
#include "feats/segmentations/prototypes_sequences_segmentations.hxx"
#include "feats/utility/na_reader.hxx"
#include "utility.hxx"

// g++ -std=c++0x SAX-EGC06.cxx -o SAXEGC06 -I../include -lstdc++ -Wall -Wno-sign-compare -O4

using namespace feats::segmentations;
using namespace feats::utility;


/**
prog de comparaison entre SAX et mes représentations symboliques de séries temporelles

Entrée: séries temporelles de m^eme taille, une par ligne

Paramètres: nb min segments, nb max segments, nb min symboles, nb max symboles

Sorties : nbSegs, NbSymboles, DPP/DSAX, DPP/DSBSR pour les nbs de segments et symboles demandés


Classes utilisées:

SaxParameters

with
inc_segments()
inc_prototypes()

gives


segments_size()
prototypes()


SBSR

il faut qqchose comme prototypes_segmentation mais pas besoin d'^etre aussi générique et l'optimisation des niveaux est un peu différente: il faut que je repasse à une dimension à partir des modèles


I will compute the n^2/2 DPP DSAX and DSBSR


**/

template<typename In>
struct SeriesComparator{
  typedef typename std::iterator_traits<In>::value_type serieIt;
  typedef typename std::iterator_traits<serieIt>::value_type value_type;
  typedef typename std::iterator_traits<serieIt>::difference_type size_type;
  SeriesComparator( In seriesBeg, In seriesEnd, size_type size):sb_(seriesBeg),se_(seriesEnd),  size_(size){}
  template<typename Out>
  Out operator()(Out o)const{
    for(In si(sb_); si!=se_; ++si){
      for( In sj(si);sj!=se_; ++sj,++o){
	*o=sqrt(dist(*si,*sj));
      }
    }
    return o;
  }
  value_type dist(serieIt s1, serieIt s2) const{
    value_type res(0.);
    for(size_type i(0); i!=size_; ++i, ++s1, ++s2){
      value_type tmp(*s1 - *s2);
      res+= tmp * tmp;
    }
    return res;
  }
  
  const In sb_, se_;
  const size_type size_;
};

template<typename Size, typename Out, typename ProtoOut>
Out saxBoundaries(Size s, Out o, ProtoOut po){
    *o=-std::numeric_limits<double>::max(); ++o;
  // normal distribution integral is not easy to compute
  switch (s){
  case 2: *o=0.; *po=-0.43; ++po; *po=.43; ++po; break;
  case 3: *o=-0.43; ++o;  *o=0.43; ++o; *po=-.67; ++po; *po=.0;++po;*po=.67;++po; break;
  case 4: *o=-.67; ++o;*o=0.; ++o;*o=.67; *po=-.84;++po; *po=-.25; ++po; *po=.25; ++po;*po=.84;++po; break; 
  case 5: *o=-0.84; ++o;*o=-0.25; ++o;*o=0.25; ++o;*o=0.84;
 *po=-0.97; ++po;*po=-0.43; ++po;  *po=0.; ++po;*po=0.43; ++po;*po=0.97; ++po;
  break;
  case 6: *o=-0.97; ++o;*o=-0.43; ++o;  *o=0.; ++o;*o=0.43; ++o;*o=0.97; 
*po=-1.07; ++po;*po=-0.57; ++po;  *po=-0.18; ++po;*po=0.18; ++po; *po=.57;++po; *po=1.07; ++po;
 break;
  case 7: *o=-1.07; ++o;*o=-0.57; ++o;  *o=-0.18; ++o;*o=0.18; ++o; *o=.57;++o; *o=1.07;
 *po=-1.15; ++po;*po=-0.67; ++po;  *po=-0.32; ++po;*po=0.; ++po;*po=0.32; ++po; *po=.67;++po; *po=1.15; ++po;
  break;
  case 8: *o=-1.15; ++o;*o=-0.67; ++o;  *o=-0.32; ++o;*o=0.; ++o;*o=0.32; ++o; *o=.67;++o; *o=1.15;
 *po=-1.22; ++po;*po=-0.76; ++po;  *po=-0.43; ++po;*po=-0.14; ++po;*po=0.14; ++po;*po=0.43; ++po; *po=.76;++po; *po=1.22; ++po;
  break;
  case 9: *o=-1.22; ++o;*o=-0.76; ++o;  *o=-0.43; ++o;*o=-0.14; ++o;*o=0.14; ++o;*o=0.43; ++o; *o=.76;++o; *o=1.22;  
 *po=-1.28; ++po;*po=-0.84; ++po;  *po=-0.52; ++po;*po=-0.25; ++po; *po=0.; ++po;*po=0.25; ++po;*po=0.52; ++po; *po=.84;++po; *po=1.28;  ++po;
break;
//  case 10: *o=-1.28; ++o;*o=-0.84; ++o;  *o=-0.52; ++o;*o=-0.25; ++o; *o=0.; ++o;*o=0.25; ++o;*o=0.52; ++o; *o=.84;++o; *o=1.28;  break;
  default: std::cerr<<"undefined nb of clusters in Sax\n";
  }
  ++o;
  return o;
}
template<typename Size, typename Out>
Out saxSizes(Size nbSegs, Size seqSize, Out o){

  for(Size remain=seqSize, currentSize=0; remain !=0; --nbSegs, remain-=currentSize,++o){
    currentSize=remain/nbSegs;
    *o=currentSize;
  }
  return o;
}
template<typename In, typename SeqSize, typename InSize, typename InProtosBoundaries, typename  InProtos, typename Out>
Out toSymboles(In sB,SeqSize seqSize, InSize siB, InProtosBoundaries pbB, InProtosBoundaries pbE, InProtos protosBegin, InProtos protosEnd, Out o, double& cost){
  typedef typename std::iterator_traits<InSize>::value_type size_type;
  typedef typename std::iterator_traits<In>::value_type value_type;
/*  
std::cerr<<"protos Boundaries:";
  std::copy(pbB, pbE, std::ostream_iterator<value_type>(std::cerr,"\t"));
  std::cerr<<'\n';
*/
  for(SeqSize iPoint(0);iPoint!=seqSize; ++o, ++siB){
    value_type tmp(0.),sum(0.),sum2(.0);
	linear0_prototype_level_model<InProtos,value_type> model(protosBegin, protosEnd);
    if(*siB){
      for(size_type i(0); i!= *siB; ++i, ++sB, ++iPoint){
	sum+=*sB;model+=*sB;
      }
      tmp=sum/(*siB);
    }
    *o=(std::distance(pbB, std::find_if(pbB, pbE, std::bind2nd(std::greater<value_type>(),tmp)))-1);
	cost+= model.cost();
//    std::cerr<<"proto="<<(std::distance(pbB, std::find_if(pbB, pbE, std::bind2nd(std::greater<value_type>(),tmp)))-1) <<" cost="<<model.cost()<<'\n';
  }
  return o;
}

template <typename InProtos, typename OutProtosB>
OutProtosB protosToBoundaries( InProtos pB, InProtos pE, OutProtosB o){
  if(pB!=pE){
    *o=-std::numeric_limits<typename std::iterator_traits<InProtos>::value_type>::max(); ++o;
  for( InProtos n(boost::next(pB)); n!=pE; ++n, ++pB,++o){
    *o=(*n + *pB)/2;
  }
  }
  return o;
}    
template<typename In, typename InSizes, typename InProtosB>
struct ProtosSeriesComparator{
  typedef typename std::iterator_traits<In>::value_type serieIt;
  typedef typename std::iterator_traits<InProtosB>::value_type value_type;
  typedef typename std::iterator_traits<InSizes>::value_type size_type;
  ProtosSeriesComparator( In seriesBeg, In seriesEnd, InSizes sizesBeg, InSizes sizesEnd, InProtosB protosBegin, InProtosB protosEnd):sb_(seriesBeg),se_(seriesEnd), sib_(sizesBeg), sie_(sizesEnd),  pb_(protosBegin),nbProtos_(std::distance(protosBegin,protosEnd)), distances_(nbProtos_,std::vector<value_type>(nbProtos_)){
    for(size_type i(0); i!=nbProtos_; ++i, ++protosBegin){
    InProtosB p1(protosBegin);
      for(size_type j(i); j!=nbProtos_; ++j, ++p1){
	if((j-i)<2) distances_[i][j]=distances_[j][i]= 0.;
	else{
	  value_type delta(*boost::next(protosBegin) - *boost::prior(p1));
	   distances_[i][j]=distances_[j][i]= delta*delta;
	}
      }
    }
	    
}
  template<typename Out>
  Out operator()(Out o)const{
    for(In si(sb_); si!=se_; ++si){
      for( In sj(si);sj!=se_; ++sj,++o){
	*o=sqrt(dist(*si,*sj));
      }
    }
    return o;
  }
  value_type dist(serieIt s1, serieIt s2) const{
    value_type res(0.);
    for(InSizes si(sib_); si!=sie_; ++si, ++s1, ++s2){
      res+= distances_[*s1][*s2]*(*si);
    }
    return res;
  }
  
  const In sb_, se_;
  const InSizes sib_, sie_;
  const InProtosB pb_;
  const size_type nbProtos_;
  std::vector<std::vector<value_type> > distances_;
};




/**

New distances mesure for symbolic representation.
We cannot use euclidian distance between reconstructions if we want to assert lower-bounding the real euclidian
distance bewteen the original time series. 
Simple apriori distance between symbols as defined in SAX does not give good results *because* the symbols are optimal hence no symbols are "wasted" between "useful" symbols, and those symbols "wasted" for representation are those essential for ensuring a useful distance.

Two kind of solutions to:
1) inserting "dummy" symbols only useful to ensure non-nul margin between symbols
2) use additional information about segments associated to symbols

I go for 2).

Usgin triangle inequality of distances, I get Dist(s1,s2)>= Dist(S1,S2)-Dist(s1,S1)-Dist(s2,S2), using modelling costs of S1 and S2 to compute Dist(s1,S1) and Dist(s2,S2), however, it is not good. However, we can do better with additional informations: Dist(s1,s2)>=max(Dist(s1,S2)-Dist(S2,s2), Dist(s2,S1)-Dist(S1,s1)). This is possible because the limited number if prototype levels allows to compute and store cost(Li) for each episod, for each prototype level Li. Even better, we do not have to choose only between the two triangle inequalities for the whole s1 and s2, we can switch at each episod. Dist(s1,s2)>=max(sqrt(sum(d²(s1,S2)))-sqrt(sum(d²(S2,s2))),sqrt(sum(d²(s2,S1)))-sqrt(sum(d²(S1,s1)))). The choices are independant so can can just maximize the partial distance function at each episod when computing distance.

 **/


template<typename DataIn, typename SeqSize, typename ProtoIn, typename SizesIn, typename OutCosts, typename OutDists, typename OutCorrectDist>
OutCosts symbDistances(DataIn seriesBeg,  DataIn seriesEnd,SeqSize seqSize, ProtoIn protosBeg,ProtoIn protosEnd, SizesIn sizesBeg, SizesIn sizesEnd, OutCosts o, OutDists oD, OutCorrectDist oCD){
  typedef typename std::iterator_traits<DataIn>::value_type serieIt;
  typedef typename std::iterator_traits<serieIt>::value_type value_type;
  typedef typename std::iterator_traits<SizesIn>::value_type size_type;

  typedef std::vector<size_type> symb_series_type;

  typedef std::vector<value_type> serie_type;
  typedef std::vector<serie_type> cont_type;
  typedef std::vector<cont_type> dists_cont_type;
  size_type nbSeries(std::distance(seriesBeg, seriesEnd)), nbSegs(std::distance(sizesBeg, sizesEnd)), nbProtos(std::distance(protosBeg, protosEnd));
  cont_type series(nbSeries, serie_type(nbSegs)), segsDist(nbSeries, serie_type(nbSegs));
  dists_cont_type segsDistsProtos(nbSeries, cont_type(nbSegs, serie_type(nbProtos)));

  //  std::vector<symb_series_type >
  linear0_prototype_level_model<ProtoIn> model(protosBeg, protosEnd);
  size_type nSerie(0);
  for(DataIn it(seriesBeg); it!=seriesEnd; ++it,++o, ++nSerie){
    SizesIn sib(sizesBeg);
    serieIt sb(*it);
    value_type serieCost(0.);
    size_type nSegment(0);
    for(size_type iPoint(0); iPoint!=seqSize; ++sib, ++nSegment){
      for(size_type i(0); i!= *sib; ++i, ++sb, ++iPoint){
	  model+=*sb;
      }
      series[nSerie][nSegment]=boost::get<0>(model());segsDist[nSerie][nSegment]=model.cost();
      serieCost+=model.cost();model-=model;
    }
    *o=sqrt(serieCost);
  }
  for(size_type nSerie(0); nSerie!=series.size(); ++nSerie){
    for(size_type nOtherSerie(nSerie); nOtherSerie!=series.size(); ++nOtherSerie, ++oD,++oCD){
      value_type dist(0.), correctDist(0.), distS1(0.), distS2(0.);
	SizesIn si(sizesBeg);
	for(size_type nSeg(0); nSeg!=series[nSerie].size(); ++nSeg, ++si){
	  value_type delta(series[nSerie][nSeg]-series[nOtherSerie][nSeg]);
	  dist+=delta*delta* (*si); correctDist +=delta*delta*(*si);distS1+=segsDist[nSerie][nSeg]; distS2+=segsDist[nOtherSerie][nSeg];
	}
	*oD=dist;*oCD=std::max(sqrt(correctDist)-sqrt(distS1)-sqrt(distS2),0.);
    }
  }
  return o;
}

template<typename RefDIn, typename SymbSymbDIn, typename SymbDIn>
typename std::iterator_traits<RefDIn>::value_type symbDistRatio(RefDIn rdB,SymbSymbDIn ssdB, SymbDIn sdB, SymbDIn sdE){
  typename std::iterator_traits<RefDIn>::value_type res(0.);
  int i(0);
  for(SymbDIn si(sdB); si!=sdE; ++si){
    for(SymbDIn sj(si); sj!=sdE; ++sj, ++i,++rdB, ++ssdB){
       typename std::iterator_traits<RefDIn>::value_type tmp((*rdB!=0.)?(sqrt(std::max(*ssdB-*si-*sj,0.)))/(*rdB):1.);
       if(tmp>1.)std::cout<<"!!! error in lower bounding!!!"<<tmp<<' '<<*rdB<<' '<<*ssdB<<' '<<*si<<' '<<*sj<<'\n';
      res+=tmp;
    }
  }
  res/=i;
  return res;
}
template<typename RefDIn, typename SymbCorrectDIn>
typename std::iterator_traits<RefDIn>::value_type symbCorrectDistRatio(RefDIn rdB, RefDIn rdE,SymbCorrectDIn scdB){
  typename std::iterator_traits<RefDIn>::value_type res(0.);
  int i(0);
  while(rdB!=rdE){
    typename std::iterator_traits<RefDIn>::value_type tmp((*rdB!=0.)?(*scdB)/(*rdB):1.);
    if(tmp>1.)std::cout<<"!!! error in lower bounding!!!"<<tmp<<' '<<*rdB<<' '<<*scdB<<'\n';
    res+=tmp;
    ++rdB; ++scdB;++i;
    }
  res/=i;
  return res;
}

template<typename DataIn, typename SeqSize, typename ProtoBIn, typename ProtosIt, typename SizesIn, typename RealDistIn>
typename std::iterator_traits<RealDistIn>::value_type distancesRatio(DataIn seriesBeg,  DataIn seriesEnd,SeqSize seqSize, ProtoBIn protosBeg,ProtoBIn protosEnd, ProtosIt realProtosBegin, ProtosIt realProtosEnd, SizesIn sizesBeg, SizesIn sizesEnd, RealDistIn realDistBeg){
  typedef typename std::iterator_traits<DataIn>::value_type serieIt;
  typedef typename std::iterator_traits<serieIt>::value_type value_type;
  typedef typename std::iterator_traits<SizesIn>::value_type size_type;

  typedef std::vector<size_type> symb_series_type;
  typedef std::vector<symb_series_type > cont_type;
  cont_type symbolic;
	double costs(0.);
  for(DataIn it(seriesBeg); it!=seriesEnd; ++it){
    symbolic.push_back(symb_series_type());
    toSymboles(*it, seqSize, sizesBeg, protosBeg, protosEnd, realProtosBegin, realProtosEnd, std::back_inserter(symbolic.back()),costs);
  }
	std::cerr<<"costs="<<costs<<'\n';
  return 0.;

}

template<template<typename, typename,typename>class ProtoModelT
	 , template<typename, typename>class Seg
	 , typename In>
void tester(const std::string& name, boost::tuple<In, In> dataBegEnd , const boost::tuple<int,int>& nbSegsMinMax, const boost::tuple<int,int>& nbProtosMinMax, int seqSize){

  typedef prototypes_sequence_segmentations<ProtoModelT, In, Seg,Seg> protos_segs_type;
  typedef typename protos_segs_type::prototypes_segmentation prototypes_segmentation ;

  SeriesComparator<In> seriesComp(boost::get<0>(dataBegEnd), boost::get<1>(dataBegEnd),seqSize);
  std::vector<typename protos_segs_type::data_type> refDists;
//  seriesComp(std::back_inserter(refDists));
 // protos_segs_type    protoSegs(dataBegEnd, seqSize);
    boost::posix_time::ptime startChrono(boost::posix_time::microsec_clock::local_time());
    for( int nbSegs(1); nbSegs!= boost::get<1>(nbSegsMinMax); ++nbSegs/*,protoSegs.free_segmentation().inc_segments()*/){
      if(nbSegs>=boost::get<0>(nbSegsMinMax)){
	for( int nbProtos(1); nbProtos != boost::get<1>(nbProtosMinMax); ++nbProtos/*,protoSegs.parameter_quantizer().inc_prototypes()*/){
	  
	if(nbProtos>=boost::get<0>(nbProtosMinMax)){
	  // output twice, calling optimize betwen the two
	  for(int nbIt(0); nbIt!=1 /*2*/; ++nbIt){
//	    if(nbIt==1){ protoSegs.optimize(); std::cout<<"Optimized";}
	    boost::posix_time::ptime stopChrono(boost::posix_time::microsec_clock::local_time());
	    std::cout<<name<<'\t'<<boost::posix_time::to_simple_string(stopChrono-startChrono)<<'\t'<<nbProtos<<'\t'<<nbSegs<<'\t';
//	    std::cout<<protoSegs.get_segmentation().cost();

	    std::vector<typename protos_segs_type::size_type> sizes;
/*
	    {
	      std::vector<typename protos_segs_type::proto_segment> res;
	      protoSegs.get_segmentation().segments(std::back_inserter(res));
	      // get segments size
	      for(int i(0); i!=res.size(); ++i)
		{sizes.push_back(std::distance(res[i].begin(), res[i].end()));}
	    }
*/
	    // get prototypes
	    std::vector<typename protos_segs_type::data_type> protos, protosBoundaries,realProtos;
	    /*
	    protoSegs.get_parameter_quantizer().prototypes(std::back_inserter(protos));
	    protosToBoundaries(protos.begin(), protos.end(), std::back_inserter(protosBoundaries));


	    std::cout<<'\t'<<distancesRatio(boost::get<0>(dataBegEnd), boost::get<1>(dataBegEnd),seqSize, protosBoundaries.begin(), protosBoundaries.end(), sizes.begin(), sizes.end(), refDists.begin());
	    std::vector<typename protos_segs_type::data_type> sDists, sCosts, correctDists;
	    symbDistances(boost::get<0>(dataBegEnd), boost::get<1>(dataBegEnd),seqSize, protos.begin(), protos.end(), sizes.begin(), sizes.end(), std::back_inserter(sCosts), std::back_inserter(sDists), std::back_inserter(correctDists));

	    std::cout<<'\t'<<symbCorrectDistRatio(refDists.begin(),refDists.end(),correctDists.begin());
	    //	    std::cout<<'\t'<<symbDistRatio(refDists.begin(),sDists.begin(),sCosts.begin(),sCosts.end());
	    */
	    saxBoundaries(nbProtos, std::back_inserter(protosBoundaries),std::back_inserter(realProtos));
	    saxSizes(nbSegs, seqSize, std::back_inserter(sizes));
	    std::cout<<'\t'<<distancesRatio(boost::get<0>(dataBegEnd), boost::get<1>(dataBegEnd),seqSize, protosBoundaries.begin(), protosBoundaries.end(), realProtos.begin(),realProtos.end(),sizes.begin(), sizes.end(), refDists.begin())<<std::endl;
	    // here compute distances and distances ratio to print result
	    startChrono+=(boost::posix_time::microsec_clock::local_time()-stopChrono); //substract time spend in op
	  }
	}
      }
    }
  }
}


 
int main(int argc, char* argv[]){
  typedef double data_type;
  bool badArgs(argc!=7);

  if(!badArgs){
    try{
      std::string const model(boost::lexical_cast<std::string>(argv[1]));
      std::string const algo(boost::lexical_cast<std::string>(argv[2]));
      int nbProtosMin(boost::lexical_cast<int>(argv[3]));
      int nbProtosMax(boost::lexical_cast<int>(argv[4]));
      int nbSegsMin(boost::lexical_cast<int>(argv[5]));
      int nbSegsMax(boost::lexical_cast<int>(argv[6]));
	     
  typedef std::vector<data_type> container_type;
  typedef container_type::difference_type size_type;
  typedef container_type::const_iterator c_iter_type;
  typedef std::vector<container_type> series_set_type;
  
  std::string lineBuffer;
  series_set_type series;

  while(std::getline(std::cin, lineBuffer)){
    series.push_back(container_type());
    std::istringstream tmp(lineBuffer);
    std::istream_iterator<na_reader<data_type> > b(tmp),e ;
    std::copy(b,e,std::back_inserter(series.back()));
  }
  std::cerr<<"reading "<<series.size()<<" series of size"<<series.back().size()<<'\n';
  int seqSize(series.back().size());  
  typedef std::vector<c_iter_type> sit_type;
  typedef sit_type::iterator set_it;
  typedef boost::tuple<set_it, set_it> sit_it_pair;
  bool shouldCenter(false), shouldReduce(false);
  if(shouldCenter){
    for( series_set_type::iterator sIt(series.begin()); sIt!=series.end(); ++sIt){
      center(sIt->begin(), sIt->end(), sIt->begin());
    }
  }
  if(shouldReduce){
    for( series_set_type::iterator sIt(series.begin()); sIt!=series.end(); ++sIt){
      reduce(sIt->begin(), sIt->end(), sIt->begin());
    }
  }
  /* must transpose the matrix for efficiency reasons*/
/*
  {
    series_set_type transposedSeries(series.front().size(),container_type(series.size()));
    for(size_type i(0); i!=series.size(); ++i){
      for(size_type j(0); j!=series.front().size(); ++j){
	transposedSeries[j][i]=series[i][j];
      }
    }
    series.swap(transposedSeries);
  }
*/
  sit_type seriesIts;
  for( series_set_type::const_iterator it(series.begin()); it != series.end(); ++it)
    { seriesIts.push_back(it->begin()); }
seqSize=series.size();
  //!!! FIXME verifier SeqSize est-ce vraiment le bon ???
	if((nbSegsMin>=1) && (nbSegsMax <= seqSize+1)){
	  sit_it_pair serieBegEnd(seriesIts.begin(), seriesIts.end());
	  boost::tuple<int, int>nbSegsMinMax(nbSegsMin,nbSegsMax);
	  boost::tuple<int, int>nbProtosMinMax(nbProtosMin,nbProtosMax);
	  if(model=="Linear0"){
	    if(algo == "TopDown"){
	      tester<linear0_prototype_level_sequence_model, binary_splits>(model+algo, serieBegEnd, nbSegsMinMax, nbProtosMinMax, seqSize);
	    }else if(algo=="Optimal"){
	      tester<linear0_prototype_level_sequence_model, optimal_splits>(model+algo, serieBegEnd, nbSegsMinMax, nbProtosMinMax, seqSize);
	    }else{std::cerr<<"Algo "<<algo<<" not implemented (yet?)\n"; badArgs=true;}
	  }else {std::cerr<<"Model "<<model<<"not implemented (yet?)\n"; badArgs=true;}
	}else{std::cerr<<"invalid nbseg interval, must be between 1 and "<<seqSize+1<<'\n';}
    }catch(boost::bad_lexical_cast &)
      {badArgs=true; }
  }
  if(badArgs){
    std::cerr<<"Usage: "<<argv[0]<<" Linear0|Linear1 TopDown|Optimal nProtosMin nProtosMax NbSegsMin nbSegsMax \n"
	     <<"Where 1 < nbSegsMin <= nbSegsMax < (nb points in input stream) is the number of segments\n"
	     <<"Where 1 < nbProtosMin <= nbProtosMax < nbSegsMin is the number of prototypes\n"
	     <<" segmentations are written to output stream in the following format:\n"
	     <<"Name\\ttime elapsed (s)\\tnb of segments\\tSSE\\tresulting segmentation\n";
  }
  return 0;

}
