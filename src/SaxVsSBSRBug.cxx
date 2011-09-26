#include <iterator>
#include <iostream>
#include <vector>
#include <functional>
#include <ext/algorithm> //std::copy_n
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

#include "utility.hxx"

// g++ -std=c++0x SaxVsSBSR.cxx -o SaxVsSBSR -I../include -lstdc++ -Wall -Wno-sign-compare

using namespace feats::segmentations;


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

template<typename Size, typename Out>
Out saxBoundaries(Size s, Out o){
  *o=-std::numeric_limits<double>::max(); ++o;
  // normal distribution integral is not easy to compute
  switch (s){
  case 2: *o=0.; break;
  case 3: *o=-0.43; ++o;  *o=0.43; ++o;  break;
  case 4: *o=-.67; ++o;*o=0.; ++o;*o=.67;break; 
  case 5: *o=-0.84; ++o;*o=-0.25; ++o;*o=0.25; ++o;*o=0.84;  break;
  case 6: *o=-0.97; ++o;*o=-0.43; ++o;  *o=0.; ++o;*o=0.43; ++o;*o=0.97;  break;
  case 7: *o=-1.07; ++o;*o=-0.57; ++o;  *o=-0.18; ++o;*o=0.18; ++o; *o=.57;++o; *o=1.07;  break;
  case 8: *o=-1.15; ++o;*o=-0.67; ++o;  *o=-0.32; ++o;*o=0.; ++o;*o=0.32; ++o; *o=.67;++o; *o=1.15;  break;
  case 9: *o=-1.22; ++o;*o=-0.76; ++o;  *o=-0.43; ++o;*o=-0.14; ++o;*o=0.14; ++o;*o=0.43; ++o; *o=.76;++o; *o=1.22;  break;
  case 10: *o=-1.28; ++o;*o=-0.84; ++o;  *o=-0.52; ++o;*o=-0.25; ++o; *o=0.; ++o;*o=0.25; ++o;*o=0.52; ++o; *o=.84;++o; *o=1.28;  break;
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

// in seriesIt (episode begin), nb of series, episode size, protoBoundaries begin & end, output of symboles number 
template<typename In, typename SeqSize, typename EpisodeSize, typename InProtosBoundaries, typename Out>
Out toSymboles(In sB,SeqSize seqSize, EpisodeSize s, InProtosBoundaries pbB, InProtosBoundaries pbE, Out o){
  typedef typename std::iterator_traits<In>::value_type seq_it_type;
  typedef typename std::iterator_traits<seq_it_type>::value_type value_type;
  std::cerr<<"protos Boundaries:";
  std::copy(pbB, pbE, std::ostream_iterator<value_type>(std::cerr,"\t"));
  std::cerr<<'\n';
  std::vector<value_type> tmp(seqSize,0.);
  for(EpisodeSize i(0); i!=s; ++i, ++sB){
    seq_it_type sit(*sB);
    for(SeqSize iS(0);iS!=seqSize;  ++iS,++sit)
      { tmp[iS]+=*sit; }
  }
  for(SeqSize iS(0);iS!=seqSize;  ++iS,++o){ 
    tmp[iS] /= s;
    *o=(std::distance(pbB, std::find_if(pbB, pbE, std::bind2nd(std::greater<value_type>(),tmp[iS])))-1);    
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

// retourne les distances entre pB et les symboles de pB à pE
template <typename InProtosB, typename InProtosB2, typename OutProtosDists>
OutProtosDists protosBoundariesToDists( InProtosB pB, InProtosB2 pE, OutProtosDists o){
  typedef typename std::iterator_traits<InProtosB>::value_type value_type;
  if(pB!=pE){
    ++pB;
    *o=0.;++o;
    if(pB!=pE){
      *o=0.;++o;
      for( InProtosB n(boost::next(pB)); n!=pE; ++n, ++o){
	value_type delta(*n - *pB);
	*o=delta*delta;
      }
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
    int i(0),j(0);
    for(In si(sb_); si!=se_; ++si,++i){
      j=i;
      for( In sj(si);sj!=se_; ++sj,++o,++j){
	std::cerr<<"dist symb entre "<<i<<" et "<<j<<'\n';
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

#ifdef OLD
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

  std::vector<symb_series_type >
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
#endif
#ifdef BUG
template<typename DataIn, typename SeqSize, typename ProtoIn, typename SizesIn, typename OutCosts, typename OutWorstDist, typename OutGlobalDist, typename OutLocalDist, typename OutEpisodesDist>
OutCosts symbDistances(DataIn seriesBeg,  DataIn seriesEnd,SeqSize seqSize, ProtoIn protosBeg,ProtoIn protosEnd, SizesIn sizesBeg, SizesIn sizesEnd, OutCosts o, OutWorstDist osd, OutGlobalDist ogd, OutLocalDist old, OutEpisodesDist oed){
  typedef typename std::iterator_traits<DataIn>::value_type serieIt;
  typedef typename std::iterator_traits<serieIt>::value_type value_type;
  typedef typename std::iterator_traits<SizesIn>::value_type size_type;

  typedef std::vector<value_type> vals_vect;
  typedef std::vector<vals_vect> vals_vect_vect;
  typedef std::vector<vals_vect_vect> vals_vect_vect_vect;
  typedef std::vector<vals_vect_vect> vals_vect_vect_vect_vect;

  typedef std::vector<size_type> symb_vect;
  typedef std::vector<symb_vect> symb_vect_vect;

  int nbProtos(std::distance(protosBeg, protosEnd)), nbPoints(std::distance(seriesBeg, seriesEnd)), nbEpisodes(std::distance(sizesBeg, sizesEnd));
  
  vals_vect_vect globalDists(nbProtos, vals_vect(nbProtos, std::numeric_limits<double>::max()));//[nbProto1][nbProto2]
  vals_vect_vect_vect localDists(seqSize, globalDists);//[ns][nbProto1][nbProto2]
  vals_vect_vect_vect episodesDists(seqSize, vals_vect_vect(nbEpisodes,vals_vect(nbProtos,std::numeric_limits<double>::max())));  // [ns][nE][nbProto2]

  vals_vect costs(seqSize,0.);
  vals_vect_vect realTotalDists(seqSize, vals_vect(seqSize, 0.));

  symb_vect_vect symbolic(nbEpisodes, symb_vect(seqSize,0));

  linear0_prototype_level_sequence_model<ProtoIn> model(protosBeg, protosEnd, seqSize);
  std::cerr<< "nbProtos:"<<nbProtos<<" nbPoints:"<< nbPoints<<" nbEpisodes:"<<nbEpisodes<<"seqSize:"<<seqSize<<'\n';
  
  DataIn it(seriesBeg); 
  size_type nE(0);
  for(SizesIn sib(sizesBeg); sib!=sizesEnd; ++sib, ++nE){
    for(size_type i(0); i!=*sib; ++i, ++it)
      {  
	model+=*it; 
	serieIt jIt(*it);
	for( size_type j(0); j!=seqSize; ++j, jIt){
	  serieIt kIt(jIt); 
	  for( size_type k(j); k!=seqSize; ++k,++kIt){
	    value_type tmp(*jIt - *kIt);
	    realTotalDists[j][k]+= tmp*tmp;
	  }
	}
      }
    for( size_type i(0); i!=seqSize; ++i){
      model.mean_costs(i, episodesDists[i][nE].begin());
      costs[i]+=model.cost(i);
    }
    model.protos_number(symbolic[nE].begin());
    model-=model;
  }

  // completer la matrice sym des distances réelles
  for( size_type j(0); j!=seqSize; ++j){
    for( size_type k(j); k!=seqSize; ++k)
      { realTotalDists[k][j]=realTotalDists[j][k]; }
  }

  std::cerr<<"aggregations\n";
  // aggregations
  for(size_type i(0); i!=seqSize; ++i){
    for(size_type j(0); j!=nbEpisodes; ++j){
      for(size_type k(0); k!=nbProtos; ++k){
	// update min dist for serie i from symbol of episode j (symbolic[j][i]) to proto k
	//	std::cerr<<"min"<<(localDists[i][symbolic[j][i]][k])<<" and"<<episodesDists[i][j][k]<<" is";
	localDists[i][symbolic[j][i]][k]=std::min(localDists[i][symbolic[j][i]][k], episodesDists[i][j][k]);
	//std::cerr<<(localDists[i][symbolic[j][i]][k])<<"\n";
      }
    }
  }
  for(size_type i(0); i!=nbProtos; ++i){
    for(size_type j(0); j!=nbProtos; ++j){
      for(size_type k(0); k!=seqSize; ++k){
	globalDists[i][j]=std::min(globalDists[i][j], localDists[k][i][j]);
      }
    }
  }
  // calculer les dists theoriques (worst) à la SAX
  vals_vect boundaries;
  protosToBoundaries( protosBeg, protosEnd, std::back_inserter(boundaries));
  vals_vect_vect worstDists(nbProtos, vals_vect(nbProtos,0.));
  size_type j(0);
  for(typename vals_vect::const_iterator i(boundaries.begin());  i!= boundaries.end(); ++i,++j){
    protosBoundariesToDists(i, boundaries.end(), worstDists[j].begin()+j);
  }
  for(size_type i(0); i!=nbProtos; ++i){
    for(size_type j(i); j!=nbProtos; ++j){
      worstDists[j][i]=worstDists[i][j];
    }
  }

  for(size_type nSerie(0); nSerie !=seqSize; ++nSerie){
    std::cerr<<"\nserie "<<nSerie<<'\n';
    SizesIn sIt(sizesBeg);
    for(size_type nSeg(0); nSeg!=nbEpisodes; ++nSeg,++sIt){
      std::cerr<<"\nsegment"<<nSeg<<" of size"<<*sIt<<" proto"<<symbolic[nSeg][nSerie];
      ProtoIn pIt(protosBeg);
      for(size_type p(0); p!=nbProtos; ++p, ++pIt){
	std::cerr<<"mean costs to"<<*pIt<<" is"<<episodesDists[nSerie][nSeg][p];
      }
    }
  }
  vals_vect_vect forDebug(nbPoints, vals_vect(seqSize));
  DataIn si(seriesBeg); 

  for(size_type i(0); i!=nbPoints; ++i, ++si){
    __gnu_cxx::copy_n(*si, seqSize, forDebug[i].begin());
  }
  for(size_type nSerie(0); nSerie!=seqSize; ++nSerie) {
    for(size_type nOtherSerie(nSerie); nOtherSerie!=seqSize; ++nOtherSerie, ++osd,++ogd, ++old, ++oed){
      {
	std::cerr<<"\nreal dist between"<<nSerie <<" and "<<nOtherSerie<<" ";
	value_type res(0.);
	for(size_type i(0); i!=nbPoints;++i){
	  value_type tmp(forDebug[i][nSerie]-forDebug[i][nOtherSerie]);
	  tmp*=tmp;
	  res+=tmp;
	  std::cerr<<" "<<forDebug[i][nSerie]<<" and "<<forDebug[i][nOtherSerie]<<":"<<tmp<<" ";
	}
	std::cerr<<"total:"<<res<<'\n';
      }

      value_type worstD(0.), globalD(0.), localD(0.), episodesD(0.);
      SizesIn si(sizesBeg);
      for(size_type nSeg(0); nSeg!=nbEpisodes; ++nSeg, ++si){
	size_type p1(symbolic[nSeg][nSerie]), p2(symbolic[nSeg][nOtherSerie]);
	worstD+=std::max(worstDists[p1][p2]*(*si), 0.);
	globalD+=std::max(std::max(globalDists[p1][p2]- episodesDists[nOtherSerie][nSeg][p2]
				   , globalDists[p2][p1]- episodesDists[nSerie][nSeg][p1])*(*si)
			  ,0.);
	localD+=std::max(std::max(localDists[nSerie][p1][p2]- episodesDists[nOtherSerie][nSeg][p2]
				  , localDists[nOtherSerie][p2][p1]- episodesDists[nSerie][nSeg][p1])*(*si)
			 ,0.);
	std::cerr<<"episode local cost"<<episodesDists[nSerie][nSeg][p2]<<"-"<<episodesDists[nOtherSerie][nSeg][p2]<<"="<<episodesDists[nSerie][nSeg][p2]- episodesDists[nOtherSerie][nSeg][p2] <<" and "<<episodesDists[nOtherSerie][nSeg][p1]<<"-"<<episodesDists[nSerie][nSeg][p1]<<"="<<episodesDists[nOtherSerie][nSeg][p1]- episodesDists[nSerie][nSeg][p1]<<" max is "<<
	  std::max(std::max(episodesDists[nSerie][nSeg][p2]- episodesDists[nOtherSerie][nSeg][p2]
			    , episodesDists[nOtherSerie][nSeg][p1]- episodesDists[nSerie][nSeg][p1])*(*si)
		   ,0.);
	    
	episodesD+=std::max(std::max(episodesDists[nSerie][nSeg][p2]- episodesDists[nOtherSerie][nSeg][p2]
				     , episodesDists[nOtherSerie][nSeg][p1]- episodesDists[nSerie][nSeg][p1])*(*si)
			    ,0.);
      }
      std::cerr<<"worstD:"<<worstD<<" globalD:"<<globalD<<" localD:"<<localD<<" episodesD:"<<episodesD<<" trueD"<<realTotalDists[nSerie][nOtherSerie]<<'\n';
      if(worstD>globalD)std::cerr<<"!!!!!!!!!!!! worstD > globalD\n";
      if(globalD>localD)std::cerr<<"!!!!!!!!!!!!  globalD > localD\n";
      if(localD>episodesD)std::cerr<<"!!!!!!!!!!!!  localD>episodesD\n";
      if(episodesD>realTotalDists[nSerie][nOtherSerie])std::cerr<<"!!!!!!!!!!!!  episodesD> true Dist\n";

      *osd=worstD;*ogd=globalD; *old= localD; *oed=episodesD;
    }
  }
  
  return o;
}
#endif
/*
  en fait les distances sont calsulées par rapport aux distances d'avec les frontières entre classes de valeurs et non les valeurs prototypiques
*/

template<typename DataIn, typename SeqSize, typename ProtoIn, typename SizesIn, typename OutCosts, typename OutWorstDist, typename OutGlobalDist, typename OutLocalDist, typename OutEpisodesDist>
OutCosts symbDistances(DataIn seriesBeg,  DataIn seriesEnd,SeqSize seqSize, ProtoIn protosBeg,ProtoIn protosEnd, SizesIn sizesBeg, SizesIn sizesEnd, OutCosts o, OutWorstDist osd, OutGlobalDist ogd, OutLocalDist old, OutEpisodesDist oed){
  typedef typename std::iterator_traits<DataIn>::value_type serieIt;
  typedef typename std::iterator_traits<serieIt>::value_type value_type;
  typedef typename std::iterator_traits<SizesIn>::value_type size_type;

  typedef std::vector<value_type> vals_vect;
  typedef std::vector<vals_vect> vals_vect_vect;
  typedef std::vector<vals_vect_vect> vals_vect_vect_vect;
  typedef std::vector<vals_vect_vect> vals_vect_vect_vect_vect;

  typedef std::vector<size_type> symb_vect;
  typedef std::vector<symb_vect> symb_vect_vect;

  int nbProtos(std::distance(protosBeg, protosEnd)), nbPoints(std::distance(seriesBeg, seriesEnd)), nbEpisodes(std::distance(sizesBeg, sizesEnd));
  
  vals_vect_vect globalDists(nbProtos, vals_vect(nbProtos, std::numeric_limits<double>::max()));//[nbProto1][nbProto2]
  vals_vect_vect_vect localDists(seqSize, globalDists);//[ns][nbProto1][nbProto2]
  vals_vect_vect_vect episodesDists(seqSize, vals_vect_vect(nbEpisodes,vals_vect(nbProtos,std::numeric_limits<double>::max())));  // [ns][nE][nbProto2]

  vals_vect costs(seqSize,0.);
  vals_vect_vect realTotalDists(seqSize, vals_vect(seqSize, 0.));

  symb_vect_vect symbolic(nbEpisodes, symb_vect(seqSize,0));
  {
    std::cerr<<"protos:";
    std::copy(protosBeg, protosEnd, std::ostream_iterator<value_type>(std::cerr," "));
    std::cerr<<"\nboundaries:";
    vals_vect boundaries;
    protosToBoundaries(protosBeg, protosEnd,std::back_inserter(boundaries));
    std::copy(boundaries.begin(), boundaries.end(), std::ostream_iterator<value_type>(std::cerr," "));
    std::cerr<<'\n';
    boundaries.erase(boundaries.begin()); // first boundaries (-Inf) is irrelevent 

    linear0_prototype_level_sequence_model<ProtoIn> modelB(boundaries.begin(), boundaries.end(), seqSize),model(protosBeg, protosEnd, seqSize);
    std::cerr<< "nbProtos:"<<nbProtos<<" nbPoints:"<< nbPoints<<" nbEpisodes:"<<nbEpisodes<<"seqSize:"<<seqSize<<'\n';
  
    DataIn it(seriesBeg); 
    size_type nE(0);
    for(SizesIn sib(sizesBeg); sib!=sizesEnd; ++sib, ++nE){
      for(size_type i(0); i!=*sib; ++i, ++it)
	{  
	  modelB+=*it; model+=*it; 
	  serieIt jIt(*it);
	  for( size_type j(0); j!=seqSize; ++j, ++jIt){
	    serieIt kIt(*it); 
	    for( size_type k(0); k!=seqSize; ++k,++kIt){
	      value_type tmp(*jIt - *kIt);
	      //	      std::cerr<<"series "<<j<<":"<<*jIt<<" and "<<k<<":"<<*kIt<<"->+ "<<tmp*tmp<<'\n';
	      realTotalDists[j][k]+= tmp*tmp;

	    }
	  }
	}
      model.protos_number(symbolic[nE].begin());
      for( size_type i(0); i!=seqSize; ++i){
	modelB.mean_costs(i, episodesDists[i][nE].begin());
	// insert cost 0. for self distance
	for( size_type k(nbProtos-1); k!=symbolic[nE][i]; --k)
	  { episodesDists[i][nE][k]=episodesDists[i][nE][k-1];}
	episodesDists[i][nE][symbolic[nE][i]]=0.;
	std::cerr<<"episodes Dists for episodes nb:"<<nE<<" of series:"<<i<<'\n';
	std::copy(episodesDists[i][nE].begin(), episodesDists[i][nE].end(), std::ostream_iterator<value_type>(std::cerr," "));
	std::cerr<<"\n";
	costs[i]+=model.cost(i);
      }
      modelB -= modelB; model-=model;
    }
  }
  std::cerr<<"aggregations\n";
  // aggregations
  for(size_type i(0); i!=seqSize; ++i){
    for(size_type j(0); j!=nbEpisodes; ++j){
      for(size_type k(0); k!=nbProtos; ++k){
	// update min dist for serie i from symbol of episode j (symbolic[j][i]) to proto k
		std::cerr<<"min"<<(localDists[i][symbolic[j][i]][k])<<" and"<<episodesDists[i][j][k]<<" is";
	localDists[i][symbolic[j][i]][k]=std::min(localDists[i][symbolic[j][i]][k], episodesDists[i][j][k]);
	std::cerr<<(localDists[i][symbolic[j][i]][k])<<"\n";
      }
    }
  }
  for(size_type i(0); i!=nbProtos; ++i){
    for(size_type j(0); j!=nbProtos; ++j){
      for(size_type k(0); k!=seqSize; ++k){
	globalDists[i][j]=std::min(globalDists[i][j], localDists[k][i][j]);
      }
    }
  }
  // calculer les dists theoriques (worst) à la SAX
  vals_vect boundaries;
  protosToBoundaries( protosBeg, protosEnd, std::back_inserter(boundaries));
  vals_vect_vect worstDists(nbProtos, vals_vect(nbProtos,0.));
  size_type j(0);
  std::cerr<<std::endl;
  for(typename vals_vect::const_iterator i(boundaries.begin());  i!= boundaries.end(); ++i,++j){
    protosBoundariesToDists(i, boundaries.end(), worstDists[j].begin()+j);
    /*    protosBoundariesToDists(i, boundaries.end(), std::ostream_iterator<value_type>(std::cerr," "));
    std::cerr<<'\n';
    std::copy(worstDists[j].begin(), worstDists[j].end(), std::ostream_iterator<value_type>(std::cerr," "));
    std::cerr<<std::endl;*/
  }
  /*
  std::cerr<<"worstDist before completion\n";
  for(size_type i(0); i!=nbProtos; ++i){
    for(size_type j(0); j!=nbProtos; ++j){
      std::cerr<<worstDists[i][j]<<' ';
    }
    std::cerr<<'\n';
  }
  */
  for(size_type i(0); i!=nbProtos; ++i){
    for(size_type j(i); j!=nbProtos; ++j){
      worstDists[nbProtos-i-1][nbProtos-j-1]=worstDists[i][j];
    }
  }
  /*
  std::cerr<<"\nafter\n";
  for(size_type i(0); i!=nbProtos; ++i){
    for(size_type j(0); j!=nbProtos; ++j){
      std::cerr<<worstDists[i][j]<<' ';
    }
    std::cerr<<'\n';
  }
  */
  /*
  for(size_type nSerie(0); nSerie !=seqSize; ++nSerie){
    std::cerr<<"\nserie "<<nSerie<<'\n';
    SizesIn sIt(sizesBeg);
    for(size_type nSeg(0); nSeg!=nbEpisodes; ++nSeg,++sIt){
      std::cerr<<"\nsegment"<<nSeg<<" of size"<<*sIt<<" proto"<<symbolic[nSeg][nSerie];
      ProtoIn pIt(protosBeg);
      for(size_type p(0); p!=nbProtos; ++p, ++pIt){
	std::cerr<<"mean costs to"<<*pIt<<" is"<<episodesDists[nSerie][nSeg][p];
      }
    }
  }
  */
  vals_vect_vect forDebug(nbPoints, vals_vect(seqSize));
  DataIn si(seriesBeg); 

  for(size_type i(0); i!=nbPoints; ++i, ++si){
    __gnu_cxx::copy_n(*si, seqSize, forDebug[i].begin());
  }
  std::cerr<<"series:\n";
  for(size_type i(0); i!=seqSize; ++i){
    for(size_type j(0); j!=nbPoints; ++j){
      std::cerr<<forDebug[j][i]<<' ';
    }
    std::cerr<<std::endl;
  }
    
  for(size_type nSerie(0); nSerie!=seqSize; ++nSerie) {
    for(size_type nOtherSerie(nSerie); nOtherSerie!=seqSize; ++nOtherSerie, ++osd,++ogd, ++old, ++oed){
      {
	std::cerr<<"\nreal dist between"<<nSerie <<" and "<<nOtherSerie<<" ";
	value_type res(0.);
	for(size_type i(0); i!=nbPoints;++i){
	  value_type tmp(forDebug[i][nSerie]-forDebug[i][nOtherSerie]);
	  tmp*=tmp;
	  res+=tmp;
	  std::cerr<<" "<<forDebug[i][nSerie]<<" and "<<forDebug[i][nOtherSerie]<<":"<<tmp<<" ";
	}
	std::cerr<<"total:"<<res<<'\n';
      }

      value_type worstD(0.), globalD(0.), localD(0.), episodesD(0.);
      SizesIn si(sizesBeg);
      size_type i(0);
      for(size_type nSeg(0); nSeg!=nbEpisodes; ++nSeg, ++si){
	size_type p1(symbolic[nSeg][nSerie]), p2(symbolic[nSeg][nOtherSerie]);
	value_type res(0.);
	for(size_type j(0); j!=*si; ++j,++i){
	  value_type tmp(forDebug[i][nSerie]-forDebug[i][nOtherSerie]);
	  res+=tmp*tmp;
	}
	std::cerr<<"episode:"<<nSeg<<" of size:"<<*si<<" and cost:"<<res<<"\n";
	worstD+=std::max(worstDists[p1][p2] *(*si), 0.);
	globalD+=std::max(std::max(globalDists[p1][p2], globalDists[p2][p1])*(*si),0.);
	localD+=std::max(std::max(localDists[nSerie][p1][p2], localDists[nOtherSerie][p2][p1])*(*si),0.);
	std::cerr<<"episode local cost"<<episodesDists[nSerie][nSeg][p2]<<" and "<<episodesDists[nOtherSerie][nSeg][p1]
		 <<" max is "<<std::max(std::max(episodesDists[nSerie][nSeg][p2], episodesDists[nOtherSerie][nSeg][p1])*(*si),0.);
	    
	episodesD+=std::max(std::max(episodesDists[nSerie][nSeg][p2], episodesDists[nOtherSerie][nSeg][p1])*(*si),0.);
      }
      std::cerr<<"worstD:"<<worstD<<" globalD:"<<globalD<<" localD:"<<localD<<" episodesD:"<<episodesD<<" trueD"<<realTotalDists[nSerie][nOtherSerie]<<'\n';
      if(worstD>globalD)std::cerr<<"!!!!!!!!!!!! worstD > globalD\n";
      if(globalD>localD)std::cerr<<"!!!!!!!!!!!!  globalD > localD\n";
      if(localD>episodesD)std::cerr<<"!!!!!!!!!!!!  localD>episodesD\n";
      if(episodesD>realTotalDists[nSerie][nOtherSerie])std::cerr<<"!!!!!!!!!!!!  episodesD> true Dist\n";

      *osd=worstD;*ogd=globalD; *old= localD; *oed=episodesD;
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

template<typename DataIn, typename SeqSize, typename ProtoBIn, typename SizesIn, typename RealDistIn>
typename std::iterator_traits<RealDistIn>::value_type distancesRatio(DataIn seriesBeg,  DataIn seriesEnd,SeqSize seqSize, ProtoBIn protosBeg,ProtoBIn protosEnd, SizesIn sizesBeg, SizesIn sizesEnd, RealDistIn realDistBeg){
  typedef typename std::iterator_traits<DataIn>::value_type serieIt;
  typedef typename std::iterator_traits<serieIt>::value_type value_type;
  typedef typename std::iterator_traits<SizesIn>::value_type size_type;

  typedef std::vector<size_type> symb_series_type;
  typedef std::vector<symb_series_type > cont_type;
  cont_type symbolic;
  //  std::cerr<<"in distancesRatio dist seriesBeg, end:"<<std::distance(seriesBeg, seriesEnd)<<'\n';
  DataIn iD(seriesBeg);
  for( SizesIn iS(sizesBeg); iS != sizesEnd; iD+=*iS, ++iS){
    symbolic.push_back(symb_series_type());
    toSymboles(iD, seqSize, *iS, protosBeg, protosEnd, std::back_inserter(symbolic.back()));
  }
  // les symboles sont tranposées aussi en sequences: re transposition
  {
    cont_type tmp(symbolic.front().size(), symb_series_type(symbolic.size()));
    for(size_type i(0); i!=symbolic.size(); ++i){
      for(size_type j(0); j!=symbolic.front().size(); ++j){
	tmp[j][i]=symbolic[i][j];
      }
    }
    symbolic.swap(tmp);
  }

  typedef std::vector<typename symb_series_type::const_iterator> sit_cont;
  sit_cont symbIt;
  for(int i(0); i!=symbolic.size(); ++i)
    { symbIt.push_back(symbolic[i].begin());}

  ProtosSeriesComparator<typename sit_cont::const_iterator, SizesIn, ProtoBIn> comp(symbIt.begin(), symbIt.end(), sizesBeg, sizesEnd, protosBeg, protosEnd);
  std::vector<value_type> dists;
  comp(std::back_inserter(dists));
  value_type res(0.);
  for(size_type n(0); n!=dists.size(); ++n, ++realDistBeg){
    std::cerr<<"serie "<<n<<" dist "<<dists[n]<<" real: "<<*realDistBeg<<'\n';
    res+=(*realDistBeg!=0.)? dists[n]/(*realDistBeg):1.;
  }
  res/=dists.size();
  return res;

}

template<template<typename, typename,typename>class ProtoModelT
	 , template<typename, typename>class Seg
	 , typename In>
void tester(const std::string& name, boost::tuple<In, In> dataBegEnd , const boost::tuple<int,int>& nbSegsMinMax, const boost::tuple<int,int>& nbProtosMinMax, int seqSize){

  typedef prototypes_sequence_segmentations<ProtoModelT, In, Seg,Seg> protos_segs_type;
  typedef typename protos_segs_type::prototypes_segmentation prototypes_segmentation ;

  SeriesComparator<In> seriesComp(boost::get<0>(dataBegEnd), boost::get<1>(dataBegEnd),seqSize);
  std::vector<typename protos_segs_type::data_type> refDists;
  seriesComp(std::back_inserter(refDists));
  protos_segs_type    protoSegs(dataBegEnd, seqSize);
  boost::posix_time::ptime startChrono(boost::posix_time::microsec_clock::local_time());
  for( int nbSegs(1); nbSegs!= boost::get<1>(nbSegsMinMax); ++nbSegs,protoSegs.free_segmentation().inc_segments()){
    if(nbSegs>=boost::get<0>(nbSegsMinMax)){
      for( int nbProtos(1); nbProtos != boost::get<1>(nbProtosMinMax); ++nbProtos,protoSegs.parameter_quantizer().inc_prototypes()){
	  
	if(nbProtos>=boost::get<0>(nbProtosMinMax)){
	  // output twice, calling optimize betwen the two
	  for(int nbIt(0); nbIt!=2; ++nbIt){
	    if(nbIt==1){ protoSegs.optimize(); std::cout<<"Optimized";}
	    boost::posix_time::ptime stopChrono(boost::posix_time::microsec_clock::local_time());
	    std::cout<<name<<'\t'<<boost::posix_time::to_simple_string(stopChrono-startChrono)<<'\t'<<nbProtos<<'\t'<<nbSegs<<'\t';
	    std::cout<<protoSegs.get_segmentation().cost();

	    std::vector<typename protos_segs_type::size_type> sizes;
	    {
	      std::vector<typename protos_segs_type::proto_segment> res;
	      protoSegs.get_segmentation().segments(std::back_inserter(res));
	      // get segments size
	      for(int i(0); i!=res.size(); ++i)
		{sizes.push_back(std::distance(res[i].begin(), res[i].end()));}
	    }
	    // get prototypes
	    std::vector<typename protos_segs_type::data_type> protos, protosBoundaries;
	    protoSegs.get_parameter_quantizer().prototypes(std::back_inserter(protos));
	    protosToBoundaries(protos.begin(), protos.end(), std::back_inserter(protosBoundaries));


	    //	    std::cout<<'\t'<<distancesRatio(boost::get<0>(dataBegEnd), boost::get<1>(dataBegEnd),seqSize, protosBoundaries.begin(), protosBoundaries.end(), sizes.begin(), sizes.end(), refDists.begin());
	    std::vector<typename protos_segs_type::data_type> sCosts, worstDists, globalDists, localDists, episodesDists;
	    symbDistances(boost::get<0>(dataBegEnd), boost::get<1>(dataBegEnd),seqSize, protos.begin(), protos.end(), sizes.begin(), sizes.end(), std::back_inserter(sCosts), std::back_inserter(worstDists), std::back_inserter(globalDists),std::back_inserter(localDists),std::back_inserter(episodesDists));

	    //	    std::cout<<'\t'<<symbCorrectDistRatio(refDists.begin(),refDists.end(),correctDists.begin());
	    //	    std::cout<<'\t'<<symbDistRatio(refDists.begin(),sDists.begin(),sCosts.begin(),sCosts.end());
	    saxBoundaries(nbProtos, protosBoundaries.begin());
	    saxSizes(nbSegs, seqSize, sizes.begin());
	    //	    std::cout<<'\t'<<distancesRatio(boost::get<0>(dataBegEnd), boost::get<1>(dataBegEnd),seqSize, protosBoundaries.begin(), protosBoundaries.end(), sizes.begin(), sizes.end(), refDists.begin())<<std::endl;
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
	std::istream_iterator<data_type> b(tmp),e ;
	std::copy(b,e,std::back_inserter(series.back()));
      }
      std::cerr<<"reading "<<series.size()<<" series\n";
      int seqSize(series.back().size());  
      typedef std::vector<c_iter_type> sit_type;
      typedef sit_type::iterator set_it;
      typedef boost::tuple<set_it, set_it> sit_it_pair;
      bool shouldCenter(true), shouldReduce(true);
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
      {
	series_set_type transposedSeries(series.front().size(),container_type(series.size()));
	for(size_type i(0); i!=series.size(); ++i){
	  for(size_type j(0); j!=series.front().size(); ++j){
	    transposedSeries[j][i]=series[i][j];
	  }
	}
	series.swap(transposedSeries);
      }

      sit_type seriesIts;
      seqSize=series.front().size();
      for( series_set_type::const_iterator it(series.begin()); it != series.end(); ++it)
	{ seriesIts.push_back(it->begin()); }

      //!!! FIXME verifier SeqSize est-ce vraiment le bon ???
      std::cerr<<"seqSize:"<<seqSize<<'\n';
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
