#include <iterator>
#include <iostream>
#include<vector>
#include <functional>

#include <numeric>
#include <string>

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
#include "feats/segmentations/Uglhybrid.hxx"

#include "feats/segmentations/utility.hxx"

#include "feats/segmentations/prototypes_sequences_models.hxx"
#include "feats/segmentations/prototypes_sequences_segmentations.hxx"


#include "feats/utility/na_reader.hxx"
#include "utility.hxx"

// g++ -std=c++0x segmentations_proto_pop.cxx -o segmentations_proto_pop -I../include -lstdc++ -Wall -Wno-sign-compare -O4

//icc segmentations_proto_pop.cxx -I ../include/ -i_dynamic -o segmentations_proto_pop_ICC -O3 -Ob2 -tpp7 -ip

using namespace feats::segmentations;
using namespace feats::utility;



/***

programme de creation de representations symboliques d'ensembles de series temporelles synchrones. Les représentations symboliques de toutes les series partagent le meme decoupage en episodes et le meme alphabet symbolique. 

Format de données en entrées
1 ere ligne (entete):
date_debut\tdate_fin\tidCapteur1\t...idCapteurN
lignes suivantes
laDateDeDebut\tlaDateDeFin\tlaValeurDuCapteur1\t...laValeurDuCapteurN

format de sortie :
1 ligne (entete)
date_debut\tdate_fin\ttaille\tidCapteur1_niveau\tidCapteur1_sse\t...idCapteurN_niveau\tidCapteurN_sse



***/

template<template<typename, typename,typename>class ProtoModelT
	 , template<typename, typename>class Seg
	 , typename In, typename Dates>
void do_it(const std::string& name, boost::tuple<In, In> dataBegEnd , int seqSize, const boost::tuple<int,int>& nbSegsMinMax, const boost::tuple<int,int>& nbProtosMinMax, const std::string& format, const Dates& dates){


  typedef prototypes_sequence_segmentations<ProtoModelT, In, Seg, hybrid_splits> protos_segs_type;
  typedef typename protos_segs_type::prototypes_segmentation prototypes_segmentation ;
  protos_segs_type    protoSegs(dataBegEnd, seqSize);
  some_sequences_output_with_dates<prototypes_segmentation, Dates> o(format,seqSize,dates);

  std::cerr<<"in do_it with"<<std::distance(boost::get<0>(dataBegEnd), boost::get<1>(dataBegEnd))<<" sequences of size "<<seqSize<<std::endl;

    boost::posix_time::ptime startChrono(boost::posix_time::microsec_clock::local_time());
    for( int nbSegs(1); nbSegs!= boost::get<1>(nbSegsMinMax); ++nbSegs,protoSegs.free_segmentation().inc_segments()){
      if(nbSegs>=boost::get<0>(nbSegsMinMax)){
	for( int nbProtos(1); nbProtos != boost::get<1>(nbProtosMinMax); ++nbProtos,protoSegs.parameter_quantizer().inc_prototypes()){
	  
	if(nbProtos>=boost::get<0>(nbProtosMinMax)){
	  protoSegs.optimize();
	    boost::posix_time::ptime stopChrono(boost::posix_time::microsec_clock::local_time());
	    std::cerr<<name<<'\t'<<boost::posix_time::to_simple_string(stopChrono-startChrono)<<'\t'<<nbProtos<<'\t'<<nbSegs<<'\t';
	    std::cerr<<protoSegs.get_segmentation().cost()<<'\n';
	    o(protoSegs.get_segmentation(), std::cout)<<std::endl;
	    startChrono+=(boost::posix_time::microsec_clock::local_time()-stopChrono); //substract time spend in op
	}
	}
      }
    }
}


/***
    Small program for testing segmentations and models. Shows how static polymorphism allows to combine segment models and segmentation algorithms.
***/

int main(int argc, char* argv[]){
  typedef double data_type;
  bool badArgs(argc!=6);
	std::cerr<<"argc:"<<argc<<'\n';
  if(!badArgs){
    try{
      int const nbDateFields(boost::lexical_cast<int>(argv[1]));
	std::cerr<<nbDateFields<<'\n';
      std::string const model(boost::lexical_cast<std::string>(argv[2]));
      std::string const algo(boost::lexical_cast<std::string>(argv[3]));
	std::cerr<<algo<<'\n';
      int nbProtosMin(boost::lexical_cast<int>(argv[4]));
      int nbProtosMax(nbProtosMin+1);
      int nbSegsMin(boost::lexical_cast<int>(argv[5]));
      int nbSegsMax(nbSegsMin+1);

	std::cerr<<"nb date fields:"<<nbDateFields<<"model:"
	<<model<<" algo:"<<algo<<"nbP:"<<nbProtosMin
	<<"nbS:"<<nbSegsMin<<std::endl;
      typedef std::vector<data_type> container_type;
	typedef container_type::const_iterator c_iter_type;
	typedef container_type::difference_type size_type;

	typedef std::vector<std::string> strings_container;

  typedef std::vector<container_type> series_set_type;
  
  std::string lineBuffer;
  series_set_type series;
  strings_container dates,ids;

  std::getline(std::cin, lineBuffer);
  std::istringstream tmp(lineBuffer);
  std::istream_iterator<std::string>b(tmp),e ;
  std::copy(b,e,std::back_inserter(ids));
  // dates != ids
  ids.erase(ids.begin(), ids.begin()+2);
   std::string date;
  while(std::getline(std::cin, lineBuffer)){
    series.push_back(container_type());
    std::istringstream tmp(lineBuffer);
    date.erase();
    for(int i(0); i!=nbDateFields; ++i){
      std::string tmpDate;
      tmp>>tmpDate;
      date.append(tmpDate);
      date.append("\t");
    }
    dates.push_back(date);
    date.erase();
    for(int i(0); i!=nbDateFields; ++i){
      std::string tmpDate;
      tmp>>tmpDate;
      date.append(tmpDate);
      date.append("\t");
    }

    // I do not use the end date, it must be the same as next start date
    std::istream_iterator<na_reader<data_type> >b(tmp),e ;
    std::copy(b,e,std::back_inserter(series.back()));
    /*
    std::cerr<<"\nseries values:";
    std::copy(series.back().begin(), series.back().end(), std::ostream_iterator<data_type>(std::cerr,"\t"));
    */
  }
  dates.push_back(date);// get the last date
  // handle empty last lines
  while(series.back().size()==0) { series.resize(series.size()-1);}
  int seqSize(series.back().size());  
  std::cerr<<"reading "<<series.back().size()<<" series of size"<<series.size()<<"\n";
  typedef std::vector<c_iter_type> sit_type;
  typedef sit_type::iterator set_it;
  typedef boost::tuple<set_it, set_it> sit_it_pair;


  sit_type seriesIts;
  for( series_set_type::const_iterator it(series.begin()); it != series.end(); ++it)
    { seriesIts.push_back(it->begin()); }

  if((nbSegsMin>=1) && (nbSegsMax <= series.size()+1)){
    std::string outputFormat("Dummy");
    std::cout<<"date_debut\tdate_fin\ttaille";
    for(strings_container::const_iterator it(ids.begin());it!=ids.end();++it){
      std::cout<<'\t'<<*it<<"_valeur\t"<<*it<<"_cout";
    }
    std::cout<<std::endl;
	  sit_it_pair serieBegEnd(seriesIts.begin(), seriesIts.end());
	  boost::tuple<int, int>nbSegsMinMax(nbSegsMin,nbSegsMax);
	  boost::tuple<int, int>nbProtosMinMax(nbProtosMin,nbProtosMax);
	  if(model=="Linear0"){
	    if(algo == "TopDown"){
	      do_it<linear0_prototype_level_sequence_model, binary_splits>(model+algo, serieBegEnd, seqSize, nbSegsMinMax, nbProtosMinMax, outputFormat, dates);
	    }else if(algo=="Hybrid"){
	      do_it<linear0_prototype_level_sequence_model, hybrid_splits>(model+algo, serieBegEnd, seqSize,nbSegsMinMax, nbProtosMinMax, outputFormat,dates);
	    }else if(algo=="Optimal"){
	      do_it<linear0_prototype_level_sequence_model, optimal_splits>(model+algo, serieBegEnd, seqSize, nbSegsMinMax, nbProtosMinMax, outputFormat,dates);	    
	    }else{std::cerr<<"Algo "<<algo<<" not implemented (yet?)\n"; badArgs=true;}
	  }else if(model=="Linear1"){
	    if(algo == "TopDown"){
	      do_it<linear0_prototype_level_sequence_model, binary_splits>(model+algo, serieBegEnd, seqSize, nbSegsMinMax, nbProtosMinMax, outputFormat,dates);
	    }else if(algo=="Hybrid"){
	      do_it<linear0_prototype_level_sequence_model, hybrid_splits>(model+algo, serieBegEnd, seqSize, nbSegsMinMax, nbProtosMinMax, outputFormat,dates);
	    }else if(algo=="Optimal"){
	      do_it<linear0_prototype_level_sequence_model, optimal_splits>(model+algo, serieBegEnd, seqSize, nbSegsMinMax, nbProtosMinMax, outputFormat,dates);
	      
	    }else{std::cerr<<"Algo "<<algo<<" not implemented (yet?)\n"; badArgs=true; }
	  }else {std::cerr<<"Model "<<model<<"not implemented (yet?)\n"; badArgs=true;}
  }else{std::cerr<<"invalid nbseg interval, must be between 1 and "<<series.size()<<'\n';}
    }catch(boost::bad_lexical_cast &)
      {badArgs=true; }
  }
  if(badArgs){
    std::cerr<<"Usage: "<<argv[0]<<" 0|1 for center 0|1 for reduce Linear0|Linear1 TopDown|Hybrid|Optimal nProtosMin nProtosMax NbSegsMin nbSegsMax Segments|Points|None\n"
	     <<"Where 1 < nbSegsMin <= nbSegsMax < (nb points in input stream) is the number of segments\n"
	     <<"Where 1 < nbProtosMin <= nbProtosMax < nbSegsMin is the number of prototypes\n"
	     <<" segmentations are written to output stream in the following format:\n"
	     <<"Name\\ttime elapsed (s)\\tnb of segments\\tSSE\\tresulting segmentation\n"
	     <<"Where resulting segmentation format depends on the last argument:\n"
	     <<" Segments -> (((model parameters)cost)(startSegment endSegment))...\n"
	     <<" Points -> reconstructed time-series\n"
	     <<" None -> no further output (for benchmarking)\n";
  }
  return 0;

}
    
