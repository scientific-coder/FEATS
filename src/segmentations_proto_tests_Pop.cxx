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

// g++ -std=c++0x segmentations_proto_tests_Pop.cxx -o segmentations_proto_tests_Pop -I../include -I/home/bernard/Code/repositories/boost-trunk -lstdc++ -Wall -Wno-sign-compare



using namespace feats::segmentations;
using namespace feats::utility;



/***
    main segmentation tester functor. creates the objets and call the functors.
***/

template<template<typename, typename,typename>class ProtoModelT
	 , template<typename, typename>class Seg
	 , typename In>
void do_it(const std::string& name, boost::tuple<In, In> dataBegEnd , int seqSize, const boost::tuple<int,int>& nbSegsMinMax, const boost::tuple<int,int>& nbProtosMinMax, const std::string& format, int segsMod, int protosMod){


  typedef prototypes_sequence_segmentations<ProtoModelT, In, Seg, binary_splits> protos_segs_type;
  typedef typename protos_segs_type::prototypes_segmentation prototypes_segmentation ;
  protos_segs_type    protoSegs(dataBegEnd, seqSize);
  some_sequences_output<prototypes_segmentation> o(format,seqSize);


    boost::posix_time::ptime startChrono(boost::posix_time::microsec_clock::local_time());
    for( int nbSegs(1); nbSegs!= boost::get<1>(nbSegsMinMax); ++nbSegs,protoSegs.free_segmentation().inc_segments()){
      if(nbSegs>=boost::get<0>(nbSegsMinMax) && ((nbSegs%segsMod)==0)){
	for( int nbProtos(1); nbProtos != boost::get<1>(nbProtosMinMax); ++nbProtos,protoSegs.parameter_quantizer().inc_prototypes()){
	  
	  if(nbProtos>=boost::get<0>(nbProtosMinMax)&& ((nbProtos %protosMod)==0)){
	  // output twice, calling optimize betwen the two
	  for(int nbIt(0); nbIt!=2; ++nbIt){ //2 for Optim
	    if(nbIt==1){ protoSegs.optimize(); std::cout<<"Optimized";}
	    boost::posix_time::ptime stopChrono(boost::posix_time::microsec_clock::local_time());
	    std::cout<<name<<'\t'<<boost::posix_time::to_simple_string(stopChrono-startChrono)<<'\t'<<nbProtos<<'\t'<<nbSegs<<'\t';
	    std::cout<<protoSegs.get_segmentation().cost()<<'\n';
	    o(protoSegs.get_segmentation(), std::cout)<<std::endl;
	    startChrono+=(boost::posix_time::microsec_clock::local_time()-stopChrono); //substract time spend in op

	  }
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
  bool badArgs(argc!=11);

  if(!badArgs){
    try{
      std::string const model(boost::lexical_cast<std::string>(argv[1]));
      std::string const algo(boost::lexical_cast<std::string>(argv[2]));
      int nbProtosMin(boost::lexical_cast<int>(argv[3]));
      int nbProtosMax(boost::lexical_cast<int>(argv[4]));
      int protosMod(boost::lexical_cast<int>(argv[5]));
      int nbSegsMin(boost::lexical_cast<int>(argv[6]));
      int nbSegsMax(boost::lexical_cast<int>(argv[7]));
      int segsMod(boost::lexical_cast<int>(argv[8]));


      std::string const outputFormat(boost::lexical_cast<std::string>(argv[9]));
      bool const mustTranspose(boost::lexical_cast<bool>(argv[10]));
      if((outputFormat=="Segments")||(outputFormat=="Points")||(outputFormat=="None")){
	     
	typedef std::vector<data_type> container_type;
	typedef container_type::const_iterator c_iter_type;
	typedef container_type::difference_type size_type;

  typedef std::vector<container_type> series_set_type;
  
  std::string lineBuffer;
  series_set_type series;

  while(std::getline(std::cin, lineBuffer)){
    series.push_back(container_type());
    std::istringstream tmp(lineBuffer);
    std::istream_iterator<na_reader<data_type> >b(tmp),e ;
    std::copy(b,e,std::back_inserter(series.back()));
  }


  typedef std::vector<c_iter_type> sit_type;
  typedef sit_type::iterator set_it;
  typedef boost::tuple<set_it, set_it> sit_it_pair;

  /* must transpose the matrix for efficiency reasons*/
  if(mustTranspose){
    series_set_type transposedSeries(series.front().size(),container_type(series.size()));
    for(size_type i(0); i!=series.size(); ++i){
      for(size_type j(0); j!=series.front().size(); ++j){
	transposedSeries[j][i]=series[i][j];
      }
    }
    series.swap(transposedSeries);
  }
  int seqSize(series.front().size());  
  std::cout<<"reading "<<series.front().size()<<" series of size"<<series.size()<<"\n";
  sit_type seriesIts;
  for( series_set_type::const_iterator it(series.begin()); it != series.end(); ++it)
    { seriesIts.push_back(it->begin()); }

  if((nbSegsMin>=1) && (nbSegsMax <= series.size()+1)){
	  sit_it_pair serieBegEnd(seriesIts.begin(), seriesIts.end());
	  boost::tuple<int, int>nbSegsMinMax(nbSegsMin,nbSegsMax);
	  boost::tuple<int, int>nbProtosMinMax(nbProtosMin,nbProtosMax);
	  if(model=="Linear0"){
	    if(algo == "TopDown"){
	      do_it<linear0_prototype_level_sequence_model, binary_splits>(model+algo, serieBegEnd, seqSize, nbSegsMinMax, nbProtosMinMax, outputFormat,segsMod, protosMod);
	    }else if(algo=="Hybrid"){
	      do_it<linear0_prototype_level_sequence_model, hybrid_splits>(model+algo, serieBegEnd, seqSize,nbSegsMinMax, nbProtosMinMax, outputFormat,segsMod, protosMod);
	    }else if(algo=="Optimal"){
	      do_it<linear0_prototype_level_sequence_model, optimal_splits>(model+algo, serieBegEnd, seqSize, nbSegsMinMax, nbProtosMinMax, outputFormat,segsMod, protosMod);	    
	    }else{std::cerr<<"Algo "<<algo<<" not implemented (yet?)\n"; badArgs=true;}
	  }else if(model=="Linear1"){
	    if(algo == "TopDown"){
	      do_it<linear0_prototype_level_sequence_model, binary_splits>(model+algo, serieBegEnd, seqSize, nbSegsMinMax, nbProtosMinMax, outputFormat,segsMod, protosMod);
	    }else if(algo=="Hybrid"){
	      do_it<linear0_prototype_level_sequence_model, hybrid_splits>(model+algo, serieBegEnd, seqSize, nbSegsMinMax, nbProtosMinMax, outputFormat,segsMod, protosMod);
	    }else if(algo=="Optimal"){
	      do_it<linear0_prototype_level_sequence_model, optimal_splits>(model+algo, serieBegEnd, seqSize, nbSegsMinMax, nbProtosMinMax, outputFormat,segsMod, protosMod);
	      
	    }else{std::cerr<<"Algo "<<algo<<" not implemented (yet?)\n"; badArgs=true; }
	  }else {std::cerr<<"Model "<<model<<"not implemented (yet?)\n"; badArgs=true;}
  }else{std::cerr<<"invalid nbseg interval, must be between 1 and "<<series.size()<<'\n';}
      }else{std::cerr<<" unknown Format\n";}
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
    
