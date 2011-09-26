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
#include "feats/segmentations/utility.hxx"
#include "feats/segmentations/prototypes_models.hxx"
#include "feats/segmentations/prototypes_segmentations.hxx"

#include "utility.hxx"

// g++ -std=c++0x segmentations_proto_tests.cxx -o segmentations_proto_tests -I../include -lstdc++  -I/home/bernard/Code/repositories/boost-trunk -Wall -Wno-sign-compare

using namespace feats::segmentations;




/***
    main segmentation tester functor. creates the objets and call the functors.
***/

template<template<typename, typename,typename>class ProtoModelT
	 , template<typename, typename>class Seg
	 , typename In>
void tester(const std::string& name, boost::tuple<In, In> dataBegEnd , const boost::tuple<int,int>& nbSegsMinMax, const boost::tuple<int,int>& nbProtosMinMax, const std::string& format){
  typedef serie_for_model<ProtoModelT<int*,boost::tuple<int,int>,int>, In> serie_converter_type;// MUST find a less ugly way to adapt to model template than using a dummy: traits ?
  serie_converter_type conv;
  typedef prototypes_segmentations<ProtoModelT, typename serie_converter_type::iterator, Seg,Seg> protos_segs_type;
  typedef typename protos_segs_type::prototypes_segmentation prototypes_segmentation ;
  some_output<prototypes_segmentation> o(format);
  protos_segs_type    protoSegs(conv(dataBegEnd));
    boost::posix_time::ptime startChrono(boost::posix_time::microsec_clock::local_time());
    for( int nbSegs(1); nbSegs!= boost::get<1>(nbSegsMinMax); ++nbSegs,protoSegs.get_free_segmentation().inc_segments()){
      if(nbSegs>=boost::get<0>(nbSegsMinMax)){
	for( int nbProtos(1); nbProtos != boost::get<1>(nbProtosMinMax); ++nbProtos,protoSegs.get_parameter_quantizer().inc_prototypes()){
	  
	if(nbProtos>=boost::get<0>(nbProtosMinMax)){
	  // output twice, calling optimize betwen the two
	  for(int nbIt(0); nbIt!=1; ++nbIt){// !!! !=2 for Optim
	    if(nbIt==1){ protoSegs.optimize(); std::cout<<"Optimized";}
	    boost::posix_time::ptime stopChrono(boost::posix_time::microsec_clock::local_time());
	    std::cout<<name<<'\t'<<boost::posix_time::to_simple_string(stopChrono-startChrono)<<'\t'<<nbProtos<<'\t';
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
  bool badArgs(argc!=8);

  if(!badArgs){
    try{
      std::string const model(boost::lexical_cast<std::string>(argv[1]));
      std::string const algo(boost::lexical_cast<std::string>(argv[2]));
      int nbProtosMin(boost::lexical_cast<int>(argv[3]));
      int nbProtosMax(boost::lexical_cast<int>(argv[4]));
      int nbSegsMin(boost::lexical_cast<int>(argv[5]));
      int nbSegsMax(boost::lexical_cast<int>(argv[6]));
      std::string const outputFormat(boost::lexical_cast<std::string>(argv[7]));
      if((outputFormat=="Segments")||(outputFormat=="Points")||(outputFormat=="None")){
	     
	typedef std::vector<data_type> container_type;
	typedef container_type::const_iterator c_iter_type;

	std::istream_iterator<data_type> b(std::cin), e;
	container_type serie(b,e);
	std::cerr<<"reading "<<serie.size()<<" elements\n";
	if((nbSegsMin>=1) && (nbSegsMax <= serie.size()+1)){
	  boost::tuple<c_iter_type, c_iter_type> serieBegEnd(serie.begin(), serie.end());
	  boost::tuple<int, int>nbSegsMinMax(nbSegsMin,nbSegsMax);
	  boost::tuple<int, int>nbProtosMinMax(nbProtosMin,nbProtosMax);
	  if(model=="Linear0"){
	    if(algo == "TopDown"){
	      tester<linear0_prototype_level_model, binary_splits>(model+algo, serieBegEnd, nbSegsMinMax, nbProtosMinMax, outputFormat);
	    }else if(algo=="Optimal"){
	      tester<linear0_prototype_level_model, optimal_splits>(model+algo, serieBegEnd, nbSegsMinMax, nbProtosMinMax, outputFormat);
	    }else{std::cerr<<"Algo "<<algo<<" not implemented (yet?)\n"; badArgs=true;}
	  }else if(model=="Linear1"){
	    if(algo == "TopDown"){
	      tester<linear1_prototype_slope_model, binary_splits>(model+algo, serieBegEnd, nbSegsMinMax, nbProtosMinMax, outputFormat);
	    }else if(algo=="Optimal"){
	      tester<linear1_prototype_slope_model, optimal_splits>(model+algo, serieBegEnd, nbSegsMinMax, nbProtosMinMax, outputFormat);
	    }else{std::cerr<<"Algo "<<algo<<" not implemented (yet?)\n"; badArgs=true; }
	  }else {std::cerr<<"Model "<<model<<"not implemented (yet?)\n"; badArgs=true;}
	}else{std::cerr<<"invalid nbseg interval, must be between 1 and "<<serie.size()+1<<'\n';}
      }else{std::cerr<<" unknown Format\n";}
    }catch(boost::bad_lexical_cast &)
      {badArgs=true; }
  }
  if(badArgs){
    std::cerr<<"Usage: "<<argv[0]<<" Linear0|Linear1 TopDown|Optimal nProtosMin nProtosMax NbSegsMin nbSegsMax Segments|Points|None\n"
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
    
