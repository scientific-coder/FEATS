#include <iterator>
#include <iostream>
#include<vector>
#include <functional>

#include <numeric>
#include <string>

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

#include "utility.hxx"

// g++ -std=c++0x segmentations_tests.cxx -o segmentations_tests -I../include -I. -lstdc++ 

using namespace feats::segmentations;

/***
    main segmentation functor, dispathing according to algorithm
    default implementation do topdown segmentation increasing the number of segments
    specialized for merges to decrease the number of segments
 ***/
template<typename Segmentation, typename OutputOp=some_output<Segmentation> > struct do_segmentation{
  typedef typename Segmentation::segment_type segment_type;
  void operator()(const std::string& name, Segmentation& seg, int nbMin, int nbMax, OutputOp op){
    boost::posix_time::ptime startChrono(boost::posix_time::microsec_clock::local_time());
    // looping over an index and not nb_segs() because overflowing the nb of segments is undefined (would happen if nbMax ==nbPoins+1)
    for( int nbSegs(seg.nb_segs()); nbSegs<nbMax;++nbSegs ){
      if(nbSegs >= nbMin){
	boost::posix_time::ptime stopChrono(boost::posix_time::microsec_clock::local_time());
	std::cout<<name<<'\t'<<boost::posix_time::to_simple_string(stopChrono-startChrono)<<'\t';
	op(seg, std::cout)<<std::endl;
	startChrono+=(boost::posix_time::microsec_clock::local_time()-stopChrono); //substract time spend in op
      }if(nbSegs < (nbMax-1))seg.inc_segments(); // can afford to overflow if nbMax=nb of points +1
    }
  }    
};

template<typename Seg, typename Nb, typename OutputOp >struct do_segmentation<merges<Seg, Nb> ,OutputOp> {
  void operator()(const std::string& name, merges<Seg,Nb>& seg, int nbMin, int nbMax, OutputOp op){
    boost::posix_time::ptime startChrono(boost::posix_time::microsec_clock::local_time());
    while(seg.nb_segs()>nbMax){seg.dec_segments(); }
    while(seg.nb_segs()>nbMin){
      seg.dec_segments();
      boost::posix_time::ptime stopChrono(boost::posix_time::microsec_clock::local_time());
      std::cout<<name<<'\t'<<boost::posix_time::to_simple_string(stopChrono-startChrono)<<'\t';
      op(seg, std::cout)<<std::endl;
      startChrono+=(boost::posix_time::microsec_clock::local_time()-stopChrono); //substract time spend in op
    }
  }    
};

/***
main segmentation tester functor. creates the objets and call the functors.
 ***/
template<typename Model, template<typename, typename>class Seg, typename In>
void tester(const std::string& name, boost::tuple<In, In> dataBegEnd, int nbMin, int nbMax, const std::string& format){
  // must dispatch three things: 
  // iterator must be changed for linear1_model -> serie_for_model template functor
  // initial trivial segmentation(1 or N segments) ->init_segments template functor
  // and order of segmetations for bottom_up_segmentation -> do_segmentation functor
  typedef Model model_type;
  typedef serie_for_model<model_type, In> serie_converter_type;
  typedef segment<model_type, typename serie_converter_type::iterator> segment_type;
  typedef Seg<segment_type, int> Segmentation;
  serie_converter_type serieConverter;
  init_segments<Segmentation> initSegConverter;
  Segmentation seg(initSegConverter(serieConverter(dataBegEnd)));
  typedef some_output<Segmentation> output_type;
  output_type ouput(format);
  do_segmentation<Segmentation> doIt;

  doIt(name, seg, nbMin, nbMax, ouput);
}


/***
    Small program for testing segmentations and models. Shows how static polymorphism allows to combine segment models and segmentation algorithms.
 ***/

int main(int argc, char* argv[]){
  typedef double data_type;
  bool badArgs(argc!=6);

  if(!badArgs){
    try{
      std::string const model(boost::lexical_cast<std::string>(argv[1]));
      std::string const algo(boost::lexical_cast<std::string>(argv[2]));
      int nbSegsMin(boost::lexical_cast<int>(argv[3]));
      int nbSegsMax(boost::lexical_cast<int>(argv[4]));
      std::string const outputFormat(boost::lexical_cast<std::string>(argv[5]));
      if((outputFormat=="Segments")||(outputFormat=="Points")||(outputFormat=="None")){
	     
	typedef std::vector<data_type> container_type;
	typedef container_type::const_iterator c_iter_type;

	std::istream_iterator<data_type> b(std::cin), e;
	container_type serie(b,e);
	std::cerr<<"reading "<<serie.size()<<" elements\n";
	if((nbSegsMin>=1) && (nbSegsMax <= serie.size()+1)){
	  boost::tuple<c_iter_type, c_iter_type> serieBegEnd(serie.begin(), serie.end());
	  if(model=="Linear0"){
	    typedef linear0_model<data_type> model_type;
	    if(algo == "TopDown"){
	      tester<model_type, binary_splits>(model+algo, serieBegEnd, nbSegsMin, nbSegsMax, outputFormat);
	    }else if(algo=="BottomUp"){
	      tester<model_type, merges>(model+algo, serieBegEnd, nbSegsMin, nbSegsMax, outputFormat);
	    }else if(algo=="Optimal"){
	      tester<model_type, optimal_splits>(model+algo, serieBegEnd, nbSegsMin, nbSegsMax, outputFormat);
	    }else{std::cerr<<"Algo "<<algo<<" not implemented (yet?)\n"; badArgs=true;}
	  }else if(model=="Linear1"){
	    typedef linear1_model<boost::tuple<data_type, data_type> > model_type;
	    if(algo == "TopDown"){
	      tester<model_type, binary_splits>(model+algo, serieBegEnd, nbSegsMin, nbSegsMax, outputFormat);
	    }else if(algo=="BottomUp"){
	      tester<model_type, merges>(model+algo, serieBegEnd, nbSegsMin, nbSegsMax, outputFormat);
	    }else if(algo=="Optimal"){
	      tester<model_type, optimal_splits>(model+algo, serieBegEnd, nbSegsMin, nbSegsMax, outputFormat);
	    }else{std::cerr<<"Algo "<<algo<<" not implemented (yet?)\n"; badArgs=true; }
	  }else {std::cerr<<"Model "<<model<<"not implemented (yet?)\n"; badArgs=true;}
	}else{std::cerr<<"invalid nbseg interval, must be between 1 and "<<serie.size()+1<<'\n';}
      }else{std::cerr<<" unknown Format\n";}
    }catch(boost::bad_lexical_cast &)
      {badArgs=true; }
  }
  if(badArgs){
    std::cerr<<"Usage: "<<argv[0]<<" Linear0|Linear1 TopDown|BottomUp|Optimal nMin nMax Segments|Points|None\n"
	     <<"Where 1 < nMin <= nMax < (nb points in input stream) is the number of segments\n"
	     <<" segmentations are written to output stream in the following format:\n"
	     <<"Name\\ttime elapsed (s)\\tnb of segments\\tSSE\\tresulting segmentation\n"
	     <<"Where resulting segmentation format depends on the last argument:\n"
	     <<" Segments -> (((model parameters)cost)(startSegment endSegment))...\n"
	     <<" Points -> reconstructed time-series\n"
	     <<" None -> no further output (for benchmarking)\n";
  }
  return 0;

}
    
