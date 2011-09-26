#include <boost/lexical_cast.hpp>
#include <iterator>
#include <iostream>
#include <sstream>
#include <vector>

#include<algorithm>
#include "utility.hxx"
#include "feats/utility/na_reader.hxx"

using namespace feats::utility;

//g++ centerReduceDownsampleFormat2SDFS.cxx -std=c++0x -O4 -I../include -I/home/bernard/Code/repositories/boost-trunk  -lstdc++ -o centerReduceDownSampleFormat2SDFS

int main(int argc, char* argv[]){
  typedef double data_type;
  bool badArgs(argc!=8);

  if(!badArgs){
    int headerSize(boost::lexical_cast<int>(argv[1]));
    int const nonNumRowsBegin(boost::lexical_cast<int>(argv[2]));
    int const nonNumRowsEnd(boost::lexical_cast<int>(argv[3]));
    bool const shouldCenter(boost::lexical_cast<bool>(argv[4]));
    bool const shouldReduce(boost::lexical_cast<bool>(argv[5]));
    int downsampleStep(boost::lexical_cast<int>(argv[6]));
    int outputLineSize(boost::lexical_cast<int>(argv[7]));

    typedef std::vector<data_type> container_type;
    typedef std::vector<std::string> strings_container;
    typedef std::vector<strings_container> strings_cont_container;
    typedef container_type::const_iterator c_iter_type;
    typedef container_type::difference_type size_type;
    
    typedef std::vector<container_type> series_set_type;
    
    std::string lineBuffer;
    series_set_type series;
    strings_cont_container beginStrings(nonNumRowsBegin), endStrings(nonNumRowsEnd);

    while(std::getline(std::cin, lineBuffer)){
      // copy header lines unchanged
      if(headerSize){
	--headerSize;
	std::cout<<lineBuffer<<std::endl;
      }else{
	series.push_back(container_type());

	std::istringstream tmp(lineBuffer);
	for(int i(0); i!= nonNumRowsBegin;++i){
	  std::string tmpString;
	  tmp>>tmpString;
	  beginStrings[i].push_back(tmpString);
	}
	for(int i(0); i!= nonNumRowsEnd;++i){
	  std::string tmpString;
	  tmp>>tmpString;
	  endStrings[i].push_back(tmpString);
	}
	/*
	  std::istream_iterator<std::string>b(tmp),e ;
	  std::copy(b,e,std::ostream_iterator<std::string>(std::cerr,"|"));
	*/
	std::istream_iterator<na_reader<data_type> >b(tmp),e ;
	std::copy(b,e,std::back_inserter(series.back()));
	//std::copy(b,e,std::ostream_iterator<data_type>(std::cerr,"|"));
      }
    }
    std::cerr<<"read done\n";
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
    std::cerr<<"transposed\n";
    /*
      downsample begin and end non numeric strings
    */
    if(downsampleStep>1){
      int target(0),sourceBegin(0), sourceEnd(downsampleStep), initSize(beginStrings.front().size());
      //      std::cerr<<"initSize"<<initSize<<'\n';
      while(sourceBegin<initSize){
	//	std::cerr<<"target:"<<target<<" sourceBegin:"<<sourceBegin<<" sourceEnd:"<<sourceEnd<<std::endl;
	for(int i(0); i!= nonNumRowsBegin; ++i){
	beginStrings[i][target]=beginStrings[i][sourceBegin];
	}
	for(int i(0); i!= nonNumRowsEnd; ++i){
	endStrings[i][target]=endStrings[i][sourceEnd-1];
	}
	++target;
	sourceBegin=sourceEnd;
	sourceEnd+=downsampleStep;
	if(sourceEnd>initSize) sourceEnd=initSize;
      }
      for(int i(0); i!= nonNumRowsBegin; ++i){
	beginStrings[i].resize(target);
      }
      for(int i(0); i!= nonNumRowsEnd; ++i){
	endStrings[i].resize(target);
      }
    }
    std::cerr<<"downsampled\n";


    for(int i(0); i!= series.size(); ++i){
      container_type& s(series[i]);
      if(downsampleStep>1){
	s.resize(std::distance(s.begin(),downsample(downsampleStep,s.begin(),s.end(),s.begin())));
      }
      if(shouldCenter){
	center(s.begin(), s.end(), s.begin());
      }
      if(shouldReduce){
	reduce(s.begin(), s.end(), s.begin());
      }
    }

    // output

    for(size_type i(0); i!=series.size(); ++i){
      for(size_type j(0); j!=series.front().size(); ++j){
	std::cout<<series[i][j];
	//	std::cerr<<"j:"<<j<<' ';
	if(((j+1)%outputLineSize)==0) {std::cout<<'\n';}
	else std::cout<<std::cout<<'\t';
      }
    }
  }
  return 0;
}
