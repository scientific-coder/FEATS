#include <iterator>
#include <iostream>
#include <iomanip>

#include <vector>
#include <functional>
#include <algorithm>// std::next_permutation
#include <cstdlib> // std::abs(int)

#include <numeric>
#include <string>
#include <sstream>
#include <boost/utility.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/progress.hpp>

#include "utility.hxx"

#include "feats/segmentations/segmentations_for_total_number.hxx"
#include "feats/segmentations/models.hxx"
#include "feats/segmentations/segment.hxx"
#include "feats/segmentations/optimal_splits.hxx"

#include "feats/misc/clusters_of_adaptive_segmentations.hxx"
#include "feats/utility/na_reader.hxx"

using namespace feats::segmentations;
using namespace feats::utility;

// g++ -std=c++0x test_clust_adapt_segPattern.cxx -I../include/ -o  test_clust_adapt_segPattern -lstdc++ -O4

//for i in .0125 .025 .0375 .05 .0625 .075 .0875 .0.1 .125; do for j in .05 .1 .125 .15 .175 .2 .225 .25; do echo $i $j |tee args.txt; R BATCH  --slave --no-restore --no-restore --no-save <synthData.R  | sed -e '1,2d' |head -n 1000 | tee `echo dat-$i $j .txt|tr -d ' ' `|./test_clust_adapt_segPattern 4 5 60 61 10 1 2>`echo errD-60-$i $j .txt|tr -d ' ' ` |tee `echo resD-60-$i $j .txt|tr -d ' ' `|fgrep miss ; done; done

//for i in  .05 .0625 .075 .0875 .1 .125 .15 .175 .2 .225 .25 ; do for j in .05 .0625 .075 .0875 .1 .125 .15 .175 .2 .225 .25; do echo$i$j; bzcat `echo dat-$i $j .txt.bz2|tr -d ' ' `|./test_clust_adapt_segPattern 4 5 30 31 500 1 2>`echo errC-30-500-$i $j .txt|tr -d ' ' `|tee `echo resC-30-500-$i $j .txt|tr -d ' ' `|fgrep miss ; done; done

// for i in res*.*.*.txt; do echo `echo $i |sed -e "s/[^0-9]*\(\.[0-9]*\)/ \1/g"|sed -e "s/.txt//g"` `fgrep miss $i| sed -e "s/missClassification//g"`; done

//for(i in 1:8){tmp<- miss[miss[,1]==miss[,2][i],];plot(tmp[,2],tmp[,3]*100,type="l",ylim=c(0,max(tmp[,3:4])),ylab="erreur de classif",xlab="perturbation sur les valeurs", main=paste("perturbation sur les instants:",tmp[1,1])); lines(tmp[,2],tmp[,4]*100,col="red")}

/*

Pb: évaluation suivant le taux de bonne classification, mais pour les kmeans, je ne prends que celui ayant la meilleure SSE :-(
Donc je ne paux pas utiliser bestKmeans...
*/



template<typename ClustType>
float missClassif( ClustType& cl){
  typedef std::vector<int> vect_type;
  vect_type clusts;
  cl.output_clusters(std::back_inserter(clusts));
  std::vector<vect_type> counts(4,vect_type(4,0));
  // count clusters
  for(int i(0),k(0); i!=4; ++i){
    for(int j(0); j!=250; ++j,++k){
      ++counts[i][clusts[k]];
    }
  }
  // testing combinations
  int combi[]={0,1,2,3};
  int currentErr, bestErr(std::numeric_limits<int>::max());
  for(int i(0); i!=24; ++i, std::next_permutation(combi,combi+4)){ // 4!=24
    currentErr=0;
    for(int j(0); j!=4; ++j){
      currentErr+= std::abs(250-counts[j][combi[j]]);
    }
    if(currentErr<bestErr) bestErr=currentErr;
  }
  return bestErr/(2.0*10); // 2: err are counted twice 10: /nbSeries * 100 (%)
}

int main(int argc, char* argv[]){
  int nbClustersMin(boost::lexical_cast<int>(argv[1]));
  int nbClustersMax(boost::lexical_cast<int>(argv[2]));
  int nbSegmentsMin(boost::lexical_cast<int>(argv[3]));
  int nbSegmentsMax(boost::lexical_cast<int>(argv[4]));

  int nbRetries(boost::lexical_cast<int>(argv[5]));

  int printSeries(boost::lexical_cast<int>(argv[6]));

  // type de selection des classifications
  // 0 suivant les SSE 1 suivant les missclassification TODO
  // je pourrai aussi sortir les valeurs moyennes et 
  // rang de la valeur min et nb de cas ou l'un est meilleur que l'autre
  int selectionType(boost::lexical_cast<int>(argv[7]));


  typedef double data_type;
  typedef std::vector<data_type> container_type;
  typedef container_type::const_iterator c_iter_type;
  typedef std::vector<container_type> series_set_type;
  
  std::string lineBuffer;
  series_set_type series;

  while(std::getline(std::cin, lineBuffer)){
    series.push_back(container_type());
    std::istringstream tmp(lineBuffer);
    std::istream_iterator<na_reader<data_type> >b(tmp),e ;
    std::copy(b,e,std::back_inserter(series.back()));
  }
  std::cerr<<"reading "<<series.size()<<" series\n";
  
  typedef std::vector<boost::tuple<c_iter_type, c_iter_type> > sit_type;
  typedef sit_type::iterator set_it;
  typedef boost::tuple<set_it, set_it> sit_it_pair;

  sit_type seriesIts;
  for( series_set_type::const_iterator it(series.begin()); it != series.end(); ++it)
    { seriesIts.push_back(boost::make_tuple(it->begin(), it->end())); }

  std::cerr.setf(std::ios::fixed);
  std::cerr<<std::setprecision(4);

  std::cout.setf(std::ios::fixed);
  std::cout<<std::setprecision(4);
  for(int i(nbClustersMin); i!=nbClustersMax; ++i){
    for(int j(std::max(nbClustersMin,nbSegmentsMin)); j!=nbSegmentsMax; ++j){
      float bestKmeans(std::numeric_limits<float>::max()), myBest(std::numeric_limits<float>::max()), sseKmeans, classifKmeans, sseMy, classifMy;
      for(int r(0); r!= nbRetries; ++r,std::cout<<"\nretrying:\n"){
	typedef feats::misc::clusters_of_adaptive_segmentations<sit_it_pair> cl_adapt_seg;
	cl_adapt_seg clustAdapSegs(boost::make_tuple(seriesIts.begin(), seriesIts.end()), i, j, nbRetries, true,true);
	std::cout<<"initialClusters\t";
	clustAdapSegs.output_clusters(std::ostream_iterator<int>(std::cout,","));
	std::cout<<"\nkmeansClusters\t";
	cl_adapt_seg::kmeans_type km(boost::make_tuple(boost::make_tuple(seriesIts.begin(), seriesIts.end()), i));
	km.compute();
	float tmpMiss(missClassif(km)), tmpCost(km.cost());
	float tmp(selectionType ? tmpMiss : tmpCost);
	if(bestKmeans>tmp){
	  bestKmeans=tmp;
	  sseKmeans=tmpCost;
	  classifKmeans=tmpMiss;
	}
	std::cout<<"\nKmeansmissClassification\t"<<tmpMiss<<"\tKmeansSSE\t"<<tmpCost;
	if(printSeries){
	  std::cout<<"\nKmeansPrototypes\t";
	  for(int k(0); k!=i; ++k, std::cout<<'\n')
	    { km.output_prototype(k,std::ostream_iterator<data_type>(std::cout,",")); }
	}
      
	clustAdapSegs.optimize();
	std::cout<<"\nnbIters\t"<<clustAdapSegs.nb_iter();
	std::cout<<"\nfinalClusters\t";
	
	clustAdapSegs.output_clusters(std::ostream_iterator<int>(std::cout,","));
	tmpMiss=missClassif(clustAdapSegs); tmpCost=clustAdapSegs.cost();
	tmp=(selectionType ? tmpMiss : tmpCost);
	if(myBest>tmp){
	  myBest=tmp;
	  sseMy=tmpCost;
	  classifMy=tmpMiss;
	}
	std::cout<<"\nMymissClassification\t"<<tmpMiss<<"\tmySSE\t"<<tmpCost;
	std::cout<<"\tbestKmeansMissClassification\t"<<classifKmeans<<"\tbestKmeansSSE\t"<<sseKmeans<<"\tbestMyMissClassification\t"<<classifMy<<"\tbestMySSE\t"<<sseMy<<std::endl;
	// early exit:
	if(myBest < std::numeric_limits<float>::min() && bestKmeans < std::numeric_limits<float>::min()){ r=nbRetries-1;}
	if(printSeries){
	  std::cout<<"\nfinalPrototypes\t";
	  for(int k(0); k!=i; ++k, std::cout<<'\n')
	    { clustAdapSegs.output_prototype(k,std::ostream_iterator<data_type>(std::cout,"\t")); }
	}
	
	
	
      }
    }
  }
  return 0;
}


      

