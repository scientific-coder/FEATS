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

#include "feats/segmentations/models.hxx"
#include "feats/segmentations/segment.hxx"
#include "feats/segmentations/optimal_splits.hxx"

#include "feats/misc/clusters_of_adaptive_segmentations.hxx"
#include "feats/segmentations/segmentations_for_total_number.hxx"
#include "feats/segmentations/same_episodes.hxx"
#include "feats/segmentations/equal_nb_of_segments.hxx"


#include "feats/utility/na_reader.hxx"

using namespace feats::segmentations;
using namespace feats::utility;

// g++ -std=c++0x test_clust_adapt_segPKDD.cxx -I../include/ -I/home/bernard/Code/repositories/boost-trunk -o  test_clust_adapt_segPKDD -lstdc++ -O4


// objectifs de l'evaluation:
// comparer les resultats d'une classif et seg croisee adaptative avec ceux d'une seg et d'une classif
// etudier les deux cas, seg globale et seg avec equirepartition des segments
// dans chaque cas sse , complexite en nb clusters et nb total de segments
// (a faire sur benchmark, synth et jeu de donnees rellees de l'INRETS
//
// Pour les differents algo a comparer, il suffit de faire differents segmentations_type dans clusters_of_adaptove_segmentations
// 1) segmentations_for_total_number : mon algo qui repartit optimalement mais sans optimisation croisee
// 1.bis) mon algo complet avec optimisation croisee
// 2) equal_nb_of_segments : algo qui repartit "equitablement"
// 3) global_segmentations : algo qui effectue une segmentation globale avec nbTot/nbClusts segments


//for i in .0125 .025 .0375 .05 .0625 .075 .0875 .0.1 .125; do for j in .05 .1 .125 .15 .175 .2 .225 .25; do echo $i $j |tee args.txt; R BATCH  --slave --no-restore --no-restore --no-save <synthData.R  | sed -e '1,2d' |head -n 1000 | tee `echo dat-$i $j .txt|tr -d ' ' `|./test_clust_adapt_segPattern 4 5 60 61 10 1 2>`echo errD-60-$i $j .txt|tr -d ' ' ` |tee `echo resD-60-$i $j .txt|tr -d ' ' `|fgrep miss ; done; done

//for i in  .05 .0625 .075 .0875 .1 .125 .15 .175 .2 .225 .25 ; do for j in .05 .0625 .075 .0875 .1 .125 .15 .175 .2 .225 .25; do echo$i$j; bzcat `echo dat-$i $j .txt.bz2|tr -d ' ' `|./test_clust_adapt_segPattern 4 5 30 31 500 1 2>`echo errC-30-500-$i $j .txt|tr -d ' ' `|tee `echo resC-30-500-$i $j .txt|tr -d ' ' `|fgrep miss ; done; done

// for i in res*.*.*.txt; do echo `echo $i |sed -e "s/[^0-9]*\(\.[0-9]*\)/ \1/g"|sed -e "s/.txt//g"` `fgrep miss $i| sed -e "s/missClassification//g"`; done

//for(i in 1:8){tmp<- miss[miss[,1]==miss[,2][i],];plot(tmp[,2],tmp[,3]*100,type="l",ylim=c(0,max(tmp[,3:4])),ylab="erreur de classif",xlab="perturbation sur les valeurs", main=paste("perturbation sur les instants:",tmp[1,1])); lines(tmp[,2],tmp[,4]*100,col="red")}


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

  typedef feats::misc::clusters_of_adaptive_segmentations<sit_it_pair,feats::segmentations::segmentations_for_total_number> cl_adapt_seg;
  typedef feats::misc::clusters_of_adaptive_segmentations<sit_it_pair,feats::segmentations::equal_nb_of_segments> cl_equal_seg;
  typedef feats::misc::clusters_of_adaptive_segmentations<sit_it_pair,feats::segmentations::same_episodes> cl_same_seg;

  for(int i(nbClustersMin); i!=nbClustersMax; ++i){
    series_set_type bestProtosMy(i,container_type(series.front().size()));
    series_set_type bestProtosMyOpt(bestProtosMy), bestProtosEqual(bestProtosMy), bestProtosEqualOpt(bestProtosMy), bestProtosSame(bestProtosMy), bestProtosSameOpt(bestProtosMy);
    for(int j(std::max(nbClustersMin,nbSegmentsMin)); j!=nbSegmentsMax; ++j){
      float bestKmeans(std::numeric_limits<float>::max()), myBest(std::numeric_limits<float>::max()), sseKmeans, classifKmeans, sseMy, classifMy;
      for(int r(0); r!= nbRetries; ++r,std::cout<<"\nretrying:\n"){
	cl_adapt_seg clustAdapSegs(boost::make_tuple(seriesIts.begin(), seriesIts.end()), i, j, nbRetries, true,true);
	keepIfBetter(clustAdapSegs, best
	cl_equal_seg clustEqualSegs(boost::make_tuple(seriesIts.begin(), seriesIts.end()), i, j, nbRetries, true,true);

	if( j % nbClustersMax ==0){
	  cl_same_seg clustAdapSegs(boost::make_tuple(seriesIts.begin(), seriesIts.end()), i, j, nbRetries, true,true);
	}
	std::cout<<"initialClusters\t";
	clustAdapSegs.output_clusters(std::ostream_iterator<int>(std::cout,","));
	clustAdapSegs.optimize();
	std::cout<<"\nnbIters\t"<<clustAdapSegs.nb_iter()<<"\tcost\t"<<clustAdapSegs.cost();
	if(printSeries){
	  std::cout<<"\nfinalPrototypes\t";
	  for(int k(0); k!=i; ++k, std::cout<<'\n')
	    { clustAdapSegs.output_prototype(k,std::ostream_iterator<data_type>(std::cout,"\t")); }
	}
      }
      // print best costs and prototypes
    }
  }
  return 0;
}


      

