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

// g++ -std=c++0x test_clust_adapt_segPKDD2.cxx -I../include/ -o  test_clust_adapt_segPKDD2 -lstdc++ -O4


// objectifs de l'evaluation:
// comparer les resultats d'une classif et seg croisee adaptative avec ceux d'une seg et d'une classif
// etudier les deux cas, seg globale et seg avec equirepartition des segments
// dans chaque cas sse , complexite en nb clusters et nb total de segments
// (a faire sur benchmark, synth et jeu de donnees rellees de l'INRETS
//
// Pour les differents algo a comparer, il suffit de faire differents segmentations_type dans clusters_of_adaptove_segmentations
// 1) segmentations_for_total_number : mon algo qui repartit optimalement
// 2) equal_nb_of_segments : algo qui repartit "equitablement"
// 3) global_segmentations : algo qui effectue une segmentation globale avec nbTot/nbClusts segments
//jdavidou@wanadoo.fr
//Jany Davidou
// 13 rue des sables - 54425 Pulnoy

//06 61 84 34 70

//for i in .0125 .025 .0375 .05 .0625 .075 .0875 .0.1 .125; do for j in .05 .1 .125 .15 .175 .2 .225 .25; do echo $i $j |tee args.txt; R BATCH  --slave --no-restore --no-restore --no-save <synthData.R  | sed -e '1,2d' |head -n 1000 | tee `echo dat-$i $j .txt|tr -d ' ' `|./test_clust_adapt_segPattern 4 5 60 61 10 1 2>`echo errD-60-$i $j .txt|tr -d ' ' ` |tee `echo resD-60-$i $j .txt|tr -d ' ' `|fgrep miss ; done; done

//for i in  .05 .0625 .075 .0875 .1 .125 .15 .175 .2 .225 .25 ; do for j in .05 .0625 .075 .0875 .1 .125 .15 .175 .2 .225 .25; do echo$i$j; bzcat `echo dat-$i $j .txt.bz2|tr -d ' ' `|./test_clust_adapt_segPattern 4 5 30 31 500 1 2>`echo errC-30-500-$i $j .txt|tr -d ' ' `|tee `echo resC-30-500-$i $j .txt|tr -d ' ' `|fgrep miss ; done; done

// for i in res*.*.*.txt; do echo `echo $i |sed -e "s/[^0-9]*\(\.[0-9]*\)/ \1/g"|sed -e "s/.txt//g"` `fgrep miss $i| sed -e "s/missClassification//g"`; done

//for(i in 1:8){tmp<- miss[miss[,1]==miss[,2][i],];plot(tmp[,2],tmp[,3]*100,type="l",ylim=c(0,max(tmp[,3:4])),ylab="erreur de classif",xlab="perturbation sur les valeurs", main=paste("perturbation sur les instants:",tmp[1,1])); lines(tmp[,2],tmp[,4]*100,col="red")}






template<typename ClustType, typename ItPairType, typename SizeType, typename CostType, typename ProtoType>
void eval( ItPairType begEnd, SizeType nbClusts, SizeType nbSegs, SizeType nbRetries
	   , CostType& cost, ProtoType& proto
	   , CostType& costOpt, ProtoType& protoOpt){
  ClustType cl(begEnd, nbClusts, nbSegs, nbRetries, true, true); // changer l'API pour remplacer random par le cas particulier de nbretries à 0
  cl.next_segmentation(); cl.next_clustering();
  {
  CostType const currCost(cl.cost());
  if(currCost<cost){
    cost=currCost;
    for (SizeType p(0); p!= nbClusts; ++p){
      cl.output_prototype(p, proto[p].begin());
    }
  }
  }
  cl.optimize();
  {
    CostType const currCost(cl.cost());
    if(currCost<costOpt){
    costOpt=currCost;
    for (SizeType p(0); p!= nbClusts; ++p){
      cl.output_prototype(p, protoOpt[p].begin());
    }
  }
  }
  return ;
}
template<typename Protos>
void printProtos(const Protos& p){
  for(int k(0); k!=p.size(); ++k)
    {std::copy(p[k].begin(), p[k].end(), std::ostream_iterator<typename Protos::value_type::value_type>(std::cout,"\t"));}
  return;
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
  int const nbPoints(series.front().size()), nbSeries(series.size());

  typedef std::vector<boost::tuple<c_iter_type, c_iter_type> > sit_type;
  typedef sit_type::iterator set_it;
  typedef boost::tuple<set_it, set_it> sit_it_pair;

  sit_type seriesIts;
  for( series_set_type::const_iterator it(series.begin()); it != series.end(); ++it)
    { seriesIts.push_back(boost::make_tuple(it->begin(), it->end())); }

  std::cerr.setf(std::ios::fixed);
  std::cerr<<std::setprecision(4);

  std::cout.setf(std::ios::fixed);  std::cout<<std::setprecision(4);

  sit_it_pair begEnd(boost::make_tuple(seriesIts.begin(), seriesIts.end()));

  series_set_type nanProtos(nbClustersMax-1, container_type(nbPoints,std::numeric_limits<data_type>::quiet_NaN()));
  series_set_type myProtos(nanProtos), myProtosOpt(nanProtos)
    , equalProtos(nanProtos), equalProtosOpt(nanProtos)
    , sameProtos(nanProtos), sameProtosOpt(nanProtos);
  for(int i(nbClustersMin); i!=nbClustersMax; ++i){
    data_type myCost(std::numeric_limits<float>::max());
    data_type myCostOpt(myCost), equalCost(myCost), equalCostOpt(myCost), sameCost(myCost), sameCostOpt(myCost);

    for(int j(std::max(nbClustersMin,nbSegmentsMin)); j!=nbSegmentsMax; ++j,  std::cout<<std::endl){
      for(int r(0); r!= nbRetries; ++r){
	if(j%i==0){
// pour ne refaire que ce qui a ete loupe
	eval<feats::misc::clusters_of_adaptive_segmentations<sit_it_pair,feats::segmentations::segmentations_for_total_number> >(begEnd, i, j, nbRetries, myCost, myProtos, myCostOpt, myProtosOpt);
	  equalCost=equalCostOpt=sameCost=sameCostOpt= std::numeric_limits<data_type>::max();

	  eval<feats::misc::clusters_of_adaptive_segmentations<sit_it_pair,feats::segmentations::equal_nb_of_segments> >(begEnd, i, j, nbRetries, equalCost, equalProtos, equalCostOpt, equalProtosOpt);
	  eval<feats::misc::clusters_of_adaptive_segmentations<sit_it_pair,feats::segmentations::same_episodes> >(begEnd, i, j, nbRetries, sameCost, sameProtos, sameCostOpt, sameProtosOpt);
	}else{
	  // mettre des na quand les segmentatiosn ne sont pas possibles
	  //  ou laisser les valeurs pre
	  equalCost=equalCostOpt=sameCost=sameCostOpt= std::numeric_limits<data_type>::quiet_NaN();
	  equalProtos=equalProtosOpt=sameProtos=sameProtosOpt= nanProtos;
	}
      }
	//mettre aussi l eval de la classification sur les données synth. mais dans un autre prog car le nb de protos n'est pas variable
      std::cout<<i<<'\t'<<j<<'\t'<<myCost<<'\t'<<myCostOpt<<'\t'<<equalCost<<'\t'<<equalCostOpt<<'\t'<<sameCost<<'\t'<<sameCostOpt<<'\t';
      if(printSeries){
	// comment gerer simplement le nb variable de protos ??
	// pour simplifier les traitement,s je mets tout sur une ligne
	// et je mets des protos NA pour completer
	//	std::cout<<"\nmyProtos:";
	printProtos(myProtos);
	// 	std::cout<<"\nmyProtosOpt:";
	printProtos(myProtosOpt); 
	//std::cout<<"\nequalProtos:";
	printProtos(equalProtos);
	//std::cout<<"\nequalProtosOpt:";
	printProtos(equalProtosOpt); 
	//std::cout<<"\nsameProtos:";
	printProtos(sameProtos);
	//std::cout<<"\nsameProtosOpt:";
	printProtos(sameProtosOpt); 
      }
    }
  }
  return 0;
}


      

