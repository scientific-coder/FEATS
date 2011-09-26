#include <iterator>
#include <iostream>
#include <vector>
#include <functional>

#include <numeric>
#include <string>

#define BOOST_NO_EXCEPTIONS

#include <boost/utility.hpp>

#include <boost/tuple/tuple_io.hpp>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>


#include "feats/segmentations/models.hxx"
#include "feats/segmentations/segment.hxx"
#include "feats/segmentations/binary_splits.hxx"
#include "feats/segmentations/merges.hxx"
#include "feats/segmentations/optimal_splits.hxx"
#include "feats/segmentations/utility.hxx"
#include "feats/segmentations/prototypes_models.hxx"
#include "feats/segmentations/prototypes_segmentations.hxx"

#include "utility.hxx"

#define OTL_ORA10G
#define OTL_STREAM_READ_ITERATOR_ON

#include "otlv4.h"


#define WIN32COMMON


#ifndef OCI_ORACLE
# include <oci.h>
#endif





// g++ caddy.cxx -fPIC -I../include -lstdc++ -Wall -Wno-sign-compare
// cl caddy.cxx  -IC:\Program Files\Microsoft Visual Studio .NET
2003\Vc7\PlatformSDK\Include -I..\include -I..\..\BOOST /LD -D_WINDLL /link 
// cl caddy.cxx  -I..\..\Sys\Include -I..\include -I
C:\Applications\Oracle10g\xdk\include -I C:\Applications\Oracle10g\OCI\include -I
D:..\..\BOOST /LD -D_WINDLL /link
/DEFAULTLIB:C:\Applications\Oracle10g\OCI\lib\MSVC/oci.lib

using namespace feats::segmentations;


/***
    main segmentation tester functor. creates the objets and call the functors.
***/

template<template<typename, typename,typename>class ProtoModelT
         , template<typename, typename>class Seg
         , typename In>
static void do_proc(const std::string& name, boost::tuple<In, In> dataBegEnd , const int nbSegs, const int nbProtos, otl_connatc& otc){
           typedef serie_for_model<ProtoModelT<int*,boost::tuple<int,int>,int>, In>
serie_converter_type;// MUST find a less ugly way to adapt to model template than using a dummy: traits ?
           serie_converter_type conv;
           typedef prototypes_segmentations<ProtoModelT, typename
serie_converter_type::iterator, Seg,Seg> protos_segs_type;
           typedef typename protos_segs_type::prototypes_segmentation
prototypes_segmentation ;
           //  some_output<prototypes_segmentation> o(format);
           protos_segs_type    protoSegs(conv(dataBegEnd));
           for( int nbSegs(1); nbSegs!= nbSegs+1; ++nbSegs)
	     {protoSegs.free_segmentation().inc_segments();}
	   for( int nbProtos(1); nbProtos !=  nbProtos+1; ++nbProtos)
	     { protoSegs.parameter_quantizer().inc_prototypes();}
  
                     // CADDY: remplacer par l'insertion des donn-Aées dans la base avec OCI-b

                     //    o(protoSegs.get_segmentation(), std::cout)<<std::endl;
            



                   }
                 }
               }
             }
           }
         }




extern "C" 
#ifdef WIN32COMMON
__declspec(dllexport)
#endif





int test(OCIExtProcContext* contextPtr){
/* recuperer les infos de connection avec OCI */
  OCIEnv     *envhp;                           /* For OCI Environment Handle */
  OCISvcCtx  *svchp;                               /* For OCI Service Handle */
  OCIError   *errhp;                                /* For OCI Error Handle  */
        
int err = OCIExtProcGetEnv(contextPtr,&envhp,&svchp,&errhp);
                        /*connection OTL */
                                                         otl_connect connectOTL;
                                                         connectOTL.rlogon(envhp, svchp); 
                                                        /*        charger les donnees par OTL */

                          otl_stream i(50, // buffer size
              "select valeur from symbole_valeur ",
                          connectOTL // connect object
             ); 
   // create select stream
 
 int tmpId;
 otl_stream_read_iterator<otl_stream,otl_exception> rs;

 rs.attach(i); // attach the iterator "rs" to the stream "i".

 std::vector<int> values;

 while(rs.next_row()){ // while not end-of-data
  rs.get(1,tmpId); values.push_back(tmpId);
 }

 rs.detach(); // detach the itertor from the stream

        return std::accumulate(values.begin(), values.end(),0);

}

extern "C" 
        int proc_sbsr(OCIExtProcContext* contextPtr
                                                        , char* idCapteur, int idCapteurLength
                                                        , char * dateDebut, int dateDebutLength
                                                        , char* dateFin, int dateFinLength
                                                        , char* representation, int representationLength // des donnees en entree
                                                        , int nbSymbols, int nbSegs
                                                        , char* argAlgo, int algoLength
                                                        , char* argModel, int modelLength){
  typedef float data_type; // ou double ?
  std::string const model(argModel);
  std::string const algo(argAlgo);
             
  typedef std::vector<data_type> container_type;
  typedef container_type::const_iterator c_iter_type;

  container_type serie;

  if((nbSegsMin>=1) && (nbSegsMax <= serie.size()+1)){
    
    
    
    int err = OCIExtProcGetEnv(contextPtr,&envhp,&svchp,&errhp);
    /*connection OTL */
    otl_connect connectOTL;
    connectOTL.rlogon(envhp, svchp); 
    /*        charger les donnees par OTL */
    otl_stream i(50, // buffer size
		 "select valeur from symbole_valeur ",
		 connectOTL // connect object
		 ); 
    // create select stream
    
    int tmpId;
    otl_stream_read_iterator<otl_stream,otl_exception> rs;
    
    rs.attach(i); // attach the iterator "rs" to the stream "i".
    
    while(rs.next_row()){ // while not end-of-data
      rs.get(1,tmpId); serie.push_back(tmpId);
    }
    
    rs.detach(); // detach the itertor from the stream
    
    std::string reprName (model+":"+Algo+":"+nbSymbols+":"+nbSegs); // il faudra recuperer les autres info caracteristiques de la serie traitee
    if((nbSegs>=1) && (nbSegs <= serie.size()+1)){
      boost::tuple<c_iter_type, c_iter_type> serieBegEnd(serie.begin(), serie.end());
      
      if(model=="Linear0"){
      if(algo == "TopDown"){
        do_sbsr<linear0_prototype_level_model, binary_splits>(reprName, serieBegEnd,nbSegs, nbSymbols, connectOTL);
      }else if(algo=="Optimal"){
        tester<linear0_prototype_level_model, optimal_splits>(reprName, serieBegEnd,nbSegs, nbSymbols, connectOTL);
      }else{std::cerr<<"Algo "<<algo<<" not implemented (yet?)\n";}
    }else if(model=="Linear1"){
      if(algo == "TopDown"){
        tester<linear1_prototype_slope_model, binary_splits>(reprName, serieBegEnd,nbSegs, nbSymbols, connectOTL);
      }else if(algo=="Optimal"){
        tester<linear1_prototype_slope_model, optimal_splits>(reprName, serieBegEnd,nbSegs, nbSymbols, connectOTL);
      }else{std::cerr<<"Algo "<<algo<<" not implemented (yet?)\n"; }
    }else {std::cerr<<"Model "<<model<<"not implemented (yet?)\n";}
                                                         
                                                         
                                                         return 0;
}


extern "C"
int sbsr(char const* argModel, int argModelLength, char*const argAlgo, int
argAlgoLength, int nbProtosMin, int nbProtosMax, int nbSegsMin, int nbSegsMax){ //
CADDY: ajouter des arguments pour indiquer sur quelles donnees on doit travailler
  typedef double data_type;

  std::string const model(argModel);
  std::string const algo(argAlgo);
             
  typedef std::vector<data_type> container_type;
  typedef container_type::const_iterator c_iter_type;

  // CADDY: remplacer par la recuperation des donnees avec OCI

  std::istream_iterator<data_type> b(std::cin), e;
  container_type serie(b,e);

  if((nbSegsMin>=1) && (nbSegsMax <= serie.size()+1)){
    boost::tuple<c_iter_type, c_iter_type> serieBegEnd(serie.begin(), serie.end());
    boost::tuple<int, int>nbSegsMinMax(nbSegsMin,nbSegsMax);
    boost::tuple<int, int>nbProtosMinMax(nbProtosMin,nbProtosMax);
    if(model=="Linear0"){
      if(algo == "TopDown"){
        tester<linear0_prototype_level_model, binary_splits>(model+algo, serieBegEnd,
nbSegsMinMax, nbProtosMinMax);
      }else if(algo=="Optimal"){
        tester<linear0_prototype_level_model, optimal_splits>(model+algo, serieBegEnd,
nbSegsMinMax, nbProtosMinMax);
      }else{std::cerr<<"Algo "<<algo<<" not implemented (yet?)\n";}
    }else if(model=="Linear1"){
      if(algo == "TopDown"){
        tester<linear1_prototype_slope_model, binary_splits>(model+algo, serieBegEnd,
nbSegsMinMax, nbProtosMinMax);
      }else if(algo=="Optimal"){
        tester<linear1_prototype_slope_model, optimal_splits>(model+algo, serieBegEnd,
nbSegsMinMax, nbProtosMinMax);
      }else{std::cerr<<"Algo "<<algo<<" not implemented (yet?)\n"; }
    }else {std::cerr<<"Model "<<model<<"not implemented (yet?)\n";}
  }else{std::cerr<<"invalid nbseg interval, must be between 1 and
"<<serie.size()+1<<'\n';}
  return 0;

}
