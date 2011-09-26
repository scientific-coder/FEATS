#include <boost/config.hpp>
#include <iostream>
//#include <fstream>
#include <sstream>
#include <string>
#include <utility> // pair
#include <boost/tokenizer.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <map>

using namespace boost;

//utilisation de tuples/pairs des coords pour les noeuds, l'opération < est définie, donc ok pour map

// g++ test_graph.cxx -o test_graph -lstdc++ -O4
// cat ~/Data/Arc_Coordinate.dat |cut -d ' ' -f 3-| sed -e "s/\([^ ]*\) \([^ ]* \)/\1 \2\\n/g"| sort -n |uniq

int
main()
{
  typedef int data_type;

  typedef std::pair<data_type,data_type> coords_t;

  typedef adjacency_list < vecS, vecS, undirectedS, coords_t , data_type  > Graph;
  Graph g;
  typedef graph_traits < Graph >::vertex_descriptor Vertex;
  typedef std::map < coords_t, Vertex > CoordsVertexMap;
  CoordsVertexMap sensors;
  

  std::string lineBuffer;

  while(std::getline(std::cin, lineBuffer)){
    std::istringstream tmp(lineBuffer);
    std::istream_iterator<data_type> b(tmp),e ;
    data_type axisId(*b); ++b;
    ++b; //drop nb sensors

    CoordsVertexMap::iterator pos;
    bool inserted;
    Vertex u,v;
    bool firstSensor(true);

    while(b!=e){
      if(firstSensor){
	data_type x(*b); ++b;
	data_type y(*b); ++b;
	coords_t c(std::make_pair(x,y));
	std::cout<<"reading sensor at"<<x<<" "<<y;
	tie(pos, inserted)=sensors.insert(std::make_pair(c,Vertex()));
	if(inserted){
	  u=add_vertex(g);
	  g[u]=c;
	  pos->second=u;
	  std::cout<<" new\n";
	}else {
	  u=pos->second;
	  std::cout<<"already present\n";
	}
	firstSensor=false;
      }else{
	u=v;
      }
      data_type x(*b); ++b;
      data_type y(*b); ++b;
      coords_t c(std::make_pair(x,y));
      std::cout<<"reading sensor at"<<x<<" "<<y;
      tie(pos, inserted)=sensors.insert(std::make_pair(c,Vertex()));
      if(inserted){
	v=add_vertex(g);
	g[v]=c;
	pos->second=v;
	std::cout<<" new\n";
      }else {
	v=pos->second;
	std::cout<<"already present\n";
      }
      graph_traits < Graph >::edge_descriptor e;
      tie(e, inserted) = add_edge(u, v, g);
      if (inserted) g[e] = axisId;
      
    }
  }
  std::pair<Graph::vertex_iterator, Graph::vertex_iterator> allSensors(vertices(g));
  std::cout<<"nb sensors:"<<std::distance(allSensors.first, allSensors.second);
  while(allSensors.first!=allSensors.second){
    std::cerr<<g[*(allSensors.first)].first<<" "<<g[*(allSensors.first)].second<<"\n";
    ++allSensors.first;
  }
  return 0;
}
