#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <random>
#include <tr1/random> // needed for std::tr1::variate_generator

//#include <chrono>
#include <boost/progress.hpp>

#include "feats/utility/kahan_sum.hxx"

//g++ -std=c++0x test-kahan-sum.cxx -I../include -o test-kahan-sum -O4 -march=native -Wall
template<typename In, typename V> std::tuple<V, V> 
sum_sumsq(In it, In end, std::tuple<V, V> res){
  while(it != end) {
    std::get<0>(res)+= *it;
    std::get<1>(res)+= (*it) * (*it);
    ++it;
  }
  return res;
}

int main(int argc, char* argv[]){
    std::tr1::variate_generator<std::mt19937
                                , std::tr1::normal_distribution<double> > 
      g(std::mt19937(), (std::tr1::normal_distribution<double>(10e3, 10e-7))) ;
    std::vector<double> d(100000000);
    std::tuple<double, double> s_sq, init_s_sq(0.,0.);
   std::generate_n(d.begin(), d.size(), g);
    std::cout<< std::fixed << std::setprecision(19);
    {
      boost::progress_timer t;
      std::cout<<std::accumulate(d.begin(), d.end(), 0.) << std::endl;
    }
    {
      boost::progress_timer t;
      std::cout<<kahan::accumulate(d.begin(), d.end(), 0.) << std::endl;
    }

    {
      boost::progress_timer t;
      s_sq= sum_sumsq(d.begin(), d.end(), init_s_sq);
      std::cout<<std::get<0>(s_sq)<<'\t'<<std::get<1>(s_sq) << std::endl;
    }
    {
      boost::progress_timer t;
      s_sq= kahan::sum_sumsq(d.begin(), d.end(), init_s_sq);
      std::cout<<std::get<0>(s_sq)<<'\t'<<std::get<1>(s_sq) << std::endl;
    }
    return 0;
}
