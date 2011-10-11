#ifndef KAHAN_SUM_HXX
#define KAHAN_SUM_HXX
#include <tuple>
namespace kahan {
  template<typename In, typename V> V accumulate(In it, In end, V res){
    for(V c(0.); it != end; ++it) {
      V const y(*it - c), t(res + y);
      c= (t - res) - y;
      res= t;
    }
    return res;
  }
  template<typename In, typename V> std::tuple<V, V> 
  sum_sumsq(In it, In end, std::tuple<V, V> res){
    for(V cs(0.), csq(0.); it != end; ++it) {
      V y(*it - cs), t(std::get<0>(res) + y);
      cs= (t - std::get<0>(res)) - y;
      std::get<0>(res)= t;
      y= ((*it) * (*it) - csq);
      t= (std::get<1>(res) + y);
      csq= (t - std::get<1>(res)) - y;
      std::get<1>(res)= t;
    }
    return res;
  }
}
#endif
