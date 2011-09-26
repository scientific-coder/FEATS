#ifndef FEATS_UTILITY_NEARLY_EQUALS_HXX
#define FEATS_UTILITY_NEARLY_EQUALS_HXX

#include <cmath>
#include <functional>

namespace feats{
  namespace utility{

    template<typename Float> struct nearly_equals:std::binary_function<Float, Float, bool>{
      nearly_equals(const Float& ratio=.001):ratio_(ratio){} // default value could use std::numeric_limits<Float>::digits
      bool operator()(const Float& f1, const Float& f2)const{ return (std::abs(f1/f2-1.)<=ratio_ )||((f1<ratio_)&&(f2<ratio_));}
    private:
      const Float ratio_;
    };
  }
}
#endif
