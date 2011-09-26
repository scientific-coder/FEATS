#ifndef FEATS_UTILITY_STATIC_CASTER_HXX
#define  FEATS_UTILITY_STATIC_CASTER_HXX

#include <functional>

namespace feats{
  namespace utility{
    template<typename From, typename To>
    struct static_caster:std::unary_function<From,To>{
      To operator()(const From& f)const {return static_cast<To>(f);}
    };

    template<typename From, typename To>
    struct static_caster<From,To*>:std::unary_function<From,To*>{
      To* operator()(const From& f)const {return new To(f);}
    };
  }
}

#endif
