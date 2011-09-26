#ifndef NA_READER_HXX
#define NA_READER_HXX

#include <string>
#include <limits>

namespace feats{
  namespace utility{
// this class is used to read numeric fp values from a character stream
// If the characters sequence NA is encontered std::numeric_limits<NumericType>::quietNaN() is returned instead of aborting.



template <typename NumericType>
class na_reader{
public:
  std::string naS;
  NumericType value; // ugly on s'en fout

  explicit na_reader(std::string naString=std::string("NA")):naS(naString){}
  operator NumericType () const{ return value;}
};

template<typename Istream, typename NumericType>
Istream& operator >>(Istream& is, na_reader<NumericType>& naR){
    std::string buffer;
    is>>buffer;
    if( buffer != naR.naS) {
      std::stringstream tmp(buffer); tmp >> naR.value;
    }
    else { naR.value=std::numeric_limits<NumericType>::quiet_NaN();}
    return is;
}
  }
}

#endif
