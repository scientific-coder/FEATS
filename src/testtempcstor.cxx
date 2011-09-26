#include <iostream>
#include <utility>

template<typename A> struct SA{
  SA(std::pair<A,A> p){r=p.first;}
  template<typename T>
  SA(std::pair<A,T> p){r=p.first; std::cerr<<"template cstor\n";}

  template<typename CA> SA(CA ca){r=ca;}
  A r;
};


template<typename A> struct SB: SA<A>{
  template<typename CA> SB(CA ca):SA<A>(ca){}
};

template<typename T>
struct S1:T{
  void f(){
    T::g();
  }
  void g(){std::cout<<"S1";}
};
template<typename T>
struct S2:S1<T>{
  void g(){std::cout<<"S2";}
};
int main(){
  /*
  typedef SB<int> t1;
  t1 test(1.5);
  
  S2<S2<S1> > s2;
  s2.f();
  */

  std::pair<int,int> p(1,1);
  std::pair<int,float> pt(1,.5);
  SA<int> s1(p), s2(pt);

  
  return 0;
}
