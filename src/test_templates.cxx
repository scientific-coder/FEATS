template<typename A, typename B = int> struct C1{};
template<template<typename A, typename B> class D> struct C2{
  typedef D<int> the_d;
};

int main(){
  C2<C1> c;
  return 0;
}
