#include <iostream>
#include <boost/lexical_cast.hpp>

#include "lea.hxx"
// g++ -std=c++0x testlea.cxx lea.cxx -O4 -I../include -I/home/bernard/Code/repositories/boost-trunk  -lstdc++ -o testlea
int main (int argc, char* argv[]){
	int nbProtos(boost::lexical_cast<int>(argv[1]));
	double tmp(1./nbProtos);
	std::cout<<-std::numeric_limits<double>::infinity()<<'\n';
	for(double d=tmp;nbProtos;--nbProtos, d+=tmp){ 
		std::cout<<stdnormal_inv(d)<<'\n';
	}
	return 0;
}
