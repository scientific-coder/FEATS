#ifndef FEATS_UTILITY_MUTABLE_PRIORITY_QUEUE_HXX
#define  FEATS_UTILITY_MUTABLE_PRIORITY_QUEUE_HXX
#include <iostream>
#include <functional>
#include <tree.h>

//For bottom_up_segmentation, I need a mutable priority queue
// My implementation as a leftish tree was not fast enough, so
// decided to make a mutable priority_queue interface using rb_tree.
// I cannot use heap algorithms because they move data and I need to keep valid pointers
// (std::priotity_queue and std::multiset are both interfaces to an rb tree in g++ stdlib)
// tree sorts reverse order of priority queue. I take rbegin() but I should be more efficient to reverse the Compare.

namespace feats{
  namespace utility{


    template<typename T>
    struct identity: public std::unary_function<T,T>{
      const T& operator()( const T& t) const {return t;}
    };

    template<typename T
	     , typename Compare=std::less<T>

	     , typename Alloc=std::allocator<T> >
    struct mutable_priority_queue:    rb_tree<T
					      , T
					      , identity<T>
					      , Compare
					      , Alloc>{
      typedef rb_tree<T, T, identity<T>, Compare, Alloc> rb_tree_type;
      typedef T value_type;
      typedef Compare compare_type;
      typedef typename rb_tree_type::size_type size_type;
      typedef typename rb_tree_type::_Link_type node_type;

      explicit mutable_priority_queue(const compare_type& c=compare_type(), const Alloc& a=Alloc()):rb_tree_type(c,a){}
      const value_type& top()const{return *(rbegin());}
      void pop(){erase(rbegin());}
      node_type push(const value_type& v){ return (node_type)insert_equal(v)._M_node;}

      // LGPL code lifted from g++ libsdc++
      template<typename UpdateOp>
      void update(node_type n, UpdateOp op){
	n=(node_type)_Rb_tree_rebalance_for_erase(n,
						  _M_header->_M_parent,
						  _M_header->_M_left,
						  _M_header->_M_right);

	op(n->_M_value_field);
	node_type __y = _M_header;
	node_type __x = _M_root();
	while (__x != 0) {
	  __y = __x;
	  __x = _M_key_compare(n->_M_value_field, _S_key(__x)) ? 
            _S_left(__x) : _S_right(__x);
	}
	node_type __z(n);

	if (__y == _M_header || __x != 0 || 
	    _M_key_compare(n->_M_value_field, _S_key(__y))) {
	  _S_left(__y) = __z;               // also makes _M_leftmost() = __z 
	  //    when __y == _M_header
	  if (__y == _M_header) {
	    _M_root() = __z;
	    _M_rightmost() = __z;
	  }
	  else if (__y == _M_leftmost())
	    _M_leftmost() = __z;   // maintain _M_leftmost() pointing to min node
	}
	else {
	  _S_right(__y) = __z;
	  if (__y == _M_rightmost())
	    _M_rightmost() = __z;  // maintain _M_rightmost() pointing to max node
	}
	_S_parent(__z) = __y;
	_S_left(__z) = 0;
	_S_right(__z) = 0;
	_Rb_tree_rebalance(__z, _M_header->_M_parent);
      }
    };

  }
}
#endif
