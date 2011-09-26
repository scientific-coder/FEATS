#ifndef LEFTISH_HEAP_HXX
#define LEFTISH_HEAP_HXX

#include<functional>
#include <numeric>


namespace feats{
  namespace utility{


template<typename T, typename Compare=std::less<T> >
struct leftish_tree{
  typedef T value_type;
  typedef Compare compare_type;
  typedef unsigned int size_type;
  
private:
  struct Node{
    Node(value_type v):v_(v),left_(0), right_(0), parent_(0),dist_(0){}
    const value_type& operator*()const {return v_;}
    value_type& operator*(){return v_;} // should not change key

    value_type const*operator->()const{return &v_;}
    value_type *operator->(){return &v_;}// should not change key

    void remove_child(Node* child){ 
      if(left_==child) left_=right_; //right_ will be the one removed anyway so we swap
      right_=0; 
      child->parent_=0;
      dist_=0;
    }
    void remove_children(){
      if(left_)left_=left_->parent_=0;
      if(right_)right_=right_->parent_=0;
      dist_=0;
    }
    value_type v_;
    Node* left_;
    Node* right_;
    Node* parent_;
    size_type dist_;
  };
  
public:
  typedef Node* node_type;

  leftish_tree(const compare_type& c=compare_type()):root_(0), c_(c){}
  bool empty()const{return root_==0;}
  
  node_type push(const value_type& v){
    node_type const newNode= new Node(v);
    root_= merge(root_,newNode);
    return newNode;
  }
  node_type top_node()const{return root_;}
  const value_type& top()const{return root_->v_;}

  void pop(){
    //assert(!empty());
    node_type const tmp(root_);
    node_type const r(tmp->right_), l(tmp->left_);
    tmp->remove_children();
    root_=merge(l,r);
    delete tmp;
  }
  template<typename UpdateOp>
  void update(node_type n, UpdateOp op){
    op(n->v_);
    update(n);
  }
  ~leftish_tree(){ delete_rec(root_);} // recursive->stack intensive

private:
  void update(node_type n){
    // not very efficient, but simple:
    // remove node, merge the 3 trees (root_, left_, right_), insert node (with new value)
    // I have to have a parent link to update parent_ child :-(
    node_type r(n->right_), l(n->left_);
    n->remove_children();
    if(n->parent_){ 
      n->parent_->remove_child(n); //now root_ is a valid tree
      root_=merge(n,merge(l,merge(root_,r)));//1.43 2.75
    }else{ // n is root_
      root_=merge(merge(l,r),n);
    }   
  }

  node_type merge(node_type p, node_type q){
    // from Donald E.Knuth The Art of Computer Programming book 3, p 150
    node_type r(0);
    size_type d=std::numeric_limits<size_type>::max();
    do{
      if(q){
	if(p){
	  if(c_(q->v_,p->v_)){
	    node_type t(p->right_); p->right_=r; if(r)r->parent_=p; r=p; p=t;    
	  }else{	
	    node_type t(q->right_); q->right_=r; if(r)r->parent_=q; r=q; q=t;
	  }
	}else{ p=q; d=(q?q->dist_:0);}
      }else{d=(p?p->dist_:0);}
    }while(d==std::numeric_limits<size_type>::max());
    while(r){
      q=r->right_;
      size_type dLeftR((r->left_) ? r->left_->dist_ : 0);
      if(dLeftR < d){
	d=dLeftR+1;
	r->right_=r->left_;
	r->left_=p;
      }else{
	++d;
	r->right_=p;
      }
      if(p)p->parent_=r;
      r->dist_=d;
      p=r;
      r=q;
    }
    if(p)p->parent_=0;
    return p;
  }

  void delete_rec(node_type n){ if(n){ delete_rec(n->left_); delete_rec(n->right_); delete n;} }

  node_type root_;

  compare_type c_;
};

  }
}
#endif
