#include <map>
#include <vector>
#include <iostream>

using std::map;
using std::vector;

struct Inode { public: double info; bool is_Term; void set(double,bool);};
inline void Inode::set(double _info,bool _terminalornot){ info = _info; is_Term = _terminalornot;}

template <typename T>
class ptree{
 private:
  void print(int,vector<T>);
  void insert(vector<T>,double,int);
  double culcchild(int,double);
 public:
  const static int root;
  map<int,Inode> node_info;
  map<int,map<int,vector<T> > > edge_label;
  int newnode;
  
  ptree(){Inode r;  r.set(0,false);  node_info[1] = r;  newnode = 1;};
  void insert(vector<T>,double);
  void pruning(double);
  void print();
  void clear();
};

template <typename T>
void ptree<T>::clear(){
  node_info.clear();
  edge_label.clear();
  Inode r;
  r.set(0,false);
  node_info[1] = r;
  newnode = 1;
}

template <typename T>
const int ptree<T>::root = 0;

template <typename T>
void ptree<T>::print(){
  std::cout << "Node : " <<node_info.size() << std::endl;
  vector<T> t;t.resize(0);
  print(root,t);
}

template <typename T>
void ptree<T>::print(int n,vector<T> pattern){
  //std::cout << "node : "<<n << std::endl;
  //std::cout << "pattern size : "<<pattern.size() << std::endl;
  //std::cout << node_info.size() << std::endl;
  vector<T> pat=pattern;
  int ps = pattern.size();
  if(node_info[n].is_Term){
    std::cout << node_info[n].info << " : ";
    for(auto it = pattern.begin(); it != pattern.end();++it){
      std::cout << *it << " ";
    }
    std::cout << std::endl;
  }
  //std::cout << n << std::endl;
  for(auto it = edge_label[n].begin();it != edge_label[n].end();++it){
    std::copy(it->second.begin(),it->second.end(),std::back_inserter(pat));
    print(it->first,pat);
    pat.resize(ps);
  }
  
}

template <typename T>
void ptree<T>::insert(vector<T> target,double info,int n){
  
  int tonode = 0;
  
  for(auto it = edge_label[n].begin();it != edge_label[n].end();++it){
    if(it->second[0] == target[0]){
      tonode = it->first;
    }
  }
  if(tonode == 0){
    Inode i;
    i.set(info,true);
    node_info[newnode]=i;
    edge_label[n][newnode] = target;
    newnode++;
    return;
  }
  unsigned int same_num=1;
  vector<T> newedge;
  vector<T> endedge;
  unsigned int lsize = edge_label[n][tonode].size();
  unsigned int tsize = target.size();
  for(unsigned int i = 1;i != target.size();++i){
    if(lsize <= i ) break;
    if(target[i]!=edge_label[n][tonode][i]){
      break;
    }
    same_num++;
  }
  if(tsize > same_num){
    for(auto i = same_num;i != tsize;++i){
      newedge.resize(tsize-same_num);
      newedge[i-same_num] = target[i];
    }
  }
  
  if(lsize > same_num){
    for(auto i = same_num;i != lsize;++i){
      endedge.resize(lsize - same_num);
      endedge[i - same_num] = edge_label[n][tonode][i];
    }
  }
  
  if(same_num == lsize){
    if(tsize == same_num){
      node_info[tonode].is_Term = true;
      node_info[tonode].info = info;
    }else{
      insert(newedge,info,tonode);
    }
  }else{
    if(tsize == same_num){
      Inode i;
      i.set(info,true);
      node_info[newnode]=i;
      edge_label[n][newnode] = target;
      edge_label[n].erase(tonode);
      edge_label[newnode][tonode] = endedge; 
      newnode++;
    }else{
      Inode i,j;
      i.set(0,false);
      node_info[newnode]=i;
      target.resize(same_num);
      edge_label[n][newnode] = target;
      edge_label[n].erase(tonode);
      edge_label[newnode][tonode] = endedge;
      newnode++;
      j.set(info,true);
      node_info[newnode]=j;
      edge_label[newnode-1][newnode] = newedge;
      newnode++;
    }
  }
}
template <typename T>
void ptree<T>::insert(vector<T> target,double info){
  insert(target,info,root);
}

template <typename T>
void ptree<T>::pruning(double th){
  culcchild(root,th);
}

template <typename T>
double ptree<T>::culcchild(int node,double th){
  double maxgain = th;
  if(node_info[node].is_Term){
    return node_info[node].info;
  }
  vector<int> rm;
  for(auto it = edge_label[node].begin();it != edge_label[node].end();++it){
    double childgain = culcchild(it->first,th);
    maxgain = std::max(maxgain,childgain);
    if(childgain < th){
      rm.push_back(it->first);
    }
  }
  for(vector<int>::iterator it = rm.begin();it != rm.end();++it){
    edge_label[node].erase(*it);
  }
  return maxgain;
}
