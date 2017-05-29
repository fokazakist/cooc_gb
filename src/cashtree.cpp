#include "gspan.h"
#include <sstream>
#include <set>
#include <strstream>

void Ctree::print(){
  std::cout << g2tracers.size()<<":"<<pattern <<std::endl;
  if(children.size()!=0){
    for(list<Ctree*>::iterator it=children.begin();it!=children.end();++it){
      for(unsigned int i=1;i<pattern.size();i++){
	std::cout <<" ";
      }
      std::cout << "-";
      (*it)->print();
    }
  }
}
void CRoot::print(){
    for(list<Ctree*>::iterator it=one_edge_graphs.begin();it!=one_edge_graphs.end();++it){
      (*it)->print();
    }
}
void Tdelete(Ctree* tree){
  if(tree->children.size()!=0){
    for(list<Ctree*>::iterator it=tree->children.begin();it!=tree->children.end();++it){
      Tdelete(*it);
    }
  }
  delete tree;
}

void Gspan::Crun(){
  if(first_flag==true){
    first_tree_make();
    first_flag=false;
    std::cout << TNnum <<std::endl;
    return;
  }
  coocsearch();
  /*
  DFSCode  dcode;
  for(list<Ctree*>::iterator it=croot->one_edge_graphs.begin();it!=croot->one_edge_graphs.end();++it){
    dcode = (*it)->pattern[(*it)->pattern.size()-1];
    cash_tsearch(croot->heap[dcode.labels],**it);
  }
  */
}

void Gspan::cash_tsearch(GraphToTracers& g2tracers,Ctree& node){

  double gain=0.0;
  double upos=0.0;
  double uneg=0.0;
  gain=-wbias;
  upos=-wbias;
  uneg=wbias;

  for(GraphToTracers::iterator it=g2tracers.begin();it!=g2tracers.end();++it){
    int gid = it->first;
    gain += 2 * corlab[gid] * weight[gid];
    if(corlab[gid]>0){
      upos += 2 * weight[gid];
	}else{
      uneg += 2 * weight[gid];
    }
  }

  if(fabs(opt_pat.gain)>std::max(upos,uneg)) {
    return;
  }

  double gain_abs = fabs(gain);
  if(gain_abs > fabs(opt_pat.gain) || (fabs(gain_abs - fabs(opt_pat.gain)) < 1e-10 && node.pattern.size() < opt_pat.size)){
    opt_pat.gain = gain;
    opt_pat.size = node.pattern.size();
    opt_pat.locsup.clear();
    for(GraphToTracers::iterator it=g2tracers.begin();it!=g2tracers.end();++it){
      opt_pat.locsup.push_back(it->first);
    }
    std::ostrstream ostrs;
    ostrs <<node.pattern;
    ostrs << std::ends;
    opt_pat.dfscode = ostrs.str();
  }
  if(node.children.size()==0){
    edge_grow(node);
    return;
  }
  DFSCode  dcode;
  Pair pkey;
  for(list<Ctree*>::iterator it=node.children.begin();it!=node.children.end();++it){
    dcode = (*it)->pattern[(*it)->pattern.size()-1];
    if(dcode.labels.z == -1){
      pkey.set(dcode.time.b,dcode.labels.y);
      cash_tsearch(node.b_heap[pkey],**it);
    }else{
      pkey.set(dcode.labels.y,dcode.labels.z);
      cash_tsearch(node.f_heap[dcode.time.a][pkey],**it);
    }
  }
}

void Gspan::first_tree_make(){
  /***   init CRoot         ***/
  croot = new CRoot;
  croot->one_edge_graphs.resize(0);

  /****  construct CRoot   ****/
  map<Triplet,GraphToTracers>& heap = croot->heap;
    p_count = 0;
  for(unsigned int gid = 0; gid < gdata.size(); ++gid){
    EdgeTracer cursor;
    Triplet wild_edge;
    Graph& g = gdata[gid];

    for(unsigned int v=0; v<g.size(); ++v){
      for(vector<Edge>::iterator e = g[v].begin(); e != g[v].end(); ++e){

	if (e->labels.x <= e->labels.z){
	  cursor.set(v,e->to,e->id,0);
	  heap[e->labels][gid].push_back(cursor);

	  if(wildcard_r>0){
	    wild_edge = e->labels;
	    wild_edge.z =999;
	    heap[wild_edge][gid].push_back(cursor);
	    wild_edge = e->labels.reverse();
	    wild_edge.z =999;
	    cursor.set(e->to,v,e->id,0);
	    heap[wild_edge][gid].push_back(cursor);
	  }
	}
      }
    }
  }
  pattern.resize(1);
  for(map<Triplet,GraphToTracers>::iterator it = croot->heap.begin(); it != croot->heap.end(); ++it){		
    pattern[0].labels = it->first;
    pattern[0].time.set(0,1);
    edge_grow(it->second,*croot);
    pattern.resize(1);
  }
  std::cout << p_count << std::endl;
}

void Gspan::coocsearch(){
  //this function is start position to search cooc pattern
  cooc_is_opt = false;
  opt_pat_cooc.dfscode.resize(0);
  opt_pat_cooc.locsup.resize(0);
  opt_pat_cooc.size_sum=0;
  opt_pat_cooc.gain=0.0;
  /* calc opt-pat and all gain,mac_gain */
  DFSCode  dcode;
  for(list<Ctree*>::iterator it=croot->one_edge_graphs.begin();it!=croot->one_edge_graphs.end();++it){
    dcode = (*it)->pattern[(*it)->pattern.size()-1];
    gcalc_tsearch(croot->heap[dcode.labels],**it);
  }
  if(need_to_cooc ==true){
    for(list<Ctree*>::iterator it=croot->one_edge_graphs.begin();it!=croot->one_edge_graphs.end();++it){
      dcode = (*it)->pattern[(*it)->pattern.size()-1];
      cooc_tsearch(croot->heap[dcode.labels],**it);
    }
  }
}

void Gspan::gcalc_tsearch(GraphToTracers& g2tracers,Ctree& node){
  //std::cout <<node.pattern<<std::endl;
  double gain=0.0;
  double upos=0.0;
  double uneg=0.0;
  gain=-wbias;
  upos=-wbias;
  uneg=wbias;

  for(GraphToTracers::iterator it=g2tracers.begin();it!=g2tracers.end();++it){
    int gid = it->first;
    gain += 2 * corlab[gid] * weight[gid];
    if(corlab[gid]>0){
      upos += 2 * weight[gid];
	}else{
      uneg += 2 * weight[gid];
    }
  }
  node.gain = gain;
  node.max_gain = std::max(upos,uneg);
  if(fabs(opt_pat.gain) > node.max_gain) {
    return;
  }
  double gain_abs = fabs(gain);

  if(gain_abs > fabs(opt_pat.gain) || (fabs(gain_abs - fabs(opt_pat.gain)) < 1e-10 && (node.pattern.size() < opt_pat.size||(node.pattern.size() == opt_pat.size && !(is_nomal) && opt_pat.new_node)))){

    opt_pat.gain = gain;
    opt_pat.size = node.pattern.size();
    opt_pat.new_node = false;
    opt_pat.locsup.clear();
    for(GraphToTracers::iterator it=g2tracers.begin();it!=g2tracers.end();++it){
      opt_pat.locsup.push_back(it->first);
    }
    std::ostrstream ostrs;
    ostrs <<node.pattern;
    ostrs << std::ends;
    opt_pat.dfscode = ostrs.str();
  }
  if(node.children.size()==0){
    edge_grow(node);
    return;
  }
  DFSCode  dcode;
  Pair pkey;
  
  for(list<Ctree*>::iterator it=node.children.begin();it!=node.children.end();++it){
    dcode = (*it)->pattern[(*it)->pattern.size()-1];
    if(dcode.labels.z == -1){
      pkey.set(dcode.time.b,dcode.labels.y);
      gcalc_tsearch(node.b_heap[pkey],**it);
    }else{
      pkey.set(dcode.labels.y,dcode.labels.z);
      gcalc_tsearch(node.f_heap[dcode.time.a][pkey],**it);
    }
  }
  
}

void Gspan::cooc_tsearch(GraphToTracers& g2tracers,Ctree& node){
  double cgain;
  if(cooc_is_opt == true){
    cgain = abs(opt_pat_cooc.gain);
  }else{
    cgain = fabs(opt_pat.gain);
  }
  if(cgain-node.max_gain>=-1e-10)  return;
  
  DFSCode  dcode;
  bool reach = false;
  for(list<Ctree*>::iterator it=croot->one_edge_graphs.begin();it!=croot->one_edge_graphs.end();++it){
    if(reach == true) break;
    dcode = (*it)->pattern[(*it)->pattern.size()-1];
    reach = cooc_tsearch(g2tracers,node,croot->heap[dcode.labels],**it);
  }

  Pair pkey;
  for(list<Ctree*>::iterator it=node.children.begin();it!=node.children.end();++it){
    dcode = (*it)->pattern[(*it)->pattern.size()-1];
    if(dcode.labels.z == -1){
      pkey.set(dcode.time.b,dcode.labels.y);
      cooc_tsearch(node.b_heap[pkey],**it);
    }else{
      pkey.set(dcode.labels.y,dcode.labels.z);
      cooc_tsearch(node.f_heap[dcode.time.a][pkey],**it);
    }
  }
}

bool Gspan::cooc_tsearch(GraphToTracers& base_g2tracers,Ctree& base,GraphToTracers& cand_g2tracers,Ctree& cand){
  if(base.pattern == cand.pattern) {
    return true;
  }
  bool ancestor = true;
  for(unsigned int i = 0;i!=cand.pattern.size();++i){
    if(cand.pattern[i] != base.pattern[i]){
      ancestor = false;
      break;
    }
  }
  double cgain;
  if(cooc_is_opt == true){
    cgain = abs(opt_pat_cooc.gain);
  }else{
    cgain = fabs(opt_pat.gain);
  }
  if(cgain-cand.max_gain>=-1e-10)  return ancestor;

  if(ancestor == false){
    double gain=0.0;
    double upos=0.0;
    double uneg=0.0;
    gain=-wbias;
    upos=-wbias;
    uneg=wbias;
    
    vector<int> loctemp;
    loctemp.resize(0);
    GraphToTracers::iterator bit=base_g2tracers.begin();
    GraphToTracers::iterator cit=cand_g2tracers.begin();
    while(!(bit==base_g2tracers.end() || cit==cand_g2tracers.end())){
      if(bit->first == cit->first){
	int gid = bit->first;
	gain += 2 * corlab[gid] * weight[gid];
	if(corlab[gid]>0){
	  upos += 2 * weight[gid];
	}else{
	  uneg += 2 * weight[gid];
	}
	++bit;
	++cit;
	//std::cout << gid <<","<<std::endl;
	loctemp.push_back(gid);
      }else if(bit->first > cit->first){
	++cit;
      }else{
	++bit;
      }
    }
    double ogain = cooc_is_opt ? opt_pat_cooc.gain : opt_pat.gain;
    if(fabs(ogain) - std::max(upos,uneg) >=-1e-10) {
      return false;
    }
    
    double gain_abs = fabs(gain);
    if(gain_abs > fabs(ogain) || (fabs(gain_abs - fabs(ogain)) < 1e-10 && (cooc_is_opt && opt_pat_cooc.size_sum > cand.pattern.size() + base.pattern.size()))){
      opt_pat_cooc.gain = gain;
      opt_pat_cooc.size_sum = cand.pattern.size() + base.pattern.size();
      opt_pat_cooc.locsup.clear();
      opt_pat_cooc.locsup=loctemp;
      opt_pat_cooc.dfscode.resize(0);
      cooc_is_opt=true;
      std::ostrstream ostrs1;
      ostrs1 <<cand.pattern;
      ostrs1 << std::ends;
      opt_pat_cooc.dfscode.push_back(ostrs1.str());
      std::ostrstream ostrs2;
      ostrs2 <<base.pattern;
      ostrs2 << std::ends;
      opt_pat_cooc.dfscode.push_back(ostrs2.str());
    }
  }


  DFSCode  dcode;
  Pair pkey;
  bool reach = false;
  for(list<Ctree*>::iterator it=cand.children.begin();it!=cand.children.end();++it){
    if(reach == true) return true;
    dcode = (*it)->pattern[(*it)->pattern.size()-1];
    if(dcode.labels.z == -1){
      pkey.set(dcode.time.b,dcode.labels.y);
      reach = cooc_tsearch(base_g2tracers,base,cand.b_heap[pkey],**it);
    }else{
      pkey.set(dcode.labels.y,dcode.labels.z);
      reach = cooc_tsearch(base_g2tracers,base,cand.f_heap[dcode.time.a][pkey],**it);
    }
  }
  return reach;
}
