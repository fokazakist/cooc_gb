#include "gspan.h"
#include <sstream>
#include <set>
#include <strstream>


void Gspan::pattern_search(){
  if(first_flag==true){
    first_tree_search();
    first_flag=false;
    //std::cout << TNnum <<std::endl;
    return;
  }
  bsp.clear();
  cash_tree_search();
  
  bsp.pruning(fabs(opt_pat.gain));
  bsp.print();
}

void Gspan::first_tree_search(){
  /***   init CRoot         ***/
  croot = new CashTreeRoot;

  /****  construct CRoot   ****/
  map<Triplet,GraphToTracers> heap;
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
  for(map<Triplet,GraphToTracers>::iterator it = heap.begin(); it != heap.end(); ++it){		
    pattern[0].labels = it->first;
    pattern[0].time.set(0,1);
    croot->one_edge_graphs[it->first] = new CashTree;
    edge_grow(it->second,*croot->one_edge_graphs[it->first]);
    pattern.resize(1);
  }
  std::cout << p_count << std::endl;
}

void Gspan::edge_grow(GraphToTracers& g2tracers, CashTree& place){
  
  for(GraphToTracers::iterator it = g2tracers.begin();it != g2tracers.end();++it){
    //std::cout << it->first << std::endl;
    place.locsup.push_back(it->first);
  }
  //place.locsup.push_back()
  TNnum++;
  int wild_flag = 0;
  if(pattern[pattern.size()-1].labels.z == 999){
    wildcard_r--;
    wild_flag = 1;
  }

  p_count++;
  place.wildcard_res=wildcard_r;
  place.leaf = false;
  if(can_prune(place)) {
    if(wild_flag == 1){
      wildcard_r++;
      wild_flag = 0;
    }
    if(pattern.size() + 1 < maxpat){
      place.leaf = true;
    }
    return;
  }
  if(pattern.size() + 1 > maxpat){
    if(wild_flag == 1){
      wildcard_r++;
      wild_flag = 0;
    }
    return;
  }
  PairSorter b_heap;
  map<int,PairSorter,greater<int> > f_heap;
  int maxtoc = scan_gspan(g2tracers,b_heap,f_heap);
  place.maxtoc = maxtoc;
  // projecting...
  DFSCode  dcode;
  for(PairSorter::iterator it = b_heap.begin(); it != b_heap.end(); ++it){	
    dcode.labels = Triplet(-1,it->first.b,-1);
    dcode.time.set(maxtoc, it->first.a);
    pattern.push_back(dcode);
    if(support(it->second) < minsup) { pattern.pop_back(); continue; }
    if (!is_min()){ pattern.pop_back(); continue; }
    place.b_child[it->first] = new CashTree;
    edge_grow(it->second,*place.b_child[it->first]);
    pattern.pop_back();
  }
	
  for(map<int,PairSorter,greater<int> >::iterator it = f_heap.begin(); it != f_heap.end(); ++it){	
    for(PairSorter::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2){		
      dcode.labels = Triplet(-1,it2->first.a,it2->first.b);
      dcode.time.set(it->first,maxtoc+1);
      pattern.push_back(dcode);
      if(support(it2->second) < minsup) { pattern.pop_back(); continue; }
      if (!is_min()){ pattern.pop_back(); continue; }
      place.f_child[it->first][it2->first] = new CashTree;
      edge_grow(it2->second,*place.f_child[it->first][it2->first]);
      pattern.pop_back();
    }
  }
  if(wild_flag == 1){
    wildcard_r++;
    wild_flag = 0;
  }
}

bool Gspan::can_prune(GraphToTracers& g2tracers,CashTree& node){

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
  node.gain=gain;
  node.max_gain = std::max(upos,uneg);
  if(fabs(opt_pat.gain) > node.max_gain ) {return true;  }
  //if(fabs(opt_pat.gain) - std::max(upos,uneg) >= -1e-10 )   return true;

  double gain_abs = fabs(gain);
  if(gain_abs > fabs(opt_pat.gain) || (fabs(gain_abs - fabs(opt_pat.gain))<1e-10 && pattern.size() < opt_pat.size)){
    opt_pat.gain = gain;
    opt_pat.size = pattern.size();
    opt_pat.new_node = true;
    opt_pat.locsup.clear();
    for(GraphToTracers::iterator it=g2tracers.begin();it!=g2tracers.end();++it){
      opt_pat.locsup.push_back(it->first);
    }
    std::ostrstream ostrs;
    ostrs <<pattern;
    ostrs << std::ends;
    opt_pat.dfscode = ostrs.str();
  }
  return false;
}

bool Gspan::can_prune(CashTree& node){

  double gain=0.0;
  double upos=0.0;
  double uneg=0.0;

  gain=-wbias;
  upos=-wbias;
  uneg=wbias;


  for(vector<int>::iterator it = node.locsup.begin();it != node.locsup.end();++it){
    int gid = *it;
    gain += 2 * corlab[gid] * weight[gid];
    if(corlab[gid]>0){
      upos += 2 * weight[gid];
	}else{
      uneg += 2 * weight[gid];
    }
  }
  node.gain=gain;
  node.max_gain = std::max(upos,uneg);
  if(fabs(opt_pat.gain) > node.max_gain ) {return true;  }
  //if(fabs(opt_pat.gain) - std::max(upos,uneg) >= -1e-10 )   return true;

  double gain_abs = fabs(gain);
  if(gain_abs > fabs(opt_pat.gain) || (fabs(gain_abs - fabs(opt_pat.gain))<1e-10 && pattern.size() < opt_pat.size)){
    opt_pat.gain = gain;
    opt_pat.size = pattern.size();
    opt_pat.new_node = true;
    opt_pat.locsup.clear();
    opt_pat.locsup=node.locsup;
    
    std::ostrstream ostrs;
    ostrs <<pattern;
    ostrs << std::ends;
    opt_pat.dfscode = ostrs.str();
  }
  return false;
}

void Gspan::cash_tree_search(){
  //and rewrite gain and max_gain of all patterns
  //and output patterns that are not searched theirs children
  pattern.resize(1);
  for(map<Triplet,CashTree*>::iterator edge = croot->one_edge_graphs.begin();edge != croot->one_edge_graphs.end();++edge){
    pattern[0].labels = edge->first;
    pattern[0].time.set(0,1);
    search_child(*edge->second);
  }
}
void Gspan::search_child(CashTree& cnode){
  //std::cout << cnode.locsup.size() << " : " << pattern << std::endl;
  
  if(can_prune(cnode)) return;
  if(cnode.leaf){
    //std::cout << "need to search below" << std::endl;
    bsp.insert(pattern,cnode.max_gain);
    std::cout << cnode.max_gain << " : " << pattern << std::endl;
    return;
  }

  DFSCode  dcode;
  for(map<Pair,CashTree*>::iterator it = cnode.b_child.begin(); it != cnode.b_child.end(); ++it){	
    dcode.labels = Triplet(-1,it->first.b,-1);
    dcode.time.set(cnode.maxtoc, it->first.a);
    pattern.push_back(dcode);
    search_child(*it->second);
    pattern.pop_back();
  }
	
  for(map<int,map<Pair,CashTree*>,greater<int> >::iterator it = cnode.f_child.begin(); it != cnode.f_child.end(); ++it){   
    for(map<Pair,CashTree*>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2){		
      dcode.labels = Triplet(-1,it2->first.a,it2->first.b);
      dcode.time.set(it->first,cnode.maxtoc+1);
      pattern.push_back(dcode);
      search_child(*it2->second);
      pattern.pop_back();
    }
  }
  
}
