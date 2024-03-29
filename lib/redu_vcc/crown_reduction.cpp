
#include <algorithm>
#include <iostream>
#include <fstream>

#include "crown_reduction.h"



void crown_reduction::reduce(redu_vcc &reduVCC, NodeID &node_v, NodeID &node_u ){

  type = "crown";

  reduVCC.buildKernel();

  branch_and_reduce_algorithm b_and_r(reduVCC.kernel_adj_list, reduVCC.remaining_nodes);
  if (b_and_r.lpCrownReduction()){
    unsigned int curr_cliqueID = reduVCC.next_cliqueID;
    reduVCC.addCrownCliques( crown_cliques, b_and_r.crown_cliques);
    unsigned int num_crown = reduVCC.next_cliqueID - curr_cliqueID;
    num_cliques += num_crown;

  }
}

void crown_reduction::reduce(redu_vcc &reduVCC, vertex_queue *queue,
                           NodeID &node_v, NodeID &node_u ){
    reduce(reduVCC, node_v, node_u);

    for (std::vector<NodeID> clique : crown_cliques) {
      for (NodeID a : clique) queue->adjust_queue(reduVCC, a);;
    }
}

void crown_reduction::unreduce(redu_vcc &reduVCC){

  for (unsigned int i = crown_cliques.size(); i > 0; i--) {
    std::vector<NodeID> &clique = crown_cliques[i-1];
    reduVCC.pop_clique(clique);
    reduVCC.addVertexSet(clique);
  }
}
