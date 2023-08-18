/******************************************************************************
 * vcc.cpp
 * *
 * Source of ReduVCC
 * Darren Strash <dstrash@hamilton.edu>
 * Louise Thompson <lmthomps@hamilton.edu>
 *****************************************************************************/
#include <argtable3.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <utility> // for std::pair
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <climits>
#include <chrono>

#include "gurobi_c.h"

#ifdef BETA
#include <boost/functional/hash.hpp>
#endif // BETA

#include "balance_configuration.h"
#include "data_structure/graph_access.h"
#include "data_structure/matrix/normal_matrix.h"
#include "data_structure/matrix/online_distance_matrix.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "mapping/mapping_algorithms.h"
#include "parse_parameters.h"
#include "partition/graph_partitioner.h"
#include "partition/partition_config.h"
#include "partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"
#include "mis/initial_mis/greedy_mis.h"

#include "mis/kernel/branch_and_reduce_algorithm.h"

#include "ccp/Chalupa/cli.h"
#include <time.h>

#include "redu_vcc/redu_vcc.h"
#include "redu_vcc/reducer.h"
#include "branch_and_reduce/b_and_r.h"

#include "sigmod_mis/Graph.h"
#include "mis/mis_config.h"

// Clique enumeration
#include "Algorithm.h"
#include "DegeneracyAlgorithm.h"

//ECC reductions
#include "ecc_reductions/ecc.hpp"
#include "ecc_reductions/cover.hpp"
#include "ecc_reductions/graph.hpp"
#include "ecc_reductions/ecc_to_vcc.hpp"


#ifdef BETA
namespace std
{
    template<> struct hash<std::pair<NodeID,NodeID>>
    {
        std::size_t operator()(std::pair<NodeID,NodeID> const & node_pair) const noexcept
        {
            std::size_t seed = 0;
            boost::hash_combine(seed, std::get<0>(node_pair));
            boost::hash_combine(seed, std::get<1>(node_pair));
            return seed;
        }
    };
};

void transform_graph(graph_access const &G,
                     graph_access       &transformed_graph,
                     std::vector<std::pair<NodeID,NodeID>> &vertex_to_edge) {

        // Step 1: acyclic orientation
        std::vector<NodeID> order;
        order.reserve(G.number_of_nodes());
        for (NodeID node = 0; node < G.number_of_nodes(); node++) {
            order.push_back(node);
        }

        // Step 2: line graph
        std::unordered_set<std::pair<NodeID,NodeID>> neighbor_hash(G.number_of_edges());
        forall_nodes(G, v) {
            forall_out_edges(G, e, v) {
                NodeID u = G.getEdgeTarget(e);
                neighbor_hash.insert(std::pair<NodeID,NodeID>(u,v));
                neighbor_hash.insert(std::pair<NodeID,NodeID>(v,u));
            } endfor
        } endfor

        vertex_to_edge.clear();
        vertex_to_edge.reserve(G.number_of_edges() / 2);
        std::unordered_map<std::pair<NodeID, NodeID>, NodeID> edge_to_vertex(G.number_of_edges() / 2);
        std::vector<std::vector<NodeID>> edge_to_edges(G.number_of_edges() / 2);

        std::vector<NodeID> const empty;

        std::cout << "Making transformed graph..." << std::endl;
        std::size_t full_num_edges = 0;
        forall_nodes(G, v) {
            forall_out_edges(G, e, v) {
                NodeID u = G.getEdgeTarget(e);
                std::pair<NodeID,NodeID> node_pair(u, v);
                if (order[v] < order[u]) {
                    node_pair = std::pair<NodeID,NodeID>(v, u);
                }
                if (edge_to_vertex.find(node_pair) == edge_to_vertex.end()) {
                    vertex_to_edge.push_back(node_pair);
                    edge_to_vertex[node_pair] = vertex_to_edge.size() - 1;
                    //edge_to_edges.push_back(empty);
                }
            } endfor

            forall_out_edges(G, e1, v) {
                NodeID u = G.getEdgeTarget(e1);
                std::pair<NodeID,NodeID> first_pair(v, u);
                if (order[u] < order[v]) {
                    first_pair = std::pair<NodeID,NodeID>(u, v);
                }

                forall_out_edges(G, e2, v) {
                    NodeID w = G.getEdgeTarget(e2);
                    if (w < u) continue;
                    if (w == u) continue;
                    std::pair<NodeID,NodeID> second_pair(v, w);
                    if (order[w] < order[v]) {
                        second_pair = std::pair<NodeID,NodeID>(w, v);
                    }

                    bool bad = std::get<0>(first_pair) == std::get<0>(second_pair);
                    std::pair<NodeID,NodeID> const bad_edge =
                        std::pair<NodeID,NodeID>(std::get<1>(first_pair), std::get<1>(second_pair));
                    bad = bad and (neighbor_hash.find(bad_edge) != neighbor_hash.end());
                    if (bad) continue;

                    // Step 3: filter triples, TODO/DS: appears to work
                    edge_to_edges[edge_to_vertex[first_pair]].push_back(edge_to_vertex[second_pair]);
                    edge_to_edges[edge_to_vertex[second_pair]].push_back(edge_to_vertex[first_pair]);
                    full_num_edges++;

                } endfor
            } endfor
        } endfor

        // Step 3: trim triples (integrated above)

        std::cout << "Sorting neighborhoods..." << std::endl;
        std::size_t num_edges = 0;
        for (std::vector<NodeID> &neighbors : edge_to_edges) {
            std::sort(neighbors.begin(), neighbors.end());
            num_edges += neighbors.size();
        }
        num_edges /= 2;

        std::cout << "Transformed Graph has " << edge_to_edges.size() << " vertices and " << num_edges << " edges" << std::endl;

        std::cout << "Constructing graph_access data structure..." << std::endl;

        transformed_graph.start_construction(edge_to_edges.size(), 2 * num_edges);
        for (NodeID v = 0; v < edge_to_edges.size(); v++) {
            NodeID shadow_node = transformed_graph.new_node();
            transformed_graph.setPartitionIndex(shadow_node, 0);
            transformed_graph.setNodeWeight(shadow_node, 1);
            for (NodeID const neighbor : edge_to_edges[v]) {
                EdgeID shadow_edge = transformed_graph.new_edge(shadow_node, neighbor);
                transformed_graph.setEdgeWeight(shadow_edge, 1);
            }
        }

        transformed_graph.finish_construction();

        //graph_io::writeGraph(transformed_graph, "transformed.graph");

        //graph_access new_transformed_graph;
        //graph_io::readGraphWeighted(new_transformed_graph, "transformed-sorted.graph");
}


void transform_is_to_clique_cover(
        graph_access &transformed_graph,
        std::vector<std::pair<NodeID, NodeID>> const &vertex_to_edge,
        std::size_t const num_vertices_original_graph,
        NodeID *solution_is,
        std::size_t const num_cliques,
        std::vector<std::vector<int>> &clique_cover) {

    std::cout << "Reconstructing clique cover..." << std::endl;
    clique_cover.clear();
    clique_cover.reserve(num_cliques);
    unordered_map<NodeID, NodeID> vertex_to_clique_id(num_cliques);
    std::vector<bool> covered(num_vertices_original_graph, false);

    forall_nodes(transformed_graph, v) {
        if (solution_is[v] == 1) {
            std::pair<NodeID,NodeID> const &edge = vertex_to_edge[v];
            //std::cout << "Edge " << std::get<0>(edge) << "," << std::get<1>(edge)
            //         << " is in the independent set" << std::endl;
            if (vertex_to_clique_id.find(std::get<0>(edge))
                    == vertex_to_clique_id.end()) {
                vertex_to_clique_id[std::get<0>(edge)] = clique_cover.size();
                clique_cover.push_back(std::vector<int>{std::get<0>(edge)});
                covered[std::get<0>(edge)] = true;
            }
            clique_cover[vertex_to_clique_id[std::get<0>(edge)]].push_back(std::get<1>(edge));
            covered[std::get<1>(edge)] = true;
        }
    } endfor

    for (NodeID v = 0; v < num_vertices_original_graph; v++) {
        if (!covered[v]) {
            //std::cout << "DS: Vertex " << v << " is uncovered...covering with singleton" << std::endl;
            clique_cover.emplace_back(std::vector<int>{v});
        }
    }
    std::cout << "clique_cover has size: " << clique_cover.size() << std::endl;
}

void run_ils(ils &ils_instance,
             PartitionConfig const &partition_config,
             graph_access &graph,
             timer &the_timer,
             std::size_t const is_offset) {
        std::cout << "Performing ILS..." << std::endl;
        MISConfig ils_config;
        ils_config.seed = partition_config.seed;
        ils_config.time_limit = partition_config.solver_time_limit;
        ils_config.force_cand = 4;
        ils_config.ils_iterations = UINT_MAX;
        ils_instance.perform_ils(ils_config, graph, ils_config.ils_iterations,
                                 the_timer.elapsed(), is_offset, partition_config.mis);
}

void run_peeling(graph_access &graph,
                 graph_access &peeled_graph,
                 std::vector<NodeID> &new_to_old_id,
                 std::size_t &cover_offset) {
    Graph mis_G;
    mis_G.read_graph(graph);
    std::vector<std::vector<NodeID>> kernel;
    cover_offset = 0;
    std::cout << "Running peeling..." << endl;
    std::vector<bool> in_initial_is;
    unsigned int res_mis = mis_G.near_linear_kernel_and_offset(kernel, new_to_old_id, in_initial_is, cover_offset);
    std::cout << "    Reduced " << graph.number_of_nodes() << " -> " << kernel.size() << " nodes" << std::endl;
    std::cout << "    Offset =  " << cover_offset << std::endl;

    std::size_t num_edges = 0;
    for (std::vector<NodeID> & neighbors : kernel) {
        std::sort(neighbors.begin(), neighbors.end());
        num_edges += neighbors.size();
    }
    num_edges /= 2;

    peeled_graph.set_partition_count(2);
    peeled_graph.start_construction(kernel.size(), 2 * num_edges);
    for (NodeID v = 0; v < kernel.size(); v++) {
        NodeID shadow_node = peeled_graph.new_node();
        peeled_graph.setPartitionIndex(shadow_node, in_initial_is[v] ? 1 : 0);
        peeled_graph.setNodeWeight(shadow_node, 1);
        for (NodeID const neighbor : kernel[v]) {
            EdgeID shadow_edge = peeled_graph.new_edge(shadow_node, neighbor);
            peeled_graph.setEdgeWeight(shadow_edge, 1);
        }
    }

    peeled_graph.finish_construction();
}
#endif // BETA

void make_graph_access(vector<vector<NodeID>> const &adjlist, graph_access & G) {
    // Create the adjacency array
    std::vector<int> xadj(adjlist.size() + 1);
    size_t m = 0;
    for (auto neighbors : adjlist) m += neighbors.size();
    std::vector<int> adjncy(m);
    unsigned int adjncy_counter = 0;
    for (unsigned int i = 0; i < adjlist.size(); ++i) {
        xadj[i] = adjncy_counter;
        for (int const neighbor : adjlist[i]) {
            if (neighbor == i) continue;
            if (neighbor == UINT_MAX) continue;
            adjncy[adjncy_counter++] = neighbor;
        }
        std::sort(std::begin(adjncy) + xadj[i], std::begin(adjncy) + adjncy_counter);
    }
    xadj[adjlist.size()] = adjncy_counter;

    // Build the graph
    G.build_from_metis(adjlist.size(), &xadj[0], &adjncy[0]);
}

int main(int argn, char **argv) {

    PartitionConfig partition_config;
    std::string graph_filename;

    bool is_graph_weighted = false;
    bool suppress_output   = false;
    bool recursive         = false;

    int ret_code = parse_parameters(argn, argv,
                                    partition_config,
                                    graph_filename,
                                    is_graph_weighted,
                                    suppress_output, recursive);

    if(ret_code) {
        return 0;
    }

    std::streambuf* backup = std::cout.rdbuf();
    std::ofstream ofs;
    ofs.open("/dev/null");
    if(suppress_output) {
        std::cout.rdbuf(ofs.rdbuf());
    }

    partition_config.LogDump(stdout);

    bool const RULE_ONE_ENABLED = true;
    bool const RULE_TWO_ENABLED = true;
    bool const RULE_THREE_ENABLED = false;
    bool const RULE_FOUR_ENABLED = false;

    // open file for reading
    std::ifstream in(graph_filename.c_str());
    if (!in) {
        std::cerr << "Error opening " << graph_filename << std::endl;
        return 1;
    }

    ECCGraph graph = ECCGraph(in);
    //std::cerr << "Done reading in graph.\n";

    timer total_timer;
    Cover cover = Cover(graph.n);

    std::string const basename = ""; //argv[1];
    //std::cerr << "basename=" << basename << std::endl;

    size_t const k = /** argc > 2 ? std::stoi(argv[2]) : **/ (graph.n * graph.n) / 4 + 1;

    size_t total_calls = 0;

    std::cout.setf( std::ios::fixed, std:: ios::floatfield );
    std::cout.precision(4);

    std::cout << "input_vertices=" << to_string(graph.n) << std::endl;
    std::cout << "input_edges=" << to_string(graph.e / 2) << std::endl;

    auto const start_ecc_reductions = std::chrono::high_resolution_clock::now();
    auto const start = start_ecc_reductions;
    apply_reductions(graph, cover, RULE_ONE_ENABLED, RULE_TWO_ENABLED, RULE_THREE_ENABLED, RULE_FOUR_ENABLED);

    //bool found_cover = compute_edge_clique_cover(graph, cover, basename, k, total_calls, RULE_ONE_ENABLED, RULE_TWO_ENABLED, RULE_THREE_ENABLED, RULE_FOUR_ENABLED);
    auto const end_ecc_reductions = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> const ecc_reduction_time = end_ecc_reductions- start_ecc_reductions;

    std::cout << "ecc_reduction_time=" << ecc_reduction_time.count() / 1000.0 << std::endl;
    std::cout << "ecc_reduction_offset=" << to_string(cover.cliques.size()) << std::endl;

    // convert to vcc
    auto const start_convert = std::chrono::high_resolution_clock::now();
    ECC2VCC converter(graph, cover);
    vector<vector<NodeID>> vcc_adjlist = converter.ecc_to_vcc();
    auto const end_convert = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> const time_to_convert = end_convert - start_convert;
    std::cout << "time_ecc_to_vcc=" << time_to_convert.count() / 1000.0 << std::endl;

    graph_access G;
    timer t;
    //graph_io::readGraphWeighted(G, graph_filename);


    if (partition_config.run_type == "Redu") {
        timer reduce_timer;
        redu_vcc reduVCC(vcc_adjlist);
        std::vector<unsigned int> iso_degree;
        iso_degree.assign(vcc_adjlist.size(), 0);
        std::vector<unsigned int> dom_degree;
        dom_degree.assign(vcc_adjlist.size(), 0);
        reducer R;
        R.exhaustive_reductions(reduVCC, iso_degree, dom_degree);

        double vcc_reduction_time = reduce_timer.elapsed();
        std::cout << "BEGIN_OUTPUT_FOR_TABLES" << std::endl;
        std::cout << "run_type=" << partition_config.run_type << std::endl;
        std::cout << "vcc_reduction_time=" << vcc_reduction_time << std::endl;

        //reduVCC.analyzeGraph(graph_filename, G, s, false /* don't check cover */);

        cout << "vcc_reduction_offset=" << to_string(R.get_cover_size_offset()) << endl;
        cout << "vcc_kernel_vertices=" << to_string(reduVCC.remaining_nodes) << endl;
        reduVCC.buildKernel();
        cout << "vcc_kernel_edges=" << to_string(reduVCC.kernel_edges / 2) << endl;
    } else if (partition_config.run_type == "ReduBnR") {
      redu_vcc reduVCC;
      branch_and_reduce B(vcc_adjlist, reduVCC, partition_config);

      vertex_queue *queue = nullptr;
      if (partition_config.redu_type == "cascading") queue = new vertex_queue(reduVCC);

      timer s;
      bool const finished = B.bandr(reduVCC, 0, queue, partition_config, total_timer);
      double bandr_time = s.elapsed();
      double total_time = total_timer.elapsed();
      //B.analyzeGraph(graph_filename, G, reduVCC, s);

      std::cout << "BEGIN_OUTPUT_FOR_TABLES" << std::endl;
      std::cout << "run_type=" << partition_config.run_type << std::endl;
      std::cout << "bnr_time=" << bandr_time << std::endl;
      std::cout << "total_time=" << total_time << std::endl;
      std::cout << "branch_count=" << to_string(B.branch_count) << std::endl;
      std::cout << "prune_count=" << to_string(B.prune_count) << endl;
      std::cout << "decompose_count=" << to_string(B.decompose_count) << std::endl;
      std::cout << "total_solution=" << to_string(cover.cliques.size() + reduVCC.clique_cover.size()) << std::endl;

      make_graph_access(vcc_adjlist, G);

      std::cout << "verified_cover=" << (reduVCC.validateCover(G) ? "passed" : "failed") << std::endl;
      std::cout << "optimal=" << (finished ? "yes" : "unknown") << std::endl;

      return 0;
    } else if (partition_config.run_type == "ReduIG") {
        timer reduce_timer;
        redu_vcc reduVCC(vcc_adjlist);
        std::vector<unsigned int> iso_degree;
        iso_degree.assign(vcc_adjlist.size(), 0);
        std::vector<unsigned int> dom_degree;
        dom_degree.assign(vcc_adjlist.size(), 0);
        reducer R;
        R.exhaustive_reductions(reduVCC, iso_degree, dom_degree);
        double const vcc_reduction_time = reduce_timer.elapsed();

        reduVCC.build_cover();
        double time_to_solution = 0.0;
        double time_without_ig = total_timer.elapsed();
        reduVCC.solveKernel(partition_config, total_timer, time_to_solution, cover.cliques.size() + R.get_cover_size_offset() /* clique cover offset */);
        //reduVCC.analyzeGraph(graph_filename, G, s, false /* don't check cover */);

        timer unwind_timer;
        R.unwindReductions(reduVCC, time_to_solution);
        double time_to_unwind = unwind_timer.elapsed();
        double total_time = total_timer.elapsed();

        std::cout << "BEGIN_OUTPUT_FOR_TABLES" << std::endl;
        std::cout << "run_type=" << partition_config.run_type << std::endl;
        //std::cout << "input_graph_vertices=" << G.number_of_nodes() << std::endl;
        //std::cout << "input_graph_edges=" << G.number_of_edges() / 2 << std::endl;
        std::cout << "vcc_reduction_time=" << vcc_reduction_time << std::endl;
        cout << "vcc_reduction_offset=" << to_string(R.get_cover_size_offset()) << endl;
        cout << "vcc_kernel_vertices=" << to_string(reduVCC.remaining_nodes) << endl;
        reduVCC.buildKernel();
        cout << "vcc_kernel_edges=" << to_string(reduVCC.kernel_edges / 2) << endl;
        std::cout << "total_time_to_best=" << time_to_solution << std::endl;
        std::cout << "ig_time=" << time_to_solution - time_without_ig << std::endl;
        std::cout << "time_to_best_without_unwind=" << total_time - time_to_unwind << std::endl;
        std::cout << "total_solution=" << to_string(cover.cliques.size() + reduVCC.clique_cover.size()) << std::endl;

        make_graph_access(vcc_adjlist, G);

        std::cout << "verified_cover=" << (reduVCC.validateCover(G) ? "passed" : "failed") << std::endl;
        std::cout << "optimal=" << (reduVCC.clique_cover.size() == partition_config.mis ? "yes" : "unknown") << std::endl;
        return 0;
    }
    else if (partition_config.run_type == "ReduILP") {
        double time_to_solution = 0.0;
        timer reduce_timer;
        redu_vcc reduVCC(vcc_adjlist);
        std::vector<unsigned int> iso_degree;
        iso_degree.assign(vcc_adjlist.size(), 0);
        std::vector<unsigned int> dom_degree;
        dom_degree.assign(vcc_adjlist.size(), 0);
        reducer R;
        R.exhaustive_reductions(reduVCC, iso_degree, dom_degree);

        double vcc_reduction_time = reduce_timer.elapsed();
        time_to_solution += vcc_reduction_time;

        //reduVCC.analyzeGraph(graph_filename, G, s, false /* don't check cover */);

        //cout << "Reduced graph offset size: " << R.get_cover_size_offset() << endl;

        timer clique_timer;

        // DEBUG: Print the reduced graph
        // cout << "Reduced graph:" << endl;
        // reduVCC.printReducedGraph();
        // cout << endl;
        // reduVCC.printAdjList();
        // cout << endl;

        // Initialize vector of vertices contained in cliques
        // Converted to a dynamic array later for use in gurobi
        size_t reduced_size = reduVCC.remaining_nodes;
        vector<vector<int>> vertex_cliques(reduced_size); // For each vertex, the cliques continaing it
        vector<vector<int>> clique_vertices(0); // For each clique, the vertices it contains

        int current_clique = 0;
        DegeneracyAlgorithm *pAlgorithm(nullptr);
        vector<list<int>> adjlist(0);

        for (int i = 0; i < reduced_size; i++) {
            list<int> current_row(reduVCC.kernel_adj_list[i].begin(), reduVCC.kernel_adj_list[i].end());
            adjlist.push_back(current_row);
        }

        // forall_nodes(G, v) {
        //     forall_out_edges(G, e, v) {
        //         NodeID u = G.getEdgeTarget(e);
        //         adjlist[u].push_back(v);
        //     } endfor
        // } endfor

        pAlgorithm = new DegeneracyAlgorithm(adjlist);

        auto print_clique = [&vertex_cliques, &current_clique](list<int> const &clique) {
            list<int>::const_iterator it = clique.begin();

            // so vertices print 1-based, even though internally they
            // are stored 0-based.
            int offset = 1; //bOneBasedVertexIds ? 1 : 0;

            if (it != clique.end()) {
                cout << *it + offset;
                ++it;
            }
            while (it != clique.end()) {
                cout << " " << *it + offset;
                ++it;
            }
            cout << endl;

        };

        auto fill_vertex_cliques = [&vertex_cliques, &clique_vertices, &current_clique](list<int> const &clique) {
            list<int>::const_iterator it = clique.begin();
            vector<int> vertices(0);
            if (it != clique.end()) {
                vertex_cliques[*it].push_back(current_clique);
                vertices.push_back(*it);
                ++it;
            }
            while (it != clique.end()) {
                vertex_cliques[*it].push_back(current_clique);
                vertices.push_back(*it);
                ++it;
            }
            clique_vertices.push_back(vertices);
            current_clique++;

        };

        // pAlgorithm->AddCallBack(print_clique);
        list<list<int>> unused;
        //cout << "Filling vertex_cliques..." << endl;
        pAlgorithm->AddCallBack(fill_vertex_cliques);
        long num_cliques = pAlgorithm->Run(unused);
        pAlgorithm->PopCallBack();


        // DEBUG: Print contents of vertex_cliques and clique_vertices
        /*
        for (int i = 0; i < reduced_size; i++) {
            cout << "clique list of node " << i << ": ";
            for (int j = 0; j < vertex_cliques[i].size(); j++) {
                cout << vertex_cliques[i][j] << " ";
            }
            cout << endl;
        }

        for (int i = 0; i < num_cliques; i++) {
            cout << "vertex list of clique " << i << ": ";
            for (int j = 0; j < clique_vertices[i].size(); j++) {
                cout << clique_vertices[i][j] << " ";
            }
            cout << endl;
        }
        */

        //cout << "num maximal cliques: " << num_cliques << endl;

        delete pAlgorithm; pAlgorithm = nullptr;

	double clique_enumeration_time = clique_timer.elapsed();
	time_to_solution += clique_enumeration_time;

	timer ilp_setup_time;

        //

        GRBenv   *env   = NULL;
        GRBmodel *model = NULL;
        int error = 0;

        // Create Environment
        error = GRBloadenv(&env, "clique_gurobi.log");
        if (error) {
            cout << "Error in GRBloadenv" << endl;
        }

        error = GRBsetdblparam(env, GRB_DBL_PAR_TIMELIMIT, partition_config.solver_time_limit);
        if (error) {
            cout << "Error setting time limit in in GRBsetdblparam" << endl;
        }

	int threads = 1;

        error = GRBsetintparam(env, "Threads", threads);
        if (error) {
            cout << "Error setting threads in GRBsetintparam" << endl;
        }

        // Create an empty model
        error = GRBnewmodel(env, &model, "clique_gurobi", 0, NULL, NULL, NULL, NULL, NULL);
        if (error) {
            cout << "Error in GRBnewmodel" << endl;
        }

        // Add variables
        double * obj = new double[num_cliques];
        fill_n(obj, num_cliques, 1);
        char * types = new char[num_cliques];
        fill_n(types, num_cliques, 'B');
        error = GRBaddvars(model, num_cliques, 0, NULL, NULL, NULL, obj, NULL, NULL, types, NULL);
        if (error) {
            cout << "Error in GRBaddvars" << endl;
        }
        delete[] types;

        // Set Model Sense to Minimization
        error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE);
        if (error) {
            cout << "Error in GRBsetintattr" << endl;
        }

        // Add constraints
        for (int i = 0; i < reduced_size; i++) {

            // Initialize a vector of 1's to pass to gurobi as the "weights" indicating that
            // the cliques above correspond to vertex i
            double * weights = new double[vertex_cliques[i].size()];
            for (int j = 0; j < vertex_cliques[i].size(); j++) {
                weights[j] = 1;
            }
            // fill_n(weights, vertex_cliques[i].size(), 1);

            // Add the constraint using the weights and the cliques containing vertex i
            error = GRBaddconstr(model, vertex_cliques[i].size(), vertex_cliques[i].data(), weights, GRB_GREATER_EQUAL, 1.0, NULL);
            if (error) {
                cout << "Error in GRBaddconstr" << endl;
                //break;
            }
            // Free the memory taken by the dynamic arrays created above
            delete[] weights;
        }

	double ilp_solver_setup_time = clique_timer.elapsed();
	time_to_solution += ilp_solver_setup_time;

	timer ilp_solver_timer;

        error = GRBoptimize(model);

        if (error) {
            cout << "Error in GRBoptimize" << endl;
        }

        // Write out
        error = GRBwrite(model, "clique_gurobi.lp");
        if (error) {
            cout << "Error in GRBWrite" << endl;
        }

        // Deallocate Memory
        delete[] obj;

        double ilp_solver_time = ilp_solver_timer.elapsed();
        time_to_solution += ilp_solver_time;

        std::cout << "BEGIN_OUTPUT_FOR_TABLES" << std::endl;
        std::cout << "run_type=" << partition_config.run_type << std::endl;
        //std::cout << "input_graph_vertices=" << G.number_of_nodes() << std::endl;
        //std::cout << "input_graph_edges=" << G.number_of_edges() / 2 << std::endl;
        //std::cout << "total_time_to_best=" << time_to_solution << std::endl;
        std::cout << "vcc_reduction_time=" << vcc_reduction_time << std::endl;
        cout << "vcc_kernel_vertices=" << to_string(reduVCC.remaining_nodes) << endl;
        reduVCC.buildKernel();
        cout << "vcc_kernel_edges=" << to_string(reduVCC.kernel_edges / 2) << endl;
        std::cout << "clique_enumeration_time=" << clique_enumeration_time << std::endl;
        std::cout << "ilp_solver_setup_time=" << ilp_solver_setup_time << std::endl;
        std::cout << "ilp_solver_time=" << ilp_solver_time << std::endl;
        std::cout << "total_time=" << total_timer.elapsed() << std::endl;

	int grb_threads = 0;
        error = GRBgetintparam(env, "Threads", &grb_threads);
        if (error) {
            cout << "Error getting threads in GRBgetintparam" << endl;
        }

        std::cout << "ilp_solver_threads=" << to_string(grb_threads) << std::endl;
        std::cout << "vcc_reduction_offset=" << to_string(R.get_cover_size_offset()) << std::endl;
        std::cout << "vcc_kernel_maximal_cliques=" << to_string(num_cliques) << std::endl;
        double primal_objval = 0;
        //double dual_objval = 0;
        char *vtype;
        error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &primal_objval);
        //error = GRBgetdblattr(dual_model, GRB_DBL_ATTR_OBJVAL, &dual_objval);
        // error = GRBgetstrattr(dual_model, EFTYPE, &vtype);
        long ilp_solution = primal_objval;
        cout << "ilp_solution_of_kernel=" << to_string(ilp_solution) << endl;
        cout << "total_solution=" << to_string(cover.cliques.size() + R.get_cover_size_offset() + ilp_solution) << endl;
        //cout << "Objective Dual value of: " << dual_objval << endl;
        // cout << "Variable type: " << vtype << endl;
        //
        int status_of_ilp_solver = -1;
        error = GRBgetintattr(model, "Status", &status_of_ilp_solver);

        std::cout << "optimal=" << (status_of_ilp_solver==GRB_OPTIMAL ? "yes" : "unknown") << endl;
        if (status_of_ilp_solver != GRB_OPTIMAL) {
            std::cout << "ilp_solver_status=" << status_of_ilp_solver << std::endl;
        }
        return 0;
    }
#ifdef BETA
    else if (partition_config.run_type == "transform") {

        timer transformation_timer;
        // Step 1: acyclic orientation
        std::vector<NodeID> order;
        order.reserve(G.number_of_nodes());
        for (NodeID node = 0; node < G.number_of_nodes(); node++) {
            order.push_back(node);
        }

        // Step 2: line graph
        std::unordered_set<std::pair<NodeID,NodeID>> neighbor_hash;
        forall_nodes(G, v) {
            forall_out_edges(G, e, v) {
                NodeID u = G.getEdgeTarget(e);
                neighbor_hash.insert(std::pair<NodeID,NodeID>(u,v));
                neighbor_hash.insert(std::pair<NodeID,NodeID>(v,u));
            } endfor
        } endfor

        std::unordered_map<std::pair<NodeID, NodeID>, NodeID> edge_to_vertex;
        std::vector<std::pair<NodeID, NodeID>> vertex_to_edge;
        std::vector<std::vector<NodeID>> edge_to_edges;

        std::vector<NodeID> const empty;

        std::cout << "Making filtered line graph..." << std::endl;
        std::size_t full_num_edges = 0;
        forall_nodes(G, v) {
            forall_out_edges(G, e, v) {
                NodeID u = G.getEdgeTarget(e);
                std::pair<NodeID,NodeID> node_pair(u, v);
                if (order[v] < order[u]) {
                    node_pair = std::pair<NodeID,NodeID>(v, u);
                }
                if (edge_to_vertex.find(node_pair) == edge_to_vertex.end()) {
                    vertex_to_edge.push_back(node_pair);
                    edge_to_vertex[node_pair] = vertex_to_edge.size() - 1;
                    edge_to_edges.push_back(empty);
                }
            } endfor

            forall_out_edges(G, e1, v) {
                NodeID u = G.getEdgeTarget(e1);
                std::pair<NodeID,NodeID> first_pair(v, u);
                if (order[u] < order[v])
                    first_pair = std::pair<NodeID,NodeID>(u, v);

                forall_out_edges(G, e2, v) {
                    NodeID w = G.getEdgeTarget(e2);
                    if (w <= u) continue;
                    std::pair<NodeID,NodeID> second_pair(v, w);
                    if (order[w] < order[v]) {
                        second_pair = std::pair<NodeID,NodeID>(w, v);
                    }

                    // Step 3: filter triples, TODO/DS: appears to work
                    if (std::get<0>(first_pair) == std::get<0>(second_pair) &&
                        (neighbor_hash.find(std::pair<NodeID,NodeID>(std::get<1>(first_pair),std::get<1>(second_pair))) != neighbor_hash.end()) &&
                        (neighbor_hash.find(std::pair<NodeID,NodeID>(std::get<1>(second_pair),std::get<1>(first_pair))) != neighbor_hash.end())) continue;


                    edge_to_edges[edge_to_vertex[first_pair]].push_back(edge_to_vertex[second_pair]);
                    edge_to_edges[edge_to_vertex[second_pair]].push_back(edge_to_vertex[first_pair]);
        full_num_edges++;

                } endfor
            } endfor
        } endfor

        // Step 3: trim triples (integrated above)

        double transformation_time = transformation_timer.elapsed();

        // Step 4: Compute MIS
        if (false) // just peeling
        {
        Graph mis_G;
        mis_G.read_graph(edge_to_edges);
        unsigned int res_mis = mis_G.degree_two_kernal_dominate_lp_and_remove_max_degree_without_contraction();

        std::size_t const cover_upper_bound = G.number_of_nodes() - res_mis;
        std::cout << "Original Graph = " << G.number_of_nodes() << ", MIS=" << res_mis << std::endl;
        std::cout << "Transform (Near Linear) Upper Bound = " << cover_upper_bound << std::endl;
        }

        if (false) // just ILS
        {
        std::size_t num_edges = 0;
        for (std::vector<NodeID> &neighbors : edge_to_edges) {
            std::sort(neighbors.begin(), neighbors.end());
            num_edges += neighbors.size();
        }
        num_edges /= 2;

        std::cout << "Graph has " << edge_to_edges.size() << " vertices and " << num_edges << " edges" << std::endl;

        std::cout << "Constructing graph_access data structure..." << std::endl;

        graph_access new_graph;
        new_graph.start_construction(edge_to_edges.size(), 2 * num_edges);
        for (NodeID v = 0; v < edge_to_edges.size(); v++) {
            NodeID shadow_node = new_graph.new_node();
            new_graph.setPartitionIndex(shadow_node, 0);
            new_graph.setNodeWeight(shadow_node, 1);
            for (NodeID const neighbor : edge_to_edges[v]) {
                EdgeID shadow_edge = new_graph.new_edge(shadow_node, neighbor);
                new_graph.setEdgeWeight(shadow_edge, 1);
            }
        }

        new_graph.finish_construction();

        //graph_io::writeGraph(new_graph, "transformed.graph");

        //graph_access new_new_graph;
        //graph_io::readGraphWeighted(new_new_graph, "transformed-sorted.graph");

        timer ils_timer;
        greedy_mis greedy;
        greedy.initial_partition(partition_config.seed, new_graph);

        std::cout << "Preparing ILS..." << std::endl;
        MISConfig ils_config;
        ils_config.seed = partition_config.seed;
        ils_config.time_limit = partition_config.solver_time_limit;
        ils_config.force_cand = 4;
        ils_config.ils_iterations = UINT_MAX;
        ils ils_instance;
        std::cout << "Performing ILS..." << std::endl;
        ils_instance.perform_ils(ils_config, new_graph, ils_config.ils_iterations);
        cout << "ILS Running time: " << ils_timer.elapsed() << std::endl;

        std::size_t const ils_cover = G.number_of_nodes() - ils_instance.get_best_solution_size();

        std::cout << "Transform (ILS) Upper Bound = " << ils_cover << std::endl;
        }

        if (false) // just branch-and-reduce
        {
        std::cout << "Building a sorted int version of adjacency list..." << std::endl;
        vector<vector<int>> adjlist(edge_to_edges.size());
        for (NodeID v = 0; v < edge_to_edges.size(); v++) {
            std::sort(edge_to_edges[v].begin(), edge_to_edges[v].end());
            adjlist[v].reserve(edge_to_edges[v].size());
            for (NodeID const neighbor : edge_to_edges[v])
                adjlist[v].push_back(neighbor);
        }
        branch_and_reduce_algorithm mis_bnr(adjlist, adjlist.size());

        timer bnr_timer;
        bool timeout = mis_bnr.solve(bnr_timer, partition_config.solver_time_limit) == -1;

        std::size_t const exact_cover = G.number_of_nodes() - mis_bnr.get_current_is_size();

        std::cout << "BNR Running time: " << bnr_timer.elapsed() << std::endl;
        std::cout << "Transform (BNR) Upper Bound = " << exact_cover << std::endl;
        std::cout << "    Upper bound is " << (timeout ? "inexact" : "exact") << endl;
        }

        if (true) // peeling + ILS
        {
        timer total_timer;
        Graph mis_G;
        mis_G.read_graph(edge_to_edges);
        std::vector<std::vector<NodeID>> kernel;
        std::size_t offset = 0;

        std::cout << "Running reducing-peeling..." << endl;
        std::vector<bool> in_initial_is;
        timer peeling_timer;
        std::vector<NodeID> new_to_old_id_unused;
        unsigned int res_mis = mis_G.near_linear_kernel_and_offset(kernel, new_to_old_id_unused, in_initial_is, offset);
        double peeling_time = peeling_timer.elapsed();
        std::cout << "Done with reducing-peeling..." << endl;
        std::cout << "    Reduced " << edge_to_edges.size() << " -> " << kernel.size() << " nodes" << std::endl;
        std::cout << "    Offset =  " << offset << std::endl;

        std::size_t num_edges = 0;
        for (std::vector<NodeID> & neighbors : kernel) {
            std::sort(neighbors.begin(), neighbors.end());
            num_edges += neighbors.size();
        }
        num_edges /= 2;

        graph_access new_graph;
        new_graph.set_partition_count(2);
        new_graph.start_construction(kernel.size(), 2 * num_edges);
        for (NodeID v = 0; v < kernel.size(); v++) {
            NodeID shadow_node = new_graph.new_node();
            new_graph.setPartitionIndex(shadow_node, in_initial_is[v] ? 1 : 0);
            new_graph.setNodeWeight(shadow_node, 1);
            for (NodeID const neighbor : kernel[v]) {
                EdgeID shadow_edge = new_graph.new_edge(shadow_node, neighbor);
                new_graph.setEdgeWeight(shadow_edge, 1);
            }
        }

        new_graph.finish_construction();

        timer ils_timer;
        MISConfig ils_config;
        ils_config.seed = partition_config.seed;
        ils_config.time_limit = partition_config.solver_time_limit;
        ils_config.force_cand = 4;
        ils_config.ils_iterations = UINT_MAX;
        ils ils_instance;
        std::size_t solution_offset = G.number_of_nodes() - offset;
        ils_instance.perform_ils(ils_config, new_graph, ils_config.ils_iterations, total_timer.elapsed(), solution_offset, partition_config.mis);
        cout << "ILS Running time: " << ils_timer.elapsed() << std::endl;

        std::size_t const near_linear_ils_cover = G.number_of_nodes() - (ils_instance.get_best_solution_size() + offset);

        std::cout << "Transform (Near Linear / ILS) Upper Bound = " << near_linear_ils_cover << std::endl;

        std::cout << "BEGIN_OUTPUT_FOR_TABLES" << std::endl;
        std::cout << "run_type=" << partition_config.run_type << std::endl;
        std::cout << "transformation_time=" << transformation_time << std::endl;
        std::cout << "transformed_graph_vertices=" << edge_to_edges.size() << std::endl;
        std::cout << "transformed_graph_edges=" << full_num_edges << std::endl;
        std::cout << "peeling_time=" << peeling_time << std::endl;
        std::cout << "reduced_transformed_graph_vertices=" << new_graph.number_of_nodes() << std::endl;
        std::cout << "reduced_transformed_graph_edges=" << new_graph.number_of_edges() / 2 << std::endl;

        std::cout << "ils_time_to_best=" << ils_instance.get_last_update_time() << std::endl;
        std::cout << "total_time_to_best=" << ils_instance.get_last_total_update_time() << std::endl;
        std::cout << "clique_cover_size=" << near_linear_ils_cover << std::endl;
        std::cout << "verified_cover=no" << std::endl;
        std::cout << "optimal=" << (near_linear_ils_cover == partition_config.mis ? "yes":"unknown") << std::endl;
        }

        return 0;
    } else if (partition_config.run_type == "Redutransform") {

        timer total_timer;

        timer vcc_reductions_timer;
        redu_vcc oldVCC(G);
        std::vector<unsigned int> iso_degree;
        iso_degree.assign(G.number_of_nodes(), 0);
        std::vector<unsigned int> dom_degree;
        dom_degree.assign(G.number_of_nodes(), 0);
        reducer R;
        R.exhaustive_reductions(oldVCC, iso_degree, dom_degree);
        oldVCC.build_cover();

        std::size_t const cover_size_offset = R.get_cover_size_offset();
        std::cout << "cover_size_offset=" << cover_size_offset << std::endl;

        oldVCC.analyzeGraph(graph_filename, G, s, false /* don't check cover */);

        double vcc_reduction_time = vcc_reductions_timer.elapsed();

        oldVCC.buildKernel();

        graph_access kernel_graph;
        kernel_graph.start_construction(oldVCC.kernel_adj_list.size(), oldVCC.kernel_edges);
        for (NodeID v = 0; v < oldVCC.kernel_adj_list.size(); v++) {
            NodeID shadow_node = kernel_graph.new_node();
            kernel_graph.setNodeWeight(shadow_node, 1);
            for (NodeID const neighbor : oldVCC.kernel_adj_list[v]) {
                EdgeID shadow_edge = kernel_graph.new_edge(shadow_node, neighbor);
                kernel_graph.setEdgeWeight(shadow_edge, 1);
            }
        }

        kernel_graph.finish_construction();

        timer transformation_timer;
        graph_access transformed_graph;
        std::vector<std::pair<NodeID, NodeID>> vertex_to_edge;
        transform_graph(kernel_graph, transformed_graph, vertex_to_edge);

        double transformation_time = transformation_timer.elapsed();
        //graph_io::writeGraph(new_graph, "transformed.graph");

        // Step 4: Compute MIS
        if (true) // Exhaustive MIS Reductions + Sigmod + ILS
        {
        std::cout << "Making copy of transformed graph..." << std::endl;
        std::vector<std::vector<int>> line_copy;
        line_copy.reserve(transformed_graph.number_of_nodes());
        forall_nodes(transformed_graph, v) {
            std::vector<int> int_neighbors(transformed_graph.getNodeDegree(v));
            forall_out_edges(transformed_graph, e, v) {
                NodeID u = transformed_graph.getEdgeTarget(e);
                int_neighbors.push_back(u);
            } endfor
            line_copy.emplace_back(int_neighbors);
        } endfor

        timer mis_reductions;

        std::cout << "Initializing branch-and-reduce..." << std::endl;
        branch_and_reduce_algorithm mis_bnr(line_copy, line_copy.size());
        std::cout << "Applying exhaustive MIS reductions..." << std::endl;
        mis_bnr.reduce();
        std::cout << "After exhaustive MIS reductions: " << line_copy.size() << " -> " << mis_bnr.number_of_nodes_remaining() << std::endl;

        std::vector<NodeID> reverse_mapping;

        double mis_reduction_time = mis_reductions.elapsed();

        graph_access full_kernel;
        mis_bnr.convert_adj_lists(full_kernel, reverse_mapping);
        std::cout << "Preparing ILS on full kernel of filtered line graph..." << std::endl;

        Graph mis_G;
        mis_G.read_graph(full_kernel);
        std::vector<std::vector<NodeID>> kernel;
        std::size_t offset = 0;
        std::cout << "Running reducing-peeling..." << endl;
        std::vector<bool> in_initial_is;

        timer peeling_timer;
        std::vector<NodeID> new_to_old_id;
        unsigned int res_mis = mis_G.near_linear_kernel_and_offset(kernel, new_to_old_id, in_initial_is, offset);
        double peeling_time = peeling_timer.elapsed();
        std::cout << "Done with reducing-peeling..." << endl;
        assert(offset == 0);

        std::cout << "    Reduced " << full_kernel.number_of_nodes() << " -> " << kernel.size() << " nodes" << std::endl;
        std::cout << "    Offset =  " << offset << std::endl;

        std::size_t num_kernel_edges = 0;
        for (std::vector<NodeID> & neighbors : kernel) {
            std::sort(neighbors.begin(), neighbors.end());
            num_kernel_edges += neighbors.size();
        }
        num_kernel_edges /= 2;

        graph_access peeled_graph;
        peeled_graph.set_partition_count(2);
        peeled_graph.start_construction(kernel.size(), 2 * num_kernel_edges);
        for (NodeID v = 0; v < kernel.size(); v++) {
            NodeID shadow_node = peeled_graph.new_node();
            peeled_graph.setPartitionIndex(shadow_node, in_initial_is[v] ? 1 : 0);
            peeled_graph.setNodeWeight(shadow_node, 1);
            for (NodeID const neighbor : kernel[v]) {
                EdgeID shadow_edge = peeled_graph.new_edge(shadow_node, neighbor);
                peeled_graph.setEdgeWeight(shadow_edge, 1);
            }
        }
        peeled_graph.finish_construction();

        //std::cout << "Writing out post-peeling graph peeled.graph" << std::endl;
        //graph_io::writeGraph(peeled_graph, "peeled.graph");


////        cout << "cover_size_offset=" << cover_size_offset << std::endl;
////        cout << "kernel_graph_number_of_nodes=" << kernel_graph.number_of_nodes() << std::endl;
////        cout << "mis_bnr.get_is_offset=" << mis_bnr.get_is_offset() << std::endl;
////        cout << "offset=" << offset << std::endl;

        timer ils_timer;
        ils ils_instance;
        std::size_t const ils_print_offset = cover_size_offset + kernel_graph.number_of_nodes() - (mis_bnr.get_is_offset() + offset /* from peeling **/);
        run_ils(ils_instance, partition_config, peeled_graph, total_timer, ils_print_offset);
        double time_for_ils = ils_timer.elapsed();
        cout << "ILS Running time: " << time_for_ils << std::endl;

        std::size_t const ils_cover = kernel_graph.number_of_nodes() - (ils_instance.get_best_solution_size() + mis_bnr.get_is_offset() + offset);

        double time_to_unwind = 0.0;
        // reconstruct is / cover
        {
            timer unwind_timer;
            std::vector<bool> full_is(line_copy.size(), false);
            for (NodeID v = 0; v < peeled_graph.number_of_nodes(); v++) {
                if (ils_instance.get_best_solution()[v] == 1)
                    full_is[reverse_mapping[new_to_old_id[v]]] = true;
            }
            std::cout << "Unwinding MIS reductions..." << std::endl;
            mis_bnr.extend_finer_is(full_is);
            NodeID *solution_is = new NodeID[line_copy.size()];
            for (NodeID v = 0; v < line_copy.size(); v++) {
                solution_is[v] = full_is[v] ? 1 : 0;
            }

            std::vector<std::vector<int>> clique_cover;
            std::cout << "Transforming IS to VCC..." << std::endl;
            transform_is_to_clique_cover(transformed_graph, vertex_to_edge,
                kernel_graph.number_of_nodes(), solution_is, ils_cover, clique_cover);
////            std::cout << "Adding " << clique_cover.size() << " kernel cliques..." << std::endl;
            oldVCC.addKernelCliques(clique_cover);
            std::cout << "Unwinding VCC reductions..." << std::endl;
            R.unwindReductions(oldVCC);
            std::cout << "Done unwinding..." << std::endl;
            oldVCC.analyzeGraph(graph_filename, G, s, false /* don't check cover */);
            time_to_unwind = unwind_timer.elapsed();
        }

        double const time_to_best = ils_instance.get_last_total_update_time() + time_to_unwind;


        cout << "Total Running time (beginning to end): " << total_timer.elapsed() << std::endl;


        std::cout << "Transform (Exhaustive+Peeling+ILS) Upper Bound = " << ils_cover + cover_size_offset << std::endl;

        std::size_t const clique_cover_size = ils_cover + cover_size_offset;

        std::cout << "BEGIN_OUTPUT_FOR_TABLES" << std::endl;
        std::cout << "run_type=" << partition_config.run_type << std::endl;
        std::cout << "vcc_reduction_time=" << vcc_reduction_time << std::endl;
        std::cout << "transformation_time=" << transformation_time << std::endl;
        std::cout << "transformed_graph_vertices=" << transformed_graph.number_of_nodes() << std::endl;
        std::cout << "transformed_graph_edges=" << transformed_graph.number_of_edges() / 2 << std::endl;
        std::cout << "mis_reduction_time=" << mis_reduction_time << std::endl;
        std::cout << "reduced_transformed_graph_vertices=" << full_kernel.number_of_nodes() << std::endl;
        std::cout << "reduced_transformed_graph_edges=" << full_kernel.number_of_edges() / 2 << std::endl;
        std::cout << "peeling_time=" << peeling_time << std::endl;
        std::cout << "ils_time_to_best=" << ils_instance.get_last_update_time() << std::endl;
        std::cout << "total_time_to_best=" << time_to_best << std::endl;
        std::cout << "time_to_best_without_unwind=" << time_to_best - time_to_unwind << std::endl;
        std::cout << "clique_cover_size=" << clique_cover_size << std::endl;
        std::cout << "verified_cover=" << (oldVCC.validateCover(G) ? "passed" : "failed") << std::endl;
        std::cout << "optimal=" << (clique_cover_size == partition_config.mis ? "yes":"unknown") << std::endl;

        return 0;
    }
    else if (partition_config.run_type == "test") {

        timer total_timer;
        timer transformation_timer;

        graph_access transformed_graph;
        std::vector<std::pair<NodeID, NodeID>> vertex_to_edge;
        transform_graph(G, transformed_graph, vertex_to_edge);

        double transformation_time = transformation_timer.elapsed();

        // Step 4: Compute MIS
        if (true) // just ILS
        {
        timer ils_timer;
        greedy_mis greedy;
        greedy.initial_partition(partition_config.seed, transformed_graph);

        ils ils_instance;
        run_ils(ils_instance, partition_config, transformed_graph, total_timer, G.number_of_nodes());

        cout << "ILS Running time: " << ils_timer.elapsed() << std::endl;

        std::size_t const ils_cover = G.number_of_nodes() - ils_instance.get_best_solution_size();

        std::cout << "Transform (ILS) Upper Bound = " << ils_cover << std::endl;

        std::vector<std::vector<int>> clique_cover;
        size_t const num_cliques = G.number_of_nodes() - ils_instance.get_best_solution_size();

        transform_is_to_clique_cover(transformed_graph,
                                     vertex_to_edge,
                                     G.number_of_nodes(),
                                     ils_instance.get_best_solution(),
                                     num_cliques,
                                     clique_cover);
        }
        return 0;
    }


#endif // BETA

}
