#ifndef LDP_BATCH_PROCESS_HXX
#define LDP_BATCH_PROCESS_HXX
#include<stdlib.h>
#include"lifted_disjoint_paths/ldp_vertex_groups.hxx"
#include<map>
#include"andres-graph/include/andres/graph/digraph.hxx"
#include "fstream"
#include <chrono>
#include "ldp_functions.hxx"
#include "ldp_directed_graph.hxx"

namespace LPMP {



class LdpBatchProcess{

public:
    LdpBatchProcess(VertexGroups<>& shiftedGroups, std::vector<std::array<std::size_t,2>>& shiftedVertexLabels_, std::size_t maxLabelSoFar_, std::size_t maxTimeForLabeled_, std::size_t minTimeToUse, std::size_t maxTimeToUse);  //including maxTimeForLabeled_ will be taken as computed
    LdpBatchProcess(VertexGroups<>& shiftedGroups, std::vector<std::array<std::size_t,2>>& shiftedVertexLabels_, std::size_t maxLabelSoFar_, std::size_t maxTimeForLabeled_);  //including maxTimeForLabeled_ will be taken as computed

    void initEdgesFromVector(const std::vector<std::array<std::size_t,2> > &edges, const std::vector<double> &costs);
    void initVertexScoreFromVector(const std::vector<std::size_t>& vertexList,const std::vector<double>& costs);
    void initFromFile(std::string filename);
    std::size_t globalIndexToLocalIndex(const std::size_t &globalIndex);
    std::size_t localIndexToGlobalIndex(const std::size_t &localIndex);
    const andres::graph::Digraph<>& getOutputGraph(){
        assert(edgesCreated);
        return outputGraph;
    }

    const std::vector<double>& getEdgeScore()const {
        assert(edgesCreated);
        return outputEdgeCosts;
    }

    const std::vector<double>& getVerticesScore()const {
        assert(edgesCreated);
        return outputVerticesScore;
    }

    void decode(const std::vector<std::vector<std::size_t> > &paths);
    std::vector<std::array<std::size_t,2>> getDecodedLabels(){
        return decodedLabels;
    }

    const std::size_t& getMaxLabelsSoFar() const{
        return maxLabelSoFar;
    }

    std::size_t getIndexToDel() const{
        return indexToDel;
    }

    void createLocalVG(LPMP::VertexGroups<>& localVG);


    const std::chrono::steady_clock::time_point& getContructorBegin()const {
        return constructorBegin;
    }

    const LdpDirectedGraph& getMyCompleteGraph()const{
        return myOutputGraph;
    }

private:
    std::size_t minValidVertex;
    std::size_t maxTimeForLabeled;
    VertexGroups<>* pvg;
    //std::vector<double> origVerticesScore;
    std::size_t numberOfVerticesInBatch;
    std::size_t minVertex;
    std::size_t maxVertex;
    std::size_t minTime;
    std::size_t maxTime;
    std::vector<std::size_t> shiftedLabels;
    std::map<std::size_t,std::map<std::size_t,double>> edgesFromLabeled;
    std::size_t maxLabelSoFar;
    std::size_t numberOfUsedLabels;
    std::vector<std::size_t> localIndexToGlobalLabel;
    andres::graph::Digraph<> outputGraph;
    LdpDirectedGraph myOutputGraph;
    std::vector<double> outputEdgeCosts;
    std::vector<double> outputVerticesScore;
    std::size_t numberOfOutputVertices;
    bool edgesCreated;
    bool labelsDecoded;
    std::vector<std::array<std::size_t,2>> decodedLabels;
    std::size_t indexToDel;
    std::chrono::steady_clock::time_point constructorBegin;

};
//maybe vertex labels as vector<arra<std::size_t,2>> globID->label
//alternativelly, also set min time and max time as constructor parameters (if not set, taken from vg)
//max label so far: maybe during decoding?
}
#endif // LDP_BATCH_PROCESS_HXX
