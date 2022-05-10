#ifndef LDP_TWO_LAYER_GRAPH_HXX
#define LDP_TWO_LAYER_GRAPH_HXX

#include "two_dimensional_variable_array.hxx"
namespace LPMP {

class LdpTwoLayerGraph{
public:
    LdpTwoLayerGraph(){}
    LdpTwoLayerGraph(const std::vector<std::array<std::size_t,2>>& edges,const std::vector<double>& inputEdgeCosts);

    struct edge {
        std::size_t head;
        double cost;
        std::size_t reverse_neighbor_index; 
    };


    void setForwardEdgeCost(std::size_t vertex,std::size_t neighborIndex,double value) {
        //todo asserts
        assert(vertex<numberOfInputs);
        assert(neighborIndex<forwardEdges[vertex].size());
        //forwardEdges[vertex][neighborIndex].cost = value;
        edge& e=forwardEdges[vertex][neighborIndex];
        e.cost=value;
        assert(e.head<numberOfOutputs);
        assert(e.reverse_neighbor_index<backwardEdges[e.head].size());
        backwardEdges[e.head][e.reverse_neighbor_index].cost=value;
        assert(backwardEdges[e.head][e.reverse_neighbor_index].head==vertex);
    }

    void updateForwardEdgeCost(std::size_t vertex,std::size_t neighborIndex,double value) {
        assert(vertex<numberOfInputs);
        assert(neighborIndex<forwardEdges[vertex].size());
        edge& e=forwardEdges[vertex][neighborIndex];
        e.cost += value;
        assert(e.head<numberOfOutputs);
        assert(e.reverse_neighbor_index<backwardEdges[e.head].size());
        backwardEdges[e.head][e.reverse_neighbor_index].cost+=value;
        assert(backwardEdges[e.head][e.reverse_neighbor_index].head==vertex);
    }

    void setBackwardEdgeCost(std::size_t vertex,std::size_t neighborIndex,double value) {
        //todo asserts
        assert(vertex<numberOfOutputs);
        assert(neighborIndex<backwardEdges[vertex].size());
        edge& e=backwardEdges[vertex][neighborIndex];
        e.cost=value;
        assert(e.head<numberOfInputs);
        assert(e.reverse_neighbor_index<forwardEdges[e.head].size());
        forwardEdges[e.head][e.reverse_neighbor_index].cost=value;
        assert(forwardEdges[e.head][e.reverse_neighbor_index].head==vertex);
    }

    void updateBackwardEdgeCost(std::size_t vertex,std::size_t neighborIndex,double value) {
        assert(vertex<numberOfOutputs);
        assert(neighborIndex<backwardEdges[vertex].size());
        edge& e=backwardEdges[vertex][neighborIndex];
        e.cost += value;
        forwardEdges[e.head][e.reverse_neighbor_index].cost+=value;
        assert(forwardEdges[e.head][e.reverse_neighbor_index].head==vertex);
    }

    double getForwardEdgeCost(std::size_t vertex,std::size_t neighborIndex) const{
        assert(vertex<numberOfInputs);
        assert(neighborIndex<forwardEdges[vertex].size());
        return forwardEdges[vertex][neighborIndex].cost;
    }

    std::size_t getForwardEdgeVertex(std::size_t vertex,std::size_t neighborIndex) const{
        assert(vertex<numberOfInputs);
        assert(neighborIndex<forwardEdges[vertex].size());
        return forwardEdges[vertex][neighborIndex].head;
    }

    double getBackwardEdgeCost(std::size_t vertex,std::size_t neighborIndex) const{
        assert(vertex<numberOfOutputs);
        assert(neighborIndex<backwardEdges[vertex].size());
        return backwardEdges[vertex][neighborIndex].cost;
    }

    std::size_t getBackwardEdgeVertex(std::size_t vertex,std::size_t neighborIndex) const{
        assert(vertex<numberOfOutputs);
        assert(neighborIndex<backwardEdges[vertex].size());
        return backwardEdges[vertex][neighborIndex].head;
    }


    const edge* forwardNeighborsBegin(std::size_t i)const {
        assert(i<numberOfInputs);
        return forwardEdges[i].begin();
    }

    const edge* forwardNeighborsEnd(std::size_t i)const {
        assert(i<numberOfInputs);
        return forwardEdges[i].end();
    }



    edge* forwardNeighborsBegin(std::size_t i){
        assert(i<numberOfInputs);
        return forwardEdges[i].begin();
    }

    edge* forwardNeighborsEnd(std::size_t i){
        assert(i<numberOfInputs);
        return forwardEdges[i].end();
    }


    const std::size_t & getNumberOfInputs()const{
        return numberOfInputs;
    }



    const edge* backwardNeighborsBegin(std::size_t i)const {
        assert(i<numberOfOutputs);
        return backwardEdges[i].begin();
    }

    const edge* backwardNeighborsEnd(std::size_t i)const {
        assert(i<numberOfOutputs);
        return backwardEdges[i].end();
    }



    edge* backwardNeighborsBegin(std::size_t i){
        assert(i<numberOfOutputs);
        return backwardEdges[i].begin();
    }

    edge* backwardNeighborsEnd(std::size_t i){
        assert(i<numberOfOutputs);
        return backwardEdges[i].end();
    }


    const std::size_t & getNumberOfOutputs()const{
        return numberOfOutputs;
    }


private:

    //two_dim_variable_array<std::pair<std::size_t,double>> forwardEdges;
    two_dim_variable_array<edge> forwardEdges;
    two_dim_variable_array<edge> backwardEdges;
    std::size_t numberOfInputs;
    std::size_t numberOfOutputs;

};

}



#endif // LDP_TWO_LAYER_GRAPH_HXX
