/*
 * disjoint-paths-init.hxx
 *
 *  Created on: Sep 10, 2018
 *      Author: fuksova
 */

#ifndef INCLUDE_DISJOINT_PATHS_DISJOINT_PATHS_INIT_HXX_
#define INCLUDE_DISJOINT_PATHS_DISJOINT_PATHS_INIT_HXX_

#include <stdexcept>
#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <andres/graph/digraph.hxx>
#include <lifted_disjoint_paths/ldp_parameters.hxx>
#include <stack>
#include <unordered_set>
#include <iterator>
#include <unordered_map>
#include <string>
#include <bitset>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <deque>
#include <queue>
#include <set>
#include <map>
#include <list>
#include "ldp_functions.hxx"
#include <utility>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include "lifted_disjoint_paths/ldp_directed_graph.hxx"
#include "config.hxx"
#include "lifted_disjoint_paths/ldp_complete_structure.hxx"
#include "lifted_disjoint_paths/ldp_vertex_groups.hxx"
#include "ldp_batch_process.hxx"
#include <chrono>
#include "ldp_interval_connection.hxx"


namespace py = pybind11;
namespace LPMP{
namespace lifted_disjoint_paths {



class LdpInstance {
public:



     LdpInstance(LdpParameters<>& configParameters,CompleteStructure<>& cs);
     LdpInstance(LdpParameters<>& configParameters,LdpBatchProcess& BP);
     LdpInstance(LdpParameters<>& configParameters,LdpIntervalConnection& IC);
     LdpInstance(LdpParameters<>& configParameters, const py::array_t<std::size_t>& baseEdges, const py::array_t<std::size_t>& liftedEdges, const  py::array_t<double>& baseCosts, const  py::array_t<double>& liftedCosts, const py::array_t<double> &verticesCosts, VertexGroups<>& pvg);

	bool isReachable(std::size_t i,std::size_t j) const{
		if(i==t_||j==s_) return false;  //Assume no path from the terminal node
		if(i==s_||j==t_) return true;
		if(reachable.size()==0) return true;
		return reachable[i].count(j)>0;
	}

	const std::unordered_set<std::size_t>& reachableFromVertex(std::size_t v)const{
		return reachable.at(v);
	}



	const std::size_t getGapLifted() const {
		return parameters.getMaxTimeLifted();
	}

    const std::size_t getGapBase() const {
        return parameters.getMaxTimeBase();
    }

    void setOutputFileName(const std::string& fileName){
        if(!fileName.empty()){
            std::size_t lastDot=std::min(fileName.find_last_of("."),fileName.size());
            outputFilePrefix=fileName.substr(0,lastDot);
            outputFileSuffix=fileName.substr(lastDot,fileName.size());
            outputFileName=fileName;
        }
        else{
            outputFilePrefix="output";
            outputFileSuffix=".txt";
            outputFileName="output.txt";
        }
    }

    const std::string& getOutputFilePrefix()const{
        return outputFilePrefix;
    }

    const std::string& getOutputFileSuffix()const{
        return outputFileSuffix;
    }

    const std::string& getOutputFileName()const{
        return outputFileName;
    }

    const std::vector<std::unordered_set<std::size_t>>* getPReachable(){
        return &reachable;
    }

	std::size_t getSourceNode() const {
		return s_;
	}

	std::size_t getTerminalNode() const {
		return t_;
	}


	const std::vector<double>& getVerticesScore() const{
		return vertexScore;
	}


	double getVertexScore(std::size_t v) const {
                assert(v<vertexScore.size());
		return vertexScore[v];
	}

    const VertexGroups<std::size_t>& getVertexGroups()const {
        return vertexGroups;
	}

	std::size_t getGroupIndex(std::size_t v)const {
		return vertexGroups.getGroupIndex(v);
	}


	std::size_t getEdgeVarIndex(std::size_t edgeIndex)const {
		return edgeIndex+numberOfVertices;
	}


	std::size_t getLiftedEdgeVarIndex(std::size_t liftedEdgeIndex)const {
		return liftedEdgeIndex+numberOfEdges+numberOfVertices;
	}

	std::size_t getVertexVarIndex(std::size_t vertexIndex)const{
		return vertexIndex;
	}

    bool existLiftedEdge(const std::size_t v,const std::size_t w)const{
        bool value=liftedStructure.at(v).isWithinBounds(w)&&liftedStructure.at(v)[w]>0;
        return value;
    }


    bool isStrongBase(std::size_t v,std::size_t w) const;

    double evaluateClustering(const std::vector<std::size_t>& labels) const;

    std::vector<std::unordered_set<std::size_t>> initReachableLdp(const LdpDirectedGraph &graph, LdpParameters<> &parameters, const VertexGroups<std::size_t> *vg=nullptr);




    LdpParameters<>& parameters;
    LPMP::VertexGroups<> vertexGroups;
	std::size_t minV=0;
	std::size_t maxV=0;

    mutable std::vector<std::size_t> sncNeighborStructure;
    mutable std::vector<std::size_t> sncBUNeighborStructure;

    mutable std::vector<char> isBSF;

    mutable std::vector<double> sncTDStructure;
    mutable std::vector<double> sncBUStructure;
    mutable std::vector<char> sncClosedVertices;
    mutable std::vector<double> sncLiftedMessages;
    mutable std::vector<char> sncVerticesInScope;


    const std::size_t& getNumberOfVertices() const{
        return numberOfVertices;
    }

    const LdpDirectedGraph& getMyGraph() const{
        return myGraph;
    }

    const LdpDirectedGraph& getMyGraphLifted() const{
        return myGraphLifted;
    }

    bool checkStrongBase(const std::size_t& v,const std::size_t& w)const;

    bool canJoin(const std::size_t& v,const std::size_t& w)const{
        if(parameters.isMustCutMissing()){
            bool value=canJoinStructure.at(v).isWithinBounds(w)&&canJoinStructure.at(v)[w]>0;
            return value;
        }
        else{
            return true;
        }

    }

    void increaseLBTime(double time) const{
        timeInSncLB+=time;
        callsOfSncLB++;
    }

    void increaseBaseMMTime(double time) const{
        timeInSncBaseMM+=time;
        callsOfSncBaseMM++;
    }

    void increaseLiftedMMTime(double time) const{
        timeInSncLiftedMM+=time;
        callsOfSncLiftedMM++;
    }

    double getAverageLBTime() const{
        return timeInSncLB/callsOfSncLB;
    }

    double getAverageBaseMMTime() const{
        return timeInSncBaseMM/callsOfSncBaseMM;
    }

    double getAverageLiftedMMTime() const{
        return timeInSncLiftedMM/callsOfSncLiftedMM;
    }



private:


    void initAdaptiveThresholds(const LdpDirectedGraph *pBaseGraph, const LdpDirectedGraph *pLiftedGraph);
    void init();

    void sparsifyBaseGraphNew(const LdpDirectedGraph& inputGraph,bool zeroCost=false);


    void sparsifyLiftedGraphNew(const LdpDirectedGraph& inputLiftedGraph);


    void initLiftedStructure();

    void initCanJoinStructure(const LdpDirectedGraph &completeGraph);

	std::size_t s_;
	std::size_t t_;

	std::vector<double> vertexScore;
    std::vector<std::unordered_set<std::size_t>> reachable;


    LdpDirectedGraph myGraph;
    LdpDirectedGraph myGraphLifted;
    const LdpDirectedGraph* pCompleteGraph;

    std::vector<std::unordered_set<std::size_t>> strongBaseEdges;

	std::vector<bool> baseEdgeLabels;

	std::size_t numberOfVertices;
	std::size_t numberOfEdges;
	std::size_t numberOfLiftedEdges;


    double negativeLiftedThreshold;
    double positiveLiftedThreshold;
    double baseThreshold;
    bool keepAllToFirst;
    bool keepAllToLast;

    std::string outputFileName;
    std::string outputFilePrefix;
    std::string outputFileSuffix;



    std::vector<ShiftedVector<char>> liftedStructure;
    std::vector<ShiftedVector<char>> canJoinStructure;

    mutable double timeInSncLB;
    mutable std::size_t callsOfSncLB;

    mutable double timeInSncBaseMM;
    mutable std::size_t callsOfSncBaseMM;

    mutable double timeInSncLiftedMM;
    mutable std::size_t callsOfSncLiftedMM;



};





}
}
#endif /* INCLUDE_DISJOINT_PATHS_DISJOINT_PATHS_INIT_HXX_ */
