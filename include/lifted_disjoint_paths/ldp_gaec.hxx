#ifndef LDP_GAEC_HXX
#define LDP_GAEC_HXX

#include<stdlib.h>
#include<map>
#include<vector>
#include<lifted_disjoint_paths/ldp_directed_graph.hxx>
#include "ldp_instance.hxx"


namespace LPMP {


template <class SNC_FACTOR>
class LdpGreedyAdditive{

public:
    LdpGreedyAdditive(std::vector<std::array<SNC_FACTOR*,2>>* _pSNCFactors, const lifted_disjoint_paths::LdpInstance * _pInstance)
    {
        pSNCFactors=_pSNCFactors;
        pInstance=_pInstance;
        pvg=&pInstance->getVertexGroups();

        numberOfVertices=pSNCFactors->size();
        t=numberOfVertices+1;
        s=numberOfVertices;

        expectedPrimal=0;

        costs=std::vector<std::vector<double>>(numberOfVertices);
        //liftedCosts=std::vector<std::map<std::size_t,double>>(numberOfVertices);
        liftedCostsFull=std::vector<ShiftedVector<double>>(numberOfVertices);
        liftedPotentialForward=std::vector<ShiftedVector<double>>(numberOfVertices);
        liftedPotentialBackward=std::vector<ShiftedVector<double>>(numberOfVertices);

        timeStamps=std::vector<std::vector<std::size_t>>(numberOfVertices);

        startingVertices=std::vector<char>(numberOfVertices,1);
        lastVertices=std::vector<std::size_t>(numberOfVertices);
        firstVertices=std::vector<std::size_t>(numberOfVertices);
      //  endingVertices=std::vector<char>(numberOfVertices,1);

        neighbors=std::vector<std::size_t>(numberOfVertices,t);

        inOutCost=pInstance->parameters.getInputCost()+pInstance->parameters.getOutputCost();


        const std::vector<std::array<SNC_FACTOR*,2>>& single_node_cut_factors_=*pSNCFactors;

        for (std::size_t i = 0; i < numberOfVertices; ++i) { //init costs, lifted costs and time stamps, costs in forward direction only
            lastVertices[i]=i;
            firstVertices[i]=i;

            auto pOutFactor=single_node_cut_factors_[i][1]->get_factor();
            costs[i]=pOutFactor->getBaseCosts();
            timeStamps[i]=std::vector<std::size_t>(costs[i].size());

            std::vector<double> lc=pOutFactor->getLiftedCosts();
            std::vector<std::size_t> liftedIDs=pOutFactor->getLiftedIDs();

            liftedCostsFull[i]=ShiftedVector<double>(getMinForwardNeighbor(i),getMaxForwardNeighbor(i),0);



            assert(lc.size()==liftedIDs.size());
            for (std::size_t j = 0; j < liftedIDs.size(); ++j) {
                std::size_t vertex2=liftedIDs[j];
                //liftedCosts[i][vertex2]=lc[j];
                liftedCostsFull[i][vertex2]=lc[j];
            }

            liftedPotentialForward[i]=ShiftedVector<double>(getMinForwardNeighbor(i),getMaxForwardNeighbor(i),0);
            pOutFactor->LowerBound(false);
            const std::vector<std::size_t>& traverseOrder=pOutFactor->getTraverseOrder();

            for (std::size_t j = 0; j < traverseOrder.size(); ++j) {
                std::size_t vertex2=traverseOrder[j];
                if(vertex2>=numberOfVertices||vertex2==i) continue;
                //liftedCosts[i][vertex2]=lc[j];
                double potential=0;
                std::size_t next=pInstance->sncNeighborStructure[vertex2];
                if(next!=t){
                    potential=pInstance->sncTDStructure[next];
                }
                liftedPotentialForward[i][vertex2]=potential;
            }



        }

        for (std::size_t i = 0; i < numberOfVertices; ++i) { //update costs with opposite direction

            auto pInFactor=single_node_cut_factors_[i][0]->get_factor();
            const std::vector<double>& baseCosts=pInFactor->getBaseCosts();
            const std::vector<std::size_t>& baseIDs=pInFactor->getBaseIDs();

            for (int j = 0; j < baseCosts.size(); ++j) {
                std::size_t vertexID=baseIDs[j];
                if(vertexID>=numberOfVertices) continue;
                auto pOutFactor=single_node_cut_factors_[vertexID][1]->get_factor();

                std::size_t index=pOutFactor->getBaseIDToOrder(i);
                assert(index<costs[vertexID].size());
                costs[vertexID][index]+=baseCosts[j];

            }

            const std::vector<double>& localLiftedCosts=pInFactor->getLiftedCosts();
            const std::vector<std::size_t>& liftedIDs=pInFactor->getLiftedIDs();

            for (int j = 0; j < liftedIDs.size(); ++j) {
                std::size_t vertexID=liftedIDs[j];
                //liftedCosts[vertexID][i]+=localLiftedCosts[j];
                liftedCostsFull[vertexID][i]+=localLiftedCosts[j];
            }



            liftedPotentialBackward[i]=ShiftedVector<double>(getMinBackwardNeighbor(i),getMaxBackwardNeighbor(i),0);
            pInFactor->LowerBound(false);
            const std::vector<std::size_t>& traverseOrder=pInFactor->getTraverseOrder();

            for (std::size_t j = 0; j < traverseOrder.size(); ++j) {
                std::size_t vertex2=traverseOrder[j];
                if(vertex2>=numberOfVertices||vertex2==i) continue;
                //liftedCosts[i][vertex2]=lc[j];
                double potential=0;
                std::size_t next=pInstance->sncNeighborStructure[vertex2];
                if(next!=t){
                    potential=pInstance->sncTDStructure[next];
                }
                liftedPotentialBackward[i][vertex2]=potential;
            }
        }





//        for (int i = 0; i < costs.size(); ++i) {
//            auto pOutFactor=single_node_cut_factors_[i][1]->get_factor();
//            const std::vector<std::size_t>& liftedIDs=pOutFactor->getLiftedIDs();
//            const std::vector<std::size_t>& baseIDs=pOutFactor->getBaseIDs();

//            for(std::size_t j=0;j<baseIDs.size();j++){  //add lifted costs to base costs
//                    costs[i][j]+=liftedCostsFull[i][baseIDs[j]];
//            }
//        }

        liftedCostsOriginal=liftedCostsFull;
        liftedPotentialForwardOrig=liftedPotentialForward;
        liftedPotentialBackwardOrig=liftedPotentialBackward;

        for (std::size_t i = 0; i < costs.size(); ++i) {
            auto pOutFactor=single_node_cut_factors_[i][1]->get_factor();
            const std::vector<std::size_t>& baseIDs=pOutFactor->getBaseIDs();

            std::vector<double> mmCheckForward=pOutFactor->getAllBaseMinMarginalsForMCF();

            assert(baseIDs.size()==costs[i].size());
            for (std::size_t j = 0; j < costs[i].size(); ++j) {
                std::size_t vertex2=baseIDs[j];

                if(vertex2>=numberOfVertices) continue;
                std::size_t timeGap=getTimeGap(i,vertex2);
                double mmCheckF=mmCheckForward[j];
                auto pInFactor=single_node_cut_factors_[vertex2][0]->get_factor();
                std::size_t index=pInFactor->getBaseIDToOrder(i);
                std::vector<double> mmCheckBackward=pInFactor->getAllBaseMinMarginalsForMCF();
                double mmCheckB=mmCheckBackward[index];

                //TODO compare the cost with base min marginals for mcf
                double cost=costs[i][j]+liftedCostsFull[i][vertex2]+liftedPotentialForward[i][vertex2]+liftedPotentialBackward[vertex2][i];
                assert(abs(cost-mmCheckB-mmCheckF)<eps);
                if(cost<inOutCost){
                    Edge e={i,j,cost,timeGap,0};
                    Q.push(e);
                }
            }
        }



    }


    void runGAEC(){

        while(!Q.empty()){
            Edge e=Q.top();
            auto pOutFactor=(*pSNCFactors)[e.firstVertex][1]->get_factor();

            std::size_t secondVertex=pOutFactor->getBaseID(e.neighborIndex);
            std::size_t timeStamp=timeStamps[e.firstVertex][e.neighborIndex];
            if(startingVertices[secondVertex]&&neighbors[e.firstVertex]==t&&timeStamp==e.timeStamp){
               // std::cout<<"connecting "<<e.firstVertex<<","<<secondVertex<<std::endl;
                std::size_t first=firstVertices[e.firstVertex];
                std::size_t last=lastVertices[secondVertex];

                startingVertices[secondVertex]=false;
                //endingVertices[e.firstVertex]=false;
                neighbors[e.firstVertex]=secondVertex;

                //assert(abs(e.cost-costs[e.firstVertex][e.neighborIndex]-liftedCostsFull[e.firstVertex][secondVertex])<eps);
                std::size_t f=first;
                double controlLiftedCost=0;
                double controlPotentialForward=0;
                double controlPotentialBackward=0;
                while(f!=secondVertex){
                    assert(f!=t);
                    if(last<=getMaxForwardNeighbor(f)) controlPotentialForward+=liftedPotentialForwardOrig[f][last];
                    std::size_t sec=secondVertex;
                    while(sec!=t){
                        if(liftedCostsOriginal[f].isWithinBounds(sec)){
                            controlLiftedCost+=liftedCostsOriginal[f][sec];
                        }
                        if(f==first){
                            if(first>=getMinBackwardNeighbor(sec)) controlPotentialBackward+=liftedPotentialBackwardOrig[sec][first];
                        }
                        sec=neighbors[sec];

                    }
                    f=neighbors[f];
                }

                double currentLiftedCost=liftedCostsFull[e.firstVertex][secondVertex];
                if(abs(controlLiftedCost-currentLiftedCost)>=eps){
                    f=first;
                    double controlLiftedCost=0;
                    while(f!=secondVertex){
                        assert(f!=t);
                        std::size_t sec=secondVertex;
                        while(sec!=t){
                            if(liftedCostsOriginal[f].isWithinBounds(sec)){
                                controlLiftedCost+=liftedCostsOriginal[f][sec];
                            }

                            std::cout<<f<<","<<sec<<":"<<liftedCostsOriginal[f][sec]<<std::endl;
                            sec=neighbors[sec];
                            //assert(sec!=t);

                        }
                        f=neighbors[f];
                    }

                }

                //TODO cost with potentials can be checked with all base min marginals for mcf
                assert(abs(controlLiftedCost-currentLiftedCost)<eps);
                assert(abs(controlPotentialForward-liftedPotentialForward[e.firstVertex][secondVertex])<eps);
                assert(abs(controlPotentialBackward-liftedPotentialBackward[secondVertex][e.firstVertex])<eps);
                controlLiftedCost+=costs[e.firstVertex][e.neighborIndex];
                controlLiftedCost+=controlPotentialForward;
                controlLiftedCost+=controlPotentialBackward;
                assert(abs(controlLiftedCost-e.cost)<eps);


                updateCostOutgoing(e.firstVertex,secondVertex);
                updateCostIncomming(e.firstVertex,secondVertex);

                firstVertices[last]=first;
                lastVertices[first]=last;
                expectedPrimal+=e.cost;
            }

            Q.pop();
        }
    }

    void finalizeResults(){

        for (std::size_t graph_node = 0; graph_node < numberOfVertices; ++graph_node) {
            (*pSNCFactors)[graph_node][0]->get_factor()->setNoBaseEdgeActive();
            (*pSNCFactors)[graph_node][1]->get_factor()->setNoBaseEdgeActive();
        }

        for (std::size_t i = 0; i < numberOfVertices; ++i) {
            if(startingVertices[i]){
               // if(neighbors[i]!=t){
                    startingFinal.push_back(i);
                    auto pInFactor=(*pSNCFactors)[i][0]->get_factor();
                    pInFactor->setBaseEdgeActiveWithID(s);
               // }
                //else{//is isolated

                //}

            }
        }

        for (int i = 0; i < numberOfVertices; ++i) {
            std::size_t neighbor=neighbors[i];
            auto pOutFactor=(*pSNCFactors)[i][1]->get_factor();
            pOutFactor->setBaseEdgeActiveWithID(neighbor);
            if(neighbor!=t){
                auto pInFactor=(*pSNCFactors)[neighbor][0]->get_factor();
                pInFactor->setBaseEdgeActiveWithID(i);
            }

        }
        if(diagnostics()) std::cout<<"expected primal "<<expectedPrimal<<std::endl;
    }
    const std::vector<std::size_t>& getStartingVertices()const{
        return startingFinal;
    }

    const std::vector<std::size_t>& getNeighbors()const{
        return neighbors;
    }


private:

    void updateCostOutgoing(std::size_t i, std::size_t j){
        assert(j>i);
       // ShiftedVector<double>& neighborsOfFirst=liftedCostsFull[i];
       // ShiftedVector<double>& neighborsOfSecond=liftedCostsFull[j];
        std::size_t last=lastVertices[j];
        std::size_t first=firstVertices[i];

        for (std::size_t l = liftedCostsFull[last].getMinVertex(); l <= liftedCostsFull[i].getMaxVertex(); ++l) {
            liftedCostsFull[last][l]+=liftedCostsFull[i][l];
            liftedPotentialBackward[l][last]=liftedPotentialBackward[l][i];
            liftedPotentialForward[last][l]+=liftedPotentialForward[i][l];
        }
        std::size_t m=std::max(liftedCostsFull[i].getMaxVertex()+1,liftedCostsFull[last].getMinVertex());
        if(m<numberOfVertices){
            for (std::size_t l = m; l<=liftedCostsFull[last].getMaxVertex() ;++l) {
                liftedPotentialBackward[l][last]=0;
            }
        }


        auto pOutFactor=(*pSNCFactors)[last][1]->get_factor();
        const std::vector<std::size_t>& baseIDs=pOutFactor->getBaseIDs();

        assert(baseIDs.size()==costs[last].size());
        for (std::size_t k = 0; k < costs[last].size(); ++k) {
            std::size_t vertex2=baseIDs[k];
            if(vertex2==t) continue;
            std::size_t timeGap=getTimeGap(last,vertex2);
            double cost=costs[last][k]+liftedCostsFull[last][vertex2]+liftedPotentialForward[last][vertex2]+liftedPotentialBackward[vertex2][last];
            timeStamps[last][k]++;
            if(cost<inOutCost){
                Edge e={last,k,cost,timeGap,timeStamps[last][k]};
                Q.push(e);
            }
        }

    }


    void updateCostIncomming(std::size_t i, std::size_t j){
        assert(j>i);

        std::size_t first=firstVertices[i];
        std::size_t minCommonVertex=getMinBackwardNeighbor(j);
        std::size_t maxCommonVertex=getMaxBackwardNeighbor(first);

//        if(i==1290&&j==1303){
//            std::cout<<"interesting join "<<std::endl;
//        }


        if(maxCommonVertex!=s){
            for (std::size_t l = minCommonVertex; l <= maxCommonVertex; ++l) {
                if(l==1277&&first==1290&&j==1303){
                    std::cout<<"interesting join, orig "<<liftedCostsFull[l][first]<<", to add "<<liftedCostsFull[l][j]<<std::endl;
                }
                liftedCostsFull[l][first]+=liftedCostsFull[l][j];
                 liftedPotentialForward[l][first]=liftedPotentialForward[l][j];

                liftedPotentialBackward[first][l]+=liftedPotentialBackward[j][l];

            }
            std::size_t m=std::min(minCommonVertex,getMaxBackwardNeighbor(first)+1);
            if(m<numberOfVertices){
                for (std::size_t l = getMinBackwardNeighbor(first); l < m; ++l) {
                    liftedPotentialForward[l][first]=0;

                }
            }

        }

        auto pInFactor=(*pSNCFactors)[first][0]->get_factor();
        const std::vector<std::size_t>& baseIDs=pInFactor->getBaseIDs();

        //assert(baseIDs.size()==costs[j].size());
        for (std::size_t k = 0; k < baseIDs.size(); ++k) {
            std::size_t vertex2=baseIDs[k];
            if(vertex2==s) continue;
            auto pOutFactor=(*pSNCFactors)[vertex2][1]->get_factor();
            std::size_t index=pOutFactor->getBaseIDToOrder(first);
            std::size_t timeGap=getTimeGap(vertex2,first);
            double cost=costs[vertex2][index]+liftedCostsFull[vertex2][first]+liftedPotentialForward[vertex2][first]+liftedPotentialBackward[first][vertex2];
            timeStamps[vertex2][index]++;
            if(cost<inOutCost){
                Edge e={vertex2,index,cost,timeGap,timeStamps[vertex2][index]};
                Q.push(e);
            }
        }
    }

    std::size_t getMinForwardNeighbor(std::size_t vertex){
        const VertexGroups<>& vg=*pvg;
        std::size_t time=vg.getGroupIndex(vertex);
        if(time==vg.getMaxTime()){
            return t;
        }
        else{
            return vg.getMinVertexInTime(time+1);
        }
    }

    std::size_t getMaxForwardNeighbor(std::size_t vertex){
         const VertexGroups<>& vg=*pvg;
        std::size_t time=vg.getGroupIndex(vertex);
        std::size_t shiftedTime=time+pInstance->parameters.getMaxTimeGapComplete();
        if(time==vg.getMaxTime()){
            return t;
        }
        else if(shiftedTime>=vg.getMaxTime()){
            return vg.getMaxVertex();
        }
        else{
            return vg.getMaxVertexInTime(shiftedTime);
        }
    }

    std::size_t getMaxBackwardNeighbor(std::size_t vertex){
         const VertexGroups<>& vg=*pvg;
        std::size_t time=vg.getGroupIndex(vertex);
        if(time==1){
            return s;
        }
        else{
            return vg.getMaxVertexInTime(time-1);
        }
    }

    std::size_t getMinBackwardNeighbor(std::size_t vertex){
         const VertexGroups<>& vg=*pvg;
        std::size_t time=vg.getGroupIndex(vertex);
        //std::size_t shiftedTime=time+pInstance->parameters.getMaxTimeGapComplete();
        if(time==1){
            return s;
        }
        else if(time<=pInstance->parameters.getMaxTimeGapComplete()){
            return 0;
        }
        else{
            return vg.getMinVertexInTime(time-pInstance->parameters.getMaxTimeGapComplete());
        }
    }



    std::size_t getTimeGap(std::size_t vertex1,std::size_t vertex2){
        std::size_t l0=pInstance->getGroupIndex(vertex1);
        std::size_t l1=pInstance->getGroupIndex(vertex2);

        assert(l1>l0);
        return l1-l0;
    }

    std::size_t getTimeForIndex(std::size_t vertex1,std::size_t index){
        assert(vertex1<numberOfVertices);
        auto pOutFactor=(*pSNCFactors)[vertex1][1]->get_factor();
        const std::vector<std::size_t>& baseIDs=pOutFactor->getBaseIDs();
        assert(index<baseIDs.size());
        return getTimeGap(vertex1,baseIDs[index]);
    }

    std::size_t getVertexToIndex(std::size_t vertex1,std::size_t index){
        assert(vertex1<numberOfVertices);
        auto pOutFactor=(*pSNCFactors)[vertex1][1]->get_factor();
        const std::vector<std::size_t>& baseIDs=pOutFactor->getBaseIDs();
        assert(index<baseIDs.size());
        return baseIDs.at(index);
    }


    struct Edge
    {
         std::size_t firstVertex;
         std::size_t neighborIndex;
         double cost;
         std::size_t timeGap;
         std::size_t timeStamp;

         bool operator <(Edge const& other) const
         {
             //return -cost/double(timeGap) < -other.cost/double(other.timeGap);
            // return -cost/double(timeGap*timeGap) < -other.cost/double(other.timeGap*other.timeGap);
              return cost > other.cost;
         }

    };

    std::size_t numberOfVertices;
    std::size_t t;
    std::size_t s;
    //std::vector<std::map<std::size_t, double>> costs;   //TODO replace with LdpDirectedGraph
   // std::vector<std::map<std::size_t, double>> liftedCosts; //TODO replace with LdpDirectedGraph

    std::vector<ShiftedVector<double>> liftedCostsFull;
    std::vector<ShiftedVector<double>> liftedPotentialForward;
    std::vector<ShiftedVector<double>> liftedPotentialBackward;

    std::vector<ShiftedVector<double>> liftedPotentialForwardOrig;
    std::vector<ShiftedVector<double>> liftedPotentialBackwardOrig;

    std::vector<ShiftedVector<double>> liftedCostsOriginal;

    std::vector<std::vector<double>> costs;   //TODO replace with LdpDirectedGraph
    //std::vector<std::vector<double>> liftedCosts; //TODO replace with LdpDirectedGraph

    std::vector<std::vector<std::size_t>> timeStamps;

    std::vector<char> startingVertices; //char for representing bool: is starting, to check if can be connected in forward direction
    std::vector<std::size_t> lastVertices;   //just for vertices that are itself firs or last in a path
    std::vector<std::size_t> firstVertices;
    std::vector<std::size_t> startingFinal;
   // std::vector<char> endingVertices;  //to check if it can be connected in a backward direction

    std::vector<std::size_t> neighbors;

    std::vector<std::array<SNC_FACTOR*,2>>* pSNCFactors;

    const lifted_disjoint_paths::LdpInstance* pInstance;

     std::priority_queue<Edge> Q;

     const VertexGroups<>* pvg;

     double inOutCost;

     double expectedPrimal;

     //std::vector<std::array<std::size_t,2>> boundariesOfTimeLayers;  //max vertex and min vertex in each time layer


};






}
#endif // LDP_GAEC_HXX
