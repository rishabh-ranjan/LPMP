#ifndef LDP_PATH_SEPARATOR_HXX
#define LDP_PATH_SEPARATOR_HXX

#include"lifted_disjoint_paths/ldp_instance.hxx"
#include"ldp_min_marginals_extractor.hxx"
#include"ldp_path_factor.hxx"
#include "ldp_functions.hxx"
//#include "stable_priority_queue.hxx"
#include "ldp_factor_queue.hxx"
#include <set>

namespace LPMP {

template <class PATH_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT>
struct LdpPathMessageInputs{

    void init(PATH_FACTOR* myPathFactor,SINGLE_NODE_CUT_FACTOR_CONT* sncFactor,std::size_t index);

    std::vector<std::size_t> edgeIndicesInPath;  //Mostly one, two for the first and the last path vertices
    std::vector<std::size_t> indicesInSnc;
    std::vector<char> isLiftedForMessage;

};



template <class PATH_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT>
inline void LdpPathMessageInputs<PATH_FACTOR,SINGLE_NODE_CUT_FACTOR_CONT>::init(PATH_FACTOR* myPathFactor,SINGLE_NODE_CUT_FACTOR_CONT* sncFactor,std::size_t index) {

    bool isOut=sncFactor->get_factor()->isNodeOutFlow();

    std::size_t numberOfEdges=myPathFactor->getNumberOfEdges();
    assert(index<numberOfEdges);

    const std::vector<std::size_t>& pathVertices=myPathFactor->getListOfVertices();
    assert(pathVertices.at(index)==sncFactor->get_factor()->nodeID);
    const std::vector<char>& liftedInfo=myPathFactor->getLiftedInfo();


    if(index==0){
        assert(isOut);
        auto * pFactorOutFactor=sncFactor->get_factor();
        indicesInSnc={0,0};
        isLiftedForMessage={liftedInfo.at(0),1};
        edgeIndicesInPath={0,numberOfEdges-1};


        //std::size_t secondPathVertexIndex;
        if(liftedInfo.at(0)){
            indicesInSnc[0]=pFactorOutFactor->getLiftedIDToOrder(pathVertices.at(1));
        }
        else{
            indicesInSnc[0]=pFactorOutFactor->getBaseIDToOrder(pathVertices.at(1));
        }
        indicesInSnc[1]=pFactorOutFactor->getLiftedIDToOrder(pathVertices.back());

    }
    else if(index==numberOfEdges-1){
        assert(!isOut);
        auto * pFactorInFactor=sncFactor->get_factor();
        indicesInSnc={0,0};
        isLiftedForMessage={liftedInfo.at(numberOfEdges-2),1};
        edgeIndicesInPath={numberOfEdges-2,numberOfEdges-1};

        if(liftedInfo.at(numberOfEdges-2)){
            indicesInSnc[0]=pFactorInFactor->getLiftedIDToOrder(pathVertices.at(numberOfEdges-2));
        }
        else{
            indicesInSnc[0]=pFactorInFactor->getBaseIDToOrder(pathVertices.at(numberOfEdges-2));
        }
        indicesInSnc[1]=pFactorInFactor->getLiftedIDToOrder(pathVertices[0]);

    }
    else if(isOut){
        auto * pFactorOutFactor=sncFactor->get_factor();
        indicesInSnc={0};
        isLiftedForMessage={liftedInfo.at(index)};
        edgeIndicesInPath={index};
        if(liftedInfo.at(index)){
            indicesInSnc[0]=pFactorOutFactor->getLiftedIDToOrder(pathVertices.at(index+1));
        }
        else{
            indicesInSnc[0]=pFactorOutFactor->getBaseIDToOrder(pathVertices.at(index+1));
        }


    }
    else {
        auto * pFactorInFactor=sncFactor->get_factor();
        indicesInSnc={0};
        isLiftedForMessage={liftedInfo.at(index-1)};
        edgeIndicesInPath={index-1};
        if(liftedInfo.at(index-1)){
            indicesInSnc[0]=pFactorInFactor->getLiftedIDToOrder(pathVertices.at(index-1));
        }
        else{
            indicesInSnc[0]=pFactorInFactor->getBaseIDToOrder(pathVertices.at(index-1));
        }
    }


}


template <class PATH_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT> //PATH_FACTOR is ldp_path_factor, SINGLE_NODE_CUT_FACTOR_CONT is the container wrapper
class ldp_path_separator {

public:
    ldp_path_separator(const lifted_disjoint_paths::LdpInstance * _pInstance, ldp_min_marginals_extractor<SINGLE_NODE_CUT_FACTOR_CONT>& _mmExtractor);


    void separatePathInequalities(std::size_t maxConstraints,double minImprovement);

    bool checkWithBlockedEdges(const PATH_FACTOR& pFactor,const std::vector<std::set<std::size_t>>& blockedBaseEdges,const std::vector<std::set<std::size_t>>& blockedLiftedEdges)const;
    void updateUsedEdges(const PATH_FACTOR& pFactor,std::vector<std::set<std::size_t>>& blockedBaseEdges,std::vector<std::map<std::size_t,std::size_t>>& usedBaseEdges,std::vector<std::set<std::size_t>>& blockedLiftedEdges,std::vector<std::map<std::size_t,std::size_t>>& usedLiftedEdges,const std::size_t& maxUsage)const;

    LdpFactorQueue<PATH_FACTOR>& getFactorQueue(){
        return factorQueue;
    }


private:

    PATH_FACTOR* createPathFactor(const std::size_t& lv1,const std::size_t& lv2,const std::size_t& bv1,const std::size_t& bv2,bool isLifted,bool isMustCut);  //lifted edge vertices and the connecting base edge vertices
    std::list<std::pair<std::size_t,bool>> findShortestPath(const std::size_t& firstVertex,const std::size_t& lastVertex);
    // void insertToQueue(const double& improvementValue,PATH_FACTOR* pPathFactor );


    const lifted_disjoint_paths::LdpInstance * pInstance;
    ldp_min_marginals_extractor<SINGLE_NODE_CUT_FACTOR_CONT>& mmExtractor;
    std::size_t numberOfVertices;
    std::size_t addedFactors;

    std::vector<std::map<std::size_t,double>> baseMM;
    std::vector<std::map<std::size_t,double>> liftedMM;
    std::vector<std::list<std::size_t>> predecessors;  //Can I use list? Maybe yes, just predecessors will not be sorted!
    std::vector<std::list<std::size_t>> descendants;
    std::vector<std::vector<std::pair<std::size_t,bool>>> usedEdges;

    std::vector<char> isInQueue;
    std::vector<std::size_t> predInQueue;
    std::vector<char> predInQueueIsLifted;

    std::size_t maxTimeGap;

    LdpFactorQueue<PATH_FACTOR> factorQueue;


};




template <class PATH_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT>
inline void ldp_path_separator<PATH_FACTOR,SINGLE_NODE_CUT_FACTOR_CONT>::updateUsedEdges(const PATH_FACTOR& pFactor,std::vector<std::set<std::size_t>>& blockedBaseEdges,std::vector<std::map<std::size_t,std::size_t>>& usedBaseEdges,std::vector<std::set<std::size_t>>& blockedLiftedEdges,std::vector<std::map<std::size_t,std::size_t>>& usedLiftedEdges,const std::size_t& maxUsage)const{
    const std::vector<std::size_t>& vertices= pFactor.getListOfVertices();
    const std::vector<char>& liftedInfo= pFactor.getLiftedInfo();

    //assert(liftedInfo.size()==vertices.size());
    for (int i = 0; i < liftedInfo.size(); ++i) {

        std::size_t vertex1=vertices[i];
        std::size_t vertex2;
        if(i<vertices.size()-1){
            vertex2=vertices[i+1];
        }
        else{
            vertex2=vertices.back();
            vertex1=vertices.front();
        }

        if(liftedInfo[i]){
            assert(vertex1<blockedLiftedEdges.size());
            assert(vertex1<usedLiftedEdges.size());
            std::size_t& currentUsage=usedLiftedEdges[vertex1][vertex2];
            assert(currentUsage<maxUsage&&blockedLiftedEdges[vertex1].count(vertex2)==0);
            currentUsage++;
            if(currentUsage==maxUsage){
                blockedLiftedEdges[vertex1].insert(vertex2);
            }
        }
        else{
            assert(vertex1<blockedBaseEdges.size());
            assert(vertex1<usedBaseEdges.size());
            std::size_t& currentUsage=usedBaseEdges[vertex1][vertex2];
            assert(currentUsage<maxUsage&&blockedBaseEdges[vertex1].count(vertex2)==0);
            currentUsage++;
            if(currentUsage==maxUsage){
                blockedBaseEdges[vertex1].insert(vertex2);
            }
        }

    }
}



template <class PATH_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT>
inline bool ldp_path_separator<PATH_FACTOR,SINGLE_NODE_CUT_FACTOR_CONT>::checkWithBlockedEdges(const PATH_FACTOR& pFactor,const std::vector<std::set<std::size_t>>& blockedBaseEdges,const std::vector<std::set<std::size_t>>& blockedLiftedEdges)const{
    const std::vector<std::size_t>& vertices= pFactor.getListOfVertices();
    const std::vector<char>& liftedInfo= pFactor.getLiftedInfo();

    bool isFree=true;


    if(!pFactor.isMustCut()){
        assert(liftedInfo.size()==vertices.size());
    }
    else{
        assert(liftedInfo.size()+1==vertices.size());
    }

    for (int i = 0; i < liftedInfo.size(); ++i) {

        std::size_t vertex1=vertices[i];
        std::size_t vertex2;
        if(i<vertices.size()-1){
            vertex2=vertices[i+1];
        }
        else{
            vertex2=vertices.back();
            vertex1=vertices.front();
        }

        if(liftedInfo[i]){
            assert(vertex1<blockedLiftedEdges.size());
            auto f=blockedLiftedEdges[vertex1].find(vertex2);
            if(f!=blockedLiftedEdges[vertex1].end()){
                isFree=false;
            }
        }
        else{
            assert(vertex1<blockedBaseEdges.size());
            auto f=blockedBaseEdges[vertex1].find(vertex2);
            if(f!=blockedBaseEdges[vertex1].end()){
                isFree=false;
            }
        }
    }

    return isFree;
}


template <class PATH_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT>
inline ldp_path_separator<PATH_FACTOR,SINGLE_NODE_CUT_FACTOR_CONT>::ldp_path_separator(const lifted_disjoint_paths::LdpInstance * _pInstance, ldp_min_marginals_extractor<SINGLE_NODE_CUT_FACTOR_CONT>& _mmExtractor):

    pInstance(_pInstance),
    mmExtractor(_mmExtractor)

{
    numberOfVertices=pInstance->getNumberOfVertices()-2;
    isInQueue=std::vector<char>(numberOfVertices);
    predInQueue=std::vector<std::size_t>(numberOfVertices,std::numeric_limits<std::size_t>::max());
    predInQueueIsLifted=std::vector<char>(numberOfVertices,2);
    maxTimeGap=std::max(pInstance->getGapLifted(),pInstance->getGapBase());
    addedFactors=0;

}





template <class PATH_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT>
inline std::list<std::pair<std::size_t,bool>> ldp_path_separator<PATH_FACTOR,SINGLE_NODE_CUT_FACTOR_CONT>::findShortestPath(const std::size_t& firstVertex,const std::size_t& lastVertex){
    assert(firstVertex!=lastVertex);

    std::list<std::pair<std::size_t,bool>> shortestPath;
    //Just BSF
    std::list<std::size_t> queue; //vertex and is lifted (between vertex and its predecessor in queue)

    // std::cout<<"finding shortest path between "<<firstVertex<<", "<<lastVertex<<std::endl;
    const auto& vg=pInstance->getVertexGroups();
    std::size_t firstIndex=vg.getGroupIndex(firstVertex);
    std::size_t lastIndex=vg.getGroupIndex(lastVertex);
    //  std::cout<<"time of last vertex "<<lastIndex<<std::endl;
    queue.push_back(firstVertex) ;  //is lifted for the first does not matter
    bool pathFound=false;
    std::vector<std::size_t> toDeleteFromQueue;
    toDeleteFromQueue.push_back(firstIndex);


    while(!queue.empty()&&!pathFound){
        std::size_t& vertex=queue.front();

        assert(vertex<usedEdges.size());
        std::vector<std::pair<std::size_t,bool>>& neighbors=usedEdges[vertex];

        // std::cout<<"vertex in queue "<<vertex<<", neighbors "<<std::endl;
        for(std::size_t i=0;i<neighbors.size();i++){
            std::pair<std::size_t,bool> p=neighbors[i];
            std::size_t neighborOfVertex=p.first;
            assert(neighborOfVertex<numberOfVertices);
            bool isLifted=p.second;
            // std::cout<<neighborOfVertex<<", "<<std::endl;
            if(neighborOfVertex==lastVertex){
                //std::cout<<"is last vertex"<<std::endl;
                assert(lastVertex<numberOfVertices);
                predInQueue[lastVertex]=vertex;
                predInQueueIsLifted[lastVertex]=isLifted;
                isInQueue[lastVertex]=1; //To be cleared
                toDeleteFromQueue.push_back(lastVertex);

                pathFound=true;
                break;
            }

            std::size_t nodeTimeIndex=vg.getGroupIndex(neighborOfVertex);
            if(nodeTimeIndex<lastIndex&&!isInQueue[neighborOfVertex]){
                //  std::cout<<"gets to queue "<<std::endl;
                queue.push_back(neighborOfVertex);
                isInQueue[neighborOfVertex]=1;
                predInQueue[neighborOfVertex]=vertex;
                predInQueueIsLifted[neighborOfVertex]=isLifted;
                toDeleteFromQueue.push_back(neighborOfVertex);
            }
        }
        queue.pop_front();
    }
    assert(pathFound);

    std::size_t currentVertex=lastVertex;
    while(currentVertex!=firstVertex){  //path contains vertex and info about edge starting in it
        assert(currentVertex<numberOfVertices);
        std::size_t newVertex=predInQueue[currentVertex];
        assert(predInQueueIsLifted.at(currentVertex)<2);
        bool isEdgeLifted=predInQueueIsLifted[currentVertex];
        shortestPath.push_front({newVertex,isEdgeLifted});
        currentVertex=newVertex;
    }

    std::size_t maxValue=std::numeric_limits<std::size_t>::max();
    for(auto& v:toDeleteFromQueue){
        assert(v<numberOfVertices);
        isInQueue[v]=0;
        predInQueue[v]=maxValue;
        predInQueueIsLifted[v]=2;
    }
    return shortestPath;
}



template  <class PATH_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT>
inline PATH_FACTOR* ldp_path_separator<PATH_FACTOR,SINGLE_NODE_CUT_FACTOR_CONT>::createPathFactor(const std::size_t& lv1, const std::size_t& lv2, const std::size_t &bv1, const std::size_t &bv2,bool isLifted,bool isMustCut){


//    if(isMustCut&&!pInstance->existLiftedEdge(lv1,lv2)){
//        std::cout<<"required lifted edge between "<<lv1<<","<<lv2<<", not found"<<std::endl;
//    }
    assert(isMustCut||pInstance->existLiftedEdge(lv1,lv2));

    std::list<std::pair<std::size_t,bool>> beginning;
    if(lv1!=bv1) beginning=findShortestPath(lv1,bv1);
    std::list<std::pair<std::size_t,bool>> ending;
    if(bv2!=lv2) ending=findShortestPath(bv2,lv2);
    if(lv1==bv1&&lv2==bv2) std::cout<<"base covering lifted for path "<<std::endl;
    std::size_t edgeLength=0;
    std::size_t verticesLength=0;
    if(!isMustCut){
        edgeLength=beginning.size()+ending.size()+2;
        verticesLength=edgeLength;
    }
    else{
        edgeLength=beginning.size()+ending.size()+1;
        verticesLength=edgeLength+1;
    }
    std::vector<std::size_t> pathVertices(verticesLength);
    std::vector<char> liftedEdgesIndices(edgeLength);


    //std::size_t numberOfVerticesInPath=pathVertices.size();

    auto iter=beginning.begin();
    auto end=beginning.end();
    std::size_t counter=0;
    // std::cout<<"path vertices"<<std::endl;
    for (;iter!=end;iter++) {
        std::size_t vertex=iter->first;
        bool isLiftedEdge=iter->second;

        assert(counter<pathVertices.size());
        pathVertices[counter]=vertex;
        liftedEdgesIndices[counter]=isLiftedEdge;
        counter++;
    }
    pathVertices[counter]=bv1;
    liftedEdgesIndices[counter]=isLifted;

    counter++;
    //Careful with the bridge edge!
    if(bv2!=lv2) assert(ending.front().first==bv2);
    iter=ending.begin();
    end=ending.end();
    for(;iter!=end;iter++){
        std::size_t vertex=iter->first;
        bool isLiftedEdge=iter->second;

        assert(counter<pathVertices.size());
        pathVertices[counter]=vertex;
        liftedEdgesIndices[counter]=isLiftedEdge;
        counter++;
    }
    assert(counter<pathVertices.size());
    pathVertices[counter]=lv2;
    if(!isMustCut){
        liftedEdgesIndices[counter]=true; //this relates to the big lifted edge connecting the first and the last path vertex
    }
    assert(counter==verticesLength-1);
    assert(pathVertices.front()<pInstance->getNumberOfVertices()-2&&pathVertices.back()<pInstance->getNumberOfVertices());


    std::vector<double> costs(edgeLength,0); //last member is the lifted edge cost
    assert(pathVertices.front()==lv1);
    assert(pathVertices.back()==lv2);

    PATH_FACTOR* pPathFactor=new PATH_FACTOR(pathVertices,costs,liftedEdgesIndices,pInstance,isMustCut);
    return pPathFactor;



}

template  <class PATH_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT>
inline void ldp_path_separator<PATH_FACTOR,SINGLE_NODE_CUT_FACTOR_CONT>::separatePathInequalities(std::size_t maxConstraints, double minImprovement){

    baseMM=mmExtractor.getBaseEdgesMinMarginals();
    liftedMM=mmExtractor.getLiftedEdgesMinMarginals();
    usedEdges= std::vector<std::vector<std::pair<std::size_t,bool>>> (numberOfVertices);

    bool removeUsedPositive=false;
    bool onlyBasePaths=false;



    assert(baseMM.size()==numberOfVertices+2);
    assert(liftedMM.size()==numberOfVertices);

    predecessors=std::vector<std::list<std::size_t>> (numberOfVertices);  //Can I use list? Maybe yes, just predecessors will not be sorted!
   //  predecessors=std::vector<std::set<std::size_t>> (numberOfVertices);  //Can I use list? Maybe yes, just predecessors will not be sorted!
    descendants=std::vector<std::list<std::size_t>> (numberOfVertices);

    addedFactors=0;

    bool useOuterPath=true;
    std::size_t numberOfOuterPaths=0;
    std::size_t numberOfStandardPaths=0;

    assert(factorQueue.isQueueEmpty());


    std::size_t constraintsCounter=0;

    std::vector<std::tuple<double,std::size_t,std::size_t,bool>> edgesToSort; //contains negative base and lifted edges: cost,vertex1,vertex2,isLifted

    for(std::size_t i=0;i<numberOfVertices;i++){
        auto iter=baseMM[i].begin();
        auto end=baseMM[i].end();
        for(;iter!=end;iter++){
            if(iter->first<numberOfVertices&&iter->second<-minImprovement){
                edgesToSort.push_back(std::tuple(iter->second,i,iter->first,false));
            }
        }
    }

    std::vector<std::map<std::size_t,double>> positiveLifted(numberOfVertices);
    for(std::size_t i=0;i<numberOfVertices;i++){
        auto iter=liftedMM[i].begin();
        auto end=liftedMM[i].end();
        for(;iter!=end;iter++){
            //if(iter->second<-eps){
            if(iter->second<-minImprovement){
               if(!onlyBasePaths) edgesToSort.push_back(std::tuple(iter->second,i,iter->first,true));
            }
            //else if(iter->second>eps){
            else if(iter->second>minImprovement){
                positiveLifted[i][iter->first]=iter->second;
            }
        }
    }

    std::stable_sort(edgesToSort.begin(),edgesToSort.end(),lifted_disjoint_paths::edgeCompare<double>);


    for(std::size_t i=0;i<numberOfVertices;i++){
        predecessors[i].push_back(i);
      //    predecessors[i].insert(i);
        descendants[i].push_back(i);

    }



    for(std::size_t i=0;i<edgesToSort.size();i++){
        std::tuple<double,std::size_t,std::size_t,bool>& edge=edgesToSort[i];
        //   std::tuple<float,std::size_t,std::size_t,bool>& edge=edgesToSort[i];
        std::size_t& vertex1=std::get<1>(edge);
        std::size_t& vertex2=std::get<2>(edge);
        bool isLifted=std::get<3>(edge);
        double edgeCost=std::get<0>(edge);


        assert(vertex1<numberOfVertices&&vertex2<numberOfVertices);
        auto iterPredV1=predecessors[vertex1].begin();   //first vertex to process is always vertex1 itself
        auto endPredV1=predecessors[vertex1].end();


        bool alreadyConnected=false;
        if(debug()){
            for(auto iter=descendants[vertex1].begin();iter!=descendants[vertex1].end();iter++){
                if(*iter==vertex2){
                    alreadyConnected=true;
                    break;
                }
            }
        }


        for (;iterPredV1!=endPredV1;iterPredV1++) {
            const std::size_t& pred=*iterPredV1;
            std::size_t l0=pInstance->getGroupIndex(pred);

            assert(pred<numberOfVertices);

            auto iterDescPred=descendants[pred].begin();  //Put descendants of V2 into descendants of pred
            auto endDescPred=descendants[pred].end();

            auto iterDescV2=descendants[vertex2].begin();
            auto endDescV2=descendants[vertex2].end();

            auto iterLifted=positiveLifted[pred].begin();

            if(pred==vertex1&&*iterDescPred==vertex2){
                alreadyConnected=true;
                if(!debug())break;
            }
            if(debug()||!alreadyConnected){
                while(iterDescV2!=endDescV2){
                    const std::size_t& descV2=*iterDescV2;
                    while(iterLifted!=positiveLifted[pred].end()&&iterLifted->first<descV2){
                        iterLifted++;
                    }
                    while(iterDescPred!=endDescPred&&*iterDescPred<descV2){
                        iterDescPred++;
                        if(pred==vertex1&&*iterDescPred==vertex2){
                            alreadyConnected=true;
                            if(!debug())break;
                        }
                    }
                    if(alreadyConnected&&!debug()) break;
                    if(iterDescPred==endDescPred||descV2<*iterDescPred){ //exists desc of v2 not contained in desc of pred

                        std::size_t l1=pInstance->getGroupIndex(descV2);

                        if(debug()) assert(l1-l0>maxTimeGap||!alreadyConnected);

                        assert(descV2<numberOfVertices);
                        if(l1-l0<=maxTimeGap){                           // std::cout<<"exists new descendant "<<std::endl;

                            if(iterLifted->first==descV2){  //Contradicting lifted edge exists!
                                if(pred!=vertex1||descV2!=vertex2){

                                    PATH_FACTOR* pPathFactor= createPathFactor(pred,descV2,vertex1,vertex2,isLifted,false); //TODO first just put to a queue (list of vertices and information if the edges are lifted) and then select the best
                                    double improvementValue=std::min(abs(edgeCost),iterLifted->second);

                                    factorQueue.insertToQueue(improvementValue,pPathFactor);
                                    constraintsCounter++;
                                    numberOfStandardPaths++;
                                }
                                if(removeUsedPositive){
                                    auto iterToDel=iterLifted;
                                    iterLifted++;
                                    positiveLifted[pred].erase(iterToDel);
                                }

                            }

                            descendants[pred].insert(iterDescPred,descV2);
                            predecessors[descV2].push_back(pred);
                            //  predecessors[descV2].insert(pred);

                        }

                    }
                    iterDescV2++;
                }

            }

        }


        if(!onlyBasePaths&&useOuterPath&&isLifted){
           //TODO
            //for desc: descendants of vertex1 lower than vertex2:
            //  for pl: positvelifted[pl] lower than vertex2
            //     if(descendants[pl.end] contain vertex2) add factor
               std::size_t l1=pInstance->getGroupIndex(vertex2);
               std::size_t l0=pInstance->getGroupIndex(vertex1);
               if(l1-l0>=4){
//                   std::sort(predecessors[vertex2].begin(),predecessors[vertex2].end()); //or use set for predecessors
                   for (auto idV1=descendants[vertex1].begin();idV1!=descendants[vertex1].end();idV1++) {
                       std::size_t descV1=*idV1;
                       std::size_t lMiddle=pInstance->getGroupIndex(descV1);
                       if(lMiddle>=l1) break;
                       auto plIt=positiveLifted[descV1].begin();
                       auto predIt=predecessors[vertex2].begin();
                       auto endPred=predecessors[vertex2].end();
                       auto descDesc=descendants[descV1].begin();
                       while(plIt!=positiveLifted[descV1].end()&&predIt!=endPred){
                           while(descDesc!=descendants[descV1].end()&&*descDesc<plIt->first){
                               descDesc++;
                           }
                           if(*predIt<plIt->first){
                               predIt++;
                           }
                           else if(*predIt>plIt->first){
                               plIt++;
                               if(plIt->second>vertex2) break;
                           }
                           else{
                               //if(descDesc==descendants[descV1].end()||*descDesc>plIt->first){ //We do not wand add positive lifted edges that have already been processed via an inner path
                                  // std::cout<<"outer path"<<std::endl;
                                   PATH_FACTOR* pPathFactor= createPathFactor(vertex1,vertex2,descV1,plIt->first,isLifted,false); //TODO first just put to a queue (list of vertices and information if the edges are lifted) and then select the best
                                   double improvementValue=std::min(abs(edgeCost),plIt->second);

                                   factorQueue.insertToQueue(improvementValue,pPathFactor);
                                   constraintsCounter++;
                                   numberOfOuterPaths++;
                               //}
                               predIt++;
                               plIt++;

                           }
                       }

                   }

               }
        }
        //}
        usedEdges[vertex1].push_back(std::pair(vertex2,isLifted));
        if(debug()){
            if(i>1){
                std::tuple<double,std::size_t,std::size_t,bool>& e=edgesToSort[i-1];

                std::size_t& v1=std::get<1>(e);
                std::size_t& v2=std::get<2>(e);
                bool il=std::get<3>(e);
                double ec=std::get<0>(e);
                if(abs(ec-edgeCost)<eps){
                    std::cout<<"same edge cost "<<v1<<" "<<v2<<", cost: "<<ec<<". is lifted "<<il<<std::endl;
                    std::cout<<"same edge cost "<<vertex1<<" "<<vertex2<<", cost: "<<edgeCost<<". is lifted "<<isLifted<<std::endl;
                    std::cout<<"equal values "<<(ec==edgeCost)<<std::endl;
                }


            }
        }

    }

    if(diagnostics()) std::cout<<"standart paths: "<<numberOfStandardPaths<<", outer paths "<<numberOfOuterPaths<<std::endl;
}



}



#endif // LDP_PATH_SEPARATOR_HXX
