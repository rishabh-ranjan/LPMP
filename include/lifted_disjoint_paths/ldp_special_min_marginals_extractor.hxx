#ifndef LDP_SPECIAL_FACTORS_EXTRACTOR_HXX
#define LDP_SPECIAL_FACTORS_EXTRACTOR_HXX

#include"lifted_disjoint_paths/ldp_instance.hxx"

namespace LPMP {

template<class CUT_FACTOR_CONT,class PATH_FACTOR_CONT>
class LdpSpecialMinMarginalsExtractor{
public:
    LdpSpecialMinMarginalsExtractor(std::vector<CUT_FACTOR_CONT*>& cutFactorCont,std::vector<PATH_FACTOR_CONT*>& pathFactorCont,const lifted_disjoint_paths::LdpInstance* pLdpInstance);
    void initMinMarginals(bool doCostUpdate);

    template<class SNC_FACTOR_CONT>
    void sendMessagesToSncFactors(std::vector<std::array<SNC_FACTOR_CONT*,2>>& SNCFactors);

    const std::vector<std::map<std::size_t,double>>& getBaseEdgesWithCosts()const {
        return baseEdgesWithCosts;
    }
    const std::vector<std::map<std::size_t,double>>& getLiftedEdgesWithCosts()const {
        return liftedEdgesWithCosts;
    }


private:
    const lifted_disjoint_paths::LdpInstance& instance;
    std::vector<CUT_FACTOR_CONT*>& cutFactors;
    std::vector<PATH_FACTOR_CONT*>& pathFactors;
    std::size_t numberOfVertices;
    std::vector<std::map<std::size_t,double>> baseEdgesWithCosts;
    std::vector<std::map<std::size_t,double>> liftedEdgesWithCosts;
};


template<class CUT_FACTOR_CONT,class PATH_FACTOR_CONT>
inline LdpSpecialMinMarginalsExtractor<CUT_FACTOR_CONT,PATH_FACTOR_CONT>::LdpSpecialMinMarginalsExtractor(std::vector<CUT_FACTOR_CONT*>& cutFactorCont,std::vector<PATH_FACTOR_CONT*>& pathFactorCont,const lifted_disjoint_paths::LdpInstance* pLdpInstance):
    instance(*pLdpInstance),
    pathFactors(pathFactorCont),
    cutFactors(cutFactorCont),
    numberOfVertices(pLdpInstance->getNumberOfVertices()-2)
{

   // std::cout<<"constructor ok"<<std::endl;
}

template<class CUT_FACTOR_CONT,class PATH_FACTOR_CONT>
inline void LdpSpecialMinMarginalsExtractor<CUT_FACTOR_CONT,PATH_FACTOR_CONT>::initMinMarginals(bool doCostUpdate){
    baseEdgesWithCosts=std::vector<std::map<std::size_t,double>>(numberOfVertices);
    liftedEdgesWithCosts=std::vector<std::map<std::size_t,double>>(numberOfVertices);


    for (int i = 0; i < pathFactors.size(); ++i) {
        auto * pFactor=pathFactors[i]->get_factor();
        std::vector<double> minMarginals=pFactor->getAllMinMarginals();
        const auto& pathVertices=pFactor->getListOfVertices();
        assert(pathVertices.size()==minMarginals.size());
        for (int j = 0; j < pathVertices.size()-1; ++j) {
            std::size_t v1=pathVertices[j];
            std::size_t v2=pathVertices[j+1];
            assert(j<minMarginals.size());
            if(pFactor->getLiftedInfo().at(j)==0){
                baseEdgesWithCosts[v1][v2]+=minMarginals[j];
            }
            else{
                liftedEdgesWithCosts[v1][v2]+=minMarginals[j];
            }
            if(doCostUpdate){
                pFactor->updateEdgeCost(j,-minMarginals[j]);
            }
        }


        liftedEdgesWithCosts[pathVertices.front()][pathVertices.back()]+=minMarginals.back();
        if(doCostUpdate){
            pFactor->updateEdgeCost(pathVertices.size()-1,-minMarginals.back());
        }

#ifndef NDEBUG
        if(doCostUpdate){
            std::vector<double> minMarginalsControl=pFactor->getAllMinMarginals();
            for (int i = 0; i < pathVertices.size()-1; ++i) {
                assert(abs(minMarginalsControl[i])<eps);
            }
        }
#endif
    }

   // std::cout<<"path factors processed"<<std::endl;



    for (int i = 0; i < cutFactors.size(); ++i) {
//        if(i==26){
//            std::cout<<i<<std::endl;
//        }
        auto * cFactor=cutFactors[i]->get_factor();
        auto allMinMarginals=cFactor->getAllMinMarginals();
        auto minMarginals=allMinMarginals.first;
        auto liftedMinMarginal=allMinMarginals.second;

        const auto & inputs=cFactor->getInputVertices();
        const auto & outputs=cFactor->getOutputVertices();

        for (int j = 0; j < inputs.size(); ++j) {
            std::size_t v1=inputs[j];
            const auto * iter=minMarginals.forwardNeighborsBegin(j);
            const auto * end=minMarginals.forwardNeighborsEnd(j);
            std::size_t outputCounter=0;
            for (;iter!=end;iter++) {
                std::size_t v2=outputs.at(iter->head);
                double delta=iter->cost;
                baseEdgesWithCosts[v1][v2]+=delta;
                if(doCostUpdate){
                    cFactor->updateCostBaseForward(j,outputCounter,-delta);
                }
                outputCounter++;
            }
        }


        liftedEdgesWithCosts[cFactor->getLiftedInputVertex()][cFactor->getLiftedOutputVertex()]+=liftedMinMarginal;
        if(doCostUpdate){
            cFactor->updateCostLifted(-liftedMinMarginal);
        }

#ifndef NDEBUG
        // if(debug()){
        auto minMarginalsControl=cFactor->getAllMinMarginals().first;
        auto liftedMinMarginalControl=cFactor->getAllMinMarginals().second;
        assert(abs(liftedMinMarginalControl)<eps);
        for (int j = 0; j < inputs.size(); ++j) {
            std::size_t v1=inputs[j];
            const auto * iter=minMarginalsControl.forwardNeighborsBegin(j);
            const auto * end=minMarginalsControl.forwardNeighborsEnd(j);
            std::size_t outputCounter=0;
            for (;iter!=end;iter++) {
                double delta=iter->cost;
                assert(abs(delta)<eps);
                outputCounter++;
            }
        }
        //}
#endif

    }
    //std::cout<<"cut factors processed"<<std::endl;

}

template<class CUT_FACTOR_CONT,class PATH_FACTOR_CONT>
template<class SNC_FACTOR_CONT>
inline void LdpSpecialMinMarginalsExtractor<CUT_FACTOR_CONT,PATH_FACTOR_CONT>::sendMessagesToSncFactors(std::vector<std::array<SNC_FACTOR_CONT *,2>>& SNCFactors)
{
    for (std::size_t i = 0; i < baseEdgesWithCosts.size(); ++i) {
        assert(i<SNCFactors.size());
        if(i>=SNCFactors.size()) std::cout<<"wrong size"<<std::endl;
        auto* pOutFactor=SNCFactors[i][1]->get_factor();
        auto iter=baseEdgesWithCosts[i].begin();
        const std::vector<std::size_t>& baseIDs= pOutFactor->getBaseIDs();
        std::size_t counter=0;
        while(iter!=baseEdgesWithCosts[i].end()){
            if(baseIDs[counter]<iter->first){
                counter++;
            }
            else if(baseIDs[counter]>iter->first){
                std::cout<<baseIDs[counter]<<","<<iter->first<<std::endl;
                throw std::runtime_error("base min marginal vertex does not exist");
            }
            else{
                assert(baseIDs[counter]==iter->first);
                std::size_t outputVertex=iter->first;
                double costUpdate=iter->second;
                auto * pInFactor=SNCFactors[outputVertex][0]->get_factor();
                pOutFactor->updateEdgeCost(0.5*costUpdate,counter,false);
                std::size_t indexForInput=pInFactor->getBaseIDToOrder(i);
                pInFactor->updateEdgeCost(0.5*costUpdate,indexForInput,false);
                iter++;
                counter++;
            }
        }
    }

    for (std::size_t i = 0; i < liftedEdgesWithCosts.size(); ++i) {
        if(i>=SNCFactors.size()) std::cout<<"wrong size"<<std::endl;
        assert(i<SNCFactors.size());
        auto* pOutFactor=SNCFactors[i][1]->get_factor();
        auto iter=liftedEdgesWithCosts[i].begin();
        const std::vector<std::size_t>& liftedIDs= pOutFactor->getLiftedIDs();
        std::size_t counter=0;

        while(iter!=liftedEdgesWithCosts[i].end()){
            if(liftedIDs[counter]<iter->first){
                counter++;
            }
            else if(liftedIDs[counter]>iter->first){
                std::cout<<liftedIDs[counter]<<","<<iter->first<<std::endl;
                throw std::runtime_error("base min marginal vertex does not exist");
            }
            else{
                assert(liftedIDs[counter]==iter->first);
                std::size_t outputVertex=iter->first;
                double costUpdate=iter->second;
                auto * pInFactor=SNCFactors[outputVertex][0]->get_factor();
                pOutFactor->updateEdgeCost(0.5*costUpdate,counter,true);
                std::size_t indexForInput=pInFactor->getLiftedIDToOrder(i);
                pInFactor->updateEdgeCost(0.5*costUpdate,indexForInput,true);
                counter++;
                iter++;
            }
        }
    }
}

}


#endif // LDP_SPECIAL_FACTORS_EXTRACTOR_HXX
