#ifndef LDP_PATH_FACTOR_HXX
#define LDP_PATH_FACTOR_HXX

#include<cstdlib>
#include<vector>
#include <config.hxx>
#include"lifted_disjoint_paths/ldp_instance.hxx"

namespace LPMP {

class ldp_path_factor{

public:
    ldp_path_factor(const std::vector<std::size_t>& _listOfVertices,const std::vector<double>& _listOfCosts,const std::vector<char>& _isLifted,const lifted_disjoint_paths::LdpInstance* pInstance,bool _mustCut=false);


    double LowerBound() const;

    double EvaluatePrimal() const;

    void setPrimal(const std::vector<std::size_t>& primalDescendants, const std::vector<std::size_t> &vertexLabels);

    const std::vector<char>& getPrimal(){
        return primalSolution;
    }

    template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(); }
    template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(listOfCosts);}

    auto export_variables() { return std::tie(listOfCosts); }

    void init_primal();

    // std::array<double,3> getAllMinMarginals();


    void updateEdgeCost(const std::size_t& edgeIndex,const double& value);

    double getMinMarginal(const std::size_t& edgeIndex,const std::vector<double> *pCosts) const;


    const std::vector<std::size_t>& getListOfVertices()const{
        return listOfVertices;
    }

    const std::vector<double>& getCosts()const{
        return listOfCosts;
    }

    const std::vector<char>& getLiftedInfo()const{
        return isLifted;
    }

    const std::size_t& getNumberOfEdges()const{
        return numberOfEdges;
    }


    const double& getPrimalLiftedCost()const{
        return primalLiftedCost;
    }
    const double& getPrimalBaseCost()const{
        return primalBaseCost;
    }

    bool isMustCut() const{
        return mustCut;
    }

    std::vector<double> getAllMinMarginals()const;

    void print() const;

private:
    double minimize(const std::vector<double> *pCosts, std::size_t edgeIndex, bool edgeLabel)const;
    double minimizeForMustCut(const std::vector<double> *pCosts, std::size_t edgeIndex, bool edgeLabel)const;

    const std::vector<std::size_t> listOfVertices;
    std::vector<double> listOfCosts;
    const std::vector<char> isLifted;
    std::vector<char> isStrongBase;
    std::size_t numberOfEdges;
    std::vector<char> primalSolution;
    mutable double primalBaseCost;
    mutable double primalLiftedCost;
    mutable std::vector<char> optSolution;
    bool mustCut;

};



class ldp_snc_path_message
{
public:
    ldp_snc_path_message(const std::vector<std::size_t>& _edgeIndexInPath,  //Mostly one, two for the first and the last path vertices
                         const std::vector<std::size_t>& _vertexIndexInSnc,
                         const std::vector<char>& _isLifted, bool debugInfo=false)
    {
        edgeIndexInPath=_edgeIndexInPath;  //Mostly one, two for the first and the last path vertices
        vertexIndexInSnc=_vertexIndexInSnc;
        isLifted=_isLifted;
        dimension=_edgeIndexInPath.size();
        debInfo=debugInfo;

        assert(dimension==1||dimension==2);
        assert(vertexIndexInSnc.size()==dimension);
        assert(isLifted.size()==dimension);
        assert(edgeIndexInPath.size()==dimension);

    }

    template<typename PATH_FACTOR>
    void RepamLeft(PATH_FACTOR& l, const double msg, const std::size_t msg_dim) const
    {

        assert(msg_dim<dimension);
        assert(vertexIndexInSnc.size()==dimension);
        assert(isLifted.size()==dimension);
        assert(edgeIndexInPath.size()==dimension);

#ifndef NDEBUG
        //if(debug()){
        std::size_t indexInPath=edgeIndexInPath[msg_dim];
        std::size_t v0=l.getListOfVertices().at(indexInPath);
        std::size_t v1;
        if(indexInPath==l.getNumberOfEdges()-1){
            v1=l.getListOfVertices().at(0);
            lastV1=v1;
            lastV2=v0;
            lastValue=msg;
        }
        else{
            v1=l.getListOfVertices().at(indexInPath+1);
            lastV1=v0;
            lastV2=v1;
            lastValue=msg;
        }
        assert(isLifted.at(msg_dim)==l.getLiftedInfo().at(indexInPath));
        //std::cout<<"Update cost of path, vertices "<<v0<<", "<<v1<<", edge index "<<edgeIndexInPath[msg_dim]<<", value "<<msg<<std::endl;
        // }
#endif

        l.updateEdgeCost(edgeIndexInPath[msg_dim],msg);


    }

    template<typename SINGLE_NODE_CUT_FACTOR>
    void RepamRight(SINGLE_NODE_CUT_FACTOR& r, const double msg, const std::size_t msg_dim) const
    {

        assert(dimension==1||dimension==2);
        assert(msg_dim<dimension);
        assert(vertexIndexInSnc.size()==dimension);
        assert(isLifted.size()==dimension);
        assert(edgeIndexInPath.size()==dimension);

        r.updateEdgeCost(msg,vertexIndexInSnc[msg_dim],isLifted[msg_dim]);
#ifndef NDEBUG
        //if(debug()){
        std::size_t centralNodeID=r.nodeID;

        std::size_t secondVertex;
        std::size_t v0;
        std::size_t v1;

        if(isLifted.at(msg_dim)){

            secondVertex=r.getLiftedIDs().at(vertexIndexInSnc[msg_dim]);

        }
        else{
            secondVertex=r.getBaseIDs().at(vertexIndexInSnc[msg_dim]);
        }
        v0=std::min(centralNodeID,secondVertex);
        v1=std::max(centralNodeID,secondVertex);
        assert(v0==lastV1);
        assert(v1==lastV2);
        if(std::abs(msg+lastValue)>=eps){
            std::cout<<"WRONG update cost "<<v0<<", "<<v1<<", left was "<<lastValue<<", right is "<<msg<<std::endl;
        }
        assert(std::abs(msg+lastValue)<eps);
        // std::cout<<"Update cost of SNC, central node "<<centralNodeID<<", second vertex "<< secondVertex<<", value "<<msg<<std::endl;
        // }
#endif
    }


    template<typename PATH_FACTOR, typename MSG_ARRAY>
    static void SendMessagesToRight(const PATH_FACTOR& l, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const double omega)
    {
        std::vector<double> allMM=l.getAllMinMarginals();

        for(auto it=msg_begin; it!=msg_end; ++it)
        {

            auto& msg = (*it).GetMessageOp();
            for (std::size_t j = 0; j < msg.dimension; ++j) {
                std::size_t index=msg.edgeIndexInPath.at(j);
                double delta=0;
                delta=0.5*allMM.at(index);
                (*it)[j] -= omega*delta;

            }
        }
    }



//TODO: Do not delete this! Make it work later!
    template<typename SINGLE_NODE_CUT_FACTOR, typename MSG_ARRAY>
    static void SendMessagesToLeft(const SINGLE_NODE_CUT_FACTOR& r, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const double omega)
    {
        //assert(omega>=0.29);
         std::vector<double> liftedCosts=r.getLiftedCosts();
         const std::vector<double>& baseCosts=r.getBaseCosts();

         std::vector<double> liftedMM=r.getAllLiftedMinMarginals();
         for (std::size_t i = 0; i < liftedCosts.size(); ++i) {
             liftedMM[i]*=0.5;
             liftedCosts[i]-=liftedMM[i];
         }
         std::vector<double> baseMM=r.getAllBaseMinMarginals(&baseCosts, &liftedCosts);

         std::vector<std::size_t> scalingLifted(liftedCosts.size());
         std::vector<std::size_t> scalingBase(baseCosts.size());


         for(auto it=msg_begin; it!=msg_end; ++it)
         {

             auto& msg = (*it).GetMessageOp();
             if(msg.dimension==0){
                 std::size_t index=msg.vertexIndexInSnc[0];
                 if(msg.isLifted.at(0)){
                     scalingLifted.at(index)++;
                 }
                 else{
                     scalingBase.at(index)++;
                 }
             }
             else{
                 for (std::size_t j = 0; j < msg.dimension; ++j) {
                     std::size_t index=msg.vertexIndexInSnc.at(j);
                     if(msg.isLifted.at(j)){
                         scalingLifted.at(index)++;
                     }
                     else{
                         scalingBase.at(index)++;
                     }
                 }
             }

         }

         for(auto it=msg_begin; it!=msg_end; ++it)
         {

             auto& msg = (*it).GetMessageOp();
             if(msg.dimension==0){
                 std::size_t index=msg.vertexIndexInSnc[0];
                 double delta=0;
                 if(msg.isLifted.at(0)){
                     delta=liftedMM.at(index)/double(scalingLifted.at(index));
                 }
                 else{
                     delta=baseMM.at(index)/double(scalingBase.at(index));
                 }
                 (*it)[0] -= omega*delta;


             }
             else{
                 for (std::size_t j = 0; j < msg.dimension; ++j) {
                     std::size_t index=msg.vertexIndexInSnc.at(j);
                     double delta=0;
                     if(msg.isLifted.at(j)){
                         delta=liftedMM.at(index)/double(scalingLifted.at(index));
                     }
                     else{
                         delta=baseMM.at(index)/double(scalingBase.at(index));
                     }
                     (*it)[j] -= omega*delta;
                 }
             }

         }
    }


    template<typename SINGLE_NODE_CUT_FACTOR, typename MSG>
    void send_message_to_left(const SINGLE_NODE_CUT_FACTOR& r, MSG& msg, const double omega = 1.0)
    {
        assert(dimension==1||dimension==2);
        assert(vertexIndexInSnc.size()==dimension);
        assert(isLifted.size()==dimension);
        assert(edgeIndexInPath.size()==dimension);
        assert(isLifted.size()==dimension);

       // std::cout<<"send to left (snc to path)"<<std::endl;
        double delta=0;
        double controlDelta=0;
        if(dimension==1){
            const std::vector<double>& baseCosts=r.getBaseCosts();
            const std::vector<double>& liftedCosts=r.getLiftedCosts();
            if(isLifted.at(0)){
                delta = r.getOneLiftedMinMarginal(vertexIndexInSnc[0],&baseCosts,&liftedCosts);

#ifndef NDEBUG
                //if(debug()){
                std::vector<double> controlLiftedCosts=liftedCosts;
                controlLiftedCosts[vertexIndexInSnc[0]]-=delta;
                controlDelta=r.getOneLiftedMinMarginal(vertexIndexInSnc[0],&baseCosts,&controlLiftedCosts);
                if(abs(controlDelta)>eps){
                    throw std::runtime_error("wrong control snc lifted min marginal in path message");
                }

#endif
                //}
            }
            else{

                delta = r.getOneBaseEdgeMinMarginal(vertexIndexInSnc[0],&baseCosts,&liftedCosts);
#ifndef NDEBUG
                //if(debug()){
                std::vector<double> controlBaseCosts=baseCosts;
                controlBaseCosts[vertexIndexInSnc[0]]-=delta;
                controlDelta=r.getOneBaseEdgeMinMarginal(vertexIndexInSnc[0],&controlBaseCosts,&liftedCosts);
                if(abs(controlDelta)>eps){
                    throw std::runtime_error("wrong control snc base min marginal in path message");
                }
                //}
#endif
            }

            msg[0] -= omega * delta;
        }
        else{
            std::vector<double> baseCosts=r.getBaseCosts();
            std::vector<double> liftedCosts=r.getLiftedCosts();
            for (std::size_t i=0;i<dimension;i++) {
                if(isLifted.at(i)){
                    // std::cout<<"lifted "<<std::endl;
                    delta = r.getOneLiftedMinMarginal(vertexIndexInSnc[i],&baseCosts,&liftedCosts);
                    liftedCosts[vertexIndexInSnc[i]]-=delta;

#ifndef NDEBUG
                    //if(debug()){
                    controlDelta=r.getOneLiftedMinMarginal(vertexIndexInSnc[i],&baseCosts,&liftedCosts);
                    if(abs(controlDelta)>eps){
                        throw std::runtime_error("control snc min marginal unequall to zero");
                    }
                    //}
#endif

                }
                else{
                    //std::cout<<"base "<<std::endl;
                    delta = r.getOneBaseEdgeMinMarginal(vertexIndexInSnc[i],&baseCosts,&liftedCosts);
                    baseCosts[vertexIndexInSnc[i]]-=delta;
#ifndef NDEBUG
                    // if(debug()){
                    controlDelta=r.getOneBaseEdgeMinMarginal(vertexIndexInSnc[i],&baseCosts,&liftedCosts);
                    if(abs(controlDelta)>eps){
                        throw std::runtime_error("control snc min marginal unequall to zero");
                    }
                    //}
#endif

                }
                msg[i] -= omega * delta;
            }
        }

    }




    template<typename PATH_FACTOR, typename MSG>
    void send_message_to_right(const PATH_FACTOR& l, MSG& msg, const double omega)
    {

        double delta=0;
        double controlDelta=0;
        if(dimension==1){
            const std::vector<double>& costs=l.getCosts();
            delta=l.getMinMarginal(edgeIndexInPath[0],&costs);
#ifndef NDEBUG
            //if(debug()){
            std::vector<double> controlCosts=costs;
            controlCosts.at(edgeIndexInPath[0])-=delta;
            controlDelta=l.getMinMarginal(edgeIndexInPath[0],&controlCosts);
            if(abs(controlDelta)>eps){
                std::cout<<"control value from path factor"<<controlDelta<<", original delta: "<<delta<<std::endl;
                throw std::runtime_error("wrong path min marginal check in path message");
            }
            //}
#endif
            msg[0] -= omega * delta;
        }
        else{
            std::vector<double> localCosts=l.getCosts();

            for (std::size_t i=0;i<dimension;i++) {
                delta=l.getMinMarginal(edgeIndexInPath[i],&localCosts);
                //  std::cout<<"sending "<<delta<<std::endl;
                assert(edgeIndexInPath[i]<localCosts.size());
                localCosts[edgeIndexInPath[i]]-=delta;
#ifndef NDEBUG
                //if(debug()){
                controlDelta=l.getMinMarginal(edgeIndexInPath[i],&localCosts);
                if(abs(controlDelta)>eps){
                    std::cout<<"control value from path factor"<<controlDelta<<", original delta: "<<delta<<std::endl;
                    throw std::runtime_error("wrong path min marginal check in path message");
                }
                //}
#endif
                msg[i] -= omega * delta;
            }
        }

    }


    template<typename SINGLE_NODE_CUT_FACTOR,typename CUT_FACTOR>
    bool check_primal_consistency(const CUT_FACTOR& l, const SINGLE_NODE_CUT_FACTOR& r) const
    {
        return true;

    }

private:
    std::vector<std::size_t> edgeIndexInPath;  //Mostly one, two for the first and the last path vertices
    std::size_t dimension;
    std::vector<std::size_t> vertexIndexInSnc;
    std::vector<char> isLifted;
    bool debInfo;

    mutable std::size_t lastV1;
    mutable std::size_t lastV2;
    mutable double lastValue;

};




}


#endif // LDP_PATH_FACTOR_HXX
