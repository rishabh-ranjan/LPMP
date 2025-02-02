#ifndef LDP_VERTEX_GROUPS_HXX
#define LDP_VERTEX_GROUPS_HXX

#include <stdexcept>
//#include "disjoint-paths/disjointPathsMethods.hxx"
#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <andres/graph/digraph.hxx>
#include <iterator>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <utility>
#include <unordered_map>
#include "ldp_file_processing_methods.hxx"


namespace LPMP {


template<class T=std::size_t>
class VertexGroups{
public:





    VertexGroups(std::unordered_map<std::size_t,std::vector<std::size_t>> groups_,std::vector<std::size_t> vToGroup_): //no shift
        vToGroup(vToGroup_)
    {
        timeShift=0;
        vertexShift=0;
        vertexShiftBack=0;
        maxVertex=vToGroup_.size()-3;
        maxTime=vToGroup_.back()-1;
        std::cout<<"vg max vertex "<<maxVertex<<std::endl;
        std::cout<<"vg max time "<<maxTime<<std::endl;
//
        std::size_t groupsSize=maxTime+2-timeShift;
        groups=std::vector<std::vector<std::size_t>>(groupsSize);
        for(auto pair:groups_){
            assert(pair.first<groupsSize);
            groups[pair.first]=pair.second;
        }
    }

    VertexGroups(){
        maxVertex=0;
        maxTime=0;
        timeShift=0;
        vertexShiftBack=0;
        vertexShift=0;
    }




    template<class PAR> VertexGroups(PAR& parameters){
        timeShift=0;
        vertexShift=0;
        vertexShiftBack=0;
        initFromFile(parameters.getTimeFileName(),parameters);

    }

    template<class PAR>
    void initFromFile(const std::string& fileName,const PAR& parameters);

    void initFromVector(const std::vector<std::size_t>& verticesInFrames,std::size_t timeShift_,std::size_t vertexShift_);

    void initFromVector(const std::vector<std::size_t>& verticesInFrames){
        initFromVector(verticesInFrames,0,0);
    }

    void print()const;


  //  std::vector<std::vector<std::size_t>> extractCorePaths(const std::vector<std::vector<std::size_t>>& paths, std::size_t cutoff,bool isFirst,bool isLast,std::size_t vertexShift) const;


    std::size_t getGroupIndex(std::size_t v) const{
        if(v>maxVertex+2){
            std::cout<<"v is "<<v<<", max vertex "<<maxVertex<<std::endl;
        }
        assert(v>=vertexShift&&v<=maxVertex+2);
        assert(v-vertexShift<vToGroup.size());
        if(v<vertexShift){
            std::cout<<"v "<<v<<", time shift "<<vertexShift<<std::endl;
        }
        if(v>maxVertex+2){
            std::cout<<"v "<<v<<", max time "<<maxVertex<<std::endl;
        }        
        return vToGroup[v-vertexShift];
    }

    const std::vector<std::size_t>& getGroupVertices(std::size_t index)const{  //e.g. first layer 201 -> timeShift =200, normally first layer = 1

        assert(index>=timeShift&&index<=maxTime+1);
        return groups.at(index-timeShift);
    }

    //ID of maximal valid vertex (i.e. without s and t)
    std::size_t getMaxVertex() const {
        return maxVertex;
    }

    //Time of the last video frame
    std::size_t getMaxTime() const {
        return maxTime;
    }

    std::size_t getMinTime() const {
        return timeShift+1;
    }

    std::size_t getMinVertex() const {
        return vertexShift;
    }

    std::size_t getVertexShiftBack() const{
        return vertexShiftBack;
    }

    std::size_t getMinVertexInTime(std::size_t time)const{
        assert(time<=maxTime);
        assert(time>timeShift);
        std::size_t i=time;
        while(i<=maxTime&&groups.at(i-timeShift).size()==0){
            i++;
        }
        if(i>maxTime) return std::numeric_limits<std::size_t>::max();
        else return groups.at(i-timeShift).front();
    }

    std::size_t getMaxVertexInTime(std::size_t time)const {
        assert(time<=maxTime);
        assert(time>timeShift);
        std::size_t i=time;
        while(i>timeShift&&groups.at(i-timeShift).size()==0){
            i--;
        }
        if(i==timeShift) return std::numeric_limits<std::size_t>::max();
        else return groups.at(i-timeShift).back();
    }


    std::vector<std::vector<std::size_t>> extractInnerPaths(const std::vector<std::vector<std::size_t> > &paths, const std::size_t minT, const std::size_t maxT, const std::size_t vShift=0) const;

    void testCorrectness(bool printToo)const ;

private:
    //std::vector<std::vector<std::size_t>> groups;
    std::vector<std::vector<std::size_t>> groups;
    std::vector<std::size_t> vToGroup;
    std::size_t maxVertex;
    std::size_t maxTime;
    std::size_t timeShift;  //If first layer is not one, but timeShift+1
    std::size_t vertexShift;  //First vertex is not zero but vertexShift

    std::size_t vertexShiftBack;


};

template<class T>
inline void VertexGroups<T>::print() const{
    for (std::size_t i = 0; i <= maxTime+1; ++i) {
        std::cout<<"layer "<<i<<": ";
        for (int j = 0; j < groups[i].size(); ++j) {
            std::cout<<groups[i][j]<<",";

        }
        std::cout<<std::endl;
    }

}

template<class T>
inline void VertexGroups<T>::testCorrectness(bool printToo) const{
    for (std::size_t i=timeShift;i<=maxTime+1;i++) {
        const std::vector<std::size_t>&vertices=getGroupVertices(i);
        if(printToo)std::cout<<"group "<<i<<": "<<std::endl;
        for(auto &v:vertices){
          if(printToo)  std::cout<<v<<",";
            std::size_t gi=getGroupIndex(v);
            assert(gi==i);
        }
       if(printToo) std::cout<<std::endl;

    }
}

template<class T>
inline void VertexGroups<T>::initFromVector(const std::vector<std::size_t>& verticesInFrames, std::size_t timeShift_,std::size_t vertexShift_){
    timeShift=timeShift_;
    vertexShift=vertexShift_;
    maxTime=verticesInFrames.size()+timeShift;
    std::size_t groupsSize=maxTime+2-timeShift;
    groups=std::vector<std::vector<std::size_t>>(groupsSize);
    std::size_t inFrameCounter=0;
    std::size_t vertexCounter=vertexShift;
    std::size_t frameCounter=timeShift+1;
    std::vector<std::size_t> verticesInGroup;
    vToGroup=std::vector<std::size_t>();
    while(frameCounter<=maxTime){
        while(frameCounter<=maxTime&&inFrameCounter==verticesInFrames.at(frameCounter-1-timeShift)){
            groups.at(frameCounter-timeShift)=verticesInGroup;
            inFrameCounter=0;
            frameCounter++;
            verticesInGroup=std::vector<std::size_t>();
        }
        if(frameCounter<=maxTime){
            verticesInGroup.push_back(vertexCounter);
            vToGroup.push_back(frameCounter);
            inFrameCounter++;
            vertexCounter++;
        }
    }

    maxVertex=vertexCounter-1;
    assert(maxVertex==vertexShift+vToGroup.size()-1);
    std::size_t s=maxVertex+1;
    std::size_t t=maxVertex+2;

    verticesInGroup=std::vector<std::size_t>();
    verticesInGroup.push_back(s);
    vToGroup.push_back(timeShift);
    groups.at(0)=verticesInGroup;

    verticesInGroup=std::vector<std::size_t>();
    verticesInGroup.push_back(t);
    vToGroup.push_back(frameCounter);
    groups.at(frameCounter)=verticesInGroup;

    for (int i = 0; i < groups.size(); ++i) {
        for (int j = 0; j < groups[i].size(); ++j) {
            assert(vToGroup.at(groups[i][j]-vertexShift)==i);
        }
    }
    //    std::cout<<"max vertex "<<maxVertex<<std::endl;
    //    std::cout<<"max time "<<maxTime<<std::endl;



}


template<class T>
template<class PAR>
inline void VertexGroups<T>::initFromFile(const std::string& fileName, const PAR &parameters){
    std::size_t lineCounter=0;
    std::vector<std::size_t> currentGroup;
    std::vector<std::string> strings;
    std::string line;
    char delim=',';


    currentGroup=std::vector<std::size_t>();
    groups.push_back(currentGroup);  //For the s layer
    currentGroup=std::vector<std::size_t>();


    std::size_t maxTimeToRead=parameters.getMaxTimeFrame();
//    std::size_t minTimeToRead=parameters.getMinTimeFrame();
//    assert(minTimeToRead<maxTimeToRead);
//    assert(minTimeToRead>0);
    timeShift=0;  //These two are not used
    vertexShift=0;

    std::size_t minTimeToRead=parameters.getMinTimeFrame();
    std::size_t timeShiftBack=minTimeToRead-1;
    vertexShiftBack=0;


    std::ifstream timeData;
    try{
        timeData.open(fileName);
        if(!timeData){
            std::cout<<"file name "<<fileName<<std::endl;
            throw std::system_error(errno, std::system_category(), "failed to open file with vertices in time layers "+fileName);
        }

        unsigned int previousTime=1;
        unsigned int time;

        bool firstVertexFound=false;
        std::size_t firstVertex=0;

        //Skip lines until the first time within scope
        while (!firstVertexFound&&std::getline(timeData, line) && !line.empty()) {
            lineCounter++;
            strings = split(line,delim);

            if (strings.size() < 2) {
                throw std::runtime_error(
                            std::string("Vertex and time frame expected"));
            }

            unsigned int v = std::stoul(strings[0]);
            time = std::stoul(strings[1]);
            // std::cout<<v<<" "<<time<<std::endl;
            if(time>=minTimeToRead){
                firstVertexFound=true;
                firstVertex=v;
            }
        }
        assert(firstVertexFound);
        if(minTimeToRead==1){
            assert(firstVertex==0);
            vertexShiftBack=0;
        }
        else{
            vertexShiftBack=firstVertex;
        }


        currentGroup.push_back(0);
        vToGroup.push_back(time-timeShiftBack);
        previousTime=time;


        while (std::getline(timeData, line) && !line.empty()) {
            lineCounter++;
            strings = split(line,delim);

            if (strings.size() < 2) {
                throw std::runtime_error(
                            std::string("Vertex and time frame expected"));
            }

            unsigned int v = std::stoul(strings[0]);
            time = std::stoul(strings[1]);
            //std::cout<<v<<" "<<time<<std::endl;

            if(time>maxTimeToRead){
                break;
            }
          // assert(v>vertexShift);
            if(vToGroup.size()!=v-vertexShiftBack){
                std::cout<<"v to gr size "<<vToGroup.size()<<std::endl;
                std::cout<<"v - shift "<<(v-vertexShiftBack)<<std::endl;
                throw std::runtime_error(
                            std::string("Wrong vertex numbering in time file"));
            }
            else{

                vToGroup.push_back(time-timeShiftBack);

                if(time==previousTime){
                    currentGroup.push_back(v-vertexShiftBack);
                }

                else{

                    groups.push_back(currentGroup);
                    currentGroup=std::vector<std::size_t>();
                    while(groups.size()<time-timeShiftBack){
                        groups.push_back(currentGroup);
                        currentGroup=std::vector<std::size_t>();
                    }
                    currentGroup.push_back(v-vertexShiftBack);
                }
                previousTime=time;


            }
        }
        groups.push_back(currentGroup);


        //assert(vToGroup.back()<=maxTime);
        maxTime=vToGroup.back();

        maxVertex=vToGroup.size()-1;

        std::size_t s=maxVertex+1;
        std::size_t t=maxVertex+2;

        vToGroup.push_back(0);  //For s
        vToGroup.push_back(maxTime+1);  //For t

        groups.at(0).push_back(s);

        currentGroup=std::vector<std::size_t>();
        currentGroup.push_back(t);
        groups.push_back(currentGroup);

        assert(maxVertex==vToGroup.size()-3);
        //testCorrectness(false);

    }
    catch (std::system_error& er) {
        std::clog << er.what() << " (" << er.code() << ")" << std::endl;

    }
    //print();

}

//template<class T>
//inline std::vector<std::vector<std::size_t>> VertexGroups<T>::extractCorePaths(const std::vector<std::vector<std::size_t>>& paths, std::size_t cutoff,bool isFirst,bool isLast,std::size_t vertexShift) const{
//    std::size_t minT;
//    std::size_t maxT;
//    if(isFirst){
//        minT=1;
//    }
//    else{
//        minT=cutoff+1;
//    }
//    if(isLast){
//        maxT=maxTime;
//    }
//    else{
//        maxT=maxTime-cutoff;
//    }
//    return extractInnerPaths(paths,minT,maxT,vertexShift);
//}


//Vertices in paths are assumed to be numbered from zero, vShift should transform them to globalIDs
template<class T>
inline std::vector<std::vector<std::size_t>> VertexGroups<T>::extractInnerPaths(const std::vector<std::vector<std::size_t>>& paths, const std::size_t minT, const std::size_t maxT,const std::size_t vShift) const{
    //maxT inclusive
    std::vector<std::vector<std::size_t>> outputPaths;
    for (int i = 0; i < paths.size(); ++i) {
        const std::vector<std::size_t>& path=paths[i];
        std::vector<std::size_t> outputPath;
        std::size_t j=0;
        std::size_t time=0;
        while(j<paths[i].size()){

            time=getGroupIndex(paths[i][j]);
            if(time<minT){
                j++;
            }
            else{
                break;
            }
        }
        while(j<paths[i].size()){

            time=getGroupIndex(paths[i][j]);
            if(time<=maxT){
                outputPath.push_back(paths[i][j]+vShift);
                j++;
            }
            else{
                break;
            }
        }

//        for(std::size_t v: path){
//            assert(v>=vShift);
//            std::size_t vertex=v-vShift;
//            std::size_t time=getGroupIndex(vertex);
//            if(time<=maxT){
//                if(time>=minT){
//                    outputPath.push_back(vertex);
//                }
//            }
//            else break;
//        }
        if(outputPath.size()>0){
            outputPaths.push_back(outputPath);
        }
    }
    return outputPaths;
}






}

#endif // LDP_VERTEX_GROUPS_HXX
