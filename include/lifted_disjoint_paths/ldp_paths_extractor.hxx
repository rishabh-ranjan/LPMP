#ifndef LDP_PATHS_EXTRACTOR_HXX
#define LDP_PATHS_EXTRACTOR_HXX

#include "ldp_vertex_groups.hxx"
#include "config.hxx"

namespace LPMP {

class LdpPathsExtractor{
public:
    LdpPathsExtractor( const VertexGroups<>& vertexGroups,const std::vector<std::vector<std::size_t>>& paths, std::size_t cutoff,bool isFirst,bool isLast,std::size_t _vertexShift);
//    const VertexGroups<>& getVertexGroups()const{
//        return vg;
//    }


    const std::size_t & getMinPathsVertex()const{
        return minPathsVertex;
    }

    const std::size_t & getMaxPathsVertex()const{
        return maxPathsVertex;
    }

    const std::size_t & getMinIntevalVertex()const{
        return minIntervalVertex;
    }


    const std::vector<std::vector<std::size_t>>& getExtractedPaths()const{
        return extractedPaths;
    }

    void printExtractedPaths()const;

    const std::size_t& pathToVertex(std::size_t vertexGlobalID) const{
        assert(vertexGlobalID>=minPathsVertex);
        assert(vertexGlobalID-minPathsVertex<vertexToPath.size());
        return vertexToPath[vertexGlobalID-minPathsVertex];
    }

    std::vector<std::size_t> getTimeLayersBeforePaths()const;
    std::vector<std::size_t> getTimeLayersAfterPaths()const;





private:
  const VertexGroups<>& vg;
  bool isFirstInterval;
  bool isLastInterval;
  std::vector<std::vector<std::size_t>> extractedPaths;
  std::size_t minPathsVertex;
  std::size_t maxPathsVertex;
  std::size_t minPathsTime;
  std::size_t maxPathsTime;
  std::vector<std::size_t> vertexToPath;
  std::size_t minIntervalVertex;
  //std::size_t lastFreeBeforePaths;  //Do not use these, they do not need to exist!
  //std::size_t firstFreeAfterPaths;

};


}


#endif // LDP_PATHS_EXTRACTOR_HXX
