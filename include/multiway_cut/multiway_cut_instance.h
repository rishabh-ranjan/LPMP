#pragma once

#include "multicut/multicut_instance.h"
#include <vector>
#include <cassert>

namespace LPMP {

    class multiway_cut_node_costs
    {
        private:
            std::size_t nr_labels_ = 0;
            std::vector<double> costs;

        public:
            std::size_t nr_labels() const { return nr_labels_; }
            std::size_t nr_nodes() const { return costs.size()/nr_labels_; }
            double operator()(const std::size_t var, const std::size_t label) const; 
            template<typename ITERATOR>
                void push_back(ITERATOR cost_begin, ITERATOR cost_end);
    };

    inline double multiway_cut_node_costs::operator()(const std::size_t var, const std::size_t label) const
    {
        assert(var < nr_nodes());
        assert(label < nr_labels());
        return costs[var*nr_labels() + label]; 
    }

    class multiway_cut_labeling : public std::vector<std::size_t> {
    };

    class multiway_cut_instance {
        public:
            multiway_cut_node_costs node_costs;
            multicut_instance edge_costs;
            std::size_t nr_nodes() const { return node_costs.nr_nodes(); }
            std::size_t nr_labels() const { return node_costs.nr_labels(); }
            std::size_t nr_edges() const { return edge_costs.no_edges(); }
            double evaluate(const multiway_cut_labeling& labeling) const;
            bool feasible(const multiway_cut_labeling& labeling) const;
            template<typename STREAM>
                void write(STREAM& s) const;
    };

    ////////////////////
    // implementation //
    ////////////////////

    template<typename ITERATOR>
        void multiway_cut_node_costs::push_back(ITERATOR cost_begin, ITERATOR cost_end)
        {
            if(nr_labels_ == 0)
                nr_labels_ = std::distance(cost_begin, cost_end);

            assert(std::distance(cost_begin, cost_end) > 0);
            assert(std::distance(cost_begin, cost_end) == nr_labels_);

            for(auto it=cost_begin; it!=cost_end; ++it)
                costs.push_back(*it); 
        }

    template<typename STREAM>
        void multiway_cut_instance::write(STREAM& s) const
        {           
            s << "ASYMMETRIC MULTIWAY CUT\n";
            s << "MULTICUT\n";
            for(const auto& e : edge_costs.edges())
                s << e[0] << " " << e[1] << " " << e.cost << "\n";

            s << "NODE COSTS\n";
            for(std::size_t i=0; i<nr_nodes(); ++i)
            {
                for(std::size_t l=0; l<nr_labels(); ++l)
                {
                    if(l > 0)
                        s << " ";
                    s << node_costs(i,l);
                }
                s << "\n";
            }
        }
}

