#include "multiway_cut/multiway_cut_instance.h"
#include "multiway_cut/multiway_cut_gaec.h"
#include <fstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

void construct_instance(LPMP::multiway_cut_instance& instance, const Eigen::Array<std::size_t, Eigen::Dynamic, 2>& edge_indices, const Eigen::Array<double, Eigen::Dynamic, 1>& edge_costs, const py::array_t<double> node_costs)
{

    for(std::size_t e=0; e<edge_costs.rows(); ++e)
    {
        const std::size_t i = edge_indices(e, 0);
        const std::size_t j = edge_indices(e, 1);
        const double w = edge_costs(e);
        instance.edge_costs.add_edge(i, j, w);
    }

    const auto _node_costs = node_costs.unchecked<2>();
    const std::size_t nr_nodes = _node_costs.shape(0);
    const std::size_t nr_labels = _node_costs.shape(1);
    std::vector<double> costs;
    for(std::size_t i=0; i<nr_nodes; ++i)
    {
        costs.clear();
        for(std::size_t l=0; l<nr_labels; ++l)
        {
            costs.push_back(_node_costs(i,l)); 
        } 
        instance.node_costs.push_back(costs.begin(), costs.end());
    } 
}

py::array_t<char> get_edge_mask(const LPMP::multiway_cut_instance& instance, const LPMP::multiway_cut_labeling& labeling)
{
    char* edge_mask = new char[instance.nr_edges()];
    for(std::size_t e=0; e<instance.nr_edges(); ++e)
        if(labeling[instance.edge_costs.edges()[e][0]] != labeling[instance.edge_costs.edges()[e][1]])
            edge_mask[e] = 1;
        else
            edge_mask[e] = 0;

    return py::array({pybind11::ssize_t(instance.nr_edges())}, edge_mask); 
} 

py::array_t<char> get_label_mask(const LPMP::multiway_cut_instance& instance, const LPMP::multiway_cut_labeling& labeling)
{
    char* label_mask = new char[instance.nr_nodes()*instance.nr_labels()];
    assert(instance.nr_nodes() == labeling.size());
    for(std::size_t i=0; i<instance.nr_nodes(); ++i)
    {
        for(std::size_t l=0; l<instance.nr_labels(); ++l)
        {
            assert(labeling[i] < instance.nr_labels());
            label_mask[i*instance.nr_labels() + l] = labeling[i] == l;
        }
    }

    return py::array({pybind11::ssize_t(instance.nr_nodes()), pybind11::ssize_t(instance.nr_labels())}, label_mask); 
} 


PYBIND11_MODULE(multiway_cut_py, m) {
    m.doc() = "python binding for LPMP asymmetric multiway cut";

    py::class_<LPMP::multiway_cut_labeling>(m, "multiway_cut_labeling")
        .def(py::init<>()); 

    py::class_<LPMP::multiway_cut_instance>(m, "multiway_cut_instance")
        .def(py::init<>())
        .def(py::init([](const Eigen::Array<std::size_t, Eigen::Dynamic, 2>& edge_indices, const Eigen::Array<double, Eigen::Dynamic, 1>& edge_costs, const py::array_t<double> node_costs) {
                    LPMP::multiway_cut_instance instance;
                    construct_instance(instance, edge_indices, edge_costs, node_costs);
                    return instance;
                    }))
        .def("evaluate", &LPMP::multiway_cut_instance::evaluate)
        .def("result_mask", [](const LPMP::multiway_cut_instance& instance, const LPMP::multiway_cut_labeling& labeling) {
                return std::make_pair(
                        get_edge_mask(instance, labeling),
                        get_label_mask(instance, labeling)
                        ); 
                })
        .def("write", [](const LPMP::multiway_cut_instance& instance, const std::string& filename)
                {
                std::ofstream f;
                f.open(filename);
                instance.write(f); 
                f.close();
                })
//        .def("read", [](LPMP::multiway_cut_instance& instance, const std::string& filename)
//                {
//                instance = LPMP::multiway_cut_parser::parse_file(filename);
//                })
                ;

        m.def("multiway_cut_gaec", [](const LPMP::multiway_cut_instance& instance) {
                return LPMP::multiway_cut_gaec(instance);
                });
}
