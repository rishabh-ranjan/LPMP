
#include "graph_matching/graph_matching.h"
#include "visitors/standard_visitor.hxx"

using namespace LPMP;
int main(int argc, char** argv) {
ProblemConstructorRoundingSolver<Solver<LP<FMC_HUNGARIAN_BP_T>,StandardTighteningVisitor>>solver(argc,argv);
auto input = LPMP::TorresaniEtAlInput::parse_file(solver.get_input_file());
solver.GetProblemConstructor().construct(input);
return solver.Solve();
}