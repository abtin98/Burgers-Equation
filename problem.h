#include "util.h"

class Problem
{
	public:
		Problem();
		void run();
	private:
		void make_grid(int n_refinements);
		void setup_system();
		void assemble_system();
		void solve();
		void output();

		Triangulation<1> triangulation;
		FE_DGQ<1> fe;
		DoFHandler<1> dof_handler;

		SparsityPattern sparsity_pattern;
		SparseMatrix<double> system_matrix;

		Vector<double> solution;
		Vector<double> system_rhs;

};
