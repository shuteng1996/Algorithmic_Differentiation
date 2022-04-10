DCO_INC_PATH=$(HOME)/Software/dco_cpp/include
DCO_LIB_PATH=$(HOME)/Software/dco_cpp/lib

all: compile_ref_para compile_draft_para


compile_ref_para:
	mpicxx -std=c++17 -O3 -DDCO_DISABLE_AUTO_WARNING -DDCO_DISABLE_AVX2_WARNING ref_parallel.cpp -o result_ref.out

compile_draft_para:
	mpicxx -std=c++17 -O3 -I $(DCO_INC_PATH) -L $(DCO_LIB_PATH) -DDCO_DISABLE_AUTO_WARNING -DDCO_DISABLE_AVX2_WARNING draft_parallel.cpp -ldcoc -o result_draft.out

compile_draft_multivar:
	g++ -std=c++17 -O3 -I $(DCO_INC_PATH) -L $(DCO_LIB_PATH) -DDCO_DISABLE_AUTO_WARNING -DDCO_DISABLE_AVX2_WARNING  draft_s_multivar_vector.cpp -ldcoc -o multi_var.out

compile_draft_demo:
	g++ -std=c++17 -O3 -I $(DCO_INC_PATH) -L $(DCO_LIB_PATH) software_deco.cpp -ldcoc -o demo.out
compile_p_draft_multivar:
	mpicxx -std=c++17 -O3 -I $(DCO_INC_PATH) -L $(DCO_LIB_PATH) -DDCO_DISABLE_AUTO_WARNING -DDCO_DISABLE_AVX2_WARNING draft_p_multivar_vector.cpp -ldcoc -o multi_var_p.out

compile_p_draft_multivar_new:
	mpicxx -std=c++17 -O3 -I $(DCO_INC_PATH) -L $(DCO_LIB_PATH) -DDCO_DISABLE_AUTO_WARNING -DDCO_DISABLE_AVX2_WARNING draft_p_multivar_vec_new3.cpp -ldcoc -o draft_p_multivar_vector_new3.out
 
run_ref: compile_ref_para
	mpiexec -n 8 ./result_ref.out

run_draft: compile_draft_para
	mpiexec -n 8 ./result_draft.out

run_multi_var_draft: compile_draft_multivar
	./multi_var.out

run_p_multi_var_draft: compile_p_draft_multivar
	mpiexec -n 4 ./multi_var_p.out

run_p_multi_var_draft_new: compile_p_draft_multivar_new
	mpiexec -n 10 ./draft_p_multivar_vector_new3.out

clean:
	rm -f result_ref.out result_draft.out

