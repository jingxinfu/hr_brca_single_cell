import pickle
import scPipe as sp
RESULT_TABLE = '/home/analysis/hr_brca_single_cell/data/table'
for POPULATION in [f'TumorRandom95_{i}' for i in range(10)]:
    with open(f'{RESULT_TABLE}/MPs/{POPULATION}/nmf_basis.pickle', "rb") as input_file:
        programs_basis = pickle.load(input_file)
    
    with open(f'{RESULT_TABLE}/MPs/{POPULATION}/nmf_coef.pickle', "rb") as input_file:
        programs_coef = pickle.load(input_file)
        
    n_programs = programs_basis[list(programs_basis.keys())[2]].shape[1] * len(programs_basis)
    print(f"Generated {n_programs:,} programs")
    
    n_top=50
    robustRPH = sp.ext.getRobustRHP(programs=programs_basis,
                             n_top=n_top,
                             intra_min = 35, inter_min = 10,intra_max = 10
                            )
    with open(f'{RESULT_TABLE}/MPs/{POPULATION}/robustRPH.pickle', 'wb') as handle:
        pickle.dump(robustRPH, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    with open(f'{RESULT_TABLE}/MPs/{POPULATION}/robustRPH.pickle', "rb") as input_file:
        robustRPH = pickle.load(input_file)
    print(f"Out of {n_programs:,} detected programs, {robustRPH.shape[1]} programs are robust.")
    print(f'FINISH {POPULATION}')