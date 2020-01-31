import dynml as dml
import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def boolean_functions():
    funcspath = "_swarmout/bnfuncs_4/gbm_bnfuncs_4.txt"
    node_to_func_str, _ = dml.read_boolean_functions(funcspath)
    bool_funcs = {k : dml.get_boolean_function(v) for (k,v) in node_to_func_str.items()}
    return bool_funcs


@pytest.fixture
def input_table():
    input_nodes = ['EGF', 'PDGF', 'WNT_Canonical', 'TIMP', 'Hyaluronan', 
                   'ECM_Migratory_Stimuli', 'DNA_Damage', 'Oxygen']
    return dml.get_state_table(pd.Series(np.nan, index=input_nodes))


def test_get_boolnet_trajectory_asynch(boolean_functions, input_table):
    #-------------------------------------------------------------------------------------------
    input_index = 1
    rseed = 13
    nsteps = 1000
    expected_traj = pd.read_pickle(path='dynml_test_obj/test_get_boolnet_trajectory_asynch.pkl')
    #-------------------------------------------------------------------------------------------
    
    input_nodes = list(input_table.columns)
    non_input_nodes = sorted(list(set(boolean_functions.keys()) - set(input_nodes)))
    node_list = input_nodes + non_input_nodes
    
    start_state = pd.Series(np.nan, index=node_list)
    start_state[input_nodes] = input_table.iloc[input_index, :]

    np.random.seed(rseed)
    start_state[non_input_nodes] = np.random.randint(low=0, high=2, size=len(non_input_nodes))
    np.random.seed(rseed+1)
    traj = dml.get_boolnet_trajectory(init_state=start_state, inputs=start_state[input_nodes], 
        bool_funcs=boolean_functions, n_steps=nsteps, synch_update=False)
    
    assert traj.equals(expected_traj)


def test_get_boolnet_trajectory_synch(boolean_functions, input_table):
    #-------------------------------------------------------------------------------------------
    input_index = 1
    rseed = 13
    nsteps = 100
    expected_traj = pd.read_pickle(path='dynml_test_obj/test_get_boolnet_trajectory_synch.pkl')
    #-------------------------------------------------------------------------------------------

    input_nodes = list(input_table.columns)
    non_input_nodes = sorted(list(set(boolean_functions.keys()) - set(input_nodes)))
    node_list = input_nodes + non_input_nodes
    
    start_state = pd.Series(np.nan, index=node_list)
    start_state[input_nodes] = input_table.iloc[input_index, :]

    np.random.seed(rseed)
    start_state[non_input_nodes] = np.random.randint(low=0, high=2, size=len(non_input_nodes))
    np.random.seed(rseed+1)
    traj = dml.get_boolnet_trajectory(init_state=start_state, inputs=start_state[input_nodes], 
        bool_funcs=boolean_functions, n_steps=nsteps, synch_update=True)
    
    assert traj.equals(expected_traj)