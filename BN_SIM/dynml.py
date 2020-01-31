# Documentation Standards:
# https://numpydoc.readthedocs.io/en/latest/format.html#overview
# https://numpydoc.readthedocs.io/en/latest/example.html#example

# Unit Testing
# https://docs.pytest.org/en/latest/

# Performance Profiling:
# https://github.com/joerick/pyinstrument

# CODE REVIEW NOTES:
# + Make sure dictionaries are not iterated over and updated at the
#   at the same time (see Fluent Python, page 92).
#-----------------------------------------------------------------------------
import errno
import itertools
import os
import pickle
import random
import re
import shutil
import subprocess
import time
import xml.etree.ElementTree

#import multiprocessing as mp
import numpy as np
import pandas as pd

from collections import OrderedDict
from collections import defaultdict
from functools import partial

# For background:
# https://stackoverflow.com/questions/8804830/python-multiprocessing-pickling-error
# https://medium.com/@jwnx/multiprocessing-serialization-in-python-with-pickle-9844f6fa1812
# http://matthewrocklin.com/blog/work/2013/12/05/Parallelism-and-Serialization
# pip install --user git+https://github.com/uqfoundation/pathos (command line on linux system)
from pathos.multiprocessing import ProcessingPool as Pool

# bnsimulation ================================================================

# =============================================================================

# file system utils -----------------------------------------------------------
def copy_dir_tree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)

def make_dir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def make_tmpdir(prefix, path = None):
    if path is None:
        path = os.getcwd()
    tmpdir_path = os.path.join(path, prefix + time.strftime("%Y%m%d_%H%M%S"))
    make_dir(tmpdir_path)

    return tmpdir_path
# -----------------------------------------------------------------------------

# boolean function-related processing -----------------------------------------
def get_bool_func_inputs(rule):
        rule = rule.replace("(", " ")
        rule = rule.replace(")", " ")
        rule = re.sub(r"\band\b", " ", rule)
        rule = re.sub(r"\bor\b", " ", rule)
        rule = re.sub(r"\bnot\b", " ", rule)
        rule = re.split("\s+", rule.strip())
        return rule

def read_boolean_functions(file_path):
    with open(file_path) as f:
        tmp = f.readlines()

    if tmp[0][0] == "#":
        tmp = tmp[1:]
    
    tmp = [tuple(re.split("\s*\*\s*=\s*", x.strip())) for x in tmp]
    bool_funcs = {k:v for (k, v) in tmp}

    bool_func_inputs = {}
    for k in bool_funcs.keys():
        operands = list(set(get_bool_func_inputs(bool_funcs[k]))) #rm dups
        if operands[0] in ['0', '1']:
            bool_func_inputs[k] = int(operands[0])
        else:
            assert(all([x in bool_funcs.keys() for x in operands]))
            bool_func_inputs[k] = operands
    
    return bool_funcs, bool_func_inputs


def write_boolean_functions(bool_funcs, file_path, from_dict = True):
    with open(file_path, 'w') as fp:
        fp.write("#BOOLEAN RULES\n")
        if from_dict:
            for node in sorted(bool_funcs.keys()):
                fp.write(node + "*=" + bool_funcs[node] + "\n")
        else:
            for rule in bool_funcs:
                fp.write("%s\n" % rule)


def get_boolean_function(f_str):
    #------------------------------------------------------------------------
    def f(x, check_inputs=False):
        eval_exp = f_str.strip()
        if eval_exp == '0':
            fx = 0
        elif eval_exp == '1':
            fx = 1
        else:
            for var in x.keys():
                val = int(x[var])
                assert((val == 0) or (val == 1))
                eval_exp = re.sub(r"\b%s\b" % var, str(val), eval_exp)

            try:
                if (check_inputs):
                    tmp = eval_exp
                    tmp = re.sub(r"\band\b", "", tmp)
                    tmp = re.sub(r"\bor\b", "", tmp)
                    tmp = re.sub(r"\bnot\b", "", tmp)
                    tmp = tmp.replace(")", "")
                    tmp = tmp.replace("(", "")
                    tmp = tmp.replace(" ", "")
                    tmp = set(tmp)
                    if not tmp.issubset({'0', '1'}):
                        raise(ValueError())

                fx = eval(eval_exp)
                assert((fx == 0) or (fx == 1))
                fx = int(fx)
            except:
                raise(ValueError("Inputs: " + str(x) +
                                 " do not match function: " + f_str + "."))

        return fx
    #------------------------------------------------------------------------
    return f



def get_boolean_functions_from_pathvisio(gpml_file, one_rule_per_node=True):
    root = xml.etree.ElementTree.parse(gpml_file).getroot()
    data_nodes = root.findall('{http://pathvisio.org/GPML/2013a}DataNode')
    shape_nodes = root.findall('{http://pathvisio.org/GPML/2013a}Shape')
    node_comment_lists = [node.findall('{http://pathvisio.org/GPML/2013a}Comment') 
                          for node in data_nodes + shape_nodes]
    
    # This is a list of rule strings, e.g. A *= B and C
    bnrules = [L[0].text.strip() for L in node_comment_lists if len(L) > 0]
    bnrules = list(set(bnrules))
    
    f_output_nodes = [rule.split("*=")[0].strip() for rule in bnrules]
    if one_rule_per_node:
        assert(len(f_output_nodes) == len(set(f_output_nodes)))
    
    f_input_nodes = []
    for rule in bnrules:
        f_def_str = rule.split("*=")[1].strip()
        if (f_def_str == '0') or (f_def_str == '1'):
            continue
        f_input_nodes.extend(get_bool_func_inputs(f_def_str))
        
    f_input_nodes = list(set(f_input_nodes))
    
    nodes_to_add = set(f_input_nodes) - set(f_output_nodes)
    for n2add in nodes_to_add:
        bnrules.append(n2add + " *= 0")
    
    # sort rules--------------------------------------------------------------------
    # The aim here is to ensure an alphabetical sort order, but with input nodes
    # (i.e., with fixed boolean values 0 or 1) placed at the top of the list.
    def get_sort_key(s):
        node, rule = [x.strip() for x in s.split("*=")]
        if (rule == '0') or (rule == '1'):
            key = rule + node
        else:
            key = node

        return key
    
    bnrules.sort(key=get_sort_key)
    # -------------------------------------------------------------------------------
    
    return bnrules


def fix_boolean_functions_for_inputs(bool_funcs, inputs):
    assert(isinstance(inputs, pd.core.series.Series))
    new_bool_funcs = bool_funcs.copy()
    for node in inputs.index:
        fixed_val = int(inputs[node])
        assert((fixed_val == 0) or (fixed_val == 1))
        #assert(node in new_bool_funcs)
        new_bool_funcs[node] = get_boolean_function(str(fixed_val))
    
    return new_bool_funcs




# -----------------------------------------------------------------------------

# wrapper functions for Albert Stable Motifs/Dynamic Modeling -----------------
def run_stable_motifs(bool_funcs, lib_path,
                      work_dir_path=None, rm_work_dir=True, timeout=300):
    orig_dir = os.getcwd()
    work_dir = make_tmpdir("stable_motifs_")

    #print(" --- at start: " + os.getcwd())
    os.chdir(work_dir)
    #print(" --- for processing" + os.getcwd())
    # -------------------------------------------------------------------------
    try:
        results = {}
        copy_dir_tree(src=lib_path, dst=work_dir)
        write_boolean_functions(bool_funcs, "boolnet.txt")
        cmd = "java -jar StableMotifs.jar boolnet.txt"
        #errcode = subprocess.call([cmd], shell=True)
        #if errcode:
        #    raise Exception("StableMotifs.jar returned with: " + errcode)
        
        subprocess.check_call([cmd], shell=True, timeout=timeout)

        # process attrators ---------------------------------------------------
        attr_tab = pd.read_table("./boolnet-QuasiAttractors.txt",
                                 na_values="X")
        results["attractors"] = attr_tab[sorted(list(bool_funcs.keys()))]
        # ---------------------------------------------------------------------
    finally:
        os.chdir(orig_dir)
        if rm_work_dir:
            shutil.rmtree(work_dir)
    # -------------------------------------------------------------------------
    
    return results


def merge_stable_motifs_results(inputs, attr_sets):
    # CHECK: input variable ordering in attractor data frames is fixed.
    attr_df_colnames = [df.columns  for _, df in attr_sets]
    ref_colnames = attr_df_colnames[0]
    assert(all([all(cnames == ref_colnames)  for cnames in attr_df_colnames]))
    
    # Concatenate attractor set data frames for different input combinations
    # and label the overall attractor set.
    attr_df = pd.concat([df  for _, df in attr_sets])
    attr_label_lists = [[(label + '_attr_' + str(i+1))  for i in df.index]  
                         for label, df in attr_sets]
    attr_labels = [label for L in attr_label_lists for label in L]
    attr_df.index = attr_labels
    
    # Re-order columns so inputs come first.
    non_inputs = list(set(attr_df.columns) - set(inputs))
    reordered_cols = list(inputs) + sorted(non_inputs)
    attr_df = attr_df[reordered_cols]
    
    return attr_df


def run_stable_motifs_for_input_combinations(inputs, bool_func_strings,
    stable_motifs_lib_path, rm_work_dir=True, timeout=300, output_file=None,
    verbose=True):
    assert(all([var in bool_func_strings.keys()  for var in inputs]))
    input_vals = get_state_table(pd.Series(np.nan, index=inputs))
    input_vals.index = input_vals.index + 1
    
    # Make nan-filled attractor state data frame for instances in which
    # the stable motif (attractor finding) code fails or times out.
    bn_nodes = sorted(bool_func_strings.keys())
    nan_attr_df = pd.DataFrame(np.full((1, len(bn_nodes)), np.nan), columns=bn_nodes)
    
    attr_sets = []
    for i in input_vals.index:
        input_label = 'in_' + str(i)
        f_dict = bool_func_strings.copy()
        for var in inputs:
            f_dict[var] = str(input_vals.loc[i, var])
        
        try:
            results = run_stable_motifs(
                bool_funcs=f_dict, 
                lib_path=stable_motifs_lib_path, 
                rm_work_dir=rm_work_dir,
                timeout=timeout
            )
        except Exception:
            _nan_attr_df = nan_attr_df.copy()
            for var in inputs:
                _nan_attr_df.loc[0, var] = input_vals.loc[i, var]
            attr_sets.append((input_label, _nan_attr_df))
        else:
            attr_sets.append((input_label, results['attractors']))
        
        if len(attr_sets) > 1:
            attr_df = merge_stable_motifs_results(inputs, attr_sets)
            if output_file:
                with open(output_file, 'wb') as pkl_file:
                    pickle.dump(attr_df, pkl_file)
            
        if verbose: 
            print(str(i) + " ", end='') 
            if (i % 40) == 0: 
                print("\n", end='')
        
    if verbose: print("\n", end='')
            
    return attr_df
        
    
#------------------------------------------------------------------------------

# attractor state-related processing ------------------------------------------
def sample_from_boolean_attractor(attr, attr_name, N, off_mean, off_sd,
                                  on_mean, on_sd):
    template = attr.values
    num_nan = sum(np.isnan(template))

    d = {}
    for i in range(N):
        template_nonan = template.copy()
        if num_nan > 0:
            template_nonan[np.isnan(template)] = np.random.randint(
                    low=0, high=2, size=num_nan)
        num_off = sum(template_nonan == 0)
        num_on  = sum(template_nonan == 1)

        vals = np.zeros(template_nonan.size)
        vals[template_nonan == 0] = np.random.normal(
                loc=off_mean, scale=off_sd, size=num_off)
        vals[template_nonan == 1] = np.random.normal(
                loc=on_mean, scale=on_sd, size=num_on)
        d[attr_name + "_" + str(i)] = pd.Series(vals, index=attr.index)

    return pd.DataFrame(d).T


def sample_from_boolean_attractors(attr_df, attr_name_prefix, N, off_mean,
                                   off_sd, on_mean, on_sd):
    frames = []
    for i in range(attr_df.shape[0]):
        df = sample_from_boolean_attractor(
            attr = attr_df.loc[i],
            attr_name=attr_name_prefix + "_" + str(i),
            N=N, off_mean=off_mean, off_sd=off_sd, on_mean=on_mean, on_sd=on_sd
        )
        frames.append(df)

    return pd.concat(frames)
# -----------------------------------------------------------------------------

def get_evidence_for_direct_relationship(a, b, bool_funcs, cor_mat,
    orig_sim_data_N, abs_cor_thresh = 0.9,
    sim_data_off_mean=5, sim_data_off_sd=0.5,
    sim_data_on_mean=10, sim_data_on_sd=0.5,
    stable_motifs_lib_path="/Users/vinodh/Dropbox/pyproj/calbert/lib/StableMotifs/dist/"):

    results = {}
    results["node_a"] = a
    results["node_b"] = b
    results["cor_ab"] = cor_mat.loc[a, b]
    # (1) find shared correlates --------------------------------------------
    tmp = np.abs(cor_mat)
    tmp = tmp.loc[[a, b], :]
    tmp = tmp.drop([a, b], axis=1)
    selector = (tmp.loc[a] > abs_cor_thresh) & (tmp.loc[b] > abs_cor_thresh)
    tmp = tmp.loc[:, selector]
    results["shared_correlates"] = tmp
    # (2) update boolean rules ----------------------------------------------
    bool_funcs_new = bool_funcs.copy()
    for node in results["shared_correlates"].columns.values:
        bool_funcs_new[node] = '0'
    # (3) find attractors ---------------------------------------------------
    stmotifs_results = run_stable_motifs(
        bool_funcs=bool_funcs_new,
        lib_path=stable_motifs_lib_path,
        rm_work_dir=True
    )
    results["attractors_post_ko"] = stmotifs_results['attractors']
    # (4) sample from updated attractor set and recompute cor(a,b) ----------
    num_attractors_post_ko = results["attractors_post_ko"].shape[0]
    if num_attractors_post_ko > 0:
        X_ko = sample_from_boolean_attractors(
            attr_df=results["attractors_post_ko"],
            attr_name_prefix = "attr",
            N = orig_sim_data_N // num_attractors_post_ko,
            off_mean = sim_data_off_mean, off_sd = sim_data_off_sd,
            on_mean  = sim_data_on_mean,  on_sd  = sim_data_on_sd
        )
        results["cor_ab_post_ko"] = X_ko.corr().loc[a, b]
        results["abs_cor_ratio"]  = np.abs(results["cor_ab_post_ko"] / results["cor_ab"])
    else:
        results["cor_ab_post_ko"] = None
        results["abs_cor_ratio"]  = None

    return results


def find_similar_pairs(S, threshold):
    idset = S.index.values
    assert(all(idset == S.index.values))

    sim_pairs = set()
    for a in idset:
        for b in idset:
            if (a != b) and (S.loc[a, b] > threshold):
                sim_pairs.add(tuple(sorted([a, b])))

    return sim_pairs


def binary_state_to_int(state):
    return int("".join([str(int(i)) for i in state]), base=2)

def binary_state_to_str(state, null_char='?'):
    return "".join([null_char if pd.isnull(i) else str(int(i))  for i in state])


def binary_state_str_list_to_df(state_str_list, column_names=None):
    nrows = len(state_str_list)
    assert(nrows > 0)
    ncols = len(state_str_list[0])
    assert(all([len(s) == ncols  for s in state_str_list]))
    if column_names is None:
        column_names = list(range(ncols))
    assert(len(column_names) == ncols)
    
    df = pd.DataFrame(np.full((nrows, ncols), 0), index=state_str_list, columns=column_names)
    for s in df.index:
        s_values = [int(i) for i in s]
        assert(all([val in [0, 1]  for val in s_values]))
        df.loc[s, :] = s_values
    
    assert(all(df.apply(binary_state_to_str, axis=1).values == df.index))
    
    return df
    
    
def get_binary_states_as_integers(df):
    assert(isinstance(df, pd.core.frame.DataFrame))
    return df.apply(axis=1, func=binary_state_to_int)


def update_state(state, node, func, as_int=True):
    assert(node in state.index)
    updated_state = state.copy()
    updated_state[node] = func(state)
    updated_state.name = binary_state_to_int(updated_state)
    if as_int:
        updated_state = updated_state.name
    
    return updated_state


# -------[(Check) Added from calbert notebook]--------------------------------------
# [update_state(s, node, func) for node, func in bool_funcs.items()]
def get_markov_matrix(state_tab, bool_funcs, synch_update=False):
    assert(not synch_update) # TO DO (synch_update)
    N = state_tab.shape[0]
    T = pd.DataFrame(np.full((N, N), 0), columns=state_tab.index,
                    index = state_tab.index)
    
    for id in T.index:
        s = state_tab.loc[id, :]
        d = [update_state(s, node, func) for node, func in bool_funcs.items()]
        assert(all([(i in T.index) for i in d]))
        uniq_s, uniq_s_count = np.unique(d, return_counts=True)
        uniq_s_freq = uniq_s_count / sum(uniq_s_count)
        T.loc[id, uniq_s] = uniq_s_freq
    
    row_sums = T.apply(axis=1, func=sum)
    assert(all([np.isclose(s, 1)  for s in row_sums]))
        
    return T



def boolnet_state_generator(init_state, bool_funcs, synch_update=False):
    node_set = init_state.index.values
    assert(all([(node in bool_funcs) for node in node_set]))
    prev_state = init_state
    
    for t in itertools.count(start=1):
        if synch_update:
            update = pd.Series(
                data=[bool_funcs[node](prev_state) for node in node_set],
                index=init_state.index, name=t) 
        else:
            update = prev_state.copy()
            rand_node = np.random.choice(node_set, 1)[0]
            update[rand_node] = bool_funcs[rand_node](prev_state)
            update.name = t
        
        yield (t, update)
        prev_state = update

        

def get_boolnet_trajectory(init_state, bool_funcs, n_steps=100, synch_update=False, inputs=None,
                           target_states=None, return_none_if_target_not_reached=True):
    assert(isinstance(init_state, pd.core.series.Series))
    if target_states is not None:
        if isinstance(target_states, pd.core.series.Series):
            target_states = pd.DataFrame(target_states).T
        assert(isinstance(target_states, pd.DataFrame))
        reached_a_target = False
    
    if inputs is not None:
        assert(isinstance(inputs, pd.core.series.Series))
        assert(all([(node in init_state.index) for node in inputs.index]))
        assert(all([(node in bool_funcs) for node in inputs.index]))
        bool_funcs = fix_boolean_functions_for_inputs(bool_funcs, inputs)
    
    node_set = init_state.index.values
    trajectory = pd.DataFrame(np.full((n_steps + 1, node_set.size), np.nan), 
                              columns=node_set, dtype=np.int8)
    trajectory.iloc[0, :] = init_state
    g = boolnet_state_generator(init_state, bool_funcs, synch_update)
    
    # Note that t starts at 1, based on boolnet_state_generator implementation.
    # This could be made a parameter, to make this setting clearer.
    #n_consecutive_equal_states = 0
    for t, state in g:
        if t > n_steps:
            break
        else:
            trajectory.iloc[t, :] = state.values
            
            if target_states is not None:
                matches = target_states.apply(axis=1, func=lambda x: all(x == state))
                if any(matches):
                    reached_a_target = True
                    trajectory = trajectory.iloc[:(t+1), :]
                    break # IF one of the target states is reached ... 
            #elif n_for_steady_state > 1:
            #    if trajectory.iloc[t-1:t+1, :].drop_duplicates().shape[0] == 1:
            #        n_consecutive_equal_states += 1
            #       if n_consecutive_equal_states == n_for_steady_state:
            #           trajectory = trajectory.iloc[:(t+1), :]
            #           break # IF the state is unchanged over n_for_steady_state steps ...
            #   else:
            #       n_consecutive_equal_states = 0

                        
    if (target_states is not None) and (not reached_a_target):
        if return_none_if_target_not_reached:
            trajectory = None
            
    return trajectory



def get_trajectories_to_attractors(N, bool_funcs, attr_states, n_steps=10000, synch_update=False, 
    inputs=None, exclude_inputs=True, states_as_str=True, keep_if_attr_not_reached=False,
    verbose=False, output_file_path=""):
    if isinstance(attr_states, pd.core.series.Series):
        attr_states = pd.DataFrame(attr_states).T
    
    if inputs is not None:
        assert(isinstance(inputs, pd.core.series.Series))
        assert(all([(node in attr_states.columns)  for node in inputs.index]))
        non_input_nodes = list(set(attr_states.columns) - set(inputs.index))
    else:
        non_input_nodes = list(attr_states.columns)
        
    d = {}
    for i in range(N):
        init_state = attr_states.iloc[0, :].copy()
        if inputs is not None:
            init_state[inputs.index] = inputs
        init_state[non_input_nodes] = np.random.randint(low=0, high=2, size=len(non_input_nodes))
        
        traj = get_boolnet_trajectory(init_state=init_state, bool_funcs=bool_funcs, 
            n_steps=n_steps, synch_update=synch_update, inputs=inputs, target_states=attr_states, 
            return_none_if_target_not_reached=True)
        
        if traj is None:
            if keep_if_attr_not_reached:
                d[binary_state_to_str(init_state)] = traj
            continue
        
        if (inputs is not None) and exclude_inputs:
            traj = traj.loc[:, non_input_nodes]
        
        if states_as_str:
            traj = traj.apply(binary_state_to_str, axis=1)
        
        d[binary_state_to_str(init_state)] = traj
        
        if output_file_path and ((i % 10) == 0):
            pickle.dump(d, open(output_file_path, 'wb'))
        
        if verbose: 
            print(str(i+1) + " ", end='')
            if ((i+1) % 30) == 0: print("\n", end='')
        
    if verbose: 
        print("\n", end='')
    
    return d



def get_initial_states_from_trajectories(d, k=1):
    initial_states = []
    for init_state, traj in d.items():
        assert(traj is not None)
        if not isinstance(traj, pd.core.series.Series):
            traj = traj.apply(binary_state_to_str, axis=1)
        
        initial_states.extend(list(traj.iloc[:k].values))
    
    return initial_states



def get_transient_states_from_trajectories(d):
    transient_states = []
    for init_state, traj in d.items():
        assert(traj is not None)
        if not isinstance(traj, pd.core.series.Series):
            traj = traj.apply(binary_state_to_str, axis=1)
        
        # Exclude last state indexed by -1 (assumed to be attractor).
        transient_states.extend(list(traj.iloc[:-1].values))
    
    return transient_states



def get_near_attractor_states_from_trajectories(d, steps_from_attr=1):
    # d is map (dictionary) from initial states to trajectories (to an attractor).
    # dd is a map from attractor states to lists of near attractor states.
    dd = defaultdict(list)
    for init_state, traj in d.items():
        assert(traj is not None)
        if not isinstance(traj, pd.core.series.Series):
            traj = traj.apply(binary_state_to_str, axis=1)
        
        attr = traj.iloc[-1]
        
        # Exclude last state indexed by -1 (assumed to be attractor).
        if steps_from_attr >= len(traj):
            dd[attr].extend(list(traj.iloc[:-1].values))
        else:
            dd[attr].extend(list(traj.iloc[-(steps_from_attr+1):-1].values))
    
    for attr in dd.keys():
        # Remove duplicates.
        dd[attr] = list(set(dd[attr]))
    
    return dd


# TO DO: Check edge cases, possibly re-write.
def sample_from_transient_states(states, N, unique=True):
    sampled_states = list(np.random.choice(states, size=N, replace=True))
    if unique:
        assert(N <= len(set(states)))
        n_states = len(states)
        sampled_fraction = N / n_states
        sampled_fraction_incr = 0.01 * sampled_fraction
        sampled_states = []

        while sampled_fraction <= 1 and (len(sampled_states) < N):
            n_sampled = int(sampled_fraction * n_states)
            sampled_states = list(set(np.random.choice(states, size=n_sampled, replace=True)))
            sampled_fraction += sampled_fraction_incr
        
    return sampled_states


                         
def get_state_counts(states_list):
    uniq_state_counts = pd.Series(states_list).value_counts()
    uniq_state_counts.sort_values(ascending=False, inplace=True)
    
    return uniq_state_counts



def get_state_frequency(states_list):
    uniq_state_counts = get_state_counts(states_list)
    state_freq = uniq_state_counts / sum(uniq_state_counts)
    
    return state_freq



def compute_state_potential(d, min_num_uniq_transient_states=1000,
                           n_traj_start=1, n_traj_end=1):
    transient_states = get_transient_states_from_trajectories(d)
    transient_state_freq  = get_state_frequency(transient_states)
    U = -np.log(transient_state_freq)
    
    attr_to_nearby_states = get_near_attractor_states_from_trajectories(d, 
        steps_from_attr=n_traj_end)
    near_attr_states = set()
    for k, v in attr_to_nearby_states.items():
        U[k] = np.min(U.loc[v]) - 1
        near_attr_states.update(v)
    
    # Note: U now has potential values for all states (including attractors).
    
    init_states = get_initial_states_from_trajectories(d, k=n_traj_start)
    sampled_transient_states = sample_from_transient_states(transient_states, 
        N=min_num_uniq_transient_states)
    near_attr_states = list(near_attr_states)
    attr_states = list(attr_to_nearby_states.keys())
    
    final_state_set = set()
    final_state_set.update(init_states, sampled_transient_states,
                           near_attr_states, attr_states)
    
    U_final_state_set = U.loc[final_state_set]
    
    return U_final_state_set



def sample_from_trajectory(traj, n_kept_at_start=1, n_kept_at_end=10, middle_frac=0.10, copy=True):
    assert(0 <= middle_frac <= 1)
    N = traj.shape[0]
    n_middle = N - n_kept_at_start - n_kept_at_end
    if n_middle <= 0:
        return traj
    
    n_from_middle = int(n_middle * middle_frac)
    middle_index_set = list(range(n_kept_at_start, N - n_kept_at_end))
    
    # Note: it seems better to select each middle index set member with some
    # probability p, rather than
    sampled_middle_index_set = random.sample(middle_index_set, n_from_middle)
    
    index_set = list(range(n_kept_at_start)) + sampled_middle_index_set
    index_set.extend(list(range(N - n_kept_at_end, N)))
    index_set.sort()
    
    sampled_traj = traj.iloc[index_set, :]
    if copy:
        sampled_traj = sampled_traj.copy()
    
    return sampled_traj



def get_reached_attractor(init_state, bool_funcs, n_steps=200,
                          num_last_states_checked=20, max_cycle_length=1,
                          synch_update=False, inputs=None):
    assert(n_steps >= (num_last_states_checked * 5))
    traj = get_boolnet_trajectory(init_state=init_state, bool_funcs=bool_funcs,
                                  n_steps=n_steps, synch_update=synch_update,
                                  inputs=inputs)
    traj = add_state_info_cols(traj)
    attr = traj.iloc[-num_last_states_checked:, :].drop_duplicates()
    attr.index = attr['STATE_ID']
    
    #attr = traj.iloc[-num_last_states_checked:, :].drop_duplicates()
    #attr.index = attr.apply(binary_state_to_str, axis=1)
    
    if (attr.shape[0] > max_cycle_length):
        attr = None
    
    return attr



def simulate_to_target_states(init_state, bool_funcs, target_states,
                              n_steps=1000, synch_update=False):
    if isinstance(target_states, pd.core.series.Series):
        target_states = pd.DataFrame(target_states).T
        
    g = boolnet_state_generator(init_state, bool_funcs, synch_update)
    matched_target = None
    
    for t, state in g:
        if t > n_steps:
            break
        else:  
            matches = target_states.apply(axis=1, func=lambda x: all(x == state))
            if any(matches):
                matched_target = matches[matches].index[0]
                break
    
    return t, matched_target



def get_target_state_arrival_freq(init_state, bool_funcs, target_states,
                                 n_trials=100, max_steps_per_trial=100,
                                 synch_update=False, input_nodes=None,
                                 randomized_nodes=None):
    assert(isinstance(init_state, (pd.Series, pd.DataFrame)))
    assert('other' not in target_states.index)
    
    if isinstance(init_state, pd.Series):
        init_state = pd.DataFrame(init_state).T
    
    if input_nodes:
        assert(all([(node in init_state.columns) for node in input_nodes]))
        assert(all([(node in bool_funcs) for node in input_nodes]))
        
    if randomized_nodes:
        assert(all([(node in init_state.columns) for node in randomized_nodes]))
        
    ar_freqs = pd.DataFrame(0, index = init_state.index,
        columns = [str(id) for id in target_states.index] + ['other'])
    
    state_bool_funcs = bool_funcs
    for id in ar_freqs.index:
        # If there are input nodes, we must fix their update functions 
        # to the appropriate constant value for each state in the
        # init_state table.
        if input_nodes:
            state_bool_funcs = bool_funcs.copy()
            for node in input_nodes:
                fixed_val = int(init_state.loc[id, node])
                assert((fixed_val == 0) or (fixed_val == 1))
                state_bool_funcs[node] = get_boolean_function(str(fixed_val))
            
        for i in range(n_trials):
            trial_init_state = init_state.loc[id, :]
            # If there are nodes to be set to random states in each trial,
            # we do so here.
            if randomized_nodes:
                trial_init_state = trial_init_state.copy()
                trial_init_state[randomized_nodes] = np.random.randint(
                    low=0, high=2, size=len(randomized_nodes))
                
            t, target_id = simulate_to_target_states(trial_init_state,
                state_bool_funcs, target_states, n_steps=max_steps_per_trial,
                synch_update=synch_update)
            target_id = str(target_id)

            if target_id in ar_freqs:
                ar_freqs.loc[id, target_id] = ar_freqs.loc[id, target_id] + 1
            else:
                ar_freqs.loc[id, 'other'] = ar_freqs.loc[id, 'other'] + 1
            
    ar_freqs = ar_freqs.apply(axis=1, func=lambda x: x / n_trials)
    
    if sum(ar_freqs['other']) == 0:
        ar_freqs = ar_freqs.drop('other', axis=1)
    
    return ar_freqs



def get_median_time_to_state(init_state, target_state, bool_funcs,
                            n_trials=100, max_steps_per_trial=100,
                            min_n_to_est_median=5,
                            synch_update=False):
    # TO DO: revisit output if state is never reached on a given trial.
    dat = pd.DataFrame({'TIME_TO_STATE' : np.full(n_trials, np.nan),
                        "STATE_REACHED" : np.full(n_trials, np.nan)}, 
                       index=np.arange(start=1, stop=n_trials+1))
    for i in dat.index:
        t, matched_target = simulate_to_target_states(
            init_state=init_state, bool_funcs=bool_funcs, target_states=target_state,
            n_steps=max_steps_per_trial, synch_update=synch_update)
        dat.loc[i, 'TIME_TO_STATE'] = t
        dat.loc[i, 'STATE_REACHED'] = matched_target
    
    median_t = np.nan
    match_selector = (dat['STATE_REACHED'] == target_state.name)
    if sum(match_selector) >= min_n_to_est_median:
        dat = dat.loc[match_selector, :]
        median_t = np.median(dat['TIME_TO_STATE'])
    
    return median_t



def get_state_table(node_states):
    nonconst_nodes = node_states[np.isnan(node_states)].index
    const_nodes = node_states.index.difference(nonconst_nodes)
    N = nonconst_nodes.size
    assert(N > 0)
    
    df = pd.DataFrame(np.full((2**N, node_states.size), np.nan), 
                              columns=node_states.index, dtype=np.int8)
    for node in const_nodes:
        df.loc[:, node] = np.int8(node_states[node])
    
    format_str = '0' + str(N) + 'b'
    
    for i in range(2**N):
        state = [np.int8(digit) for digit in format(i, format_str)]
        df.loc[i, nonconst_nodes] = state
        
    df.index = get_binary_states_as_integers(df)
    
    return df


def add_state_info_cols(df, attractors=None):
    df = df.copy()
    state_ids = df.apply(axis=1, func=binary_state_to_int)
    
    if attractors is not None:
        attr_ids = get_binary_states_as_integers(attractors).values
        is_attr = lambda x: 'ATTRACTOR' if binary_state_to_int(x) in attr_ids else ""
        state_types = df.apply(axis=1, func=is_attr)
        df['STATE_ID'] = state_ids
        df['STATE_TYPE'] = state_types
    else:
        df['STATE_ID'] = state_ids

    return df

# -------[(Check) Added from calbert notebook (END)]--------------------------------------


# -------[Processing Boolean design matrices (obs x vars)]--------------------------------
def get_non_na_bool_df(df, col1, col2):
    return df.loc[:, [col1, col2]].dropna(axis=0, how='any').astype('bool')



def get_pairwise_true_counts(table, variables):
    N = len(variables)
    counts_df = pd.DataFrame(np.full((N, N), 0), index=variables, columns=variables)
    for a1 in variables:
        for a2 in variables:
            if a1 == a2:
                counts_df.loc[a1, a2] = sum(table[a1].dropna().astype('bool'))
            else:
                df = get_non_na_bool_df(table, a1, a2)
                counts_df.loc[a1, a2] = sum(df[a1] & df[a2])
    
    return counts_df



def split_boolean_obs_table_by_inputs(table, input_vars, 
        set_null_in_table_inputs_to_none=True):
    assert(table.loc[:, input_vars].notnull().values.all())
    table = table.copy()
    input_id_table = table.loc[:, input_vars].copy().astype('bool').astype('int')
    assert(set(input_id_table.values.flatten()) == {0, 1})

    input_id_table = add_state_info_cols(input_id_table)
    assert('__ID__' not in table.columns)
    table['__ID__'] = input_id_table['STATE_ID']
    cols = table.columns.tolist()
    table = table[[cols[-1]] + cols[:-1]] 

    d = {k:v for k, v in table.groupby('__ID__')}
    input_ids = sorted(d.keys())
    odict = OrderedDict()
    for i in input_ids:
        odict[i] = d[i].drop(labels=['__ID__'], axis=1)
    
    if set_null_in_table_inputs_to_none:
        for i in odict.keys():
            if odict[i].isnull().values.any():
                odict[i] = None
            else:
                odict[i] = odict[i].astype('bool').astype('int')
                
    return odict


# -----------------------------------------------------------------------------------------
def get_boolnet_node_expression(X, node_to_genes):
    n_nodes = len(node_to_genes)
    n_samples = X.shape[1]
    nodeX = pd.DataFrame(np.full((n_nodes, n_samples), np.nan), 
                         index=node_to_genes.keys(), columns=X.columns)
    for node, genes in node_to_genes.items():
        genes = list(set(genes).intersection(X.index))
        if len(genes) == 0:
            continue
        elif len(genes) == 1:
            nodeX.loc[node, :] = X.loc[genes, :].values
        else:
            nodeX.loc[node, :] = np.mean(X.loc[genes, :], axis=0)
    
    return nodeX



def binarize_boolnet_node_expression(X, on_exp_pctl=50, dropna=True):
    binX = pd.DataFrame(np.full(X.shape, np.nan), index=X.index, columns=X.columns)
    for node in binX.index:
        X_node = X.loc[node, :]
        if any(np.isnan(X_node)):
            continue
        threshold = np.percentile(X_node, on_exp_pctl)
        binX.loc[node, :] = (X_node > threshold).astype(int)
    
    if dropna:
        binX = binX.dropna().astype(int)
    
    return binX


def read_node_to_gene_symbols_file(path):
    with open(path, 'r+') as f:
        tmp = [s.strip() for s in f.read().splitlines()]

    node_to_genes = OrderedDict()
    for s in tmp:
        node, genes = [part.strip() for part in s.split('=')]
        node_to_genes[node] = genes.split()
    
    return node_to_genes


def get_uniq_state_count_tab(df, input_nodes):
    df = df.copy()
    input_nodes = list(input_nodes)

    row_state_strings = df.apply(binary_state_to_str, axis=1)
    assert(row_state_strings.index is df.index)
    assert('__attr_str__' not in df.columns)
    df.insert(loc=0, column='__attr_str__', value=row_state_strings)

    row_input_ids = df.loc[:, input_nodes].apply(func=binary_state_to_int, axis=1)
    assert(row_input_ids.index is df.index)
    assert('__input_id__' not in df.columns)
    df.insert(loc=1, column='__input_id__', value=row_input_ids)

    row_state_counts = row_state_strings.value_counts()
    df.drop_duplicates(inplace=True)
    df.index = df.__attr_str__
    df.drop(labels='__attr_str__', axis=1, inplace=True)

    assert('__count__' not in df.columns)
    df.insert(loc=1, column='__count__', value=[row_state_counts[attr_str] for attr_str in df.index])
    
    return df


def get_simulation_based_input_attractor_tab(inputs, nodes, bool_funcs, 
    n_steps=10000, n_for_steady_state=100, n_init_conds=3, synch_update=False, 
    outfile_root=None, verbose=False):
    nodes = list(nodes)
    # Inputs can be fully specified in a data frame, or given as a list of nodes.
    # In the latter case, a data frame w/all possible Boolean combinations will 
    # be generated.
    if isinstance(inputs, pd.DataFrame):
        input_nodes = list(inputs.columns)
        input_tab = inputs
    else:
        input_nodes = list(inputs)
        input_tab = get_state_table(pd.Series(np.nan, index=input_nodes))
        
    assert(all([node in nodes  for node in input_nodes]))
    assert(all([node in bool_funcs  for node in nodes]))
    
    non_input_nodes = sorted(list(set(nodes) - set(input_nodes)))
    nodes = input_nodes + non_input_nodes
    
    outfile_path = None
    if outfile_root:
        tstamp = time.strftime("%Y%m%d_%H%M%S")
        ihash = abs(hash(str(input_tab) + str(bool_funcs) + 
            str(n_steps) + str(n_for_steady_state) + str(n_init_conds) + 
            str(synch_update) + str(outfile_root)))
        outfile_path = outfile_root + '_time_' + str(tstamp) + '_hash_' + str(ihash) + '.pkl'
    
    n_inputs = input_tab.shape[0]
    max_n_attr = n_inputs * n_init_conds
    full_attr_mat = np.full((max_n_attr, len(nodes)), np.nan)
    attr_labels = [None] * max_n_attr
    k = 0
        
    for i in range(n_inputs):
        start_state = pd.Series(np.nan, index=nodes)
        start_state[input_nodes] = input_tab.iloc[i, :]

        if verbose:
            print(str(i+1) + " ", end='')
            if (i+1 % 40) == 0: print("\n", end='')

        for j in range(n_init_conds):
            start_state[non_input_nodes] = np.random.randint(low=0, high=2, 
                size=len(non_input_nodes))
            traj = get_boolnet_trajectory(init_state=start_state, inputs=start_state[input_nodes], 
                bool_funcs=bool_funcs, n_steps=n_steps, synch_update=synch_update)
            attr = traj.tail(n=n_for_steady_state).drop_duplicates()
            if attr.shape[0] == 1:
                attr_labels[k] = 'input_' + str(i+1) + '_ic_' + str(j+1)
                full_attr_mat[k, :] = attr
                k += 1
         
        if outfile_path and ((i+1 % 100) == 0):
            full_attr_df = pd.DataFrame(full_attr_mat[:k, :], index=attr_labels[:k], columns=nodes)
            input_attr_tab = get_uniq_state_count_tab(full_attr_df, input_nodes)
            #input_attr_tab.to_hdf(path_or_buf=outfile_path, key='input_attr_tab', mode='w')
            input_attr_tab.to_pickle(path=outfile_path)
    
    full_attr_df = pd.DataFrame(full_attr_mat[:k, :], index=attr_labels[:k], columns=nodes)
    input_attr_tab = get_uniq_state_count_tab(full_attr_df, input_nodes)
    
    return input_attr_tab


def apply_func_to_df_split_in_parallel(df, func, n_cores=None, **keywords):
    if (n_cores is None) or (n_cores < 2):
        n_cores = os.cpu_count()
        
    if keywords:
        # partial() returns a single (DataFrame) argument function, 
        # setting all other arguments according to keyword:argument 
        # pairs in keywords dictionary.
        func = partial(func, **keywords)
    
    df_split = np.array_split(df, n_cores)
    #pool = mp.Pool(n_cores)
    # Using: from pathos.multiprocessing import ProcessingPool as Pool
    pool = Pool(n_cores)
    output_df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    
    return output_df


def is_attractor_state(state, inputs, bool_funcs, n_steps=1000, synch_update=False):
    traj = get_boolnet_trajectory(init_state=state, inputs=inputs, bool_funcs=bool_funcs, 
        n_steps=n_steps, synch_update=False)
    
    return traj.drop_duplicates().shape[0] == 1