import random
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from random import sample
from networkx.algorithms.cycles import find_cycle
from networkx.drawing.nx_agraph import write_dot, graphviz_layout
from admixture_graph import *
from scipy.special import comb
from scipy.stats import beta

def plot_graph(g, title):
    plt.title(title)
    pos = graphviz_layout(g, prog='dot')
    nx.draw(g, pos, with_labels=True, arrows=True)
    plt.show()

def draw_random_branch(graph):
    """ draw a random branch """
    return sample(graph.edges, 1)[0]

def find_target_subgraph(graph, exclude_node):
    """ return the subgraph that does not contain exclude_node """
    graph_unadmixed = graph.copy()
    for e in list(graph_unadmixed.edges):
        if graph_unadmixed.is_admixture_edge(e):
            graph_unadmixed.remove_edge(*e)
    graphs = graph_unadmixed.to_undirected()
    graphs = list(nx.connected_component_subgraphs(graphs))
    if exclude_node in graphs[0].nodes:
        graph = graphs[1]
    elif exclude_node in graphs[1].nodes:
        graph = graphs[0]
    else:
        raiseValueError('The specified node not in either graph!')
    return graph

def get_pruned_node_candidate(admixture_graph):
    """ Return a list of nodes to be drawn from for pruning """
    ag = admixture_graph.copy()
    nodes = ag.nodes()
    root = [n for n in nodes if ag.in_degree(n) == 0]
    children_root =  [n for n in nodes if ag.parent_merge_node(n) in root]
    admixture_nodes = [n for n in nodes if ag.nodes[n]['type'] == 'admixture']
    exclusion = root + children_root + admixture_nodes
    candidate_nodes = [n for n in nodes if n not in exclusion]
    return candidate_nodes

# SPR
def SPR(admixture_graph, force_draw_node = None, print_output = False, locate_cycle = True):
    """ Subtree pruning and regrafting
    A node is drawn and the non-admixed parent of this node is detached, and
    The broken edge is reconnected, and finally the detached node is inserted
    to a branch. 
    """
    # random.seed(101)

    # TO DO: assert t_event is in the node attribute
    ag = admixture_graph.copy()
    current_attributes = dict(ag.nodes(data=True))
    original_edges = ag.edges
    if print_output:
        plot_graph(ag, 'original graph')
        print(f'events: {ag.get_events()}')
    # detach a subtree and reconnect the broken part
    candidate_pruned_nodes = get_pruned_node_candidate(ag)
    if print_output:
        print(f'candidate_pruned_nodes: {candidate_pruned_nodes}')
    if force_draw_node:
        pruned_node = force_draw_node
    else:
        pruned_node = sample(candidate_pruned_nodes, 1)[0]
    parent_pruned = ag.parent_merge_node(pruned_node)
    if print_output:
        print(f'pruned: {pruned_node}; its parent {parent_pruned}')
    
    reconnect_node_parent = next(ag.predecessors(parent_pruned))
    reconnect_node_child = ag.successors(parent_pruned)
    reconnect_node_child = [r for r in reconnect_node_child if r != pruned_node][0]
    lower_bound_node = ag.lower_bound_node(pruned_node) # for lower bound
    lower_bound_time = ag.get_event_time(lower_bound_node)
    if print_output:    
        print(f'lower bound node : {lower_bound_node}')

    t0 = ag.get_event_time(parent_pruned) - ag.get_event_time(lower_bound_node)
    ag.remove_node(parent_pruned)
    ag.add_edge(reconnect_node_parent, reconnect_node_child)
    if print_output:
        print(f'reconnet edge({reconnect_node_parent},{reconnect_node_child})')
    reconnected_nodes = (reconnect_node_parent, reconnect_node_child)

    if print_output:
        plot_graph(ag, 'detached graph')  
    # find a brach to insert the detached subtree
    # 1) pruned subtree must be younger than the inserted branch 
    # 2) must draw a different branch
    # 3) no cycle formed
    subgraph = find_target_subgraph(ag, pruned_node)
    cycle = True
    while cycle:
        ag_holder = ag.copy() 
        branch_to_attach = draw_random_branch(subgraph)
        branch_to_attach = ('RT', 'D')
        condition_event_time = True
        condition_same_branch = True
        while condition_event_time | condition_same_branch:
            branch_to_attach = draw_random_branch(subgraph)
            if branch_to_attach not in original_edges:
                branch_to_attach = branch_to_attach[::-1]
            condition_event_time = (ag_holder.get_event_time(branch_to_attach[0]) < lower_bound_time)
            # print(f'condition_event_time:{condition_event_time}')
            condition_same_branch = (branch_to_attach == reconnected_nodes)
            # print(f'condition_same_branch: {condition_same_branch}')
        # print(f'branch to attach: {branch_to_attach}')
        ag_holder.remove_edge(*branch_to_attach)
        ag_holder.add_edge(branch_to_attach[0],parent_pruned)
        ag_holder.add_edge(parent_pruned, branch_to_attach[1])
        # debug
        ag_holder.add_edge(parent_pruned, lower_bound_node)
        try:
            c = find_cycle(ag_holder)
            if locate_cycle:
                print(f'cycle:{c}; prune node {pruned_node}; branch to attach: {branch_to_attach}')
                #plot_graph(ag_holder, 'cycle')  
            cycle = True
            # return ag_holder
        except nx.exception.NetworkXNoCycle:
            cycle = False
        print(cycle)
    ag = ag_holder
    t1 = ag.get_event_time(branch_to_attach[0]) - ag.get_event_time(branch_to_attach[1]) 
    if print_output:
        print(f'add edge: {(branch_to_attach[0],parent_pruned)},{(parent_pruned, branch_to_attach[1])}, {(parent_pruned, lower_bound_node)}')
    nx.set_node_attributes(ag, {parent_pruned: current_attributes[parent_pruned]})
    
    # set random event time
    lower = max(ag.get_event_time(branch_to_attach[1]), ag.get_event_time(lower_bound_node))
    upper = ag.get_event_time(branch_to_attach[0])
    t = np.random.uniform(lower, upper)

    ag.set_event_time({parent_pruned: t})
    if print_output:
        print(f'lower_upper = ({lower},{upper})')
        print(f'ag.set_event_time({parent_pruned}:{t})')
        print(f'detach node {pruned_node}; regraft to branch {branch_to_attach}')
        plot_graph(ag, 'graph after SPR')  
    return ag, t0/t1


# admixture proposal
def admixture_edge_proposal(admixture_graph, force_draw_node = None, print_output = False):
    """
    Note an admixture edge "A --> B". 
    A has 1 parent and 2 children, and B has 2 parents and 1 child.
    Selecting every possible pair of edges, possible in this context means
    there is overlap in their times.
    --------------------------------------------------------------------
    return proposed graph, q(x | x') / q(x' | x)
    """
    # random.seed(101)
    ag = admixture_graph.copy()
    node_types = ag.get_event_type()
    admixture_edges_dict = ag.get_admixture_proportions()
    if print_output:
        plot_graph(ag, 'before')

    for e in admixture_edges_dict:
        t0 = ag.get_event_time(e[0])
        num_branch_t0 = len(ag.get_branches_at_time(t0))
        reconnect_source_branch = next(ag.predecessors(e[0])), next(ag.successors(e[0]))
        reconnect_target_branch = next(ag.predecessors(e[1])), next(ag.successors(e[1]))
        # remove two nodes and their edges
        for n in e:
            ag.remove_node(n)
        ag.add_edge(*reconnect_source_branch)
        ag.add_edge(*reconnect_target_branch)
        # attach branch    
        t_upper = ag.get_event_time(ag.root())
        t_lower = 0
        t1 = np.random.uniform(t_lower, t_upper)
        candidate_branches = ag.get_branches_at_time(t1)
        
        num_branch_t1 = len(candidate_branches)
        new_source_branch, new_target_branch = sample(candidate_branches, 2)
        if print_output:
            print(f'new source {new_source_branch}; new target {new_target_branch}')
        ag.add_edge(new_source_branch[0], e[0])
        ag.add_edge(e[0], new_source_branch[1])
        ag.add_edge(new_target_branch[0], e[1])
        ag.add_edge(e[1], new_target_branch[1])
        ag.add_edge(*e)
        ag.remove_edge(*new_source_branch)
        ag.remove_edge(*new_target_branch)
        ag.set_event_time({e[0]: t1, e[1]: t1})
        ag.set_admixture_proportion({e: admixture_edges_dict[e]})
    ag.set_event_type(node_types)
    # proposal q(x | x') / q(x' | x)
    q = (np.log(comb(num_branch_t0, 2, exact=True)) - np.log(comb(num_branch_t1, 2, exact=True)))
    if print_output:
        plot_graph(ag, 'after')

    return ag, q


def proposal_rejection(trace):
    trace['time'].append(trace['time'][-1])
    trace['loglik'].append(trace['time'][-1])
    trace['admixture'].append(trace['admixture'][-1])
    trace['topology'].append(trace['topology'][-1])

def initialize_mcmc( admixture_graph, sfs, rep = 100, log_acceptance = True,):
    pass


from collections import OrderedDict

def reject_proposal(trace):
    trace['theta'].append(trace['theta'][-1])
    trace['loglik'].append(trace['loglik'][-1])
    trace['topology'].append(trace['topology'][-1])
    trace['acceptance'].append('reject')

def accept_proposal(trace, theta, loglik, topology):
    trace['theta'].append(theta)
    trace['loglik'].append(loglik)
    trace['topology'].append(topology)
    trace['acceptance'].append('accept')

def event_time_proposal_randomwalk(graph):
    theta_current = OrderedDict(graph.get_events())
    admixture_edges = graph.get_admixture_edges()
    admixture_events = graph.get_admixture_events()
    admixture_nodes = [e[0] for e in admixture_events]
    proposal_sigmas = [(t[1]/20 + 1000) for t in graph.get_events()]

    times = dict()
    for i, e in enumerate(theta_current):
        mu = theta_current[e]
        times[e] = max(np.random.normal(mu, proposal_sigmas[i], 1)[0], 50)
    times_graph = times
    times_demo = times.copy()
    for e in admixture_edges:
        times_demo[f'{e[0]}_{e[1]}'] = times_demo[e[0]]
        times_graph[e[1]] = times_graph[e[0]]
        del times_demo[e[0]]
        del times_demo[e[1]]
    return times_graph, times_demo


def event_time_proposal_multiplicative(graph, rate = 10):
    theta_current = OrderedDict(graph.get_events())
    admixture_edges = graph.get_admixture_edges()
    admixture_events = graph.get_admixture_events()
    admixture_nodes = [e[0] for e in admixture_events]
    
    scale = np.random.exponential(rate)
    shrink = np.random.binomial(1, 0.5)
    if shrink:
        scale = 1/scale

    times = dict()
    for i, e in enumerate(theta_current):
        mu = theta_current[e]
        times[e] = mu * scale
    times_graph = times
    times_demo = times.copy()
    for e in admixture_edges:
        times_demo[f'{e[0]}_{e[1]}'] = times_demo[e[0]]
        times_graph[e[1]] = times_graph[e[0]]
        del times_demo[e[0]]
        del times_demo[e[1]]
    return times_graph, times_demo, scale


def admixture_proportion_proposal(graph):
    '''
    return proposed admixture fraction 
    and "log q(x'|x) - log q(x|x')"
    '''
    # beta proposal
    admixture_proportions = graph.get_admixture_proportions()
    proportions = dict()
    qs = 0
    for i, e in enumerate(admixture_proportions):
        mu = admixture_proportions[e]
        admixture_edge_name = f'{e[0]}_{e[1]}_proportion'
        a = mu * 30
        b = (1 - mu) * 30
        p = np.random.beta(a, b, 1)[0]
        proportions[admixture_edge_name] = p
        q_forward = beta.logpdf(p, a, b)
        a = p * 30
        b = (1 - p) * 30
        q_backward = beta.logpdf(mu, a, b)
        qs += q_backward - q_forward
    return proportions, qs


def update_theta(model, theta):
    model.set_params(theta)