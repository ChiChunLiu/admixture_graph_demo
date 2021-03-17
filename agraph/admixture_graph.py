import networkx as nx
from random import sample
import momi

class admixture_graph(nx.DiGraph):
    """
    Admixture graph class derived from the DiGraph class.
    A pulse admixture event can be represented by two nodes
    with the direction from node a (source) to node b (recepient).
    A split event can be represented by a parent with two children.
    ----------------------------------------------------------------
    Details:
    1. nodes:
    All nodes have attributes "t_event" and "type", which record the
    event time and event type. Event type can be either "merge" or
    "admixture".

    2. edges:
    each admixture edge has an attribute "admixture_proportion", 
    representing the fraction of contribution from the source 
    population.
    """
    def __init__(self):
        super().__init__(self)

    def __is_inner(self, x):
        """ check if a node x is an inner node """
        self.out_degree(x)==2 and self.in_degree(x)==1

    def is_admixture_edge(self, x):
        """ Check if an edge is an admixture edge """
        ap = nx.get_edge_attributes(self, 'admixture_proportion')
        if x in ap:
            if ap[x] > 0:
                return True
        else:
            return False

    def inner_nodes(self):
        """ return a list of inner nodes (no root nodes) """
        nodes = self.nodes()
        inner = [x for x in nodes if self.__is_inner(x)]
        return inner

    def random_inner_node(self):
        """ draw a random innder node """
        nodes = self.nodes()
        inner = [x for x in nodes if self.__is_inner(x)]
        return sample(inner, 1)[0]

    def random_nonroot_node(self):
        """ draw a random nonroot node """
        nodes = self.nodes()
        inner = [x for x in nodes if self.in_degree(x) > 0]
        return sample(inner, 1)[0]

    def random_branch(self):
        """ draw a random branch """
        return sample(self.edges, 1)[0]
    
    def root(self):
        nodes = self.nodes()
        root = [x for x in nodes if self.in_degree(x) == 0][0]
        return root


    def get_events(self):
        """ Retrieve a list of tuples with elements
        (node, event_time) 
        """
        events = []
        for n, t in list(self.nodes(data=True)):
            if t['t_event'] > 0:
                events.append((n, t['t_event']))
        return sorted(events, key=lambda t: t[1])

    def get_leaf_nodes(self):
        """ retrieve a list of leave nodes """
        nodes = self.nodes()
        leaves = [x for x in nodes if self.out_degree(x) == 0]
        return leaves 

    def __set_leaves_attribute_at_leaves(self):
        """ set the leaf attributes as themselves """
        leaves = self.get_leaf_nodes()
        leaves_dict = dict()
        for i, n in enumerate(leaves):
            leaves_dict[n] = leaves[i]
        nx.set_node_attributes(self, leaves_dict, 'leaves')

    def set_leaves_attribute(self):
        """ recursively find the leaves under a node """
        self.__set_leaves_attribute_at_leaves()
        events = self.get_events()
        for e in events:
            leaves_dict = nx.get_node_attributes(self, 'leaves')
            children = list(self.successors(e[0])) # get left and right child
            if len(children) == 1:
                child_leaves = leaves_dict[children[0]] # the only child
            elif self.nodes[e[0]]['type'] == 'admixture':
                # set leaves as non-admixed child
                if self.nodes[children[0]]['type'] == 'admixture':
                    child_leaves = leaves_dict[children[1]]
                else:
                    child_leaves = leaves_dict[children[0]]
            else:
                child_leaves = leaves_dict[children[1]] # leaves of the right child
            nx.set_node_attributes(self, {e[0]: child_leaves}, 'leaves')

    def set_event_time(self, event_dict):
        """ Set user specified event times """
        nx.set_node_attributes(self, event_dict, 't_event')

    def set_event_type(self, type_dict):
        """ Set user specified event types"""
        nx.set_node_attributes(self, type_dict, 'type')
    
    def get_event_type(self):
        return nx.get_node_attributes(self, 'type')
    
    def set_admixture_proportion(self, admixture_dict):
        """ Set admixture proportion as an edge attribute """
        nx.set_edge_attributes(self, admixture_dict, 'admixture_proportion')

    def get_admixture_events(self):
        """ retrieve a list of tuples (node, event_time) """
        events = []
        for n, t in list(self.nodes(data=True)):
            if t['type'] == 'admixture':
                events.append((n, t['t_event']))
        return sorted(events, key=lambda t: t[1])
    
    def get_admixture_proportions(self):
        """ Return a dictionary with edges as keys and admixture proportion
        as value 
        """
        admixture_proportions = nx.get_edge_attributes(self, 'admixture_proportion')
        admixture_edges = dict()
        for e in admixture_proportions:
            if admixture_proportions[e] > 0:
                admixture_edges[e] = admixture_proportions[e]
        return admixture_edges
    
    def get_admixture_edges(self):
        return list(self.get_admixture_proportions().keys())

    def draw_random_branch(graph):
        return sample(graph.edges, 1)[0]


    def get_event_time(self, node):
        """ Return event time of a node """
        return self.nodes[node]['t_event']

    def is_event_order_feasible(self):
        """ Check if the current event time configuration
        conforms to the parent-children relationships
        """
        events = self.get_events()
        for e in events:
            e_time = e[1]
            children = self.successors(e[0])
            for c in children:
                if e_time < self.get_event_time(c):
                    return False
        return True

    def to_demography(self, print_events = False):
        """ convert a networkx graph into a momi demography object """
        model = momi.DemographicModel(N_e=1e4, gen_time=25, 
                                      muts_per_gen=1.25e-8)
        leaves = self.get_leaf_nodes()
        for l in leaves:
            model.add_leaf(l, t=0)

        self.set_leaves_attribute()
        events = self.get_events()
        leaves_dict = nx.get_node_attributes(self, 'leaves')
        admixture_proportions = self.get_admixture_proportions()
        for e in events:
            n, t = e[0], e[1]
            if self.nodes[n]['type'] == 'merge':
                model.add_time_param(n)
                children = list(self.successors(n))
                model.move_lineages(leaves_dict[children[0]], leaves_dict[children[1]], t = n)
                model.set_params({n: t})
                if print_events:
                    print(f'move from {leaves_dict[children[0]]} to ' +
                          f'{leaves_dict[children[1]]} at t = {t:.2f}')
            elif self.nodes[n]['type'] == 'admixture':
                pass
                #model.set_params({n: t})
            else:
                raise ValueError('event can only be admixture or merge')

        for e in admixture_proportions:
            p = admixture_proportions[e]
            t = self.get_event_time(e[0])
            admixed_edge_name = f'{e[0]}_{e[1]}'
            model.add_time_param(admixed_edge_name)
            model.add_pulse_param(f'{admixed_edge_name}_proportion')
            model.move_lineages(leaves_dict[e[1]], 
               leaves_dict[e[0]], t = admixed_edge_name, p = f'{admixed_edge_name}_proportion')
            model.set_params({admixed_edge_name: t, f'{admixed_edge_name}_proportion': p})
            if print_events:
                print(f'move from {leaves_dict[e[1]]} to ' +
                      f'{leaves_dict[e[0]]} at t = {t:.2f}' +
                      f' and proportion = {p}')
        return model

    def parent_merge_node(self, node):
        """ Get the node with a deeper event time that is not an admixture node.
        This traces the non-admixture lineages recursively.
        """
        parents = list(self.predecessors(node))
        if len(parents) == 1:
            if self.nodes[parents[0]]['type'] == 'admixture':
                return self.parent_merge_node(parents[0])
            else:
                return parents[0]
        elif len(parents) == 2:
            # follow the one not an admixture edge
            if self.is_admixture_edge((parents[0], node)):
                return self.parent_merge_node(parents[1])
            elif self.nodes[parents[0]]['type'] == 'admixture':
                return self.parent_merge_node(parents[0])
            else:
                return parents[0]

    def lower_bound_node(self, node):
        parents = list(self.predecessors(node))
        if len(parents) == 1:
            if self.nodes[parents[0]]['type'] == 'admixture':
                return parents[0]
            else:
                return node
        elif len(parents) == 2:
            # follow the one not an admixture edge
            if self.is_admixture_edge((parents[0], node)):
                return self.lower_bound_node(parents[1])
            elif self.nodes[parents[0]]['type'] == 'admixture':
                return self.lower_bound_node(parents[0])
            else:
                return parents[0]
    
    def get_branches_at_time(self, t):
        #edges = [e for e in self.edges if e not in self.get_admixture_edges()]
        g = self.to_unadmixed_tree()
        return [e for e in g.edges if (self.get_event_time(e[0]) > t) and (self.get_event_time(e[1]) < t)]


    def to_unadmixed_tree(self):
        ag = self.copy()
        admixture_edges_dict = ag.get_admixture_proportions()
        for e in admixture_edges_dict:
                reconnect_source_branch = next(ag.predecessors(e[0])), next(ag.successors(e[0]))
                reconnect_target_branch = next(ag.predecessors(e[1])), next(ag.successors(e[1]))
                # remove two nodes and their edges
                for n in e:
                    ag.remove_node(n)
                ag.add_edge(*reconnect_source_branch)
                ag.add_edge(*reconnect_target_branch)
        return ag

    def __tree_to_newick(self, root):
        """
        modified from stackoverflow:
        46444454/save-networkx-tree-in-newick-format
        """
        if len(self[root]) == 0:
            return '(' + root + ')'

        subgs = []
        for child in self[root]:
            if len(self[child]) > 0:
                subgs.append(self.__tree_to_newick(root=child))
            else:
                subgs.append(child)
        return "(" + ','.join(subgs) + ")"

    def to_newick(self):
        g = self.to_unadmixed_tree()
        newick = {'tree': '()', 'admixture':[]}
        # unadmixed tree
        root = self.root()
        newick['tree'] = g.__tree_to_newick(root = root)
        # admixture subtree
        admixture_edges_dict = self.get_admixture_proportions()
        for e in admixture_edges_dict:
            for admixture_node in e:
                child = next(self.successors(admixture_node))
                #print(child)
                newick['admixture'].append(g.__tree_to_newick(root = child))
        return newick

# functions for generating example graphs
def generate_example_graph():
    topology = admixture_graph()
    topology.add_edge('RT', 'G')
    topology.add_edge('RT', '2')
    topology.add_edge('2','3')
    topology.add_edge('2','4')
    topology.add_edge('3','A')
    topology.add_edge('3','5')
    topology.add_edge('5','B')
    topology.add_edge('4','7')
    topology.add_edge('4','C')
    topology.add_edge('7','8')
    topology.add_edge('8','D')
    topology.add_edge('8','10')
    topology.add_edge('10','E')
    topology.add_edge('10','F')
    topology.add_edge('7','5') # admixture  edge
    return topology

def generate_easy_example_graph():
    topology = admixture_graph()
    topology.add_edge('RT', '1')
    topology.add_edge('RT', 'D')
    topology.add_edge('1','2')
    topology.add_edge('1','3')
    topology.add_edge('2','A')
    topology.add_edge('2','4')
    topology.add_edge('3','C')
    topology.add_edge('4','B')
    topology.add_edge('3','4') # admixture  edge
    return topology

