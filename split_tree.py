################################################################################
#
# split_tree.py
#
# Dependencies: networkx, matplotlib.pyplot
#
#
# This module can be used to:
#
#   * Compute the original graph corresponding to a given distance-hereditary or
#     three-leaf power split tree string (typically this string will have been
#     computed by our distance-hereditary or three-leaf power samplers, which
#     are implemented in Maple)
#
#   * Draw the original graph to the screen or export a drawing of the graph
#     to a file
#
#
# Example:
#
#   >>> from split_tree import *
#   >>> split_tree = string_to_split_tree("SR(K(SC(Z, SX(K(SC(Z, Z), Z), \
#                                         K(Z, Z, Z))), Z), Z, Z)")
#   >>> graph = split_tree_to_graph(split_tree)
#   >>> draw_graph(graph)
#   >>> export_graph(graph, "drawing.pdf")
#   >>> export_graph(graph, "drawing.png")
#   >>> export_graph(graph, "drawing") # If no format specified, defaults to png
#
# For the semantics of the input to string_to_split_tree, see the README
#
################################################################################

import matplotlib.pyplot as plt
import networkx as nx

# Every node receives a distinct id value
current_id = 0


class Node(object):
    """
    Represent a node in a split tree
    """
    def __init__(self):
        global current_id
        self.id = current_id
        current_id += 1
        self.visited = False


class StarNode(Node):
    """
    Extend Node; represent a star node in a split tree

    Attributes:
    center -- Node connected to the center of the star
    extremities -- list of Nodes connected to the extremities of the star
    """
    def __init__(self, center=None, extremities=None):
        super(StarNode, self).__init__()
        self.center = center
        if extremities is None:
            extremities = []
        self.extremities = extremities


class CliqueNode(Node):
    """
    Extend Node; represent a clique node in a split tree

    Attributes:
    extremities -- list of Nodes connected to the vertices of the star
    """
    def __init__(self, extremities=None):
        super(CliqueNode, self).__init__()
        if extremities is None:
            extremities = []
        self.extremities = extremities


class Leaf(Node):
    """
    Extend Node; represent a leaf node in a split tree

    Attributes:
    neighbor -- Node connected to this leaf (or None if the tree has just this
                one leaf)
    """
    def __init__(self, neighbor=None):
        super(Leaf, self).__init__()
        self.neighbor = neighbor


def draw_graph(graph):
    """
    Draw graph on the screen

    Arguments:
    graph -- list of adjacencies representing the graph (typically obtained from
             split_tree_to_graph)
    """
    G = graph_to_pyplot(graph)
    plt.show()
    return G


def export_graph(graph, filename):
    """
    Export a pdf drawing of graph to filename

    Arguments:
    graph -- list of adjacencies representing the graph (typically obtained from
             split_tree_to_graph)
    filename -- file path string (may provide a format extension from the list:
                bmp, eps, gif, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg,
                svgz, tif, tiff -- otherwise, png is used by default)
    """
    G = graph_to_pyplot(graph)
    plt.savefig(filename, bbox_inches="tight", transparent=True, pad_inches=0)
    plt.close()
    return G


def graph_to_pyplot(graph):
    """
    Draw graph networkx to pyplot

    Arguments:
    graph -- list of adjacencies representing the graph (typically obtained from
             split_tree_to_graph)
    """
    G = graph_to_networkx_G(graph)

    plt.close()
    pos=nx.spring_layout(G, k=0.3, iterations=50)
    nx.draw(G, pos=pos, node_size=20, node_color="k", edge_color="k",
            node_shape="o", scale=200, alpha=0.5)
    return G


def graph_to_networkx_G(graph):
    G = nx.Graph()
    if not graph:
        G.add_node(1)
    for (id1, id2) in graph:
        G.add_edge(id1, id2)
    return G


def split_tree_to_graph(tree):
    """
    Convert tree, a split tree, into its original graph (represented as a list
    of adjacencies on the nodes 1, ..., n)

    Arguments:
    tree -- Node object representing the split tree to convert
            (typically obtained from string_to_split_tree)
    """
    G = nx.Graph()

    if is_leaf(tree) and tree.neighbor is None:
        return []

    leafs = get_leafs(tree)
    adjacencies = set()
    for leaf in leafs:
        adjacent_leafs = find_adjacent_leafs(leaf)
        for adjacent_leaf in adjacent_leafs:
            if (leaf.id != adjacent_leaf.id and
                (adjacent_leaf.id, leaf.id) not in adjacencies):
                    adjacencies.add((leaf.id, adjacent_leaf.id))

    adjacencies_list = list(adjacencies)
    min_id = min(list(sum(adjacencies_list, ())))
    return map(lambda (x, y) : (x - min_id + 1, y - min_id + 1),
               adjacencies_list)


def string_to_split_tree(s, parent_node=None):
    """
    Convert the string s to a split tree, where parent_node is the previously
    visited node (or None if this is the first node we have visited)

    Syntax for s:
    Z -- a leaf node
    KR -- a clique node root (can only be used as the root of the tree)
    SR -- a star node root (can only be used as the root of the tree)
    K -- a clique that has been entered from another node
    SX -- a star that has been entered from another node at one of its
          extremities
    SC -- a star that has been entered from another node at its center
    e -- an edge that connects two nodes

    For example:

    KR(Z, Z, Z, Z) -- a 4-clique whose children are all leaves

    Z(K(Z, Z, Z)) -- a leaf attached to a 4-clique whose other three children
                     are leaves (equivalent to the previous exampmle)

    SR(Z, SX(Z, Z), K(SC(Z, Z), Z)) -- a star whose children are a leaf, the
                                       extremity of a star whose center and
                                       other extremity are leaves, and a
                                       3-clique whose other children are a leaf
                                       and the center of a star whose two
                                       extremities are leaves

    e(SC(Z, Z), SC(Z, Z)) -- two stars, each of which has two leaves as its
                             extremities, connected at their centers by an edge

    Arguments:
    s -- string following the syntax described above
    parent_node -- Node object
    """
    s = "".join(s.split())
    if s == "Z":
        return Leaf(parent_node)
    else:
        pieces = s.split("(", 1);
        root = pieces[0];
        neighbors = pieces[1][:-1]
        if root == "Z":
            l = Leaf()
            l.neighbor = string_to_split_tree(neighbors, l)
            return l
        elif root == "KR":
            c = CliqueNode()
            for subtree in custom_split(neighbors):
                c.extremities.append(string_to_split_tree(subtree, c))
            return c
        elif root == "SR":
            s = StarNode()
            first = True
            for subtree in custom_split(neighbors):
                if first:
                    s.center = string_to_split_tree(subtree, s)
                    first = False
                else:
                    s.extremities.append(string_to_split_tree(subtree, s))
            return s
        elif root == "e":
            subtrees = custom_split(neighbors)
            subtree0pieces = subtrees[0].split("(", 1)
            root0 = subtree0pieces[0]
            neighbors0 = subtree0pieces[1][:-1]
            s0 = StarNode()
            if root0 == "SX":
                s0.extremities.append(string_to_split_tree(subtrees[1], s0))
                first = True
                for subtree in custom_split(neighbors0):
                    if first:
                        s0.center = string_to_split_tree(subtree, s0)
                        first = False
                    else:
                        s0.extremities.append(string_to_split_tree(subtree, s0))
                return s0
            elif root0 == "SC":
                s0.center = string_to_split_tree(subtrees[1], s0)
                for subtree in custom_split(neighbors0):
                    s0.extremities.append(string_to_split_tree(subtree, s0))
                return s0
            else:
                raise Exception("???")
        elif root == "K":
            c = CliqueNode()
            c.extremities.append(parent_node)
            for subtree in custom_split(neighbors):
                c.extremities.append(string_to_split_tree(subtree, c))
            return c
        elif root == "SC":
            s = StarNode()
            s.center = parent_node
            for subtree in custom_split(neighbors):
                s.extremities.append(string_to_split_tree(subtree, s))
            return s
        elif root == "SX":
            s = StarNode()
            s.extremities.append(parent_node)
            first = True
            for subtree in custom_split(neighbors):
                if first:
                    s.center = string_to_split_tree(subtree, s)
                    first = False
                else:
                    s.extremities.append(string_to_split_tree(subtree, s))
            return s
        else:
            raise Exception("???")


def custom_split(s):
    """
    Split s by the character ',', but only for ',' characters
    that are not inside any pair of parentheses

    Arguments:
    s -- string
    """
    parts = []
    current = []
    level = 0
    for c in (s + ","):
        if c == "," and level == 0:
            parts.append("".join(current))
            current = []
        else:
            if c == "(":
                level += 1
            elif c == ")":
                level -= 1
            current.append(c)
    return parts


def is_of_type(node, type_string):
    """
    Determine whether or not node is of type type_string

    Arguments:
    node -- Node object
    type_string -- "Leaf", "CliqueNode", or "StarNode"
    """
    return type(node).__name__ == type_string


is_leaf = lambda node: is_of_type(node, "Leaf")

is_clique = lambda node: is_of_type(node, "CliqueNode")

is_star = lambda node: is_of_type(node, "StarNode")


def find_adjacent_leafs(leaf):
    """
    Return a list of all Leafs that are accessible from leaf

    l' is accessible from l if there is a path through the graph
    from l to l' that uses at most one edge from each internal node label

    Arguments:
    leaf -- Leaf object
    """
    nodes = get_nodes(leaf)

    leafs = []
    to_visit = [leaf]

    while len(to_visit) > 0:
        node = to_visit[0]
        del to_visit[0]
        node.visited = True

        if is_leaf(node):
            leafs.append(node)
            if not node.neighbor.visited:
                to_visit.insert(0, node.neighbor)
        elif is_clique(node):
            for extremity in node.extremities:
                if not extremity.visited:
                    to_visit.insert(0, extremity)
        elif is_star(node):
            if node.center.visited:
                for extremity in node.extremities:
                    if not extremity.visited:
                        to_visit.insert(0, extremity)
            else:
                to_visit.insert(0, node.center)

    for node in nodes:
        node.visited = False

    return leafs


def get_leafs(tree):
    """
    Return a list of all Leafs in tree

    Arguments:
    tree -- Node object
    """
    return filter(is_leaf, get_nodes(tree))


def get_nodes(tree):
    """
    Return a list of all Nodes in tree

    Arguments:
    tree -- Node object
    """
    nodes = []
    to_visit = [tree]

    while len(to_visit) > 0:
        node = to_visit[0]
        del to_visit[0]
        node.visited = True

        nodes.append(node)

        if is_leaf(node):
            if not node.neighbor.visited:
                to_visit.insert(0, node.neighbor)
        elif is_clique(node):
            for extremity in node.extremities:
                if not extremity.visited:
                    to_visit.insert(0, extremity)
        elif is_star(node):
            if not node.center.visited:
                to_visit.insert(0, node.center)
            for extremity in node.extremities:
                if not extremity.visited:
                    to_visit.insert(0, extremity)
        else:
            raise Exception("???")

    for node in nodes:
        node.visited = False

    return nodes


### Tests ###

def test1():
    G1 = graph_to_networkx_G(split_tree_to_graph(string_to_split_tree("Z")))
    G2 = nx.Graph()
    G2.add_node(0)
    assert nx.is_isomorphic(G1, G2)


def test2():
    G1 = graph_to_networkx_G(split_tree_to_graph(string_to_split_tree(
        "Z(K(Z, Z))")))
    G2 = nx.Graph()
    G2.add_edges_from([(0, 1), (1, 2), (2, 0)])
    assert nx.is_isomorphic(G1, G2)


def test3():
    G1 = graph_to_networkx_G(split_tree_to_graph(string_to_split_tree(
        "KR(Z, Z, Z)")))
    G2 = nx.Graph()
    G2.add_edges_from([(0, 1), (1, 2), (2, 0)])
    assert nx.is_isomorphic(G1, G2)


def test4():
    G1 = graph_to_networkx_G(split_tree_to_graph(string_to_split_tree(
        "Z(SX(Z, Z))")))
    G2 = nx.Graph()
    G2.add_edges_from([(0, 1), (1, 2)])
    assert nx.is_isomorphic(G1, G2)


def test5():
    G1 = graph_to_networkx_G(split_tree_to_graph(string_to_split_tree(
        "Z(SC(Z, Z))")))
    G2 = nx.Graph()
    G2.add_edges_from([(0, 1), (1, 2)])
    assert nx.is_isomorphic(G1, G2)


def test6():
    G1 = graph_to_networkx_G(split_tree_to_graph(string_to_split_tree(
        "Z(SX(Z, Z, Z))")))
    G2 = nx.Graph()
    G2.add_edges_from([(0, 1), (1, 2), (1, 3)])
    assert nx.is_isomorphic(G1, G2)


def test7():
    G1 = graph_to_networkx_G(split_tree_to_graph(string_to_split_tree(
        "SR(Z, Z, Z)")))
    G2 = nx.Graph()
    G2.add_edges_from([(0, 1), (1, 2)])
    assert nx.is_isomorphic(G1, G2)


def test8():
    G1 = graph_to_networkx_G(split_tree_to_graph(string_to_split_tree(
        "Z(SC(Z, SX(Z, Z)))")))
    G2 = nx.Graph()
    G2.add_edges_from([(0, 1), (1, 2), (2, 3)])
    assert nx.is_isomorphic(G1, G2)


def test9():
    G1 = graph_to_networkx_G(split_tree_to_graph(string_to_split_tree(
        "Z(SX(SC(Z, Z), Z))")))
    G2 = nx.Graph()
    G2.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0)])
    assert nx.is_isomorphic(G1, G2)


def test10():
    G1 = graph_to_networkx_G(split_tree_to_graph(string_to_split_tree(
        "KR(SC(Z, Z), Z, Z)")))
    G2 = nx.Graph()
    G2.add_edges_from([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)])
    assert nx.is_isomorphic(G1, G2)


def test11():
    G1 = graph_to_networkx_G(split_tree_to_graph(string_to_split_tree(
        "e(SX(Z, Z), SX(Z, Z))")))
    G2 = nx.Graph()
    G2.add_edges_from([(0, 1), (1, 2), (2, 3)])
    assert nx.is_isomorphic(G1, G2)


def test12():
    G1 = graph_to_networkx_G(split_tree_to_graph(string_to_split_tree(
        "e(SC(Z, Z), SC(Z, Z))")))
    G2 = nx.Graph()
    G2.add_edges_from([(0, 2), (0, 3), (1, 2), (1, 3)])
    assert nx.is_isomorphic(G1, G2)


def test13():
    G1 = graph_to_networkx_G(split_tree_to_graph(string_to_split_tree(
        "e(SX(Z, Z, Z), SX(Z, Z, Z))")))
    G2 = nx.Graph()
    G2.add_edges_from([(0, 1), (0, 2), (0, 3), (3, 4), (3, 5)])
    assert nx.is_isomorphic(G1, G2)


def test14():
    G1 = graph_to_networkx_G(split_tree_to_graph(string_to_split_tree(
        "KR(Z, Z, SX(Z, Z), SX(Z, Z))")))
    G2 = nx.Graph()
    G2.add_edges_from([(0, 1), (0, 2), (0, 4), (1, 2), (1, 4), (2, 3), (2, 4),
        (4, 5)])
    assert nx.is_isomorphic(G1, G2)


def test15():
    G1 = graph_to_networkx_G(split_tree_to_graph(string_to_split_tree(
        "SR(K(SC(Z, SX(K(SC(Z, Z), Z), K(Z, Z, Z))), Z), Z, Z)")))
    G2 = nx.Graph()
    G2.add_edges_from([(0, 2), (0, 3), (0, 7), (0, 8), (0, 9), (1, 2), (1, 3),
        (1, 7), (1, 8), (1, 9), (2, 3), (2, 7), (2, 8), (2, 9), (4, 5), (4, 6),
        (4, 7), (4, 8), (4, 9), (5, 6), (5, 7), (5, 8), (5, 9), (6, 7), (6, 8),
        (6, 9), (7, 8), (7, 9)])
    assert nx.is_isomorphic(G1, G2)


def test():
    test1()
    test2()
    test3()
    test4()
    test5()
    test6()
    test7()
    test8()
    test9()
    test10()
    test11()
    test12()
    test13()
    test14()
    test15()

    print "All tests passed!"

test()
