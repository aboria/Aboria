import numpy as np

def compute_tags(sources,m):
    minp = np.zeros(3,len(sources))
    maxp = np.ones(3,len(sources))
    tags = np.zeros((1,len(sources)))
    for i in range(m):
        mid = 0.5*(maxp-minp)
        bit = (source >= mid)
        minp[bit] = mid
        maxp[!bit] = mid
        tags = (tags << 1) | bit

def child_tag_mask(active_nodes,child,level,max_level):
    return active_nodes | (child << (max_level-level)*2)

EMPTY = 0
LEAF = 1
NODE = 2

def classify_children(lower,upper,level,max_level,thresh):
    if level == max_level:
        child_node_kind = np.array(len(lower),LEAF)
    else:
        child_node_kind = np.array(len(lower),NODE)

    n = upper-lower
    empty = n == 0
    child_node_kind[n==0] = EMPTY
    child_node_kind[n<thresh] = LEAF


def build_tree(sources,m):
    N = len(sources)

    # tag and sort source indices
    tags = compute_tags(sources,m)
    sorted_indicies = np.argsort(tags)

    active_nodes = np.zeros(1)
    nodes = np.array([],dtype=int)
    leaves = np.array([],dtype=int)
    for i in range(m):
        children = np.zeros(4*len(active_nodes))
        for j in range(8):
            children[j::8] = compute_tag_mask(active_nodes,j,i,m)
        lower_bounds = np.searchsorted(tags,children,side='left',sorter=sorted_indicies)
        upper_bounds = np.searchsorted(tags,children,side='right',sorter=sorted_indicies)

        empty,leaf,node = classify_children(lower_bounds,upper_bounds,i,m,thresh)

        nodes_on_this_level = np.cumsum(node,dtype=int)
        leafs_on_this_level = np.cumsum(leaf,dtype=int)

        #create child nodes
        new_nodes = np.array(child_node_kind)
        new_nodes[empty] = get_empty_id()
        new_nodes[leaf] = get_leaf_id(len(leaves) + leafs_on_this_level[leaf])
        new_nodes[node] = len(nodes) + 8*nodes_on_this_level[node]
        np.concatenate(nodes,new_nodes)

        #create leaves
        np.concatenate(leaves,new_leaves)

        #activate nodes for next level






if __name__ == "__main__":

    N = 10
    sources = np.random.rand(3,N)
    weights = np.random.rand(1,N)
    targets = sources

    tree = build_tree(sources,m)

    field_fmm = fmm(tree,sources,weights,targets)

    field_direct = direct(sources,weights,targets)

