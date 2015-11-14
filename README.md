# unlabeled-graph-samplers
## Unbiased random samplers for distance-hereditary and three-leaf power graphs

### Split-tree string

#### Semantics

Symbol | Meaning
:----: | :------
**Z** | A leaf node
**KR** | A clique root (can only appear as the root of the tree)
**SR** | A star root (can only appear as the root of the tree)
**K** | A clique that has been entered from another node
**SX** | A star that has been entered from another node at one of its extremities
**SC** | A star that has been entered from another node at its center
e(**A**, **B**) | An edge that connects nodes **A** and **B**
**A**(**B<sub>1</sub>**, ..., **B<sub>k</sub>**) | **B<sub>1</sub>**, ..., **B<sub>k</sub>** are neighbors of **A** (if **A** is **SR** or **SX**, then **B<sub>1</sub>** is connected to the center of **A**)

#### Examples

Symbol | Meaning
:----: | :------
**Z**(**K**(**Z**, **Z**)) | A leaf node connected to a clique that has two other leaves as neighbors
**KR**(**Z**, **Z**, **Z**) | A clique with three leaves as neighbors (same as previous)
**Z**(**SC**(**Z**, **Z**)) | A leaf connected to the center of a star that has two leaves as its extremities
**Z**(**SX**(**Z**, **Z**)) | A leaf connected to an extremity of a star that has a leaf as its center and a leaf as its other extremity (same as previous)
e(**SC**(**Z**, **Z**), **SC**(**Z**, **Z**)) | An edge joining two star nodes at their centers, each of which has two leaves as its extremities
**SR**(**K**(**Z**, **Z**), **Z**, **Z**, **Z**) | A star with three leaves as its extremities and whose center is connected to a clique with two other leaves as neighbors

---

### split_tree.py

Dependencies: networkx, matplotlib.pyplot


This package can be used to:

  * Compute the original graph corresponding to a given distance-hereditary or
    three-leaf power split tree string (typically this string will have been
    computed by the distance-hereditary or three-leaf power samplers
    described above)

  * Draw the original graph to the screen or export a drawing of the graph
    to a file


#### Example

```python
from split_tree import *
split_tree = string_to_split_tree("SR(K(SC(Z, SX(K(SC(Z, Z), Z), \
								  K(Z, Z, Z))), Z), Z, Z)")
graph = split_tree_to_graph(split_tree)
draw_graph(graph)
export_graph(graph, "drawing.pdf")
export_graph(graph, "drawing.png")
export_graph(graph, "drawing") # If no format is specified, defaults to png
```
