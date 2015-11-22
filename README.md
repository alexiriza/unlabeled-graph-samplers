# unlabeled-graph-samplers
## Unbiased random samplers for distance-hereditary and three-leaf power graphs <br> (Paper: http://arxiv.org/abs/1511.06037)

### Samplers (Maple)

#### Dependencies

  * **combstruct** package
  * **proba_laws.mpl** (standard probability laws - by Carine Pivoteau)

#### Usage

##### distance_hereditary_sampler.mpl

  ```
  gen_dh_cycle_pointed(z)
  ```
  For a parameter 0 < *z* <= 0.137935, this function randomly samples a
  (split-tree of a) distance-hereditary graph. This is a Boltzmann sampler for
  cycle-pointed three-leaf power graphs, hence an unbiased sampler for
  three-leaf power graphs.

  By unbiased, we mean that for a fixed value of *z*, any two graphs of the same
  size (size = number of leaves in the split tree = number of vertices in the
  graph) are equally likely to be drawn.

##### three_leaf_power_sampler.mpl

  ```
  gen_3lp_cycle_pointed(z)
  ```
  For a parameter 0 < *z* <= 0.259845, this function randomly samples a
  (split-tree of a) three-leaf power graph. This is a Boltzmann sampler for
  cycle-pointed three-leaf power graphs, hence an unbiased sampler for
  three-leaf power graphs.

##### sampler_tools.mpl

  ```
  radius(poly)
  ```
  Returns an estimate of the radius of convergence of the formal power series
  given by *poly* (see the first two .mpl files for the various power series
  that are used).

  ```
  num_Z(split_tree_string)
  num_cliques(split_tree_string)
  num_stars(split_tree_string)
  ```
  Return the number of leaf, clique, and star nodes (respectively) in the
  specified *split_tree_string*.

#### Examples

```
read "distance_hereditary_sampler.mpl":

radius(dhpoly);
>> 0.1379358986434849

distance_hereditary_sample:=gen_dh_cycle_pointed(0.137935);
>> distance_herediary_sample:=SR(K(SC(SX(SC(Z, Z, SX(Z, SX(K(Z, SX(Z, Z)),
							  SX(SC(Z, SX(Z, Z, Z)), Z)))), SX(Z, Z)),
							  SX(SC(SX(Z, Z), Z), Z)), Z), Z, Z)

num_Z(distance_hereditary_sample);
>> 20
num_cliques(distance_hereditary_sample);
>> 2
num_stars(distance_hereditary_sample);
>> 14
```

```
read "three_leaf_power_sampler.mpl":

radius(threelppoly);
>> 0.259845809752822

three_leaf_power_sample:=gen_3lp_cycle_pointed(0.259845);
>> three_leaf_power_sample:=KR(SX(Z, SX(Z, SX(Z, SX(Z, Z, SX(Z, SX(Z, SX(Z,
							SX(K(Z, Z, Z, Z), Z), SX(Z, SX(Z, SX(Z, SX(Z,
							SX(Z, SX(Z, SX(Z, SX(K(Z, Z, Z, Z), K(Z, Z))), Z),
							K(Z, Z)), Z)), SX(Z, Z, Z)))))))), SX(Z, SX(Z, Z),
							K(Z, Z, Z))), K(Z, Z)), Z, Z)

num_Z(three_leaf_power_sample);
>> 43
num_cliques(three_leaf_power_sample);
>> 7
num_stars(three_leaf_power_sample);
>> 19
```

#### Notes

  * The namespaces of **distance_hereditary_sampler.mpl** and
    **three_leaf_power_sampler.mpl** overlap, so one should not read them into
    the same Maple worksheet.

  * **sampler_tools.mpl** is read by each of the first two .mpl files, so it is
    not necessary to read it as well.

---

### Split-tree string representation

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

### Parsing module (Python)

#### Dependencies

  * **networkx**
  * **matplotlib.pyplot**


#### Usage

##### split_tree.py

  ```
  string_to_split_tree(split_tree_string)
  ```
  Returns an internal split_tree representation corresponding to the specified
  distance-hereditary or three-leaf power *split_tree_string* (typically this
  string will have been returned by the distance-hereditary or three-leaf
  power samplers described above).

  ```
  split_tree_to_graph(split_tree)
  ```
  Returns the graph, expressed as a list of adjacencies over a set {1, ..., n}
  of vertices, corresponding to the specified internal *split_tree*
  representation (typically computed by the previous function).

  ```
  draw_graph(graph)
  ```
  Draws the specified *graph* to the screen.

  ```
  export_graph(graph, file_path)
  ```
  Exports a drawing of the *graph* to the specified *file_path*. An format
  extension may be provided from the list: *bmp, eps, gif, jpeg, jpg, pdf, pgf,
  png, ps, raw, rgba, svg, svgz, tif, tiff*; otherwise, *png* is used by
  default.


#### Example

```python
from split_tree import *
split_tree = string_to_split_tree("SR(K(SC(Z, SX(K(SC(Z, Z), Z), \
								  K(Z, Z, Z))), Z), Z, Z)")
graph = split_tree_to_graph(split_tree)
draw_graph(graph) # Draws graph to the screen
export_graph(graph, "drawing.pdf") # Saves drawing of graph in pdf format
export_graph(graph, "drawing.png") # Saves drawing of graph in png format
export_graph(graph, "drawing") # If no format is specified, defaults to png
```
