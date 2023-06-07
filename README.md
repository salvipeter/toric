# Toric patch evaluator

[Toric patch](https://doi.org/10.1023/A:1015289823859) evaluation based on my [geometry library](https://github.com/salvipeter/libgeom).

The `.trp` files have the following structure:
```
<# of sides>
<1st vertex u v>
...
<# of CPs>
<1st CP u v x y z w>
...
```
... where `u` and `v` are domain (lattice) coordinates, `x`, `y`, `z` are 3D coordinates, and `w` is the rational weight.
