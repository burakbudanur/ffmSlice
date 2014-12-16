## ffmSlice

This repository holds the Matlab/Octave scripts which were used in the
in-slice integration of the Kuramoto-Sivashinsky system in 
[Budanur et al. (2014)](http://arxiv.org/abs/1405.1096). 

## Quick start

In order to test the code on your system, start a Matlab or Octave 
instance in `ffmSice/scripts/` and run `ks_statesp_red` or 
`ks_statesp_full`. When runs, `ks_statesp_red` and `ks_statesp_full`
produces Matlab/Octave-equivalents of Figure 2 b and Figure 2 a of 
[Budanur et al. (2014)](http://arxiv.org/abs/1405.1096) respectively.  

## Tested on

Octave 3.8.1, and Matlab R2013b

## File descriptions

`ffmSlice/scripts/`

* `Lg.m` : SO(2) generator

* `LieEl.m` : SO(2) element

* `SliceCond.m` : Slice condition

* `slicep.m` : Slice fixing point (0,1,0,0,...) 

* `vel.m` : Velocity function in real representation

* `vksred.m` : Reduced velocity function in real representation

* `gradV.m` : KS stability matrix in real representation

* `gradVred.m` : KS reduced stability matrix

* `ksETDRK4.m` : KS solver in full state space

* `ksETDRK4red.m` : KS solver in the first-mode slice

* `ks_statesp_red.m` : Computes the unstable manifold of TW_1 on the 
first Fourie mode slice, and plots it on a state space projection along 
with TW_2, and RPO_{33.50} (arXiv: 1405.1096, Fig 2.b).  

* `ks_statesp_full.m` : Same with `ks_statesp_red.m`, but trajectories 
are computed in the full state space, in addition, group orbits of TW_1 
and TW_2 are drawn (arXiv: 1405.1096, Fig 2.a).

`ffmSlice/docs/`

* docSlice.tex: Latex document which explains implementation details of 
the in-slice integration of Kuramoto-Sivashinsky system. 
`pdflatex docSlice`

* def.tex: Definitions used in `docSlice.tex`

## License

MIT license 

Copyright (c) 2014 Nazmi Burak Budanur, 
[http://cns.physics.gatech.edu/~burak](http://cns.physics.gatech.edu/~burak)
