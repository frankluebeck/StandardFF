The GAP 4 package `StandardFF'
==============================

This package contains an implementation of *standard* generators of
finite fields and of cyclic subgroups in the multiplicative groups
of finite fields, as described in:

     Frank Lübeck, Standard generators of finite fields and their 
     cyclic subgroups, Journal of Symbolic Computation, Vol 117,
     2023, p. 51-67
     (https://doi.org/10.1016/j.jsc.2022.11.001)

A preprint version of the article is:

     https://arxiv.org/abs/2107.02257


Installation and usage
======================

The package is redistributed with GAP (4.12 and later).

If you want to install a newer version:

Copy the package into a subdirectory below $GAPHOME/pkg.
No further installation is needed.

Load the package within your GAP session with
   LoadPackage("StandardFF");

Regenerate documentation inside the package directory with
   gap makedocrel.g

Access documentation in GAP session with

   ?StandardFF:?

