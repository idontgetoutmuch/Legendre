% The Precession of the Perihelion of Mercury via Legendre Polynomials
% Dominic Steinitz
% 25th August 2013

Introduction
------------

The planet Mercury has a highly elliptical orbit with a perihelion of
about 0.31 AU and an aphelion of about 0.47 AU. This ellipse is not
stationary but itself rotates about the Sun, a phenomenon known as the
precession of the perihelion. A calculation carried out using
Newtonian mechanics gives a value at variance with observation. The
deficit is explained using General Relativity.

This article calculates the precession in Haskell using Newtonian
methods.  Over a long enough period, the gravitational effect of each
outer planet on Mercury can be considered to be the same as a ring
with the same mass as the planet; in other words we assume that the
mass of each planet has been smeared out over its orbit.

More specifically, we model the mass of the ring as being
totally concentrated on one particular value of $\theta = \pi / 2$ and
one particular value of $r = a$ with total mass $M$.

$$
\begin{aligned}
M &= \int_0^{2\pi} \int_0^\pi \int_0^\infty K\, \delta(\theta - \pi / 2)\, \delta(r - a)\, r^2\sin\theta\, \mathrm{d} r\, \mathrm{d} \theta\, \mathrm{d} \phi \\
  &= 2\pi \int_0^\pi \int_0^\infty K\, \delta(\theta - \pi / 2)\, \delta(r - a)\, r^2\sin\theta\, \mathrm{d} \theta\, \mathrm{d} r \\
  &= 2\pi \int_0^\infty K\, \delta(r - a)\, r^2\, \mathrm{d} r \\
  &= 2\pi K a^2
\end{aligned}
$$

where $\delta$ is the Dirac delta function. Thus the density of our ring is

$$
\rho(r, \theta) = \frac{\delta(\theta - \pi / 2) \delta(r - a)}{2\pi a^2}
$$

Axially Symmetric Mass Distributions
------------------------------------

We consider axially
symmetric mass distributions


$$
P_{2n}(0) = \frac{(-1)^n(2n)!}{2^{2n}(n!)^2}
$$

> module Main (main) where
>
> import Initial
>
> phi i =   sunPotential
>         + innerPotential
>         + outerPotential
>   where
>     sunPotential = negate $ gConst * (massesOuter!!5) / (radius $ initQsOuter!!i)
>     innerPotential = undefined
>     outerPotential = undefined
>     radius [x, y, z] = sqrt $ x^2 + y^2 + z^2
>
> legendre0 n
>   | odd n = 0
>   | even n = (negate 1)^n * 1
>              where
>                m = n `div` 2
>
> main = putStrLn "Hello"
