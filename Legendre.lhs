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

Just to give a flavour of the Haskell, we will have to calculate
values of the infinite series of [Legendre
Polynomials](http://en.wikipedia.org/wiki/Legendre_polynomials
"Wikipedia definition") evaluated at 0. We have

$$
\begin{aligned}
P_{2n}(0)   &= \frac{(-1)^n(2n)!}{2^{2n}(n!)^2} \\
P_{2n+1}(0) &= 0
\end{aligned}
$$

Since we are dealing with infinite series we will want to define this
co-recursively. We could use the [Stream
package](http://hackage.haskell.org/package/Stream "Hackage") but let
us stay with lists.

> module Main (main) where
>
> import Data.List
>
> import Initial
>
> legendre0s :: [Rational]
> legendre0s = interleave legendre0Evens legendre0Odds
>   where
>     legendre0Evens = 1 : zipWith f [1..] legendre0Evens
>       where f n p = negate $ p * (2 * n * (2 * n - 1)) / (2^2 * n^2)
>     legendre0Odds = 0 : legendre0Odds
>
> interleave :: [a] -> [a] -> [a]
> interleave = curry $ unfoldr g
>   where
>     g ([],  _) = Nothing
>     g (x:xs, ys) = Just (x, (ys, xs))

This article calculates the precession in Haskell using Newtonian
methods.  Over a long enough period, the gravitational effect of each
outer planet on Mercury can be considered to be the same as a ring
with the same mass as the planet; in other words we assume that the
mass of each planet has been smeared out over its orbit. Probably one
can model Saturn's rings using this technique but that is certainly
the subject of a different blog post.

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

> phi i =   sunPotential
>         + innerPotential
>         + outerPotential
>   where
>     sunPotential = negate $ gConst * (massesOuter!!5) / (radius $ initQsOuter!!i)
>     innerPotential = undefined
>     outerPotential = undefined
>     radius [x, y, z] = sqrt $ x^2 + y^2 + z^2
>
> main = putStrLn "Hello"
