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
[co-recursively](http://en.wikipedia.org/wiki/Corecursion "Wikipedia
definition"). We could use the [Stream
package](http://hackage.haskell.org/package/Stream "Hackage") but let
us stay with lists.

> module Legendre (
>     legendre0s
>   , main) where
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

    [ghci]
    take 10 $ legendre0s

This article calculates the precession in Haskell using Newtonian
methods.  Over a long enough period, the gravitational effect of each
outer planet on Mercury can be considered to be the same as a ring
with the same mass as the planet; in other words we assume that the
mass of each planet has been smeared out over its orbit. Probably one
can model Saturn's rings using this technique but that is certainly
the subject of a different blog post.

More specifically, we model the mass of the ring as being
totally concentrated on one particular value of $\phi = \pi / 2$ and
one particular value of $r = a$ with total mass $M$.

$$
\begin{aligned}
M &= \int_0^{2\pi} \int_0^\pi \int_0^\infty K\, \delta(\phi - \pi / 2)\, \delta(r - a)\, r^2\sin\phi\, \mathrm{d} r\, \mathrm{d} \phi\, \mathrm{d} \theta \\
  &= 2\pi \int_0^\pi \int_0^\infty K\, \delta(\phi - \pi / 2)\, \delta(r - a)\, r^2\sin\phi\, \mathrm{d} \phi\, \mathrm{d} r \\
  &= 2\pi \int_0^\infty K\, \delta(r - a)\, r^2\, \mathrm{d} r \\
  &= 2\pi K a^2
\end{aligned}
$$

where $\delta$ is the Dirac delta function. Thus the density of our ring is

$$
\rho(r, \phi) = M \frac{\delta(\phi - \pi / 2) \delta(r - a)}{2\pi a^2}
$$

Acknowledgement
---------------

This blog follows the exposition given in [@Fitz:Newtonian:Dynamics]
concretized for the precession of the perihelion of Mercury with some
of the elisions expanded. More details on Legendre Polynomials can be
found in [@Bowles:Legendre:Polynomials].

Axially Symmetric Mass Distributions
------------------------------------

We consider axially symmetric mass distributions in spherical polar
co-ordinates $(r, \phi, \theta)$ where $r$ runs from $0$ to $\infty$,
$\phi$ (the polar angle) runs from $0$ to $\pi$ and $\theta$ (the
azimuthal angle) runs from $0$ to $2\pi$.

For clarity we give their conversion to cartesian co-ordinates.

$$
\begin{align}
x &= r\sin\phi\cos\theta \\
y &= r\sin\phi\sin\theta \\
z &= r\cos\phi
\end{align}
$$

The volume element in spherical polar co-ordinates is given by
$r^2\sin\phi\,\mathrm{d} r\,\mathrm{d} \phi\,\mathrm{d} \theta$.

The gravitational potential given by $N$ masses each of mass $m_i$ and
at position $\boldsymbol{r}_i$ is:

$$
\Phi(\boldsymbol{r}) = -G\sum_{i=1}^N\frac{m_i}{\|\boldsymbol{r}_i - \boldsymbol{r}\|}
$$

If instead of point masses, we have a mass distribution $\rho(\boldsymbol{r})$ then

$$
\Phi(\boldsymbol{r}) = -G\int_{\mathbb{R}^3}\frac{\rho(\boldsymbol{r}')}{\|\boldsymbol{r}' - \boldsymbol{r}\|}\, \mathrm{d} V
$$

where $\mathrm{d} V$ is the volume element.

If the mass distribution is axially symmetric then so will the
potential. In spherical polar co-ordinates:

$$
\begin{align}
\Phi(r, \phi) &= -G\int_0^{2\pi} \int_0^\pi \int_0^\infty \frac{\rho(r', \phi')}{\|\boldsymbol{r}' - \boldsymbol{r}\|}\, r'^2\sin\phi'\, \mathrm{d} r\, \mathrm{d} \phi'\, \mathrm{d} \theta' \\
              &= -2\pi G\int_0^\pi \int_0^\infty \rho(r', \phi') \langle\|\boldsymbol{r}' - \boldsymbol{r}\|^{-1}\rangle\, r'^2\sin\phi'\, \mathrm{d} r\, \mathrm{d} \phi' \\
\end{align}
$$

where $\langle\ldots\rangle$ denotes the average over the azimuthal angle.

> theta i =   sunPotential
>         + innerPotential
>         + outerPotential
>   where
>     sunPotential = negate $ gConst * (massesOuter!!5) / (radius $ initQsOuter!!i)
>     innerPotential = undefined
>     outerPotential = undefined
>     radius [x, y, z] = sqrt $ x^2 + y^2 + z^2
>
> main = putStrLn "Hello"
