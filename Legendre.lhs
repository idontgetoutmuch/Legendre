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
deficit is explained using General Relativity although we do not apply
the relativistic correction in this article.

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

$$
\|\boldsymbol{r} - \boldsymbol{r}'\|^{-1} = (r^2 - 2\boldsymbol{r}\cdot\boldsymbol{r}' + r'^2)^{-1/2}
$$

Expanding the middle term on the right hand size and noting that $\theta = 0$:

$$
\begin{align}
\boldsymbol{r}\cdot\boldsymbol{r}' &= r\sin\phi\cos\theta r'\sin\phi'\cos\theta' +
                                      r\sin\phi\sin\theta r'\sin\phi'\sin\theta' +
                                      r\cos\phi r'\cos\phi' \\
                                   &= r\sin\phi r'\sin\phi'\cos\theta' +
                                      r\cos\phi r'\cos\phi' \\
                                   &= rr'(\sin\phi\sin\phi'\cos\theta' + \cos\phi\cos\phi')
\end{align}
$$

Writing $F = \sin\phi\sin\phi'\cos\theta' + \cos\phi\cos\phi'$ and noting that

$$
\frac{1}{\sqrt{1 - 2xt + t^2}} = \sum_{n=0}^\infty t^n P_n(x)
$$

where $P_n$ are the [Legendre
Polynomials](http://en.wikipedia.org/wiki/Legendre_polynomials
"Wikipedia definition") we see that when $r' < r$

$$
\|\boldsymbol{r} - \boldsymbol{r}'\|^{-1} = \frac{1}{r}\sum_{n=0}^\infty{\bigg(\frac{r'}{r}}\bigg)^n P_n(F)
$$

Applying the [Spherical Harmonic Addition
Theorem](http://mathworld.wolfram.com/SphericalHarmonicAdditionTheorem.html
"Wolfram MathWorld definition") (or see [@arfken]) we
obtain

$$
\langle\|\boldsymbol{r} - \boldsymbol{r}'\|^{-1}\rangle = \frac{1}{r}\sum_{n=0}^\infty{\bigg(\frac{r'}{r}}\bigg)^n P_n(\cos\phi) P_n(\cos\phi')
$$

Similarly when $r < r'$ we obtain

$$
\langle\|\boldsymbol{r} - \boldsymbol{r}'\|^{-1}\rangle = \frac{1}{r'}\sum_{n=0}^\infty{\bigg(\frac{r}{r'}}\bigg)^n P_n(\cos\phi) P_n(\cos\phi')
$$

Substituting into the equation for the potential for axially symmetric
mass distributions gives us

$$
\begin{align}
\Phi(r, \phi) &= -2\pi G\int_0^\pi \int_0^\infty \rho(r', \phi') \langle\|\boldsymbol{r}' - \boldsymbol{r}\|^{-1}\rangle\, r'^2\sin\phi'\, \mathrm{d} r\, \mathrm{d} \phi' \\
             &=  -2\pi G\int_0^\pi \int_0^r \rho(r', \phi')\frac{1}{r}\sum_{n=0}^\infty{\bigg(\frac{r'}{r}}\bigg)^n P_n(\cos\phi) P_n(\cos\phi')\, r'^2\sin\phi'\, \mathrm{d} r\, \mathrm{d} \phi'
\\
             &\phantom{=}  -2\pi G\int_0^\pi \int_r^\infty \rho(r', \phi')\frac{1}{r'}\sum_{n=0}^\infty{\bigg(\frac{r}{r'}}\bigg)^n P_n(\cos\phi) P_n(\cos\phi')\, r'^2\sin\phi'\, \mathrm{d} r\, \mathrm{d} \phi'
\\
             &= \sum_{n=0}^\infty \Phi_n(r) P_n(\cos\phi)
\end{align}
$$

where

$$
\begin{align}
\Phi_n(r) &= -\frac{2\pi G}{r^{n+1}}\int_0^r\int_0^\pi r'^{n+2}\rho(r', \phi')P_n(\cos\phi')\sin\phi'\,\mathrm{d}r'\,\mathrm{d}\phi' \\
          &\phantom{=} -2\pi G r^n\int_r^\infty\int_0^\pi r'^{1-n}\rho(r', \phi')P_n(\cos\phi')\sin\phi'\,\mathrm{d}r'\,\mathrm{d}\phi'
\end{align}
$$

Note that the first integral has limits $0$ to $r$ and the second has limits $r$ to $\infty$.

It is well known that the Legendre Polynomials form an orthogonal and
complete set for continuous functions. Indeed

$$
\int_{-1}^1 P_n(x)P_m(x)\,\mathrm{d}x = \frac{2\delta_{nm}}{2n + 1}
$$

Thus we can write

$$
\rho(r, \phi) = \sum_{n=o}^\infty \rho_n(r)P_n(\cos\phi)
$$

Using the orthogonality condition we have

$$
\rho_n(r) = (n + 1/2)\int_0^\pi \rho(r, \phi) P_n(\cos\phi) \sin\phi\,\mathrm{d}\phi
$$

Hence

$$
\begin{align}
\Phi_n(r) &= -\frac{2\pi G}{(n + 1/2)r^{n+1}}\int_0^r r'^{n+2}\rho_n(r')\,\mathrm{d}r' \\
          &\phantom{=} -\frac{2\pi G r^n}{n + 1/2}\int_r^\infty r'^{1-n}\rho_(r')\,\mathrm{d}r'\,\mathrm{d}r'
\end{align}
$$

Gravitational Potential of a Ring
---------------------------------

We now substitute in the axially symmetric density of a ring

$$
\begin{align}
\rho_n(r) &= (n + 1/2)\int_0^\pi \rho(r, \phi) P_n(\cos\phi) \sin\phi\,\mathrm{d}\phi \\
          &=(n + 1/2)\int_0^\pi M \frac{\delta(\phi - \pi / 2) \delta(r - a)}{2\pi a^2} P_n(\cos\phi) \sin\phi\,\mathrm{d}\phi \\
          &= (n + 1/2) M \frac{\delta(r - a)}{2\pi a^2} P_n(0)
\end{align}
$$

Substituting again

$$
\begin{align}
\Phi_n(r) &= -\frac{2\pi G}{(n + 1/2)r^{n+1}}\int_0^r r'^{n+2}\rho_n(r')\,\mathrm{d}r' \\
          &\phantom{=} -\frac{2\pi G r^n}{n + 1/2}\int_r^\infty r'^{1-n}\rho_(r')\,\mathrm{d}r'\,\mathrm{d}r' \\
          &= -\frac{2\pi G}{(n + 1/2)r^{n+1}}\int_0^r r'^{n+2} (n + 1/2) M \frac{\delta(r' - a)}{2\pi a^2} P_n(0) \,\mathrm{d}r' \\
          &\phantom{=} -\frac{2\pi G r^n}{n + 1/2}\int_r^\infty r'^{1-n} (n + 1/2) M \frac{\delta(r' - a)}{2\pi a^2} P_n(0) \,\mathrm{d}r'
\end{align}
$$


Thus for $a < r$

$$
\begin{align}
\Phi_n(r) &= -\frac{2\pi G}{r^{n+1}} a^{n+2} M \frac{1}{2\pi a^2} P_n(0) \\
          &= -\frac{G M P_n(0)}{a}\bigg(\frac{a}{r}\bigg)^{n+1}
\end{align}
$$


And for $r < a$

$$
\begin{align}
\Phi_n(r) &= -2\pi G r^n a^{1-n} M \frac{1}{2\pi a^2} P_n(0) \\
          &= -\frac{G M P_n(0)}{a} \bigg(\frac{r}{a}\bigg)^n
\end{align}
$$

Thus at $\phi = \pi / 2$ and $r < a$ we have

$$
\Phi(r) \equiv \Phi(r, \pi / 2) = -\frac{G M}{a} \sum_{n=0}^\infty P_n^2(0) \bigg(\frac{r}{a}\bigg)^n
$$

and for $r > a$

$$
\Phi(r) = -\frac{G M}{a} \sum_{n=0}^\infty P_n^2(0) \bigg(\frac{a}{r}\bigg)^{n+1}
$$

Let $M$ be the mass of the Sun then the potential due to all the Sun
and all planets at a distance $r$ (excluding the planet positioned at
$r$) is

$$
\Phi(r) = -\frac{GM}{r} - \sum_{n=0}^\infty P_n^2(0) \Bigg[\sum_{a_i < r}\frac{G m_i}{a_i}\bigg(\frac{a_i}{r}\bigg)^{n+1} + \sum_{a_i > r}\frac{G m_i}{a_i}\bigg(\frac{r}{a_i}\bigg)^n\Bigg]
$$

Apsidal Angles
--------------

An [apsis](http://en.wikipedia.org/wiki/Apsis "Wikipedia definition") is
the closest and furthest point that a planet reaches i.e. the
perihelion and the aphelion. Without the perturbing influence of the
outer planets the angle between these points, the apsidal angle, would
be $\pi$. In presence of the outer planets this is no longer the case.

Writing down the Lagrangian for a single planet we have

$$
\mathbb{L} = \frac{1}{2}m(\dot{r}^2 + r^2\dot{\theta}^2) + \Phi(\boldsymbol{r})
$$

where $\Phi$ is the total potential due to the Sun and the other planets
(as calculated above).  $\theta$ is ignorable so we have a conserved
quantity $mr^2\dot{\theta}$. We write $h$ for $r^2\dot{\theta}$ which
is also conserved.

Applying Lagrange's equation for $r$ we have

$$
\frac{\mathrm{d}}{\mathrm{d} t}\bigg(\frac{\partial \mathbb{L}}{\partial \dot{r}}\bigg) - \frac{\partial \mathbb{L}}{\partial r} = m\ddot{r} - mr\dot{\theta}^2  + \frac{\partial \Phi}{\partial r} = 0
$$

Thus the radial equation of motion is

$$
\ddot{r} - \frac{h^2}{r^3} =
       -\frac{GM}{r^2}
       - \sum_{n=0}^\infty P_n^2(0) \Bigg[
         \sum_{a_i < r}\frac{G m_i}{a_i^2}(n+1)\bigg(\frac{a_i}{r}\bigg)^{n+2}
       - \sum_{a_i > r}\frac{G m_i}{a_i^2}n\bigg(\frac{r}{a_i}\bigg)^{n-1}\Bigg]
$$

To make further progress let us take just one term for a ring outside the planet of consideration and use the trick given in [@brown:SpaceTime].

$$
\begin{align}
\frac{A}{r_{\mathrm{ap}}^2} + \frac{B}{r_{\mathrm{ap}}^3} &=
\frac{M}{r_{\mathrm{ap}}^2} -
\sum_{n=0}^\infty P_n^2(0) \Bigg[\frac{G m}{a^2}n\bigg(\frac{r_\mathrm{ap}}{a}\bigg)^{n-1}\Bigg] \\
\frac{A}{r_{\mathrm{peri}}^2} + \frac{B}{r_{\mathrm{peri}}^3} &=
\frac{M}{r_{\mathrm{peri}}^2} -
\sum_{n=0}^\infty P_n^2(0) \Bigg[\frac{G m}{a^2}n\bigg(\frac{r_\mathrm{peri}}{a}\bigg)^{n-1}\Bigg]
\end{align}
$$

Defining $g$ by writing

$$
\frac{A}{r^2} + \frac{B}{r^3} = \frac{g(r)}{r^3}
$$

we have

$$
\begin{align}
Ar_\mathrm{peri} + B &= g(r_\mathrm{peri}) \\
Ar_\mathrm{ap} + B &= g(r_\mathrm{ap})
\end{align}
$$

giving

$$
A = \frac{g(r_\mathrm{peri}) - g(r_\mathrm{ap})}{r_\mathrm{peri} - r_\mathrm{ap}}
$$

Using the Taylor approximation

$$
\begin{align}
g(r_{\mathrm{peri}}) &\approx g(r_\mathrm{mean}) + \frac{r_{\mathrm{peri}} - r_{\mathrm{ap}}}{2} g'(r_\mathrm{mean}) \\
g(r_{\mathrm{ap}}) &\approx g(r_\mathrm{mean}) - \frac{r_{\mathrm{peri}} - r_{\mathrm{ap}}}{2} g'(r_\mathrm{mean})
\end{align}
$$

Thus

$$
A \approx g'(r_\mathrm{mean})
$$

Then since

$$
g(r) = Mr -
\sum_{n=0}^\infty P_n^2(0) \Bigg[G m a n\bigg(\frac{r_\mathrm{ap}}{a}\bigg)^{n+2}\Bigg]
$$

We have

$$
A = g'(r_\mathrm{mean}) = M -
\sum_{n=0}^\infty P_n^2(0) \Bigg[G m n(n+2)\bigg(\frac{r_\mathrm{peri}}{a}\bigg)^{n+1}\Bigg]
$$

It is a nuisance to be continually writing $r_\mathrm{peri}$. From now on this is denoted by $r$. Using

$$
B = r^3 g(r) - r A
$$

We obtain

$$
\begin{align}
B &= Mr -
\sum_{n=0}^\infty P_n^2(0) \Bigg[G m a n\bigg(\frac{r}{a}\bigg)^{n+2}\Bigg] \\
  &\phantom{=} -r\Bigg(M -
\sum_{n=0}^\infty P_n^2(0) \Bigg[G m n(n+2)\bigg(\frac{r}{a}\bigg)^{n+1}\Bigg]\Bigg) \\
  &= \sum_{n=0}^\infty P_n^2(0) \Bigg[G m a n(n+1)\bigg(\frac{r}{a}\bigg)^{n+2}\Bigg] \\
\end{align}
$$

We can therefore re-write the radial equation of motion approximately as

$$
\begin{align}
\ddot{r} - \frac{h^2}{r^3}     &= -\frac{A}{r^2} - \frac{B}{r^3} \\
\ddot{r} - \frac{h^2 - B}{r^3} &= -\frac{A}{r^2}
\end{align}
$$

Now let us re-write the equation of motion as a relation between $r$ and $\theta$.

$$
\dot{r} = \frac{\mathrm{d} \theta}{\mathrm{d} t}\frac{\mathrm{d} t}{\mathrm{d} \theta}\frac{\mathrm{d} r}{\mathrm{d} t} = \frac{h}{r^2}\frac{\mathrm{d} r}{\mathrm{d} \theta}
$$

$$
\begin{align}
\ddot{r} &= \frac{h}{r^2}\frac{\mathrm{d}}{\mathrm{d} t}\bigg(r^{-2}\frac{\mathrm{d} r}{\mathrm{d} \theta}\bigg) \\
         &= \frac{h}{r^2}\frac{\mathrm{d} \theta}{\mathrm{d} t}\frac{\mathrm{d} t}{\mathrm{d} \theta}\frac{\mathrm{d}}{\mathrm{d} t}\bigg(r^{-2}\frac{\mathrm{d} r}{\mathrm{d} \theta}\bigg) \\
         &= \frac{h}{r^2}\frac{\mathrm{d}}{\mathrm{d} \theta}\bigg(r^{-2}\frac{\mathrm{d} r}{\mathrm{d} \theta}\bigg)
\end{align}
$$

Thus we have

$$
\frac{\mathrm{d}}{\mathrm{d} \theta}\bigg(r^{-2}\frac{\mathrm{d} r}{\mathrm{d} \theta}\bigg) - \frac{(h^2 - B)}{r} = -\frac{A}{r^2}
$$

Letting $u = 1 /r$ we can re-write this as

$$
\frac{1}{(1 - B / h^2)} \frac{\mathrm{d}^2 u}{\mathrm{d} \theta^2} + u = \frac{A}{h^2(1 - B / h^2)}
$$

This is the equation for simple harmonic motion with $\omega^2 = 1 - B
/ h^2$ and since for a circular orbit $h^2 = GMr$ we can write

$$
\begin{align}
\omega &= \sqrt{1 - B / h^2} \approx 1 - \frac{1}{2}\frac{B}{h^2} \\
       &= 1 -  \frac{1}{2}\frac{m}{M}\sum_{n=0}^\infty P_n^2(0) n(n+1)\bigg(\frac{r}{a}\bigg)^{n+1} \\
\end{align}
$$

and therefore the change in radians per revolution is

$$
\Delta \theta = |2\pi (\omega - 1)| = \pi\frac{m}{M}\sum_{n=0}^\infty P_n^2(0) n(n+1)\bigg(\frac{r}{a}\bigg)^{n+1}
$$

To convert this to arc-seconds per century we apply a conversion factor

$$
414.9 \frac{360}{2\pi} 3600
$$

where 414.9 is the number of orbits of Mercury per century.

> earthPerihelion :: Double
> earthPerihelion = 1.470983e11
>
> earthAphelion   :: Double
> earthAphelion   = 1.520982e11
>
> earthMajRad :: Double
> earthMajRad = (earthPerihelion + earthAphelion) / 2
>
> venusMass = 4.8676e24
> venusMajRad = 108208000e3
>
> mercuryMajRad = 57909100e3
>
> marsAphelion = 249209300e3
> marsPerihelion = 206669000e3
> marsMajRad = (marsAphelion + marsPerihelion) / 2
> marsMass = 6.4185e23
>
> jupiterPerihelion :: Double
> jupiterPerihelion = 7.405736e11
>
> jupiterAphelion   :: Double
> jupiterAphelion   = 8.165208e11
>
> jupiterMajRad :: Double
> jupiterMajRad = (jupiterPerihelion + jupiterAphelion) / 2
>
> conv x = x * 414.9 * (360 / (2 * pi)) * 3600
>
> deltaThetas majRad mass =
>   zipWith (*)
>   (perturbations $ mercuryMajRad / majRad)
>   (repeat $ pi * mass / sunMass)
>
> coeffs :: [Rational]
> coeffs = zipWith (*) [3,5..] alphas
>
> alphas :: [Rational]
> alphas = zipWith (*)
>          (map (^2) $ drop 1 $ filter (/= 0) $ legendre0s)
>          [2,4..]
>
> poly :: Num a => a -> [a]
> poly x = map (x^) [3,5..]
>
> perturbationsR :: Rational -> [Rational]
> perturbationsR x = zipWith (*) coeffs (poly x)
>
> perturbations :: Fractional a => a -> [a]
> perturbations x = zipWith (*) (map fromRational coeffs) (poly x)

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

    [ghci]
    sum $ map conv $ take 10 $ deltaThetas jupiterMajRad jupiterMass

Appendix
--------

Note the lectures by Fitzpatrick [@Fitz:Newtonian:Dynamics] use a
different approximation for the apsidal angle

$$
\psi = \pi\bigg(3 + \frac{r \mathrm{d} F / \mathrm{d} r}{F}\bigg)^{-1/2}
$$

We do not derive this here but note that the expansion and
approximation are not entirely straightforward and are given here for
completenes. Note that the result derived this way is identical to the
result obtained in the main body of the article.

The radial force is given by $F(r) = -\mathrm{d}\Phi(r) / \mathrm{d} r$

$$
F(r) = -\frac{GM}{r^2}
       - \sum_{n=0}^\infty P_n^2(0) \Bigg[
         \sum_{a_i < r}\frac{G m_i}{a_i^2}(n+1)\bigg(\frac{a_i}{r}\bigg)^{n+2}
       - \sum_{a_i > r}\frac{G m_i}{a_i^2}n\bigg(\frac{r}{a_i}\bigg)^{n-1}\Bigg]
$$

We also have

$$
r\frac{\mathrm{d} F}{\mathrm{d} r} =
  2\frac{GM}{r^2}
  + \sum_{n=0}^\infty P_n^2(0) \Bigg[
    \sum_{a_i < r}\frac{G m_i}{a_i^2}(n+1)(n+2)\bigg(\frac{a_i}{r}\bigg)^{n+2}
  + \sum_{a_i > r}\frac{G m_i}{a_i^2}n(n-1)\bigg(\frac{r}{a_i}\bigg)^{n-1}\Bigg]
$$

Thus

$$
2F(r) + r\frac{\mathrm{d} F}{\mathrm{d} r} =
    \sum_{n=0}^\infty P_n^2(0) \Bigg[
    \sum_{a_i < r}\frac{G m_i}{a_i^2}n(n+1)\bigg(\frac{a_i}{r}\bigg)^{n+2}
  + \sum_{a_i > r}\frac{G m_i}{a_i^2}n(n+1)\bigg(\frac{r}{a_i}\bigg)^{n-1}\Bigg]
$$

Re-arranging

$$
\bigg(3 + \frac{r \mathrm{d} F / \mathrm{d} r}{F}\bigg)^{-1/2} =
\bigg(1 + 2 + \frac{r \mathrm{d} F / \mathrm{d} r}{F}\bigg)^{-1/2}
$$

we note that the last two terms can be re-written with a numerator of

$$
\sum_{n=0}^\infty P_n^2(0) \Bigg[
\sum_{a_i < r}\frac{G m_i}{a_i^2}n(n+1)\bigg(\frac{a_i}{r}\bigg)^{n+2}
+ \sum_{a_i > r}\frac{G m_i}{a_i^2}n(n+1)\bigg(\frac{r}{a_i}\bigg)^{n-1}\Bigg]
$$

and a denominator which is dominated by the $-GM / r^2$. Thus

$$
\begin{align}
2 + \frac{r \mathrm{d} F / \mathrm{d} r}{F} &\approx
-\sum_{n=0}^\infty P_n^2(0) \Bigg[
\sum_{a_i < r}\frac{m_i r^2}{M a_i^2}n(n+1)\bigg(\frac{a_i}{r}\bigg)^{n+2}
+ \sum_{a_i > r}\frac{m_i r^2}{M a_i^2}n(n+1)\bigg(\frac{r}{a_i}\bigg)^{n-1}\Bigg] \\
&=
-\sum_{n=0}^\infty P_n^2(0) n(n+1)\Bigg[
  \sum_{a_i < r}\frac{m_i}{M}\bigg(\frac{a_i}{r}\bigg)^n
+ \sum_{a_i > r}\frac{m_i}{M}\bigg(\frac{r}{a_i}\bigg)^{n+1}\Bigg]
\end{align}
$$

Since this term is $\ll 1$ we can expand the term of interest further

$$
\begin{align}
\bigg(1 + 2 + \frac{r \mathrm{d} F / \mathrm{d} r}{F}\bigg)^{-1/2} &\approx
\Bigg(1
-\sum_{n=0}^\infty P_n^2(0) n(n+1)\Bigg[
  \sum_{a_i < r}\frac{m_i}{M}\bigg(\frac{a_i}{r}\bigg)^n
+ \sum_{a_i > r}\frac{m_i}{M}\bigg(\frac{r}{a_i}\bigg)^{n+1}\Bigg]
\Bigg)^{-1/2} \\
&=
1 + \frac{1}{2}
  \sum_{n=0}^\infty P_n^2(0) n(n+1)\Bigg[
  \sum_{a_i < r}\frac{m_i}{M}\bigg(\frac{a_i}{r}\bigg)^n
+ \sum_{a_i > r}\frac{m_i}{M}\bigg(\frac{r}{a_i}\bigg)^{n+1}\Bigg]
\end{align}
$$

Bibliography
------------
