% The Precession of the Perihelion of Mercury via Legendre Polynomials
% Dominic Steinitz
% 25th August 2013

We model the mass of the ring as being totally concentrated on one
particular value of $\theta = \pi / 2$ and one particular value of $r
= a$ with total mass $M$.

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
