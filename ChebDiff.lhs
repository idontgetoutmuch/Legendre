% Chebyshev Approximations
% Dominic Steinitz
% 13th November 2014

---
bibliography: Legendre.bib
---

Introduction
============

It is well known that one can represent a sufficiently well-behaved
function on the interval $[-\pi, \pi]$ by its Fourier expansion. In
the case of a symmetric function, $f(-x) = f(x)$, this becomes

$$
f(\theta) = \frac{a_0}{2} + \sum_{n=1}^\infty a_n \cos{n\theta}
$$

If we make a change of variable, $x = \cos{\theta}$ we obtain

$$
g(x) = \frac{a_0}{2} + \sum_{n=1}^\infty a_n T_n(x)
$$

where $T_n(x) = \cos{n\theta}$ is the $n$-th Chebyshev polynomial and
$g$ is an aribtrary (sufficiently well-behaved) function on $[-1, 1]$.

Here are a few examples of Chebyshev polynomials.

```{.dia width='800'}
import ChebDiff
import ChebDiag
import Data.Vector ( toList )

dia = diag (\i x -> toList $ chebUnfold i x)
````

We would like to approximate functions by a truncated Chebyshev
expansion, for example, in order to solve differential
equations. However, the coefficients are given by

$$
\begin{aligned}
a_n &= \frac{1}{\pi}\int_{-\pi}^{\pi} f(\theta)\cos{(n\theta)}\,\mathrm{d}\theta \\
    &= \frac{1}{\pi}\int_{-\pi}^{\pi} g(cos{(\theta)})\cos{(n\theta)}\,\mathrm{d}\theta \\
    &= \frac{1}{\pi}\int_{-1}^{1} \frac{g(x)T_n(x)}{\sqrt{1 - x^2}}\,\mathrm{d}x
\end{aligned}
$$

So it seems if we have to evaluate (potentially many) integrals in
order to derive the approximation then we have not gained very much.

Fortunately it turns out by using some of the many properties of
Chebyshev polynomials and the fast Fourier transform (FFT) we can
avoid integration and calcuate the approximation very efficiently.

Haskell Preamble
----------------

> {-# OPTIONS_GHC -Wall                      #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults    #-}
> {-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
> {-# OPTIONS_GHC -fno-warn-missing-methods  #-}
> {-# OPTIONS_GHC -fno-warn-orphans          #-}

> module ChebDiff where

> import Prelude hiding ( length, sum, zipWith, zipWith3,
>                         map, (++), reverse, drop, replicate )
> -- import qualified Prelude as P
> import Data.Complex
> import Data.Vector hiding ( tail )
> import qualified Data.Vector as V
> import Numeric.FFT

> -- import Text.Printf

Chebyshev Polynomial Properties
-------------------------------

The Chebyshev polynomials satisify a discrete orthogonality condition

$$
\langle T_j, T_k\rangle \triangleq
\sum_{l=0}^n \cos j\hat{x}_l \cos k\hat{x}_l =
\begin{cases}
0           & \text{if } j \neq k \\
\frac{n}{2} & \text{if } j = k \neq 0 \\
n           & \text{if } j = k = 0
\end{cases}
$$

where

$$
\hat{x}_l = \frac{\pi(j + 1/2)}{n+1}
$$

are the Chebyshev zeros of $T_{n+1}$ (they have other names but this seems safest).

The case $j = k = 0$ is obvious.

In the case $j \neq k$ and writing

$$
\cos(j\hat{x}_i) = \frac{1}{2}\big(z^{j\hat{x}_i} + z^{-j\hat{x}_i}\big) \\
$$

we have

$$
\begin{aligned}
\sum_{l=0}^n \big(z^{j\hat{x}_l} + z^{-j\hat{x}_l}\big)\big(z^{k\hat{x}_l} + z^{-k\hat{x}_l}\big) & =
\sum_{l=0}^n z^{(j+k)\hat{x}_l} + z^{-(j+k)\hat{x}_l} + z^{(j-k)\hat{x}_l} + z^{-(j-k)\hat{x}_l} \\
& = \sum_{l=0}^n z^{(j+k)(hl + h/2)} + z^{-(j+k)(hl + h/2)} + z^{(j-k)(hl + h/2)} + z^{-(j-k)(hl + h/2)} \\
& = z^{(j+k)h/2}\frac{z^{(j+k)(n+1)h} - 1}{z^{(j+k)h} - 1} + z^{-(j+k)h/2}\frac{z^{-(j+k)(n+1)h} - 1}{z^{-(j+k)h} - 1} \\
& \quad + z^{(j+k)h/2}\frac{z^{(j+k)(n+1)h} - 1}{z^{(j+k)h} - 1} + z^{-(j+k)h/2}\frac{z^{-(j+k)(n+1)h} - 1}{z^{-(j+k)h} - 1} \\
& = \frac{z^{(j+k)\pi} - 1}{z^{(j+k)h/2} - z^{-(j+k)h/2}} + \frac{z^{-(j+k)\pi} - 1}{z^{-(j+k)h/2} - z^{(j+k)h/2}} \\
& \quad + \frac{z^{(j+k)\pi} - 1}{z^{(j+k)h/2} - z^{-(j+k)h/2}} + \frac{z^{-(j+k)\pi} - 1}{z^{-(j+k)h/2} - z^{(j+k)h/2}}
\end{aligned}
$$

All the denominators on the RHS are imaginary and all the numerators
are real thus the RHS is imaginary. Since the LHS is real, both must
be 0.

A similar argument applies when $j = k \neq 0$ except that the last
two terms in the RHS of the first line are 1. Summing these and
equating real parts gives $\frac{2n}{4} = \frac{n}{2}$ as required.

We can approximate any function by its expansion in Chebyshev polynomials

$$
f(x) \approx \sum_0^n c_i T_i(x)
$$

We can substitute in a Chebyshev zero of $T_{n+1}$, $\hat{x}_k$

$$
\sum_{k=0}^{n} f(\hat{x}_k) T_j(\hat{x}_k) =
\sum_{i=0}^n c_i \sum_{k=0}^{n} \Bigg[T_i(\hat{x}_k)T_j(\hat{x}_k)\Bigg]
$$

and then apply the discrete orthogonality condition to obtain

$$
\begin{aligned}
c_0 &= \frac{1}{n}\sum_{k=0}^{n} f(\hat{x}_k) \\
c_i &= \frac{2}{n}\sum_{k=0}^{n} f(\hat{x}_k)T_i(\hat{x}_k)
\end{aligned}
$$

Now we can write the approximation as

$$
f(x) \approx \underbrace{\frac{1}{n}\sum_{k=0}^{n} f(\hat{x}_k)}_{c_0} +
\sum_{i=0}^n\underbrace{\frac{2}{n}\Bigg(\sum_{k=0}^{n}T_i(\hat{x}_k)f(\hat{x}_k)\Bigg)}_{c_i}T_i(x)
$$

> chebPolFit :: Vector Double -> (Double -> Double) -> Vector Double
> chebPolFit zs f = (sum ys / n') `cons` (unfoldrN (n - 1) chebUnfoldAux (replicate n 1, zs))
>   where
>     n = length zs
>     n' = fromIntegral n
>     ys = map f zs
>     chebUnfoldAux (t1, t2) =
>       Just ((sum (t2 .* ys)) * 2 / n',
>             (t2, zipWith3 (\x a b -> 2 * x * b - a) zs t1 t2))

It is tempting to calculate the Chebyshev points (of the first kind) directly.

> chebPoints' :: Floating a => Int -> Vector a
> chebPoints' n = map f (enumFromN 0 (n + 1))
>   where
>     n' = fromIntegral n
>     f k = cos (pi * k' / n')
>       where
>         k' = fromIntegral k

But this has some unfortunate consequences

    [ghci]
    chebPoints' 2

It is better to use a method which produces symmetric values.

> chebPoints :: Floating a => Int -> Vector a
> chebPoints n = generate (n + 1) f
>   where
>     f i = sin (pi * fromIntegral (2 * i - n) / fromIntegral (2 * n))

    [ghci]
    chebPoints 2

> bigN :: Int
> bigN = 10

> x :: Vector (Complex Double)
> x = map realToFrac $ chebPoints bigN

> x' :: Vector (Complex Double)
> x' = generate
>      (bigN + 1)
>      (\i -> realToFrac $ sin (pi * fromIntegral (2 * i - bigN) / fromIntegral (2 * bigN)))

> v :: Vector (Complex Double)
> v = map (\x -> exp x * sin ( 5 * x)) x

> bigV :: Vector (Complex Double)
> bigV = v ++ (reverse $ slice 1 (bigN - 1) v)

> ii, jj, kk :: Vector (Complex Double)
> ii = generate bigN fromIntegral
> jj = cons 0 (generate (bigN - 1) (\i -> fromIntegral (i + 1 - bigN)))
> kk = ii ++ jj

> i1 :: Complex Double
> i1 = 0.0 :+ 1.0

> test :: IO ()
> test = do
>   bigU <- fft bigV
>   let bigU'  = map (\z -> (realPart z :+ 0)) bigU
>       preInv = zipWith (*) (replicate (2 * bigN) i1) (zipWith (*) kk bigU')
>   bigW <- ifft preInv
>   let bigW' = map (\z -> (realPart z :+ 0)) bigW
>   -- putStrLn $ show $ bigU'
>   -- putStrLn $ show preInv
>   putStrLn $ show bigW'

> chebPoly :: Int -> Double -> Double
> chebPoly 0 _ = 1
> chebPoly 1 x = x
> chebPoly n x = 2 * x * chebPoly (n - 1) x - chebPoly (n - 2) x

> chebUnfold :: Int -> Double -> Vector Double
> chebUnfold n x = unfoldrN n chebUnfoldAux (1, x)
>   where
>     chebUnfoldAux (a, b) = Just (a, (b, 2 * x * b - a))


> chebPolVal :: Vector Double -> Vector Double -> Vector Double
> chebPolVal cs xs = u .- xs .* ujp1
>   where
>     (u, ujp1) = us!(n-3)
>     us        = unfoldr chebPolValAux (initU, initUjp1, ds)
>
>     initUjp1 = replicate m (cs!(n-1))
>     initU    = map ((cs!(n-2))+) $ map (2*cs!(n-1)*) xs
>     ds       = V.init (V.init cs)
>
>     n        = length cs
>     m        = length xs
>
>     chebPolValAux (_, _, c)
>       | V.null c = Nothing
>     chebPolValAux (ujp1, ujp2, c)
>       | otherwise = Just ((u, ujp1), (u, ujp1, V.init c))
>       where
>         u = map ((V.last c)+) $
>             zipWith (-) (map (2*) $ zipWith (*) xs ujp1) ujp2

> witchOfAgnesi :: Double -> Double
> witchOfAgnesi x = 1 / (1+25 * x^2)


> chebZeros :: Int -> Vector Double
> chebZeros n = map f (enumFromN 0 n)
>   where
>     n' = fromIntegral (n - 1)
>     f k = cos (pi * (2 * k' + 1) / (2 * n' + 2))
>       where
>         k' = fromIntegral k


> chebFft :: Vector a -> (a -> Complex Double) -> IO (Vector Double)
> chebFft ps f = do
>   fcs <- fft $ map f ps'
>   return $ map realPart $ zipWith (*) ns (V.take l fcs)
>   where
>     ns = recip2l `cons` (replicate (l - 2) (1 / l')) `snoc` recip2l
>     recip2l = 1 / (2 * l')
>     l = length ps
>     l' = fromIntegral l - 1
>     ps' = ps ++ reverse (slice 1 (l - 2) ps)

> testZeros :: Int -> Vector Double
> testZeros n = map (\z -> (chebUnfold n z)!(n - 1)) (chebZeros (n - 1))

    [ghci]
    testZeros 5

```{.dia width='800'}
import ChebDiff
import ChebDiag
import Data.Vector ( toList )

dia = diagChebFit (toList x1) (toList y) witchOfAgnesi
````

> infixr 8 ^*
> (^*) :: Num a => a -> Vector a -> Vector a
> s ^* v = map (* s) v

> infixl 7 .*
> (.*) :: Num a => Vector a -> Vector a -> Vector a
> (.*) = zipWith (*)

> infixl 6 .-
> (.-) :: Num a => Vector a -> Vector a -> Vector a
> (.-) = zipWith (-)



To Be Moved
===========

> x1, c, y :: Vector Double
> x1 = generate 201 (\n -> (fromIntegral n - 100) / 100)
> c = chebPolFit (chebZeros 11) witchOfAgnesi
> y = chebPolVal c x1


https://hackage.haskell.org/package/polynomial

http://www.chebfun.org/docs/guide/guide08.html

