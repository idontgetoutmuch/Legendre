$$
p(x) = \sum_{n=0}^N a_n T_n(x)
$$

$$
P(\theta) = \sum_{n=0}^N a_n \cos n\theta
$$

where $x = \cos\theta$.

FIXME: Also $z + z^{-1}$.

$$
u_j = f(\cos \pi j / N)
$$

For $j = 0,1, \ldots, N$

$$
U_{j+1} = u_j
$$

For $j = 1, 2, \ldots, N - 1 $

$$
U_{2N - j + 1} = u_j
$$

FFT

$$
\hat{U_k} = \frac{\pi}{N}\sum_{j = 1}^{2N} e^{ik\theta}U_j
$$

> module Main where

> import Prelude hiding ( length, sum, zipWith, map, (++), reverse, drop, replicate )
> import Data.Complex
> import Data.Vector
> import Numeric.FFT

> bigN = 10

> x :: Vector (Complex Double)
> x = generate (bigN + 1) (\i -> realToFrac $ cos (pi * fromIntegral i / fromIntegral bigN))

> v :: Vector (Complex Double)
> v = map (\x -> exp x * sin ( 5 * x)) x

> bigV = v ++ (reverse $ slice 1 (bigN - 1) v)

> ii = generate bigN fromIntegral
> jj = cons 0 (generate (bigN - 1) (\i -> fromIntegral (i + 1 - bigN)))
> kk = ii ++ jj
> i1 = 0.0 :+ 1.0

> main :: IO ()
> main = do
>   bigU <- fft bigV
>   let bigU'  = map (\z -> (realPart z :+ 0)) bigU
>       preInv = zipWith (*) (replicate (2 * bigN) i1) (zipWith (*) kk bigU')
>   bigW <- ifft preInv
>   let bigW' = map (\z -> (realPart z :+ 0)) bigW
>   -- putStrLn $ show $ bigU'
>   -- putStrLn $ show preInv
>   putStrLn $ show bigW'
