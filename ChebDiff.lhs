
It is well known that one can represent a sufficiently well-behaved
function on the interval $[-\pi, \pi]$ by its Fourier expansion. In
the case of a symmetric function, $f(-x) = f(x)$, this becomes

$$
f(\theta) = \sum_{n=0}^\infty a_n \cos{n\theta}
$$

If we make a change of variable, $x = \cos{\theta}$ we obtain

$$
g(x) = \sum_{n=0}^\infty a_n T_n(x)
$$

where $T_n(x) = \cos{n\theta}$ is the $n$-th Chebyshev polynomial and
$g$ is an aribtrary (sufficiently well-behaved) function on $[-1, 1]$.

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

> {-# OPTIONS_GHC -Wall                      #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults    #-}
> {-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
> {-# OPTIONS_GHC -fno-warn-missing-methods  #-}
> {-# OPTIONS_GHC -fno-warn-orphans          #-}

> module Main where

> import Prelude hiding ( length, sum, zipWith, map, (++), reverse, drop, replicate )
> import qualified Prelude as P
> import Data.Complex
> import Data.Vector hiding ( tail )
> import qualified Data.Vector as V
> import Numeric.FFT
> import Debug.Trace

> import qualified Control.Lens as L
> import qualified Graphics.Rendering.Chart as C
> import Graphics.Rendering.Chart.Backend.Diagrams
> import Diagrams.Backend.Cairo.CmdLine
> import Diagrams.Prelude hiding ( render, Renderable )
> import Data.Default.Class

> import Diagrams.Prelude
> import Diagrams.Backend.CmdLine
> import Diagrams.Backend.Cairo.CmdLine

> import Text.Printf

> import System.IO.Unsafe

> bigN = 10

> x :: Vector (Complex Double)
> x = generate
>     (bigN + 1)
>     (\i -> realToFrac $ cos (pi * fromIntegral i / fromIntegral bigN))

> v :: Vector (Complex Double)
> v = map (\x -> exp x * sin ( 5 * x)) x

> bigV = v ++ (reverse $ slice 1 (bigN - 1) v)

> ii, jj, kk :: Vector (Complex Double)
> ii = generate bigN fromIntegral
> jj = cons 0 (generate (bigN - 1) (\i -> fromIntegral (i + 1 - bigN)))
> kk = ii ++ jj
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

> chebPoly :: (Num b, Num a, Eq a) => a -> b -> b
> chebPoly 0 _ = 1
> chebPoly 1 x = x
> chebPoly n x = 2 * x * chebPoly (n - 1) x - chebPoly (n - 2) x

> chebUnfold :: Int -> Double -> Vector Double
> chebUnfold n x = unfoldrN n chebUnfoldAux (1, x)
>   where
>     chebUnfoldAux (a, b) = Just (a, (b, 2 * x * b - a))

> chart = C.toRenderable layout
>   where
>     sinusoid n c = C.plot_lines_values .~ [[ (x, (chebUnfold n x)!(n - 1)) | x <- [-1.0,(-0.99)..1.0]]]
>                  $ C.plot_lines_style  . C.line_color .~ opaque c
>                  $ C.plot_lines_title .~ ("n = " P.++ show n)
>                  $ def
> 
>     layout = C.layout_title .~ "ChebyShev Polynomials"
>            $ C.layout_plots .~ [C.toPlot (sinusoid 3 blue),
>                                 C.toPlot (sinusoid 5 green),
>                                 C.toPlot (sinusoid 7 red)
>                                ]
>            $ def

> denv :: DEnv
> denv = unsafePerformIO $ defaultEnv C.vectorAlignmentFns 500 500

> diag :: Diagram B R2
> diag =
>   fst $ runBackend denv (C.render chart (500, 500))

> main = do
>   displayHeader "diagrams/ChebyShev.svg" diag
>   putStrLn "Hello"
> 
>

> displayHeader :: FilePath -> Diagram B R2 -> IO ()
> displayHeader fn =
>   mainRender ( DiagramOpts (Just 900) (Just 700) fn
>              , DiagramLoopOpts False Nothing 0
>              )
