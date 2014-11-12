
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

> module ChebDiff where

> import Prelude hiding ( length, sum, zipWith, zipWith3,
>                         map, (++), reverse, drop, replicate )
> import qualified Prelude as P
> import Data.Complex
> import Data.Vector hiding ( tail )
> import qualified Data.Vector as V
> import Numeric.FFT

> import qualified Graphics.Rendering.Chart as C
> import Graphics.Rendering.Chart.Backend.Diagrams
> import Diagrams.Backend.Cairo.CmdLine
> import Diagrams.Prelude hiding ( render, Renderable, trace )
> import Data.Default.Class

> import Diagrams.Backend.CmdLine

> import System.IO.Unsafe


> bigN :: Int
> bigN = 10

> x :: Vector (Complex Double)
> x = generate
>     (bigN + 1)
>     (\i -> realToFrac $ cos (pi * fromIntegral i / fromIntegral bigN))

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

> chebPolFit :: Vector Double -> (Double -> Double) -> Vector Double
> chebPolFit zs f = (sum ys / n') `cons` (unfoldrN (n - 1) chebUnfoldAux (replicate n 1, zs))
>   where
>     n = length zs
>     n' = fromIntegral n
>     ys = map f zs
>     chebUnfoldAux (t1, t2) =
>       Just ((sum (t2 .* ys)) * 2 / n',
>             (t2, zipWith3 (\x a b -> 2 * x * b - a) zs t1 t2))

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

> f :: Double -> Double
> f x = 1 / (1+25 * x^2)

```{.dia width='800'}
import ChebDiff

dia = diag
````

> chebZeros :: Int -> Vector Double
> chebZeros n = map f (enumFromN 0 n)
>   where
>     n' = fromIntegral (n - 1)
>     f k = cos (pi * (2 * k' + 1) / (2 * n' + 2))
>       where
>         k' = fromIntegral k

> testZeros :: Int -> Vector Double
> testZeros n = map (\z -> (chebUnfold n z)!(n - 1)) (chebZeros (n - 1))

    [ghci]
    testZeros 5

```{.dia width='800'}
import ChebDiff

dia = diagChebFit
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

> chartChebFit :: C.Renderable ()
> chartChebFit = C.toRenderable layout
>   where
>     c = chebPolFit (chebZeros 11) f
>     x = generate 201 (\n -> (fromIntegral n - 100) / 100)
>     y = chebPolVal c x
>     z = toList $ V.zip x y
>     fit = C.plot_lines_values .~ [z]
>           $ C.plot_lines_style  . C.line_color .~ opaque blue
>           $ C.plot_lines_title .~ ("Interpolation")
>           $ def
>     org = C.plot_lines_values .~ [[ (x, f x) | x <- [-1.0,(-0.99)..1.0]]]
>           $ C.plot_lines_style  . C.line_color .~ opaque red
>           $ C.plot_lines_title .~ ("Original")
>           $ def
>
>     layout = C.layout_title .~ "Chebyshev Interpolation"
>            $ C.layout_plots .~ [C.toPlot fit,
>                                 C.toPlot org
>                                ]
>            $ def

> diagChebFit :: Diagram B R2
> diagChebFit =
>   fst $ runBackend denv (C.render chartChebFit (500, 500))

> chart :: C.Renderable ()
> chart = C.toRenderable layout
>   where
>     cheby n c = C.plot_lines_values .~ [[ (x, (chebUnfold n x)!(n - 1)) | x <- [-1.0,(-0.99)..1.0]]]
>                  $ C.plot_lines_style  . C.line_color .~ opaque c
>                  $ C.plot_lines_title .~ ("n = " P.++ show n)
>                  $ def
>
>     layout = C.layout_title .~ "Chebyshev Polynomials"
>            $ C.layout_plots .~ [C.toPlot (cheby 3 blue),
>                                 C.toPlot (cheby 5 green),
>                                 C.toPlot (cheby 7 red)
>                                ]
>            $ def

> denv :: DEnv
> denv = unsafePerformIO $ defaultEnv C.vectorAlignmentFns 500 500

> diag :: Diagram B R2
> diag =
>   fst $ runBackend denv (C.render chart (500, 500))

> main :: IO ()
> main = do
>   displayHeader "diagrams/Chebyshev.svg" diag
>   displayHeader "diagrams/ChebyInterp.svg" diagChebFit
>   putStrLn "Hello"

> displayHeader :: FilePath -> Diagram B R2 -> IO ()
> displayHeader fn =
>   mainRender ( DiagramOpts (Just 900) (Just 700) fn
>              , DiagramLoopOpts False Nothing 0
>              )
