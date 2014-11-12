
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

> chebPolVal :: Vector Double -> Vector Double -> Vector (Vector Double)
> chebPolVal cs xs =
>   unfoldr chebPolValAux (initUjp2, initUjp1, initUjp1, V.init (V.init cs), 2)
>   where
>     initUjp1 = replicate m (cs!(n-1))
>     initUjp2 = map (+cs!(n-2)) $ map (2*cs!(n-1)*) xs
>     n = length cs
>     m = length xs
>     chebPolValAux (_, _, _, c, _) | V.null c =
>       Nothing
>     chebPolValAux (u', ujp1', ujp2', c, n) | otherwise =
>       Just (u, (u, u', ujp1, V.init c, n - 1))
>       where
>         ujp2 = ujp1'
>         ujp1 = u'
>         u = map ((V.last c)+) $ foo
>         foo = zipWith (-) (map (2*) $ zipWith (*) xs ujp1) ujp2

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

> y :: Int -> (Double -> Double) -> Vector Double
> y n f = map f (chebZeros n)

> c0 :: Int -> (Double -> Double) -> Double
> c0 n f = 1 / (fromIntegral n + 1) *
>          sum (y (n + 1) f .* (map ((!0) . (chebUnfold 3)) (chebZeros (n + 1))))

> c1 :: Int -> (Double -> Double) -> Double
> c1 n f = 2 / (fromIntegral n + 1) *
>          sum (y (n + 1) f .* (map ((!1) . (chebUnfold 3)) (chebZeros (n + 1))))

> c2 :: Int -> (Double -> Double) -> Double
> c2 n f = 2 / (fromIntegral n + 1) *
>          sum (y (n + 1) f .* (map ((!2) . (chebUnfold 3)) (chebZeros (n + 1))))

> infixl 7 ^*
> (^*) :: Num a => a -> Vector a -> Vector a
> s ^* v = map (* s) v

> infixl 7 .*
> (.*) :: Num a => Vector a -> Vector a -> Vector a
> (.*) = zipWith (*)

> chart :: C.Renderable ()
> chart = C.toRenderable layout
>   where
>     sinusoid n c = C.plot_lines_values .~ [[ (x, (chebUnfold n x)!(n - 1)) | x <- [-1.0,(-0.99)..1.0]]]
>                  $ C.plot_lines_style  . C.line_color .~ opaque c
>                  $ C.plot_lines_title .~ ("n = " P.++ show n)
>                  $ def
> 
>     layout = C.layout_title .~ "Chebyshev Polynomials"
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

> main :: IO ()
> main = do
>   displayHeader "diagrams/Chebyshev.svg" diag
>   putStrLn "Hello"
> 
>

> displayHeader :: FilePath -> Diagram B R2 -> IO ()
> displayHeader fn =
>   mainRender ( DiagramOpts (Just 900) (Just 700) fn
>              , DiagramLoopOpts False Nothing 0
>              )
