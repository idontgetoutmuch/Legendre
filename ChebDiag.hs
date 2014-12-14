{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

module ChebDiag where

import qualified Graphics.Rendering.Chart as C
import Graphics.Rendering.Chart.Backend.Diagrams
import Diagrams.Prelude hiding ( render, Renderable, trace )
import Data.Default.Class

import Diagrams.Backend.Cairo.CmdLine
import Diagrams.Backend.CmdLine

import System.IO.Unsafe


displayHeader :: FilePath -> Diagram B R2 -> IO ()
displayHeader fn =
  mainRender ( DiagramOpts (Just 900) (Just 700) fn
             , DiagramLoopOpts False Nothing 0
             )

chartChebFit :: [Double] -> [Double] -> (Double -> Double) -> C.Renderable ()
chartChebFit x y f = C.toRenderable layout
  where
    z = zip x y
    fit = C.plot_lines_values .~ [z]
          $ C.plot_lines_style  . C.line_color .~ opaque blue
          $ C.plot_lines_title .~ ("Interpolation")
          $ def
    org = C.plot_lines_values .~ [[ (x, f x) | x <- [-1.0,(-0.99)..1.0]]]
          $ C.plot_lines_style  . C.line_color .~ opaque red
          $ C.plot_lines_title .~ ("Original")
          $ def

    layout = C.layout_title .~ "Chebyshev Interpolation"
           $ C.layout_plots .~ [C.toPlot fit,
                                C.toPlot org
                               ]
           $ def

diagChebFit :: [Double] -> [Double] -> (Double -> Double) -> Diagram B R2
diagChebFit x y f =
  fst $ runBackend denv (C.render (chartChebFit x y f) (500, 500))

chart :: (Int -> Double -> [Double]) -> C.Renderable ()
chart chebUnfold = C.toRenderable layout
  where
    cheby n c = C.plot_lines_values .~ [[ (x, (chebUnfold n x)!!(n - 1)) | x <- [-1.0,(-0.99)..1.0]]]
                 $ C.plot_lines_style  . C.line_color .~ opaque c
                 $ C.plot_lines_title .~ ("n = " ++ show n)
                 $ def

    layout = C.layout_title .~ "Chebyshev Polynomials"
           $ C.layout_y_axis . C.laxis_generate .~ C.scaledAxis def (-1,1)
           $ C.layout_x_axis . C.laxis_generate .~ C.scaledAxis def (-1,1)
           $ C.layout_plots .~ [C.toPlot (cheby 3 blue),
                                C.toPlot (cheby 5 green),
                                C.toPlot (cheby 7 red)
                               ]
           $ def

denv :: DEnv
denv = unsafePerformIO $ defaultEnv C.vectorAlignmentFns 500 500

diag :: (Int -> Double -> [Double]) -> Diagram B R2
diag chebUnfold =
  fst $ runBackend denv (C.render (chart chebUnfold) (500, 500))
