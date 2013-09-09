{-# OPTIONS_GHC -Wall                     #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
{-# OPTIONS_GHC -fno-warn-type-defaults   #-}

{-# LANGUAGE NoMonomorphismRestriction    #-}
{-# LANGUAGE FlexibleContexts             #-}
{-# LANGUAGE ScopedTypeVariables          #-}

module Initial (
    gConst
  , massesOuter
  , initQsOuter
  , sunMass
  , earthMass
  , jupiterMass
  ) where

import Data.List.Split

gConst :: Double
gConst = 6.67384e-11

sunMass, jupiterMass, earthMass :: Double
sunMass     = 1.9889e30
jupiterMass = 1.8986e27
earthMass   = 5.9722e24

massesOuter :: [Double]
massesOuter =
  [ 9.54786104043e-4
  , 2.85583733151e-4
  , 4.37273164546e-5
  , 5.17759138449e-5
  , 1.0 / 1.3e8
  , 1.00000597682
  ]

initQsOuter :: [[Double]]
initQsOuter = chunksOf 3 xs
  where
    xs = [  -3.5023653
         ,  -3.8169847
         ,  -1.5507963
         ,   9.0755314
         ,  -3.0458353
         ,  -1.6483708
         ,   8.3101420
         , -16.2901086
         ,  -7.2521278
         ,  11.4707666
         , -25.7294829
         , -10.8169456
         , -15.5387357
         , -25.2225594
         ,  -3.1902382
         ,   0.0
         ,   0.0
         ,   0.0
         ]
