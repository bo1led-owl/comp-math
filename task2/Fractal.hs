{-# OPTIONS_GHC -Wno-name-shadowing #-}
{-# OPTIONS_GHC -Wno-unused-do-bind #-}

module Main where

import Codec.Picture
import Data.Complex
import Data.List
import GHC.Float

type C = Complex Double

output :: String
output = "task2/fractal.png"

resolution :: Int
resolution = 1000

gridSize :: Double
gridSize = 2

eps :: Double
eps = 1e-5

closeTo :: C -> C -> Bool
closeTo u v = magnitude (u - v) <= eps

newton :: (C -> C) -> (C -> C) -> C -> C
newton f df x0 = newton' f df x0 cutoff
  where
    cutoff = 75

newton' :: (C -> C) -> (C -> C) -> C -> Int -> C
newton' f df x0 cutoff
  | cutoff <= 0 = x0
  | f x0 == 0 = x0
  | x1 `closeTo` x0 = x1
  | otherwise = newton' f df x1 (cutoff - 1)
  where
    x1 = x0 - f x0 / df x0

f :: C -> C
f z = z * z * z - 1

df :: C -> C
df z = 3 * (z * z)

regions :: [(PixelRGB8, C)]
regions =
  [ PixelRGB8 83 105 170,
    PixelRGB8 68 180 128,
    PixelRGB8 236 113 141
  ]
    `zip` [ 1,
            (-0.5) :+ (sqrt 3 / 2),
            (-0.5) :+ (-(sqrt 3 / 2))
          ]

renderer :: Int -> Int -> PixelRGB8
renderer x y = maybe (PixelRGB8 0 0 0) fst (find (closeTo root . snd) regions)
  where
    root = newton f df point
    point =
      ((-gridSize) + gridStep * int2Double x)
        :+ ((-gridSize) + gridStep * int2Double y)
    gridStep = gridSize * 2 / int2Double resolution

main :: IO ()
main = writePng output $ generateImage renderer resolution resolution
