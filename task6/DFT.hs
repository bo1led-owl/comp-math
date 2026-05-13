module DFT (dft, idft) where

import Data.Complex

dft :: (RealFloat a) => a -> Int -> [Complex a] -> [Complex a]
dft _ _ [] = []
dft period n fks = f <$> [0 .. (n - 1)]
  where
    grid = [fromIntegral k * (period / fromIntegral n) | k <- [0 .. (n - 1)]]
    f j =
      (period :+ 0)
        * sum (zipWith (\tk fk -> fk * exp (0 :+ (omega j * tk))) grid fks)
        / fromIntegral n
    omega j = 2 * pi * fromIntegral j / period

idft :: (RealFloat a) => a -> Int -> [Complex a] -> [Complex a]
idft _ _ [] = []
idft period n fjs = f <$> [0 .. (n - 1)]
  where
    grid = [2 * pi * fromIntegral j / period | j <- [0 .. (n - 1)]]
    f k =
      sum (zipWith (\omegaj fj -> fj * exp (0 :+ (-(omegaj * t k)))) grid fjs)
        / (period :+ 0)
    t k = period * fromIntegral k / fromIntegral n
