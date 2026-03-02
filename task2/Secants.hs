module Secants where

import Counted

secants :: (Fractional a, Ord a) => (a -> a) -> a -> a -> a -> Counted a
secants f delta x0 x1
  | abs (x0 - x1) <= delta = pure x2
  | otherwise = step $ secants f delta x1 x2
  where
    x2 = x1 - f x1 * (x1 - x0) / (f x1 - f x0)
