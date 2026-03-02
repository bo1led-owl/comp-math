module Newton where

import Counted

newton :: (Fractional a, Ord a) => (a -> a) -> (a -> a) -> a -> a -> Counted a
newton f df delta x0
  | f x0 == 0 = pure x0
  | abs (x1 - x0) <= delta = pure x1
  | otherwise = step $ newton f df delta x1
  where
    x1 = x0 - f x0 / df x0
