module Bisection where

import Counted
import Utils

bisection :: (Fractional a, Ord a) => (a -> a) -> a -> (a, a) -> Counted a
bisection f delta (a, b)
  | abs (b - a) <= 2 * delta = pure c
  | f c == 0 = pure c
  | signum (f a) == signum (f c) = step $ bisection f delta (c, b)
  | otherwise = step $ bisection f delta (a, c)
  where
    c = midpoint (a, b)
