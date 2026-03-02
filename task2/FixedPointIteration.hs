module FixedPointIteration where
import Counted

fixedPointIteration :: (Fractional a, Ord a) => (a -> a) -> a -> a -> a -> Counted a
fixedPointIteration f lambda delta x0
  | abs (x1 - x0) <= delta = pure x1
  | otherwise = step $ fixedPointIteration f lambda delta x1
  where
    x1 = x0 - lambda * f x0
