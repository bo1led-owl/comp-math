module Utils where

midpoint :: (Fractional a) => (a, a) -> a
midpoint (a, b) = (a + b) / 2

square :: (Num a) => a -> a
square x = x * x
