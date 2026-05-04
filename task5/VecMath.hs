module VecMath where

import Vec

magnitudeSquared :: (Num a) => Vec n a -> a
magnitudeSquared = sum . fmap (\x -> x * x)

magnitude :: (Floating a) => Vec n a -> a
magnitude = sqrt . sum . fmap (\x -> x * x)

vadd :: (Num a) => Vec n a -> Vec n a -> Vec n a
vadd = zipWithVec (+)

vsub :: (Num a) => Vec n a -> Vec n a -> Vec n a
vsub = zipWithVec (-)

infixl 5 `vadd`

infixl 5 `vsub`

scale :: (Num a) => a -> Vec n a -> Vec n a
scale x = fmap (x *)
