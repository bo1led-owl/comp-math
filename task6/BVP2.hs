module BVP2
  ( BVP2 (..),
    BCond (..),
    solveBVP2,
  )
where

import Safe
import Thomas
import Utils

data BVP2 a = BVP2 (a -> a) (a -> a) (a -> a) (a -> a)

data BCond a = BCond !a !a !a !a

makeGrid :: (Num a, Ord a) => a -> a -> a -> [a]
makeGrid l r tau = takeWhile (< r) (iterate (+ tau) l)

solveBVP2 :: (Fractional a, Ord a, Show a) => BVP2 a -> BCond a -> BCond a -> Int -> [(a, a)]
solveBVP2 (BVP2 g h s f) (BCond l a b d) (BCond r a' b' d') gridSize =
  let tau = (r - l) / fromIntegral gridSize
      midGrid = makeGrid (l + tau) (r - tau) tau
      subdiag =
        ((\x -> g (x - tau / 2) / square tau - h x / (2 * tau)) <$> midGrid)
          ++ [a' / tau]
      superdiag =
        (a / tau)
          : ((\x -> g (x + tau / 2) / square tau + h x / (2 * tau)) <$> midGrid)
      diag =
        (a / tau - b)
          : zipWith3 (\an cn sn -> an + cn - sn) (initSafe subdiag) (tailSafe superdiag) (s <$> midGrid)
          ++ [a' / tau + b']
      rhs = d : (f <$> midGrid) ++ [-d']
   in midGrid `zip` thomas subdiag (negate <$> diag) superdiag rhs
