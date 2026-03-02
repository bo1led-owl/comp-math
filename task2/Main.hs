module Main where

import Bisection
import Counted
import Data.Bifunctor
import Data.Foldable
import FixedPointIteration
import Newton
import Secants
import System.Environment
import Utils

f :: (Floating a) => a -> a
f x = tan x - x

df :: (Floating a) => a -> a
df x = recip (square (cos x)) - 1

phiInv :: (Floating a) => a -> a
phiInv = atan

eps :: (Floating a) => a
eps = 1e-5

tangentPeriod :: (Floating a) => Integer -> (a, a)
tangentPeriod k = (mid - (pi / 2), mid + (pi / 2))
  where
    mid = fromInteger k * pi

tangentPeriodDamped :: (Floating a) => Integer -> (a, a)
tangentPeriodDamped k = bimap (+ eps) (+ (-eps)) (tangentPeriod k)

methods :: Integer -> [(String, Counted Double)]
methods k =
  [ ( "bisection",
      bisection f eps dampedPeriod
    ),
    ( "fixed-point iteration",
      fixedPointIteration f (0.3 ^ (abs k + 1)) eps initialGuess
    ),
    ( "newton",
      newton f df eps initialGuess
    ),
    ( "secants",
      secants f eps (initialGuess - eps) (initialGuess + eps)
    )
  ]
  where
    mid = midpoint period
    period = tangentPeriod k
    dampedPeriod = tangentPeriodDamped k
    initialGuess = phiInv mid + mid

showResults :: (Show a) => String -> Counted a -> String
showResults name computation =
  name ++ ": " ++ show res ++ " in " ++ show iterations ++ " steps"
  where
    (res, iterations) = runCounted computation

main :: IO ()
main = do
  (rawPeriodK : _) <- getArgs
  let k = read rawPeriodK
  traverse_ (putStrLn . uncurry showResults) (methods k)
