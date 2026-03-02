module Counted where

import Control.Monad.State

type Counted = State Int

runCounted :: Counted a -> (a, Int)
runCounted c = runState c 0

step :: Counted a -> Counted a
step x = modify (+ 1) *> x
