{-# LANGUAGE DataKinds #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Ode
  ( OdeSystem (..),
    eulerExplicit,
    eulerImplicit,
    rungeKuttaGeneric,
    rungeKutta4,
  )
where

import Control.Monad.State
import PeanoNat
import Vec
import VecMath

-- | dy/dy = f(t, y)
-- | y(t0) = y0
data OdeSystem (n :: PeanoNat) a
  = OdeSystem
      -- | f: R^(n+1) -> R^n
      (a -> Vec n a -> Vec n a)
      -- | t0
      a
      -- | y0
      (Vec n a)

rungeKuttaGeneric ::
  (Traversable f, Num a) =>
  -- | Algorithm step
  ((a -> Vec n a -> Vec n a) -> a -> State (a, Vec n a) (Vec n a)) ->
  -- | System to solve
  OdeSystem n a ->
  -- | Grid
  f a ->
  f (Vec n a)
rungeKuttaGeneric step (OdeSystem f t0 y0) grid = evalState (traverse (step f) grid) (t0, y0)

eulerExplicit :: forall f a n. (Traversable f, Num a) => OdeSystem n a -> f a -> f (Vec n a)
eulerExplicit = rungeKuttaGeneric step
  where
    step :: (Num a) => (a -> Vec n a -> Vec n a) -> a -> State (a, Vec n a) (Vec n a)
    step f tn1 = do
      (tn, yn) <- get
      let hn = tn1 - tn
      let yn1 = yn `vadd` scale hn (f tn yn)
      put (tn1, yn1)
      pure yn1

eulerImplicit ::
  forall f a n.
  (Traversable f, Floating a, Ord a) =>
  ((Vec n a -> Vec n a) -> Vec n a -> Vec n a) -> OdeSystem n a -> f a -> f (Vec n a)
eulerImplicit eqSystemSolver = rungeKuttaGeneric step
  where
    step :: (Num a) => (a -> Vec n a -> Vec n a) -> a -> State (a, Vec n a) (Vec n a)
    step f tn1 = do
      (tn, yn) <- get
      let hn = tn1 - tn
      let yn1 = eqSystemSolver (\y -> yn `vadd` scale hn (f tn1 y) `vsub` y) yn
      put (tn1, yn1)
      pure yn1

rungeKutta4 :: forall f a n. (Traversable f, Fractional a) => OdeSystem n a -> f a -> f (Vec n a)
rungeKutta4 = rungeKuttaGeneric step
  where
    step :: (Num a) => (a -> Vec n a -> Vec n a) -> a -> State (a, Vec n a) (Vec n a)
    step f tn1 = do
      (tn, yn) <- get
      let hn = tn1 - tn
      let k1 = f tn yn
      let k2 = f (tn + hn / 2) (yn `vadd` scale (hn / 2) k1)
      let k3 = f (tn + hn / 2) (yn `vadd` scale (hn / 2) k2)
      let k4 = f tn (yn `vadd` scale hn k3)
      let yn1 = yn `vadd` scale (hn / 6) (k1 `vadd` scale 2 k2 `vadd` scale 2 k3 `vadd` k4)
      put (tn1, yn1)
      pure yn1
