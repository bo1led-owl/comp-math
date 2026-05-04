{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# OPTIONS_GHC -Wno-name-shadowing #-}

import Async
import Control.Monad
import Data.Foldable
import Graphics.Rendering.Chart.Backend.Diagrams
import Graphics.Rendering.Chart.Easy
import Numeric.IEEE (epsilon)
import Ode
import PeanoNat
import Safe
import Vec
import VecMath

grid :: (Fractional a) => a -> a -> Integer -> [a]
grid a b n = [a + fromInteger k * recip (fromInteger n) * (b - a) | k <- [1 .. n]]

data GridLimits a = GridLimits a a Integer Int Integer

pairwise :: (a -> a -> b) -> [a] -> [b]
pairwise f xs = zipWith f xs (tailSafe xs)

orderOfConvergence ::
  -- | Solver to check
  (OdeSystem n Double -> [Double] -> [Vec n Double]) ->
  -- | Test system
  OdeSystem n Double ->
  -- | Answer to the given system
  (Double -> Vec n Double) ->
  -- | Grid limits
  GridLimits Double ->
  -- | List of orders of convergence
  [Double]
orderOfConvergence solver system answer (GridLimits a b initialGridSize steps factor) =
  pairwise (\x y -> logBase (fromInteger factor) (x / y)) maxErrors
  where
    maxError (grid, solution) = maximum (zipWith (\t y -> magnitude (answer t `vsub` y)) grid solution)
    maxErrors = maxError . (\g -> (g, solver system g)) <$> grids
    grids = grid a b <$> gridSizes
    gridSizes = take steps (iterate (* factor) initialGridSize)

order :: (OdeSystem (ToPeano 2) Double -> [Double] -> [Vec (ToPeano 2) Double]) -> [Double]
order solver = orderOfConvergence solver testSystem testAnswer gridLimits
  where
    -- cases(
    --   du/dt = v,
    --   dv/dt = -u,
    --   u(0) = 1,
    --   v(0) = 0,
    -- )
    --
    -- u(t) = cos t, v(t) = -sin t
    testSystem :: OdeSystem (ToPeano 2) Double
    testSystem = OdeSystem f 0 (1 `Cons` 0 `Cons` Nil)
      where
        f :: Double -> Vec (ToPeano 2) Double -> Vec (ToPeano 2) Double
        f _ (u `Cons` v `Cons` Nil) = v `Cons` (-u) `Cons` Nil

    testAnswer :: Double -> Vec (ToPeano 2) Double
    testAnswer t = cos t `Cons` -sin t `Cons` Nil

    gridLimits :: GridLimits Double
    gridLimits =
      GridLimits
        0 -- a
        pi -- b
        8 -- initialGridSize
        6 -- steps
        2 -- factor

fixedPointIteration :: (Floating a, Ord a) => a -> (Vec n a -> Vec n a) -> Vec n a -> Vec n a
fixedPointIteration eps f y =
  if magnitude (y `vsub` y') < eps
    then y
    else fixedPointIteration eps f y'
  where
    y' = f y `vadd` y

printOrders :: Async ()
printOrders =
  for_
    [ ("euler's explicit", eulerExplicit),
      ("euler's implicit", eulerImplicit (fixedPointIteration epsilon)),
      ("runge-kutta 4", rungeKutta4)
    ]
    ( \(name, solver) -> async_ $ do
        let orders = order solver
        putStrLn (name ++ " order of convergence:\t" ++ show (orders))
    )

-- solveStiffOde :: Vec (ToPeano 2) Double -> Double -> Double -> Double -> Async ()
-- solveStiffOde phi t0 a b =
--   for_
--     [ ("stiff, euler's explicit", eulerExplicit),
--       ("stiff, euler's implicit", eulerImplicit (newtons 4 epsilon))
--     ]
--     ( \(name, solver) -> async_ $ do
--         let g = grid a b gridSize
--         let solution = solver system g
--         putStrLn
--           ( name
--               ++ " max error magnitude:\t"
--               ++ show (maximum $ fmap (\(t, s) -> magnitude (exact t `vsub` s)) (g `zip` solution))
--           )
--     )
--   where
--     exact t = 2 * exp (-t) - exp (-(1000 * t)) `Cons` exp (-(1000 * t)) - exp (-t) `Cons` Nil
--     gridSize = 100
--     system =
--       OdeSystem
--         ( \_ (u `Cons` v `Cons` Nil) ->
--             998 * u + 1998 * v `Cons` -999 * u - 1999 * v `Cons` Nil
--         )
--         t0
--         phi

drawHunterVictim :: Double -> Double -> Double -> Async ()
drawHunterVictim t0 l r = async_ $ do
  renderableToFile def "task5/hunter-victim.svg" (toRenderable layout)
  where
    toTuple :: Vec (ToPeano 2) a -> (a, a)
    toTuple (x `Cons` y `Cons` Nil) = (x, y)

    gridSize = 100

    system :: Vec (ToPeano 2) Double -> OdeSystem (ToPeano 2) Double
    system =
      OdeSystem
        ( \_ (x `Cons` y `Cons` _) ->
            10 * x - 2 * x * y `Cons` 2 * x * y - 10 * y `Cons` Nil
        )
        t0

    solution phi = rungeKutta4 (system phi) (grid l r gridSize)

    colors =
      [ opaque blue,
        opaque green,
        opaque red,
        opaque violet
      ] ::
        [AlphaColour Double]

    nPhases = 6

    phase phi color =
      plot_lines_style . line_color .~ color $
        plot_lines_values .~ [fmap toTuple (solution phi)] $
          def
    phis =
      [(let x = (7 / nPhases) * n in x `Cons` x `Cons` Nil) | n <- [1 .. nPhases - 1]]
        ++ [(let x = (7 / nPhases) * n in x / 2 `Cons` x `Cons` Nil) | n <- [1 .. nPhases - 1]]
        ++ [(let x = (7 / nPhases) * n in x `Cons` x / 2 `Cons` Nil) | n <- [1 .. nPhases - 1]]

    phisPlot =
      plot_points_style .~ filledCircles 2 (opaque red) $
        plot_points_values .~ fmap toTuple phis $
          def

    layout =
      layout_plots
        .~ (fmap (toPlot . uncurry phase) (phis `zip` join (repeat colors)) ++ [toPlot phisPlot])
        $ def

main :: IO ()
main = runAsync_ $ do
  printOrders
  -- solveStiffOde (1 `Cons` 0 `Cons` Nil) 0 0 0.1
  drawHunterVictim 0 0 1
