import BVP2
import DFT
import Data.Complex
import Data.Foldable
import Safe
import Utils

pairwise :: (a -> a -> b) -> [a] -> [b]
pairwise f xs = zipWith f xs (tailSafe xs)

bvp2 :: BVP2 Double
bvp2 = BVP2 (const 1) (const 0) (const 0) sin

exact :: Double -> Double
exact x = -(sin x)

conditions :: [(String, BCond Double, BCond Double)]
conditions =
  -- BCond x A B D
  -- Au'(x) + Bu(x) = D
  [ ("function left-right", BCond 0 0 1 0, BCond pi 0 1 0),
    ("function left, derivative right", BCond 0 0 1 0, BCond pi 1 0 1),
    ("derivative left, function right", BCond 0 1 0 (-1), BCond pi 0 1 0)
  ]

main :: IO ()
main = do
  printOrdersForThomas
  printErrorsForDft

printOrdersForThomas :: IO ()
printOrdersForThomas = do
  putStrLn "thomas"
  for_
    conditions
    ( \(name, cl, cr) -> do
        putStrLn $ name ++ ":"
        putStrLn "n\torder"
        traverse_
          (\(n, order) -> putStrLn $ show n ++ "\t" ++ show order)
          (ordersOfConvergence (solveBVP2 bvp2 cl cr) exact)
        putStrLn ""
    )

ordersOfConvergence :: (Int -> [(Double, Double)]) -> (Double -> Double) -> [(Int, Double)]
ordersOfConvergence solver solution =
  let gridSizes = [initialSize * factor ^ i | i <- [0 .. (nSteps - 1)]]
      errors = fmap (\(x, y) -> abs (solution x - y)) . solver <$> gridSizes
      maxErrors = maximum <$> errors
      orders = pairwise (\prevErr curErr -> logBase (fromIntegral factor) (prevErr / curErr)) maxErrors
   in tailSafe gridSizes `zip` orders
  where
    factor = 2
    nSteps :: Int
    nSteps = 5
    initialSize = 16

solveDft :: (RealFloat a) => (a -> a) -> Int -> [(a, a)]
solveDft rhs n =
  let fs = dft period n ((:+ 0) . rhs <$> grid)
      ys = 0 : tailSafe (zipWith yk fs [0 .. (n - 1)])
   in grid `zip` (realPart <$> idft period n ys)
  where
    period = 2 * pi
    yk fk k
      | k <= n `div` 2 = -(fk / fromIntegral (square k))
      | otherwise = -(fk / fromIntegral (square (n - k)))
    grid = [fromIntegral i * (period / fromIntegral n) | i <- [0 .. (n - 1)]]

printErrorsForDft :: IO ()
printErrorsForDft = do
  putStrLn "DFT"
  putStrLn "n\tmax error"
  for_
    [initialSize * factor ^ i | i <- [0 .. (nSteps - 1)]]
    ( \n ->
        let solution = solveDft sin n
            errors = (\(x, y) -> abs (exact x - y)) <$> solution
            maxError = maximum errors
         in putStrLn $ show n ++ "\t" ++ show maxError
    )
  where
    initialSize = 16
    factor = 2
    nSteps :: Int
    nSteps = 5
