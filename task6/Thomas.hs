module Thomas (thomas) where

thomas :: (Fractional a, Show a) => [a] -> [a] -> [a] -> [a] -> [a]
thomas a b c d =
  let aExt = 0 : a
      cExt = c ++ [0]
   in substitute (makeAlphaBeta aExt b cExt d)

substitute :: (Fractional a) => [(a, a)] -> [a]
substitute [] = undefined
substitute [(_, beta)] = [beta]
substitute ((alpha, beta) : rest) =
  let restRes = substitute rest
   in beta - alpha * head restRes : restRes

makeAlphaBeta :: (Fractional a) => [a] -> [a] -> [a] -> [a] -> [(a, a)]
makeAlphaBeta a b c d = makeAlphaBeta' a b c d 0 0

makeAlphaBeta' :: (Fractional a) => [a] -> [a] -> [a] -> [a] -> a -> a -> [(a, a)]
makeAlphaBeta' (a : as) (b : bs) (c : cs) (d : ds) prevAlpha prevBeta =
  let denom = b - a * prevAlpha
      curAlpha = c / denom
      curBeta = (d - a * prevBeta) / denom
   in (curAlpha, curBeta) : makeAlphaBeta' as bs cs ds curAlpha curBeta
makeAlphaBeta' _ _ _ _ _ _ = []
