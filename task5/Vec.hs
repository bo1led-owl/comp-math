{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UndecidableInstances #-}
{-# OPTIONS_GHC -Wno-name-shadowing #-}

module Vec
  ( Vec (..),
    SqMat,
    Mat,
    zipWithVec,
    singleMaps,
    append,
    reverseVec,
  )
where

import PeanoNat

data Vec (n :: PeanoNat) a where
  Nil :: Vec Z a
  Cons :: a -> Vec n a -> Vec (S n) a

infixr 5 `Cons`

type Mat (n :: PeanoNat) (m :: PeanoNat) a = Vec n (Vec m a)

type SqMat (n :: PeanoNat) a = Mat n n a

instance Functor (Vec n) where
  fmap _ Nil = Nil
  fmap f (Cons x xs) = Cons (f x) (fmap f xs)

instance Foldable (Vec n) where
  foldr _ acc Nil = acc
  foldr f acc (Cons x xs) = f x (foldr f acc xs)
  foldl _ acc Nil = acc
  foldl f acc (Cons x xs) = foldl f (f acc x) xs

instance Traversable (Vec Z) where
  traverse _ Nil = pure Nil

instance (Traversable (Vec n)) => Traversable (Vec (S n)) where
  traverse f (Cons x xs) = Cons <$> f x <*> traverse f xs

instance (Show a) => Show (Vec n a) where
  show v = "{" ++ foldl (\s x -> if null s then show x else s ++ ", " ++ show x) "" v ++ "}"

zipWithVec :: (a -> b -> c) -> Vec n a -> Vec n b -> Vec n c
zipWithVec _ Nil Nil = Nil
zipWithVec f (Cons x xs) (Cons y ys) = Cons (f x y) (zipWithVec f xs ys)

zipWithVecLossy :: (a -> b -> c) -> Vec n a -> Vec m b -> Vec (Min n m) c
zipWithVecLossy _ Nil _ = Nil
zipWithVecLossy _ _ Nil = Nil
zipWithVecLossy f (x `Cons` xs) (y `Cons` ys) = f x y `Cons` zipWithVecLossy f xs ys

singleMaps :: (a -> a) -> Vec n a -> Vec n (Vec n a)
singleMaps _ Nil = Nil
singleMaps f (x `Cons` xs) = (f x `Cons` xs) `Cons` fmap (x `Cons`) (singleMaps f xs)

getIth :: Int -> Vec n a -> a
getIth _ Nil = error "empty vec"
getIth 0 (x `Cons` _) = x
getIth i (_ `Cons` xs) = getIth (i - 1) xs

setIth :: Int -> a -> Vec n a -> Vec n a
setIth _ _ Nil = error "empty vec"
setIth 0 y (_ `Cons` xs) = y `Cons` xs
setIth i y (x `Cons` xs) = x `Cons` setIth (i - 1) y xs

getRow :: Int -> Mat m n a -> Vec n a
getRow = getIth

getCol :: Int -> Mat m n a -> Vec m a
getCol _ Nil = error "empty matrix"
getCol i rows = fmap (getIth i) rows

getElem :: Int -> Int -> Mat m n a -> a
getElem i j = getIth j . getIth i

setRow :: Int -> Vec n a -> Mat m n a -> Mat m n a
setRow = setIth

setCol :: Int -> Vec m a -> Mat m n a -> Mat m n a
setCol i = zipWithVec (setIth i)

setElem :: Int -> Int -> a -> Mat m n a -> Mat m n a
setElem i j x mat = setRow i (setIth j x (getRow i mat)) mat

botRightSubMatrix :: SqMat (S n) a -> SqMat n a
botRightSubMatrix (_ `Cons` rows) = fmap (\(_ `Cons` xs) -> xs) rows

append :: a -> Vec m a -> Vec (S m) a
append y Nil = y `Cons` Nil
append y (z `Cons` zs) = z `Cons` append y zs

reverseVec :: Vec n a -> Vec n a
reverseVec Nil = Nil
reverseVec (x `Cons` xs) = append x (reverseVec xs)
