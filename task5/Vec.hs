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
