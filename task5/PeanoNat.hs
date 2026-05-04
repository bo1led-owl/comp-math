{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE UndecidableInstances #-}

module PeanoNat (PeanoNat (..), ToPeano, Add, Min, KnownNat (..)) where

import GHC.TypeLits hiding (KnownNat)

data PeanoNat = Z | S PeanoNat

type family ToPeano (n :: Natural) :: PeanoNat where
  ToPeano 0 = Z
  ToPeano n = S (ToPeano (n - 1))

type family Add (n :: PeanoNat) (m :: PeanoNat) :: PeanoNat where
  Add Z m = m
  Add (S n) m = S (Add n m)

type family Min (n :: PeanoNat) (m :: PeanoNat) :: PeanoNat where
  Min Z _ = Z
  Min _ Z = Z
  Min (S n) (S m) = S (Min n m)

class KnownNat (n :: PeanoNat) where
  toInt :: Int

instance KnownNat Z where
  toInt = 0

instance (KnownNat n) => KnownNat (S n) where
  toInt = 1 + toInt @n
