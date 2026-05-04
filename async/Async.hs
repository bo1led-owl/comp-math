{-# LANGUAGE GADTs #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# OPTIONS_GHC -Wno-unused-do-bind #-}

module Async
  ( Async,
    Future,
    runAsync,
    runAsync_,
    async,
    async_,
    await,
    mapAsync,
    mapAsync_,
    forAsync,
    forAsync_,
  )
where

import Control.Concurrent
import Control.Exception
import Control.Monad
import Control.Monad.IO.Class
import Data.Foldable

newtype Async a = Async (ThreadPool -> IO a)

instance Functor Async where
  fmap f (Async c) = Async (fmap f . c)

instance Applicative Async where
  pure = Async . const . pure
  (Async f) <*> (Async g) = Async $ \t -> f t <*> g t

instance Monad Async where
  (Async c) >>= f = Async $ \t -> do
    x <- c t
    let (Async c2) = f x
    c2 t

instance MonadIO Async where
  liftIO = Async . const

data Tasks where
  Nil :: Tasks
  Cons :: Future a -> Tasks -> Tasks

newtype Future a = Future (MVar a)

awaitIO :: Future a -> IO a
awaitIO (Future mvar) = readMVar mvar

await :: Future a -> Async a
await future = liftIO $ awaitIO future

newtype ThreadPool = ThreadPool (MVar Tasks)

threadPool :: IO ThreadPool
threadPool = ThreadPool <$> newMVar Nil

stopThreadPool :: ThreadPool -> IO ()
stopThreadPool tp@(ThreadPool threads) = do
  tasks <- liftIO $ takeMVar threads
  case tasks of
    Nil -> pure ()
    Cons t ts -> do
      liftIO $ putMVar threads ts
      awaitIO t
      stopThreadPool tp

async :: IO a -> Async (Future a)
async task =
  Async $ \(ThreadPool tasks) -> liftIO $ do
    future@(Future mvar) <- Future <$> newEmptyMVar
    modifyMVar_ tasks (pure . (future `Cons`))
    forkIO (task >>= putMVar mvar)
    pure future

async_ :: IO a -> Async ()
async_ = void . async . void

runAsync :: Async a -> IO a
runAsync (Async comp) = do
  tp <- threadPool
  (res :: Either SomeException a) <- try $ comp tp
  stopThreadPool tp
  either throw pure res

runAsync_ :: Async a -> IO ()
runAsync_ = void . runAsync

mapAsync :: (Traversable t) => (a -> IO b) -> t a -> IO (t b)
mapAsync f items = runAsync (traverse (\x -> async (f x) >>= await) items)

mapAsync_ :: (Traversable t) => (a -> IO b) -> t a -> IO ()
mapAsync_ f items = runAsync (traverse_ (\x -> async (f x) >>= await) items)

forAsync :: (Traversable t) => t a -> (a -> IO b) -> IO (t b)
forAsync = flip mapAsync

forAsync_ :: (Traversable t) => t a -> (a -> IO b) -> IO ()
forAsync_ = flip mapAsync_
