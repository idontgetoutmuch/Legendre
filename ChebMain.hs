import ChebDiag

main :: IO ()
main = do
  displayHeader "diagrams/Chebyshev.svg" (diag undefined)
  displayHeader "diagrams/ChebyInterp.svg" (diagChebFit undefined undefined undefined)
  putStrLn "Hello"

