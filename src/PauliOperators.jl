PauliX = sparse([0.0 1.0; 1.0 0.0])
PauliY = sparse([0.0 -im; im 0.0])
PauliZ = sparse([1.0 0.0; 0.0 -1.0])
density_operator = sparse([0.0 0.0; 0.0 1.0])
PauliOperators = [PauliX, PauliY, PauliZ]

Sigma_minus = sparse([0.0 1.0; 0.0 0.0])
Sigma_plus = sparse([0.0 0.0; 1.0 0.0])
