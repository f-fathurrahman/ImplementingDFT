using LinearAlgebra

Haux = rand(5,5)
Haux = 0.5*(Haux + Haux')

# Eigenvalue decomposition
λ, Urot = eigen(Haux)

# Check
D = diagm(0 => λ)
# Urot' should be equal to Urot^(-1) or inv(Urot)
display(Urot * D * Urot' - Haux); println()
