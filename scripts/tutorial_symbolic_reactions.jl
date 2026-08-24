using ChemistryLab
using DynamicQuantities
import SymPy: symbols, Sym, N, subs, factor, simplify

## construction of a Reaction by a balance calculation with symbolic numbers
a, b, g = symbols("a b g", real = true) ;
C3S, H, CH = [CemSpecies(s) for s in split("C3S H CH")] ;
CSH = CemSpecies(Dict(:C => a, :S => one(a), :H => g))
r = Reaction([CSH, C3S, H, CH]; equal_sign = '→')
## application of a function to stoichimetric coefficients (here simplify)
r = apply(simplify, Reaction([C3S, H], [CH, CSH]; equal_sign = '→'))
SM = StoichMatrix([C3S], [CSH, H, CH]; involve_all_atoms = true) ; pprint(SM)
pprint(apply(x -> simplify.(x), SM))

# Chen & Brouwers
CSH = CemSpecies(Dict(:C => a, :S => one(a), :A => b, :H => g))
HT = CemSpecies("M₅AH₁₃")
HG = CemSpecies("C₆AFS₂H₈")
AFt = CemSpecies("C₆S̄₃H₃₂")
ST = CemSpecies("C₂ASH₈")
AH = CemSpecies("C₄AH₁₃")
CSM = CanonicalStoichMatrix([CSH, HT, HG, AFt, ST, AH]); pprint(CSM)
A, ox = CSM.A, CSM.primaries
A = typeof(a).(A[1:(end - 1), 1:end]) # end-1 to remove the line corresponding to water H
oxides = (CemSpecies ∘ string).(ox[1:(end - 1)])
hydrates = [CSH, HT, HG, AFt, ST, AH]
pprint(A, symbol.(oxides), symbol.(hydrates))
pprint(inv(A), symbol.(hydrates), symbol.(oxides))
Mhyd = ustrip.(us"g/mol", getproperty.(hydrates, :M))
Mox = ustrip.(us"g/mol", getproperty.(oxides, :M))
B = Mox .* A .* inv.(Mhyd)'
# or directly
mCSM = mass_matrix(CSM)
B, ox = mCSM.A, mCSM.primaries
B = B[1:(end - 1), 1:end] # to remove the H line
pprint(B, "m_" .* symbol.(oxides), "m_" .* symbol.(hydrates))
pprint(subs.(inv(B), Ref(Dict(a => 1.8, b => 1, g => 4))), "m_" .* symbol.(hydrates), "m_" .* symbol.(oxides))

# Alkane combustion with SymPy
n = symbols("n", real = true) ; CₙH₂ₙ₊₂ = Species(:C => n, :H => 2n + 2) ;
O₂, H₂O, CO₂ = Species.(split("O₂ H₂O CO₂")) ;
r = Reaction([CₙH₂ₙ₊₂, O₂], [H₂O, CO₂])
apply(factor, r)
pprint(2r)
for vn in 1:9
    print("n=$vn ⇒ "); println(colored(apply(subs, r, n => vn)))
end
println(Reaction([CₙH₂ₙ₊₂, O₂], [H₂O, CO₂]; side = :products))
println(Reaction([CₙH₂ₙ₊₂, O₂], [H₂O, CO₂]; side = :reactants))
@show r[O₂]
