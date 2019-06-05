###############################################################################
#
#   FPSyllable/FPGroupElem/FPGroup definition
#
###############################################################################

struct FPSyllable <: Syllable
   id::Symbol
   pow::Int
end

FPGroupElem = GroupWord{FPSyllable}

mutable struct FPGroup <: AbstractFPGroup
   gens::Vector{FPSyllable}
   rels::Dict{FPGroupElem, FPGroupElem}

   function FPGroup(gens::Vector{T}, rels::Dict{FPGroupElem, FPGroupElem}) where {T<:Syllable}
      G = new(gens)
      G.rels = Dict(G(k) => G(v) for (k,v) in rels)
      return G
   end
end

export FPGroupElem, FPGroup

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

parent_type(::Type{FPGroupElem}) = FPGroup

elem_type(::FPGroup) = FPGroupElem

###############################################################################
#
#   FPSyllable constructors
#
###############################################################################

FPSyllable(s::Symbol) = FPSyllable(s, 1)
FPSyllable(s::String) = FPSyllable(Symbol(s))

convert(::Type{FPSyllable}, s::FreeSyllable) = FPSyllable(s.id, s.pow)

FPGroup(gens::Vector{FPSyllable}) = FPGroup(gens, Dict{FPGroupElem, FPGroupElem}())

FPGroup(a::Vector{String}) = FPGroup([FPSyllable(i) for i in a])

FPGroup(n::Int, symbol::String="f") = FPGroup(["$symbol$i" for i in 1:n])
FPGroup(H::FreeGroup) = FPGroup([FPSyllable(s) for s in H.gens])

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function (G::FPGroup)()
   id = FPGroupElem(FPSyllable[])
   id.parent = G
   return id
end

function (G::FPGroup)(w::GWord)
   if length(w) == 0
      return G()
   end

   if eltype(w.symbols) == FreeSyllable
      w = FPGroupElem(w.symbols)
   end

   if eltype(w.symbols) == FPSyllable
      for s in w.symbols
         i = findfirst(g -> g.id == s.id, G.gens)
         i == nothing && throw(DomainError(
            "Symbol $s does not belong to $G."))
         s.pow % G.gens[i].pow == 0 || throw(DomainError(
         "Symbol $s doesn't belong to $G."))
      end
   end
   w.parent = G
   return reduce!(w)
end

(G::FPGroup)(s::FPSyllable) = G(FPGroupElem(s))

###############################################################################
#
#   Basic manipulation
#
###############################################################################

hash(s::FPSyllable, h::UInt) = hash(s.id, hash(s.pow, hash(FPSyllable, h)))

change_pow(s::FPSyllable, n::Int) = FPSyllable(s.id, n)

length(s::FPSyllable) = abs(s.pow)

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, G::FPGroup)
   print(io, "FPgroup on $(length(G.gens)) generators ")
   strrels = join(G.rels, ", ")
   if length(strrels) > 300
      print(io, "⟨ ", join(G.gens, ", "), " | $(length(G.rels)) relation(s) ⟩.")
   else
      print(io, "⟨ ", join(G.gens, ", "), " | ", join(G.rels, ", "), " ⟩.")
   end
end

###############################################################################
#
#   Comparison
#
###############################################################################

###############################################################################
#
#   Inversion
#
###############################################################################

inv(s::FPSyllable) = change_pow(s, -s.pow)

###############################################################################
#
#   Binary operations
#
###############################################################################

(*)(W::FPGroupElem, Z::FPGroupElem) = r_multiply(W, Z.symbols)
(*)(W::FPGroupElem, s::FPSyllable) = r_multiply(W, [s])
(*)(s::FPSyllable, W::FPGroupElem) = l_multiply(W, [s])

function reduce!(W::FPGroupElem)
    if length(W) < 2
        deleteat!(W.symbols, findall(x -> x.pow == 0, W.symbols))
    else
        reduced = false
        while !reduced
            reduced = freereduce!(W) || replace_all!(W, parent(W).rels)
        end
    end

    W.savedhash = hash(W.symbols, hash(typeof(W)))
    W.modified = false
    return W
end

###############################################################################
#
#   Misc
#
###############################################################################

function add_rels!(G::FPGroup, newrels::Dict{FPGroupElem,FPGroupElem})
   for w in keys(newrels)
      if !(w in keys(G.rels))
         G.rels[w] = G(newrels[w])
      end
   end
end

function /(G::FPGroup, newrels::Vector{FPGroupElem})
   for r in rels
      parent(r) == G || throw(DomainError(
      "Can not form quotient group: $r is not an element of $G"))
   end
   H = deepcopy(G)
   newrels = Dict(H(r) => H() for r in newrels)
   add_rels!(H, newrels)
   return H
end

function /(G::FreeGroup, rels::Vector{FreeGroupElem})
   for r in rels
      parent(r) == G || throw(DomainError(
         "Can not form quotient group: $r is not an element of $G"))
   end
   H = FPGroup(G)
   H.rels = Dict(H(rel) => H() for rel in unique(rels))
   return H
end
