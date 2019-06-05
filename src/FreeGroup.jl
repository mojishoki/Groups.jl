###############################################################################
#
#   FreeSyllable/FreeGroupElem/FreeGroup definition
#
###############################################################################

struct FreeSyllable <: Syllable
   id::Symbol
   pow::Int
end

FreeGroupElem = GroupWord{FreeSyllable}

mutable struct FreeGroup <: AbstractFPGroup
   gens::Vector{FreeSyllable}

   function FreeGroup(gens::Vector{T}) where {T<:Syllable}
      G = new(gens)
      G.gens = gens
      return G
   end
end

export FreeGroupElem, FreeGroup

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

elem_type(::Type{FreeGroup}) = FreeGroupElem

parent_type(::Type{FreeGroupElem}) = FreeGroup

###############################################################################
#
#   FreeSyllable constructors
#
###############################################################################

FreeSyllable(s::Symbol) = FreeSyllable(s,1)
FreeSyllable(s::String) = FreeSyllable(Symbol(s))

FreeGroup(n::Int, symbol::String="f") = FreeGroup([Symbol(symbol,i) for i in 1:n])

FreeGroup(a::Vector) = FreeGroup(FreeSyllable.(a))

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function (G::FreeGroup)()
   id = FreeGroupElem(FreeSyllable[])
   id.parent = G
   return id
end

function (G::FreeGroup)(w::GroupWord{FreeSyllable})
   if length(w) > 0
      for s in w.symbols
         i = findfirst(g -> g.id == s.id, G.gens)
         i == 0 && throw(DomainError(
            "Symbol $s does not belong to $G."))
         s.pow % G.gens[i].pow == 0 || throw(DomainError(
            "Symbol $s doesn't belong to $G."))
      end
   end
   w.parent = G
   return w
end

(G::FreeGroup)(s::FreeSyllable) = G(FreeGroupElem(s))

###############################################################################
#
#   Basic manipulation
#
###############################################################################

hash(s::FreeSyllable, h::UInt) = hash(s.id, hash(s.pow, hash(FreeSyllable, h)))

change_pow(s::FreeSyllable, n::Int) = FreeSyllable(s.id, n)

length(s::FreeSyllable) = abs(s.pow)

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, G::FreeGroup)
   print(io, "Free group on $(length(G.gens)) generators: ")
   join(io, G.gens, ", ")
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

inv(s::FreeSyllable) = change_pow(s, -s.pow)
