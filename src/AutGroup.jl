###############################################################################
#
#   AutSymbol/ AutGroup / Automorphism
#
###############################################################################

struct RTransvect
    i::Int8
    j::Int8
end

struct LTransvect
    i::Int8
    j::Int8
end

struct FlipAut
    i::Int8
end

struct PermAut
    perm::Generic.perm{Int8}
end

struct Identity end

struct AutSymbol <: GSymbol
   id::Symbol
   pow::Int8
   fn::Union{LTransvect, RTransvect, PermAut, FlipAut, Identity}
end

mutable struct AutGroup{N} <: AbstractFPGroup
   objectGroup::FreeGroup
   gens::Vector{AutSymbol}
end

mutable struct Automorphism{N} <: GWord{AutSymbol}
    symbols::Vector{AutSymbol}
    modified::Bool
    savedhash::UInt
    parent::AutGroup{N}

    function Automorphism{N}(f::Vector{AutSymbol}) where {N}
        return new{N}(f, true, zero(UInt))
    end
end

export Automorphism, AutGroup, Aut, SAut

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

elem_type(::AutGroup{N}) where N = Automorphism{N}

parent_type(::Automorphism{N}) where N = AutGroup{N}

###############################################################################
#
#   AutSymbol defining functions
#
###############################################################################

function (ϱ::RTransvect)(v, pow::Integer=1)
    @inbounds Groups.r_multiply!(v[ϱ.i], (v[ϱ.j]^pow).symbols, reduced=false)
    return v
end

function (λ::LTransvect)(v, pow::Integer=1)
    @inbounds Groups.l_multiply!(v[λ.i], (v[λ.j]^pow).symbols, reduced=false)
    return v
end

function (σ::PermAut)(v, pow::Integer=1)
   w = deepcopy(v)
   if pow == 1
       @inbounds for k in eachindex(v)
           v[k].symbols = w[σ.perm.d[k]].symbols
       end
   else
       s = (σ.perm^pow).d
       @inbounds for k in eachindex(v)
           v[k].symbols = w[s[k]].symbols
       end
   end
   return v
end

function (ɛ::FlipAut)(v, pow::Integer=1)
   @inbounds if isodd(pow)
       v[ɛ.i].symbols = inv(v[ɛ.i]).symbols
   end
   return v
end

(::Identity)(v, pow::Integer=1) = v

# taken from ValidatedNumerics, under under the MIT "Expat" License:
# https://github.com/JuliaIntervals/ValidatedNumerics.jl/blob/master/LICENSE.md
function subscriptify(n::Integer)
    subscript_0 = Int(0x2080) # Char(0x2080) -> subscript 0
    @assert 0 <= n <= 9
    return Char(subscript_0 + n)
    # return [Char(subscript_0 + i) for i in reverse(digits(n))])
end

function id_autsymbol()
   return AutSymbol(Symbol("(id)"), 0, Identity())
end

function rmul_autsymbol(i::Integer, j::Integer; pow::Integer=1)
    id = Symbol("ϱ", subscriptify(i), subscriptify(j))
    return AutSymbol(id, pow, RTransvect(i, j))
end

function lmul_autsymbol(i::Integer, j::Integer; pow::Integer=1)
    id = Symbol("λ", subscriptify(i), subscriptify(j))
    return AutSymbol(id, pow, LTransvect(i, j))
end

function flip_autsymbol(i::Integer; pow::Integer=1)
    if iseven(pow)
       return id_autsymbol()
    else
        id = Symbol("ɛ", subscriptify(i))
        return AutSymbol(id, 1, FlipAut(i))
    end
end

function perm_autsymbol(p::Generic.perm{I}; pow::Integer=one(I)) where I<:Integer
    if pow != 1
        p = p^pow
    end
    for i in eachindex(p.d)
        if p.d[i] != i
            id = Symbol("σ", [subscriptify(i) for i in p.d]...)
            return AutSymbol(id, 1, PermAut(p))
        end
    end
    return id_autsymbol()
end

function perm_autsymbol(a::Vector{T}) where T<:Integer
   return perm_autsymbol(perm(Vector{Int8}(a), false))
end

function domain(G::AutGroup{N}) where N
    F = G.objectGroup
    gg = gens(F)
    return ntuple(i->gg[i], Val(N))
end

###############################################################################
#
#   AutGroup / Automorphism constructors
#
###############################################################################

function AutGroup(G::FreeGroup; special=false)
   S = AutSymbol[]
   n = length(gens(G))
   n == 0 && return AutGroup{n}(G, S)

   indexing = [[i,j] for i in 1:n for j in 1:n if i≠j]

   rmuls = [rmul_autsymbol(i,j) for (i,j) in indexing]
   lmuls = [lmul_autsymbol(i,j) for (i,j) in indexing]

   append!(S, [rmuls; lmuls])

   if !special
      flips = [flip_autsymbol(i) for i in 1:n]
      syms = [perm_autsymbol(p) for p in PermutationGroup(n)][2:end]

      append!(S, [flips; syms])

   end
   return AutGroup{n}(G, S)
end

Aut(G::Group) = AutGroup(G)
SAut(G::Group) = AutGroup(G, special=true)

###############################################################################
#
#   Types call overloads
#
###############################################################################

Automorphism{N}(s::AutSymbol) where N = Automorphism{N}(AutSymbol[s])

function (G::AutGroup{N})() where N
   id = Automorphism{N}(id_autsymbol())
   id.parent = G
   return id
end

function (G::AutGroup{N})(f::AutSymbol) where N
   g = Automorphism{N}([f])
   g.parent = G
   return g
end

function (G::AutGroup{N})(g::Automorphism{N}) where N
   g.parent = G
   return g
end

###############################################################################
#
#   Functional call overloads for evaluation of AutSymbol and Automorphism
#
###############################################################################

function (f::AutSymbol)(v::NTuple{N, T}) where {N, T}
   if f.pow != 0
      v = f.fn(v, f.pow)::NTuple{N, T}
   end
   return v
end

function (F::Automorphism{N})(v::NTuple{N, T}) where {N, T}
    for f in reverse(F.symbols)
        v = f(v)::NTuple{N, T}
    end
    return v
end

###############################################################################
#
#   Comparison
#
###############################################################################

const HASHINGCONST = 0x7d28276b01874b19 # hash(Automorphism)

hash(s::AutSymbol, h::UInt) = hash(s.id, hash(s.pow, hash(:AutSymbol, h)))

function hash(g::Automorphism, h::UInt)
    if g.modified
        g_im = reduce!.(g(domain(parent(g))))
        g.savedhash = hash(g_im, hash(typeof(g), hash(parent(g), HASHINGCONST)))
        g.modified = false
        end
    return xor(g.savedhash, h)
end

function (==)(g::Automorphism{N}, h::Automorphism{N}) where N
    parent(g) == parent(h) || return false

    if !g.modified && !h.modified
        if g.savedhash != h.savedhash
            return false
        end
    end

    # expensive:
    g_im = reduce!.(g(domain(parent(g))))
    h_im = reduce!.(h(domain(parent(h))))
    # cheap:
    g.savedhash = hash(g_im, hash(typeof(g), hash(parent(g), HASHINGCONST)))
    g.modified = false
    h.savedhash = hash(h_im, hash(typeof(h), hash(parent(h), HASHINGCONST)))
    h.modified = false

    return g_im == h_im
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function change_pow(s::AutSymbol, n::Integer)
    if n == zero(n)
        return id_autsymbol()
    end
    symbol = s.fn
    if symbol isa FlipAut
        return flip_autsymbol(symbol.i, pow=n)
    elseif symbol isa PermAut
        return perm_autsymbol(symbol.perm, pow=n)
    elseif symbol isa RTransvect
        return rmul_autsymbol(symbol.i, symbol.j, pow=n)
    elseif symbol isa LTransvect
        return lmul_autsymbol(symbol.i, symbol.j, pow=n)
    elseif symbol isa Identity
        return s
    else
        warn("Changing power of an unknown type of symbol! $s")
        return AutSymbol(s.id, n, s.fn)
    end
end

length(s::AutSymbol) = abs(s.pow)

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, G::AutGroup)
   print(io, "Automorphism Group of $(G.objectGroup)\n")
   print(io, "Generated by $(join(G.gens, ","))")
end

###############################################################################
#
#   Binary operators
#
###############################################################################

###############################################################################
#
#   Inversion
#
###############################################################################

inv(f::AutSymbol) = change_pow(f, -f.pow)

###############################################################################
#
#   Misc
#
###############################################################################

function getperm(s::AutSymbol)
    if s.pow != 1
        @warn("Power for perm_symbol should be never 0!")
        return s.fn.perm^s.pow
    else
        return s.fn.perm
    end
end

function simplifyperms!(W::Automorphism{N}) where N
    reduced = true
    to_delete = Int[]
    for i in 1:length(W.symbols)-1
        if W.symbols[i].pow == 0
            continue
        elseif W.symbols[i].fn isa PermAut && W.symbols[i+1].fn isa PermAut
            reduced = false
            c = W.symbols[i]
            n = W.symbols[i+1]
            W.symbols[i+1] = perm_autsymbol(getperm(c)*getperm(n))
            push!(to_delete, i)
        end
    end
    deleteat!(W.symbols, to_delete)
    deleteids!(W)
    return reduced
end

function reduce!(W::Automorphism)
    if length(W) == 0
        return W
    elseif length(W.symbols) == 1
        deleteids!(W)
    else
        reduced = false
        while !reduced
            reduced = simplifyperms!(W) && freereduce!(W)
        end
    end

    W.modified = true

    return W
end

function linear_repr(A::Automorphism{N}, hom=matrix_repr) where N
    return reduce(*, linear_repr.(A.symbols, N, hom), init=hom(Identity(),N,1))
end

linear_repr(a::AutSymbol, n::Int, hom) = hom(a.fn, n, a.pow)

function matrix_repr(a::Union{RTransvect, LTransvect}, n::Int, pow)
    x = Matrix{Int}(I, n, n)
    x[a.i,a.j] = pow
    return x
end

function matrix_repr(a::FlipAut, n::Int, pow)
    x = Matrix{Int}(I, n, n)
    x[a.i,a.i] = -1^pow
    return x
end

matrix_repr(a::PermAut, n::Int, pow) = Matrix{Int}(I, n, n)[(a.perm^pow).d, :]

matrix_repr(a::Identity, n::Int, pow) = Matrix{Int}(I, n, n)
