@testset "Groups.FreeSymbols" begin
   s = Groups.FreeSymbol("s")
   t = Groups.FreeSymbol("t")
   @testset "defines" begin
      @test isa(Groups.FreeSymbol("aaaaaaaaaaaaaaaa"), Groups.GSymbol)
      @test Groups.FreeSymbol("abc").pow == 1
      @test isa(s, Groups.FreeSymbol)
      @test isa(t, Groups.FreeSymbol)
   end
   @testset "eltary functions" begin
      @test length(s) == 1
      @test Groups.change_pow(s, 0) == Groups.change_pow(t, 0)
      @test length(Groups.change_pow(s, 0)) == 0
      @test inv(s).pow == -1
      @test Groups.FreeSymbol("s", 3) == Groups.change_pow(s, 3)
      @test Groups.FreeSymbol("s", 3) != Groups.FreeSymbol("t", 3)
      @test Groups.change_pow(inv(s), -3) == inv(Groups.change_pow(s, 3))
   end
   @testset "powers" begin
      s⁴ = Groups.change_pow(s,4)
      @test s⁴.pow == 4
      @test Groups.change_pow(s, 4) == Groups.FreeSymbol("s", 4)
   end
end

@testset "FreeGroupElems" begin
   s = Groups.FreeSymbol("s")
   t = Groups.FreeSymbol("t", -2)
   @testset "defines" begin
      @test isa(Groups.GWord(s), Groups.GWord)
      @test isa(Groups.GWord(s), FreeGroupElem)
      @test isa(FreeGroupElem(s), Groups.GWord)
      @test isa(convert(FreeGroupElem, s), Groups.GWord)
      @test isa(convert(FreeGroupElem, s), FreeGroupElem)
      @test isa(Vector{FreeGroupElem}([s,t]), Vector{FreeGroupElem})
      @test length(FreeGroupElem(s)) == 1
      @test length(FreeGroupElem(t)) == 2
   end

   @testset "eltary functions" begin
      @test_skip (s*s).symbols == (s^2).symbols
      @test_skip Vector{Groups.GWord{Groups.FreeSymbol}}([s,t]) ==
         Vector{FreeGroupElem}([s,t])
      @test_skip Vector{Groups.GWord}([s,t]) ==
         [Groups.GWord(s), Groups.GWord(t)]
      @test_skip hash([t^1,s^1]) == hash([t^2*inv(t),s*inv(s)*s])
   end
end

@testset "FreeGroup" begin
   @test isa(FreeGroup(["s", "t"]), Nemo.Group)
   G = FreeGroup(["s", "t"])

   @testset "elements constructors" begin
      @test isa(G(), FreeGroupElem)
      @test eltype(G.gens) == Groups.FreeSymbol
      @test length(G.gens) == 2
      @test_skip eltype(G.rels) == FreeGroupElem
      @test_skip length(G.rels) == 0
      @test eltype(generators(G)) == FreeGroupElem
      @test length(generators(G)) == 2
   end

   s, t = generators(G)

   @testset "internal arithmetic" begin
      t_symb = Groups.FreeSymbol("t")
      tt = deepcopy(t)
      @test string(Groups.r_multiply!(tt,[inv(t_symb)]; reduced=true)) ==
         "(id)"
      tt = deepcopy(t)
      @test string(Groups.r_multiply!(tt,[inv(t_symb)]; reduced=false)) ==
         "t*t^-1"
      tt = deepcopy(t)
      @test string(Groups.l_multiply!(tt,[inv(t_symb)]; reduced=true)) ==
         "(id)"
      tt = deepcopy(t)
      @test string(Groups.l_multiply!(tt,[inv(t_symb)]; reduced=false)) ==
         "t^-1*t"
   end

   @testset "reductions" begin
      @test length(G().symbols) == 1
      @test length((G()*G()).symbols) == 0
      @test G() == G()*G()
      w = deepcopy(s)
      push!(w.symbols, (s^-1).symbols[1])
      @test Groups.reduce!(w) == parent(w)()
      o = (t*s)^3
      @test o == t*s*t*s*t*s
      p = (t*s)^-3
      @test p == s^-1*t^-1*s^-1*t^-1*s^-1*t^-1
      @test o*p == parent(o*p)()
      w = FreeGroupElem([o.symbols..., p.symbols...])
      w.parent = G
      @test Groups.reduce!(w).symbols ==Vector{Groups.FreeSymbol}([])
   end

   @testset "binary/inv operations" begin
      @test parent(s) == G
      @test parent(s) === parent(deepcopy(s))
      @test isa(s*t, FreeGroupElem)
      @test parent(s*t) == parent(s^2)
      @test s*s == s^2
      @test inv(s*s) == inv(s^2)
      @test inv(s)^2 == inv(s^2)
      @test inv(s)*inv(s) == inv(s^2)
      @test inv(s*t) == inv(t)*inv(s)
      w = s*t*s^-1
      @test inv(w) == s*t^-1*s^-1
      @test (t*s*t^-1)^10 == t*s^10*t^-1
      @test (t*s*t^-1)^-10 == t*s^-10*t^-1
   end

   @testset "replacements" begin
      a = Groups.FreeSymbol("a")
      b = Groups.FreeSymbol("b")
      @test Groups.is_subsymbol(a, Groups.change_pow(a,2)) == true
      @test Groups.is_subsymbol(a, Groups.change_pow(a,-2)) == false
      @test Groups.is_subsymbol(b, Groups.change_pow(a,-2)) == false
      @test Groups.is_subsymbol(inv(b), Groups.change_pow(b,-2)) == true
      c = s*t*s^-1*t^-1
      @test findfirst(c, s^-1*t^-1) == 3
      @test findnext(c*s^-1, s^-1*t^-1,3) == 3
      @test findnext(c*s^-1*t^-1, s^-1*t^-1,4) == 5
      @test findfirst(c*t, c) == 0
      w = s*t*s^-1
      subst = Dict{FreeGroupElem, FreeGroupElem}(w => s^1, s*t^-1 => t^4)
      @test Groups.replace(c, 1, s*t, G()) == s^-1*t^-1
      @test Groups.replace(c, 1, w, subst[w]) == s*t^-1
      @test Groups.replace(s*c*t^-1, 1, w, subst[w]) == s^2*t^-2
      @test Groups.replace(t*c*t, 2, w, subst[w]) == t*s
      @test Groups.replace_all!(s*c*s*c*s, subst) == s*t^4*s*t^4*s
   end
end