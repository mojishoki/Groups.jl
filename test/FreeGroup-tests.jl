
@testset "Groups.FreeSyllables" begin
   s = Groups.FreeSyllable(:s)
   t = Groups.FreeSyllable(:t)

   @testset "constructors" begin
      @test isa(Groups.FreeSyllable(:aaaaaaaaaaaaaaaa), Groups.Syllable)
      @test Groups.FreeSyllable(:abc).pow == 1
      @test isa(s, Groups.FreeSyllable)
      @test isa(t, Groups.FreeSyllable)
   end
   @testset "eltary functions" begin
      @test length(s) == 1
      @test Groups.change_pow(s, 0) == Groups.change_pow(t, 0)
      @test length(Groups.change_pow(s, 0)) == 0
      @test inv(s).pow == -1
      @test Groups.FreeSyllable(:s, 3) == Groups.change_pow(s, 3)
      @test Groups.FreeSyllable(:s, 3) != Groups.FreeSyllable(:t, 3)
      @test Groups.change_pow(inv(s), -3) == inv(Groups.change_pow(s, 3))
   end
   @testset "powers" begin
      s⁴ = Groups.change_pow(s,4)
      @test s⁴.pow == 4
      @test Groups.change_pow(s, 4) == Groups.FreeSyllable(:s, 4)
   end
end

@testset "FreeGroupSyllables manipulation" begin
   s = Groups.FreeSyllable("s")
   t = Groups.FreeSyllable(:t, -2)

   @test isa(Groups.GroupWord(s), Groups.GWord{Groups.FreeSyllable})
   @test isa(Groups.GroupWord(s), FreeGroupElem)
   @test isa(FreeGroupElem(s), Groups.GWord)
   @test isa(convert(FreeGroupElem, s), Groups.GWord)
   @test isa(convert(FreeGroupElem, s), FreeGroupElem)
   @test isa(Vector{FreeGroupElem}([s,t]), Vector{FreeGroupElem})
   @test length(FreeGroupElem(s)) == 1
   @test length(FreeGroupElem(t)) == 2

end

@testset "FreeGroup" begin
   @test isa(FreeGroup(["s", "t"]), AbstractAlgebra.Group)
   G = FreeGroup(["s", "t"])

   @testset "elements constructors" begin
      @test isa(G(), FreeGroupElem)
      @test eltype(G.gens) == Groups.FreeSyllable
      @test length(G.gens) == 2
      @test eltype(gens(G)) == FreeGroupElem
      @test length(gens(G)) == 2
   end

   s, t = gens(G)

   @testset "internal arithmetic" begin

      @test Vector{Groups.FreeGroupElem}([s,t]) == [Groups.GroupWord(s), Groups.GroupWord(t)]
      @test (s*s).symbols == (s^2).symbols
      @test hash([t^1,s^1]) == hash([t^2*inv(t),s*inv(s)*s])

      t_symb = Groups.FreeSyllable(:t)
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
      @test length(G().symbols) == 0
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
      @test Groups.reduce!(w).symbols ==Vector{Groups.FreeSyllable}([])
   end

   @testset "Group operations" begin
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
      a = Groups.FreeSyllable(:a)
      b = Groups.FreeSyllable(:b)
      @test Groups.issubsymbol(a, Groups.change_pow(a,2)) == true
      @test Groups.issubsymbol(a, Groups.change_pow(a,-2)) == false
      @test Groups.issubsymbol(b, Groups.change_pow(a,-2)) == false
      @test Groups.issubsymbol(inv(b), Groups.change_pow(b,-2)) == true

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
      @test Groups.replace_all(s*c*s*c*s, subst) == s*t^4*s*t^4*s

      G = FreeGroup(["x", "y"])
      x,y = gens(G)

      @test Groups.replace(x*y^9, 2, y^2, y) == x*y^8
      @test Groups.replace(x^3, 1, x^2, y) == x*y
      @test Groups.replace(y*x^3*y, 2, x^2, y) == y*x*y^2
   end
end
