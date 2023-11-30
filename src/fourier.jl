export
  vertices_array,
  fourier_transform

function vertices_array(v...)
  """-
  Compute an array in which each row is an element of v in projective coordinates.

  ## Arguments

  - v: list with generating points of a polytope.

  ## Example
  
  julia> vertices_array([0, 0], [1, 0], [0, 1])
  3x3 Matrix{Int64}:
   1  0  0
   1  1  0
   1  0  1
  """
  v0 = permutedims(hcat(v...))
  v1 = hcat(ones(size(v0, 1)), v0)
  return v1
end

function tangent_cone(j, P)
  """-
  Compute generating vectors of the tangent cone at the j-th vertice of P.

  ## Arguments

  - j: index of a vertex.
  - P: polytope.

  ## Example
  
  julia> tangent_cone(1, P)
  2-element Vector{Vector{Polymake.LibPolymake.RationalAllocated}}:
   [1, 0]
   [0, 1]
  """
  vertices = eachslice(P.VERTICES[:, 2:end], dims = 1)
  v = vertices[j]
  adjacent = polytope.vertex_figure(P, j - 1).VERTEX_LABELS
  cone = [vertices[parse(Int, l) + 1] - v for l in adjacent]
  return collect(cone)
end

function triangulate(c)
  """-
  Triangulate a cone given by generators c.

  ## Arguments

  - c: list with generators of a cone.

  ## Example

  julia> triangulate([[-1, 0, 1], [1, 0, 1], [0, -1, 1], [0, 1, 1]])
  2-element Vector{Vector{Int64}}:
   [0, 1, 2]
   [0, 1, 3]
  """
  object_cone = polytope.Cone(INPUT_RAYS = permutedims(hcat(c...)))
  triangulation = collect.(object_cone.TRIANGULATION.FACETS)
  return triangulation
end

function triangulation_determinant(w, t)
  """-
  Compute the determinant of the lattice with generators w in a triangle t.

  ## Arguments

  - w: generators of the lattice.
  - t: list of lists indices of w containing the triangulation.

  ## Example

  julia> w = [[-1, 0, 1], [1, 0, 1], [0, -1, 1], [0, 1, 1]];
  julia> triangulation_determinant(w, triangulate(w)[1])
  2.0
  """
  return abs(LinearAlgebra.det( Float64.(stack(w[1 .+ t])) ))
end

function fourier_transform(P)
  """-
  Compute the Fourier transform of polytope P based on Brion's formula.

  ## Arguments

  - P: polytope.

  ## Example

  julia> P = polytope.Polytope(POINTS = vertices_array([0, 0], [1, 0], [0, 1]));
  julia> fourier_transform(P)([2.5, 2.0])
  0.04052 - 2.48164e-17im
  """
  tangent_cones = tangent_cone.(range(1, P.N_VERTICES), Ref(P))
  triangulation = triangulate.(tangent_cones)
  vertices = eachslice(P.VERTICES[:, 2:end], dims = 1)
  dets = [triangulation_determinant.(Ref(c), t) for (c, t) in zip(tangent_cones, triangulation)]
  d = length(P.VERTICES[1, :]) - 1
  factor = 1.0 / ( 2.0 * pi * 1im ) ^ d
  exponentials(xi) = [exp(-2 * pi * 1im * LinearAlgebra.dot(v, xi)) for v in vertices]
  prods(tg, tr, xi) = [prod(LinearAlgebra.dot.(tg[1 .+ t], Ref(xi))) for t in tr]
  total_prods(xi) = prods.(tangent_cones, triangulation, Ref(xi))
  division(xi) = ((x, y) -> x ./ y).( dets, total_prods(xi) )
  function FP(xi)
    return factor * sum(exponentials(xi) .* sum.(division(xi)))
  end
  return FP
end
