export
  cube,
  cross_polytope

function cube(n::Int)
  """-
  Create the n-dimensional cube.

  ## Arguments

  - n: dimension.

  ## Example

  julia> cube(2)
  type: Polytope<Rational>

  POINTS
    1  -1  -1
    1   1  -1
    1  -1   1
    1   1   1
  """
  v = Array{Array}(undef, 2 ^ n)
  v = collect.(Iterators.product([[-1, 1] for j in range(1, n)]...))
  return polytope.Polytope(POINTS = vertices_array(v...))
end

function cross_polytope(n::Int)
  """-
  Create the n-dimensional cross-polytope.

  ## Arguments

  - n: dimension.

  ## Example

  julia> cross_polytope(2)
  type: Polytope<Rational>

  POINTS
    1   1   0
    1   0   1
    1  -1   0
    1   0  -1
  """
  v = Array{Array}(undef, 2 * n)
  v = [zeros(n) for _ in range(1, 2 * n)]
  for j in range(1, n)
    v[j][j] = 1.0
    v[n + j][j] = -1.0
  end
  return polytope.Polytope(POINTS = vertices_array(v...))
end
