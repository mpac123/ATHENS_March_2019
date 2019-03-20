using Plots; gr()

### Representation of the mesh
# assignment 5: mesh representation [p,e]
function mesh(x::Array{T,1}) where T <: Float64
    p = Array{Float64}[]
    e = Array{Integer}[]
    for i = 1:(length(x) - 1)
        push!(p, [x[i], x[i+1]])
        push!(e, [i, i+1])
    end
    p, e
end

### Computation of Element-by-Element Contributions
# assignment 6: funcion f(X)
f1(x) = sin(pi*x)
f2(x) = 2*sin(pi*x) + 4 * pi * (x-1) - pi*pi*(x-1)*(x-1)*sin(pi*x)
u(x) = 1/pi/pi*sin(pi*x)+x

# assignment 7: generate 2x2 element matrix A_ei
function element_matrix(xi, xi_1)
    c = 1/(xi_1 - xi)
    A = Array{Float64}[]
    push!(A, [c, -c])
    push!(A, [-c, c])
    A 
end

# asignment 8: generate 2x1 element vector f_ei
function element_vector(xi, xi_1, func) 
    c = (xi_1 - xi) / 2
    f = []
    push!(f, c * func(xi))
    push!(f, c * func(xi_1))
    f
end

### Assembly of Global Linear System ###
# assignment 9: assemble all the element matrices A_ei into the global matrix
# (n+1)x(n+1) A
function global_matrix(p_e)
    p = collect(p_e[1]); e=collect(p_e[2])
    n = length(p)
    A = zeros(n+1, n+1)
    for i = 1:n
        A_ei = element_matrix(p[i][1], p[i][2])
        for j = 1:2, k= 1:2
            A[e[i][j], e[i][k]] = A[e[i][j], e[i][k]] + A_ei[j][k]
        end
    end
    A
end

# assignment 10: assemble all the element vectors f_ei into the global
# matrix (n + 1) x 1 vector f
function global_vector(p_e, func)
    p = p_e[1]; e=p_e[2]
    n = length(p)
    f = zeros(n+1, 1)
    for i = 1:n
        f_ei = element_vector(p[i][1], p[i][2], func)
        for j = 1:2, k = 1:2
            f[e[i][j]] = f[e[i][j]] + f_ei[j]
        end
    end
    f
end

### Trearment of Boundary Conditions ###
# assignment 11: modify the first equation of the linear system
# A_u^h = f in such a way that the finite element solution u^h(x) satisfies the 
# Dirichlet boundary conditions at x = 0
function dirichlet_BC(p_e, func) 
    A = collect(global_matrix(p_e))
    f = collect(global_vector(p_e, func))
    A[1,1] = 1.0
    A[1,2] = 0.0
    f[1] = 0
    A, f
end

data = sort(rand(100))
mesh(data)
func = (x) -> f1(x)
element_vector(3,5,func)

global_matrix(mesh(data))
global_vector(mesh(data), func)
dirichlet_BC(mesh(data), func)

### Solving the Global Linear System ###
# assignment 12: visualise matrix A using Gadfly
# Gadfly.spy(dirichlet_BC(mesh(collect(1:0.1:10)), func)[1])

# assignment 13: compute the finite element solution
function FEM_solution(p_e, func)
    A, f = dirichlet_BC(p_e, func)
    u_h = A \ f
end

x = data
y = FEM_solution(mesh(data), func)
plotly()
plot(x,y)
plot!(u, 0, 1)
