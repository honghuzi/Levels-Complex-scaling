using PyPlot
#gr()
Potential(x) = ( 0.5x^2 - 0.8) * exp(-0.19x^2)

ns = 800;
dx = 0.05
x = linspace(-ns/2, ns/2, ns)*dx;
V = map(Potential, x);
plot(x, V, ".-")

θ = 0.2
q = e^(im*θ)

H = spzeros(Complex, ns, ns)
a11 = 5/4dx^2 * q^2;
a12 = -2/3dx^2 * q^2;
a13 = 1/24dx^2 * q^2;

for i in 1:2
    H[i, i] = a11 + V[i]
    H[i+1, i] = a12
    H[i+2, i] = a13
end
H[1, 2] = H[ns, ns-1] = a12
for i in 3:ns-2
    H[i, i] = a11 + V[i]
    H[i+1, i] = H[i-1, i] = a12
    H[i+2, i] = H[i-2, i] = a13
end

for i in ns-1:ns
    H[i, i] = a11 + V[i]
    H[i-1, i] = a12
    H[i-2, i] = a13
end

E, λ = eigs(H)

# plot(real(E), imag(E))
# show()
